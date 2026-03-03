# ── Test set inference for BacDive multihead model ───────────────────────────
# Combined data generation + inference in a single pass (no intermediate RDS)
# Usage: conda run -n tf Rscript eval_test.R

library(deepG)
library(keras)
library(tensorflow)
library(dplyr)

# ── Inline sample weight functions (avoids sourcing data_prep.R) ─────────────
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg)) else normalizePath(getwd())
proj_dir <- if (dir.exists(script_path)) script_path else dirname(script_path)
source(file.path(proj_dir, "R", "config_paths.R"))
cfg <- load_bacdive_config(
  proj_dir,
  required = c("count_list_rds", "fasta_dir", "bd_splits_csv", "bd_labels_csv", "checkpoints_dir")
)

count_list <- readRDS(cfg$count_list_rds)
cw_list <- list()
for (i in names(count_list)) {
  v <- count_list[[i]]
  cw_list[[i]] <- sum(v) / (v * length(v))
}

set_sw_cellsize <- function(y, bs, na_enc = -9999) {
  sw <- matrix(1, nrow = bs)
  sw[apply(y == na_enc, 1, any)] <- 0
  sw
}
set_sw_class <- function(y, bs, na_enc = -9999, cw) {
  sw <- matrix(1, nrow = bs)
  na_row <- which(apply(y == na_enc, 1, any) | rowSums(y) == 0)
  sw[na_row, ] <- 0
  non_na <- setdiff(1:bs, na_row)
  sw[non_na, ] <- cw[apply(y[non_na, , drop = FALSE], 1, which.max)]
  sw
}
set_sw_binary <- function(y, bs, na_enc = -9999, cw) {
  sw <- matrix(1, nrow = bs)
  na_row <- which(apply(y == na_enc, 1, any))
  sw[na_row, ] <- 0
  non_na <- setdiff(1:bs, na_row)
  sw[non_na, ] <- cw[ifelse(y[non_na, 1] == 0, 1, 2)]
  sw
}
set_sw_patho <- function(y, bs, cw, default_weight = 1e-12) {
  cw[1] <- cw[1] * default_weight
  matrix(cw[y[, 1] + 1], ncol = 1)
}
set_sw_biosafety <- function(y, bs, na_enc = -9999, cw) {
  sw <- matrix(1, nrow = bs)
  na_row <- which(apply(y == na_enc, 1, any))
  sw[na_row, ] <- 0
  non_na <- setdiff(1:bs, na_row)
  sw[non_na, ] <- cw[y[non_na, 1]]
  sw
}

set_sw <- function(y_list, na_enc = -9999, cw_list) {
  bs <- nrow(y_list[[1]])
  list(
    set_sw_cellsize(y_list[[1]], bs, na_enc),
    set_sw_class(y_list[[2]], bs, na_enc, cw_list$cellshape),
    set_sw_class(y_list[[3]], bs, na_enc, cw_list$flagellum),
    set_sw_class(y_list[[4]], bs, na_enc, cw_list$gram),
    set_sw_binary(y_list[[5]], bs, na_enc, cw_list$motility),
    set_sw_biosafety(y_list[[6]], bs, na_enc, cw_list$biosafety),
    set_sw_patho(y_list[[7]], bs, cw_list$pathogenicity_human),
    set_sw_patho(y_list[[8]], bs, cw_list$pathogenicity_animal),
    set_sw_patho(y_list[[9]], bs, cw_list$pathogenicity_plant),
    set_sw_class(y_list[[10]], bs, na_enc, cw_list$oxygen_growth),
    set_sw_class(y_list[[11]], bs, na_enc, cw_list$oxygen_facultative),
    set_sw_class(y_list[[12]], bs, na_enc, cw_list$oxygen_obligate),
    set_sw_binary(y_list[[13]], bs, na_enc, cw_list$oxygen_microaerophile),
    set_sw_binary(y_list[[14]], bs, na_enc, cw_list$spore)
  )
}

# ── Settings (matching validation setup) ─────────────────────────────────────
ttv <- "test"
maxlen <- 2000000
maxlen_fragment <- 10000
samples_per_target <- maxlen / maxlen_fragment  # 200
num_batches <- 800
batch_size <- 12
step <- maxlen
max_samples <- floor(batch_size / 6)  # 2
seed <- 1 * 77
vocabulary <- c("A", "C", "G", "T")

target_split <- list(
  c("rgStart_len", "rgEnd_len", "rgStart_wid", "rgEnd_wid"),
  c("is_cell_shape_rod.shaped", "is_cell_shape_other",
    "is_cell_shape_coccus.shaped", "is_cell_shape_vibrio.shaped",
    "is_cell_shape_filament.shaped", "is_cell_shape_sphere.shaped",
    "is_cell_shape_ovoid.shaped", "is_cell_shape_pleomorphic.shaped",
    "is_cell_shape_spiral.shaped", "is_cell_shape_curved.shaped",
    "is_cell_shape_oval.shaped"),
  c("is_flagellum_arrangement_monotrichous",
    "is_flagellum_arrangement_monotrichous_polar",
    "is_flagellum_arrangement_polar",
    "is_flagellum_arrangement_peritrichous",
    "is_flagellum_arrangement_lophotrichous",
    "is_flagellum_arrangement_gliding"),
  c("is_gram_stain_positive", "is_gram_stain_variable",
    "is_gram_stain_negative"),
  c("is_motile"), c("biosafety_level"),
  c("pathogenicity_human"), c("pathogenicity_animal"),
  c("pathogenicity_plant"),
  c("aerobe", "anaerobe"),
  c("facultative.aerobe", "facultative.anaerobe"),
  c("obligate.aerobe", "obligate.anaerobe"),
  c("microaerophile"), c("ability_spore")
)

# Named version for output columns
target_split_named <- list(
  cellsize = target_split[[1]], cellshape = target_split[[2]],
  flagellum = target_split[[3]], gram = target_split[[4]],
  motility = target_split[[5]], biosafety = target_split[[6]],
  pathogenicity_human = target_split[[7]],
  pathogenicity_animal = target_split[[8]],
  pathogenicity_plant = target_split[[9]],
  oxygen_growth = target_split[[10]],
  oxygen_facultative = target_split[[11]],
  oxygen_obligate = target_split[[12]],
  oxygen_microaerophile = target_split[[13]],
  spore = target_split[[14]]
)

# ── Load test FASTA file list ────────────────────────────────────────────────
all_files <- list.files(cfg$fasta_dir, full.names = TRUE)
all_files <- all_files[file.info(all_files)$size != 0]
all_csv <- read.csv(cfg$bd_splits_csv)
type_csv <- all_csv %>% filter(type == ttv) %>% select(file) %>% unlist()
fasta_files <- as.list(all_files[basename(all_files) %in% type_csv])
cat("Test FASTA files found:", length(fasta_files), "\n"); flush(stdout())

target_from_csv <- cfg$bd_labels_csv

# ── Create data generator ────────────────────────────────────────────────────
cat("Creating generator...\n"); flush(stdout())
t0 <- Sys.time()
gen <- get_generator(
  path = fasta_files, train_type = "label_csv", batch_size = batch_size,
  maxlen = maxlen, step = step, shuffle_file_order = FALSE,
  vocabulary = vocabulary, seed = seed, shuffle_input = TRUE,
  format = "fasta", reverse_complement = TRUE, ambiguous_nuc = "zero",
  padding = TRUE, target_from_csv = target_from_csv,
  max_samples = max_samples, concat_seq = "", target_len = 1,
  sample_by_file_size = FALSE, random_sampling = TRUE,
  reverse_complement_encoding = FALSE, read_data = FALSE,
  target_split = target_split, output_format = "target_right",
  val = FALSE, return_int = FALSE, delete_used_files = FALSE)
cat("Generator created in", round(difftime(Sys.time(), t0, units = "secs"), 1),
    "s\n"); flush(stdout())

# ── Load model ───────────────────────────────────────────────────────────────
cat("Loading model...\n"); flush(stdout())
t0 <- Sys.time()
id <- 28; ep_index <- 23
cp <- file.path(cfg$checkpoints_dir, paste0("bacdive_attn_lstm_sum_", id))
run_name <- basename(cp)
cp <- deepG:::get_cp(cp, ep_index = ep_index)
model <- load_cp(cp, compile = FALSE, mirrored_strategy = TRUE)
ep <- cp %>% basename() %>% stringr::str_extract("\\d+") %>% as.integer()
cat("Model loaded in", round(difftime(Sys.time(), t0, units = "secs"), 1),
    "s\n"); flush(stdout())

cat("Model:", run_name, "Epoch:", ep, "\n")
cat("Model inputs:", length(model$inputs), "\n")
for (i in seq_along(model$inputs))
  cat("  Input", i, "shape:",
      paste(as.integer(model$inputs[[i]]$shape), collapse = ", "), "\n")
flush(stdout())

# ── Reshape helper ───────────────────────────────────────────────────────────
to_time_dist <- function(x, spt) {
  xd <- dim(x)
  keras::k_eval(keras::k_reshape(x, c(xd[1], spt, xd[2] / spt, xd[3])))
}

# ── Generate + Predict in streaming fashion ──────────────────────────────────
set.seed(4)
pred_list <- vector("list", num_batches)
true_list <- vector("list", num_batches)
sw_list   <- vector("list", num_batches)

cat("\nStarting", num_batches, "batches of generation + inference...\n")
flush(stdout())
t0 <- Sys.time()

for (i in 1:num_batches) {
  z <- gen()
  x_raw <- z[[1]]
  y <- z[[2]]

  sw <- set_sw(y_list = y, na_enc = -9999, cw_list = cw_list)

  # Reshape for model
  x_td <- to_time_dist(x_raw, samples_per_target)
  if (length(model$inputs) == 2) {
    x_input <- list(x_raw, x_td)
  } else {
    x_input <- x_td
  }

  y_pred <- model$predict(x_input, verbose = 0L)

  pred_list[[i]] <- do.call(cbind, y_pred)
  true_list[[i]] <- do.call(cbind, y)
  sw_list[[i]]   <- do.call(cbind, sw)

  if (i %% 10 == 0 || i <= 3) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    eta <- elapsed / i * (num_batches - i)
    cat(sprintf("  Batch %d/%d (%.1f min elapsed, ~%.1f min remaining)\n",
                i, num_batches, elapsed, eta))
    flush(stdout())
  }
}

elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("\nDone! %d batches in %.1f minutes\n", num_batches, elapsed))
flush(stdout())

# ── Combine results ──────────────────────────────────────────────────────────
pred_df <- do.call(rbind, pred_list) %>% as.data.frame()
names(pred_df) <- paste0(unlist(target_split_named), "_pred")

true_df <- do.call(rbind, true_list) %>% as.data.frame()
names(true_df) <- unlist(target_split_named)

sw_df <- do.call(rbind, sw_list) %>% as.data.frame()
names(sw_df) <- names(target_split_named)

cat("Total test samples:", nrow(pred_df), "\n")
cat("Pred columns:", ncol(pred_df), "\n")

# ── Save ─────────────────────────────────────────────────────────────────────
rds_path <- cfg$bacdive_test_rds
l_all <- list(
  l = list(pred_df = pred_df, true_df = true_df, sw_df = sw_df),
  ep = ep,
  run_name = run_name
)
saveRDS(l_all, rds_path)
cat("Saved to:", rds_path, "\n")
