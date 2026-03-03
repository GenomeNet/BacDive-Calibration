# Smoke test: 3 batches of generation + inference
# Usage: conda run -n tf Rscript smoke_test.R

library(deepG)
library(keras)
library(tensorflow)
library(dplyr)

# ── Inline set_sw (avoid sourcing data_prep.R which writes to restricted paths)
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

set_sw_cellsize <- function(y, batch_size, na_enc = -9999) {
  sw <- matrix(1, nrow = batch_size)
  na_col <- apply(y == na_enc, 1, any)
  sw[na_col] <- 0
  sw
}
set_sw_class <- function(y, batch_size, na_enc = -9999, cw) {
  sw <- matrix(1, nrow = batch_size)
  na_row <- apply(y == na_enc, 1, any) | rowSums(y) == 0
  na_row <- which(na_row)
  sw[na_row, ] <- 0
  non_na_row <- setdiff(1:batch_size, na_row)
  index <- apply(y[non_na_row, , drop = FALSE], 1, which.max)
  sw[non_na_row, ] <- cw[index]
  sw
}
set_sw_binary <- function(y, batch_size, na_enc = -9999, cw) {
  sw <- matrix(1, nrow = batch_size)
  na_row <- which(apply(y == na_enc, 1, any))
  sw[na_row, ] <- 0
  non_na_row <- setdiff(1:batch_size, na_row)
  index <- ifelse(y[non_na_row, 1] == 0, 1, 2)
  sw[non_na_row, ] <- cw[index]
  sw
}
set_sw_patho <- function(y, batch_size, cw, default_weight = 1e-12) {
  cw[1] <- cw[1] * default_weight
  sw <- cw[y[, 1] + 1]
  matrix(sw, ncol = 1)
}
set_sw_biosafety <- function(y, batch_size, na_enc = -9999, cw) {
  sw <- matrix(1, nrow = batch_size)
  na_row <- which(apply(y == na_enc, 1, any))
  sw[na_row, ] <- 0
  non_na_row <- setdiff(1:batch_size, na_row)
  sw[non_na_row, ] <- cw[y[non_na_row, 1]]
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

# ── Test FASTA files ─────────────────────────────────────────────────────────
ttv <- "test"
maxlen <- 2000000
samples_per_target <- 200
batch_size <- 12
seed <- 77

target_split <- list(
  c("rgStart_len","rgEnd_len","rgStart_wid","rgEnd_wid"),
  c("is_cell_shape_rod.shaped","is_cell_shape_other","is_cell_shape_coccus.shaped",
    "is_cell_shape_vibrio.shaped","is_cell_shape_filament.shaped","is_cell_shape_sphere.shaped",
    "is_cell_shape_ovoid.shaped","is_cell_shape_pleomorphic.shaped","is_cell_shape_spiral.shaped",
    "is_cell_shape_curved.shaped","is_cell_shape_oval.shaped"),
  c("is_flagellum_arrangement_monotrichous","is_flagellum_arrangement_monotrichous_polar",
    "is_flagellum_arrangement_polar","is_flagellum_arrangement_peritrichous",
    "is_flagellum_arrangement_lophotrichous","is_flagellum_arrangement_gliding"),
  c("is_gram_stain_positive","is_gram_stain_variable","is_gram_stain_negative"),
  c("is_motile"), c("biosafety_level"),
  c("pathogenicity_human"), c("pathogenicity_animal"), c("pathogenicity_plant"),
  c("aerobe","anaerobe"), c("facultative.aerobe","facultative.anaerobe"),
  c("obligate.aerobe","obligate.anaerobe"), c("microaerophile"), c("ability_spore")
)

all_files <- list.files(cfg$fasta_dir, full.names = TRUE)
all_files <- all_files[file.info(all_files)$size != 0]
all_csv <- read.csv(cfg$bd_splits_csv)
type_csv <- all_csv %>% filter(type == ttv) %>% select(file) %>% unlist()
fasta_files <- as.list(all_files[basename(all_files) %in% type_csv])
cat("Test FASTA files:", length(fasta_files), "\n")

target_from_csv <- cfg$bd_labels_csv

cat("Creating generator...\n"); t0 <- Sys.time()
gen <- get_generator(
  path = fasta_files, train_type = "label_csv", batch_size = batch_size,
  maxlen = maxlen, step = maxlen, shuffle_file_order = FALSE,
  vocabulary = c("A","C","G","T"), seed = seed, shuffle_input = TRUE,
  format = "fasta", reverse_complement = TRUE, ambiguous_nuc = "zero",
  padding = TRUE, target_from_csv = target_from_csv,
  max_samples = floor(batch_size / 6), concat_seq = "", target_len = 1,
  sample_by_file_size = FALSE, random_sampling = TRUE,
  reverse_complement_encoding = FALSE, read_data = FALSE,
  target_split = target_split, output_format = "target_right",
  val = FALSE, return_int = FALSE, delete_used_files = FALSE)
cat("Generator created in", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")

cat("Loading model...\n"); t0 <- Sys.time()
cp <- file.path(cfg$checkpoints_dir, "bacdive_attn_lstm_sum_28")
cp <- deepG:::get_cp(cp, ep_index = 23)
model <- load_cp(cp, compile = FALSE, mirrored_strategy = TRUE)
cat("Model loaded in", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")

cat("Model inputs:", length(model$inputs), "\n")
for (i in seq_along(model$inputs))
  cat("  Input", i, "shape:", paste(as.integer(model$inputs[[i]]$shape), collapse = ", "), "\n")

to_time_dist <- function(x, spt) {
  xd <- dim(x)
  keras::k_eval(keras::k_reshape(x, c(xd[1], spt, xd[2] / spt, xd[3])))
}

cat("\n--- Testing 3 batches ---\n")
for (i in 1:3) {
  t1 <- Sys.time()
  z <- gen()
  x_raw <- z[[1]]
  y <- z[[2]]
  t_gen <- difftime(Sys.time(), t1, units = "secs")

  sw <- set_sw(y_list = y, na_enc = -9999, cw_list = cw_list)

  t1 <- Sys.time()
  x_td <- to_time_dist(x_raw, samples_per_target)
  if (length(model$inputs) == 2) {
    x_input <- list(x_raw, x_td)
  } else {
    x_input <- x_td
  }
  y_pred <- model$predict(x_input, verbose = 0L)
  t_pred <- difftime(Sys.time(), t1, units = "secs")

  cat(sprintf("Batch %d: x=%s, gen=%.1fs, pred=%.1fs, pred_len=%d\n",
      i, paste(dim(x_raw), collapse = "x"), t_gen, t_pred, length(y_pred)))
}
cat("\nSmoke test PASSED\n")
