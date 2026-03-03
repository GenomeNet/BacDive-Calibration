# Build calibrated deployment CSVs with label-preserving confidence calibration.
# Usage: Rscript R/deploy_calibrated_csv.R

library(dplyr)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg)) else normalizePath(getwd())
script_dir <- if (dir.exists(script_path)) script_path else dirname(script_path)
proj_dir <- if (basename(script_dir) == "R") dirname(script_dir) else script_dir

source(file.path(proj_dir, "R", "config_paths.R"))
source(file.path(proj_dir, "R", "calibration_common.R"))

cfg <- load_bacdive_config(
  proj_dir,
  required = c("bd_pred_csv", "bd_labels_csv", "bd_splits_csv")
)

deployment_dir <- file.path(proj_dir, "data", "deployment")
dir.create(deployment_dir, recursive = TRUE, showWarnings = FALSE)

raw_full_path <- Sys.getenv("BACDIVE_DEPLOY_RAW_FULL_CSV", unset = cfg$bd_pred_csv)

calibration_rds <- Sys.getenv(
  "BACDIVE_DEPLOY_CALIBRATION_RDS",
  unset = file.path(proj_dir, "output", "calibration_logo.rds")
)

date_tag <- format(Sys.Date(), "%Y%m%d")
out_full_cal <- file.path(deployment_dir, paste0("bd_calibrated_full_", date_tag, ".csv"))
out_integrate_cal <- file.path(deployment_dir, paste0("bd_calibrated_integrate_", date_tag, ".csv"))
out_ece <- file.path(deployment_dir, paste0("bd_calibration_ece_", date_tag, ".csv"))
param_source <- Sys.getenv("BACDIVE_DEPLOY_PARAM_SOURCE", unset = "lofo_median")

clip01 <- function(p, eps = 1e-9) pmax(pmin(p, 1 - eps), eps)

softmax_t_rows <- function(p_mat, T = 1) {
  z <- log(clip01(p_mat))
  softmax_t(z, T)
}

pick_best_method <- function(ece_none, ece_temp = NA_real_, ece_platt = NA_real_,
                             ece_dirichlet = NA_real_, tol = 1e-4) {
  cands <- c(none = ece_none, temp = ece_temp, platt = ece_platt, dirichlet = ece_dirichlet)
  cands <- cands[is.finite(cands)]
  if (!("none" %in% names(cands))) return("none")
  best <- names(cands)[which.min(cands)][1]
  if ((cands[["none"]] - cands[[best]]) <= tol) return("none")
  best
}

infer_calibration_fit_scope <- function(tbl) {
  nms <- names(tbl)
  if (all(c("ECE_lofo_temp", "median_lofo_T", "global_T") %in% nms)) return("valtest_lofo")
  if (all(c("ECE_tc_temp", "tc_T") %in% nms) || any(startsWith(nms, "tc_"))) return("test_7030_tc")
  if (all(c("ECE_val_temp", "val_T") %in% nms) || any(startsWith(nms, "val_"))) return("validation_to_test")
  if (all(c("ECE_temp", "T", "a", "b") %in% nms)) return("single_split_7030")
  "unknown"
}

choose_method_from_lofo <- function(row, tol = 1e-4) {
  type <- as.character(row$type[1])

  get_num <- function(col) {
    if (!(col %in% names(row))) return(NA_real_)
    v <- suppressWarnings(as.numeric(row[[col]]))
    if (length(v) == 0) return(NA_real_)
    v[1]
  }

  e_none <- get_num("ECE_uncal")
  e_temp <- get_num("ECE_lofo_temp")
  e_platt <- get_num("ECE_lofo_platt")

  if (identical(type, "linear_platt")) {
    if (is.finite(e_platt) && (!is.finite(e_none) || (e_none - e_platt) > tol)) return("platt_linear")
    return("none")
  }
  if (type %in% c("binary_sigmoid", "softmax_2class")) {
    return(pick_best_method(e_none, e_temp, e_platt, tol = tol))
  }
  if (startsWith(type, "multiclass_")) {
    # Deployment currently supports temperature for multiclass.
    return(pick_best_method(e_none, e_temp, tol = tol))
  }
  "none"
}

to_bool <- function(x) {
  if (is.null(x)) return(NA)
  lx <- tolower(as.character(x))
  lx[lx == ""] <- NA_character_
  ifelse(lx %in% c("true", "1", "t", "yes", "y"), TRUE,
         ifelse(lx %in% c("false", "0", "f", "no", "n"), FALSE, NA))
}

calibration_obj <- readRDS(calibration_rds)
summary_tbl <- calibration_obj$summary
cal_fit_scope <- infer_calibration_fit_scope(summary_tbl)
if (!("quality" %in% names(summary_tbl))) {
  summary_tbl$quality <- "PASS"
}
summary_tbl$target <- as.character(summary_tbl$target)

get_row <- function(target) {
  x <- summary_tbl[summary_tbl$target == target, , drop = FALSE]
  if (nrow(x) == 0) return(NULL)
  x[1, , drop = FALSE]
}

pick_param <- function(row, which, source = "lofo_median") {
  col <- if (source == "lofo_median") paste0("median_lofo_", which) else paste0("global_", which)
  suppressWarnings(as.numeric(row[[col]][1]))
}

binary_heads <- list(
  motility = list(target = "motility", val_col = "is_motile_val",
                  prob_col = "is_motile_prob", label_col = "is_motile",
                  true_col = "is_motile_true"),
  oxygen_microaerophile = list(target = "oxygen_microaerophile",
                               val_col = "oxygen_microaerophile_val",
                               prob_col = "oxygen_microaerophile_prob",
                               label_col = "oxygen_microaerophile",
                               true_col = "microaerophile"),
  spore = list(target = "spore", val_col = "ability_spore_val",
               prob_col = "ability_spore_prob", label_col = "ability_spore",
               true_col = "ability_spore_true")
)

softmax2_heads <- list(
  oxygen_growth = list(target = "oxygen_growth", val_col = "oxygen_growth_val",
                       prob_col = "oxygen_growth_prob", label_col = "oxygen_growth",
                       labels = c("aerobe", "anaerobe"),
                       true_cols = c("aerobe", "anaerobe")),
  oxygen_facultative = list(target = "oxygen_facultative", val_col = "oxygen_facultative_val",
                            prob_col = "oxygen_facultative_prob", label_col = "oxygen_facultative",
                            labels = c("facultative.aerobe", "facultative.anaerobe"),
                            true_cols = c("facultative.aerobe", "facultative.anaerobe")),
  oxygen_obligate = list(target = "oxygen_obligate", val_col = "oxygen_obligate_val",
                         prob_col = "oxygen_obligate_prob", label_col = "oxygen_obligate",
                         labels = c("obligate.aerobe", "obligate.anaerobe"),
                         true_cols = c("obligate.aerobe", "obligate.anaerobe"))
)

multiclass_heads <- list(
  gram = list(target = "gram", val_col = "gram_val", prob_col = "gram_prob",
              label_col = "gram", labels = c("positive", "variable", "negative"),
              true_cols = c("is_gram_stain_positive", "is_gram_stain_variable",
                            "is_gram_stain_negative")),
  cellshape = list(target = "cellshape", val_col = "cellshape_val", prob_col = "cellshape_prob",
                   label_col = "cellshape",
                   labels = c("rod.shaped", "other", "coccus.shaped", "vibrio.shaped",
                              "filament.shaped", "sphere.shaped", "ovoid.shaped",
                              "pleomorphic.shaped", "spiral.shaped", "curved.shaped",
                              "oval.shaped"),
                   true_cols = c("is_cell_shape_rod.shaped", "is_cell_shape_other",
                                 "is_cell_shape_coccus.shaped", "is_cell_shape_vibrio.shaped",
                                 "is_cell_shape_filament.shaped", "is_cell_shape_sphere.shaped",
                                 "is_cell_shape_ovoid.shaped", "is_cell_shape_pleomorphic.shaped",
                                 "is_cell_shape_spiral.shaped", "is_cell_shape_curved.shaped",
                                 "is_cell_shape_oval.shaped")),
  flagellum = list(target = "flagellum", val_col = "flagellum_val", prob_col = "flagellum_prob",
                   label_col = "flagellum",
                   labels = c("monotrichous", "monotrichous_polar", "polar",
                              "peritrichous", "lophotrichous", "gliding"),
                   true_cols = c("is_flagellum_arrangement_monotrichous",
                                 "is_flagellum_arrangement_monotrichous_polar",
                                 "is_flagellum_arrangement_polar",
                                 "is_flagellum_arrangement_peritrichous",
                                 "is_flagellum_arrangement_lophotrichous",
                                 "is_flagellum_arrangement_gliding"))
)

patho_heads <- list(
  pathogenicity_human = list(target = "pathogenicity_human",
                             val_col = "pathogenicity_human_val",
                             prob_col = "pathogenicity_human_prob",
                             label_col = "pathogenicity_human",
                             true_col = "pathogenicity_human_true",
                             threshold = 1.4),
  pathogenicity_animal = list(target = "pathogenicity_animal",
                              val_col = "pathogenicity_animal_val",
                              prob_col = "pathogenicity_animal_prob",
                              label_col = "pathogenicity_animal",
                              true_col = "pathogenicity_animal_true",
                              threshold = 1.4),
  pathogenicity_plant = list(target = "pathogenicity_plant",
                             val_col = "pathogenicity_plant_val",
                             prob_col = "pathogenicity_plant_prob",
                             label_col = "pathogenicity_plant",
                             true_col = "pathogenicity_plant_true",
                             threshold = 1.0)
)

cat("Loading raw full predictions:", raw_full_path, "\n")
raw_full <- read.csv(raw_full_path, check.names = FALSE)

cat("Building val+test evaluation frame for deployment summary...\n")
labels_true <- read.csv(cfg$bd_labels_csv, check.names = FALSE)
splits <- read.csv(cfg$bd_splits_csv, check.names = FALSE)

raw_eval <- raw_full
raw_eval$genome <- basename(raw_eval$file)
labels_true$genome <- labels_true$file
splits$genome <- splits$file

eval_df <- merge(raw_eval, labels_true, by = "genome", suffixes = c("_pred", "_true"))
eval_df <- merge(eval_df, splits[, c("genome", "type")], by = "genome")
eval_df <- eval_df %>% filter(type %in% c("validation", "test"))
tax_tbl <- calibration_obj$genome_taxonomy
tax_tbl$genome <- basename(as.character(tax_tbl$file))
keep_genomes_tax <- unique(tax_tbl$genome[!is.na(tax_tbl$Phylum)])
n_before <- nrow(eval_df)
eval_df <- eval_df %>% filter(genome %in% keep_genomes_tax)
taxonomy_filter <- "phylum_non_na"
cat("Applied taxonomy filter (Phylum non-NA):", n_before, "->", nrow(eval_df), "\n")
cat("Evaluation genomes (val+test):", nrow(eval_df), "\n")

rows <- list()
selected_methods <- list()

for (name in names(binary_heads)) {
  h <- binary_heads[[name]]
  srow <- get_row(h$target)
  if (is.null(srow)) stop("Missing calibration row for target: ", h$target, call. = FALSE)

  p_raw <- as.numeric(eval_df[[h$val_col]])
  y <- as.numeric(eval_df[[h$true_col]])
  ok <- !is.na(y) & y != -9999 & is.finite(p_raw)
  p_raw <- clip01(p_raw[ok]); y <- y[ok]
  pred_raw <- as.numeric(p_raw > 0.5)
  correct <- as.numeric(pred_raw == y)

  T <- pick_param(srow, "T", param_source)
  a <- pick_param(srow, "a", param_source)
  b <- pick_param(srow, "b", param_source)
  q <- as.character(srow$quality)

  p_temp <- if (is.finite(T) && T > 0) clip01(sigmoid(logit(p_raw) / T)) else rep(NA_real_, length(p_raw))
  p_platt <- if (is.finite(a) && is.finite(b)) clip01(sigmoid(a * logit(p_raw) + b)) else rep(NA_real_, length(p_raw))

  e_none <- ece(ifelse(pred_raw == 1, p_raw, 1 - p_raw), correct)
  e_temp <- if (all(is.na(p_temp))) NA_real_ else ece(ifelse(pred_raw == 1, p_temp, 1 - p_temp), correct)
  e_platt <- if (all(is.na(p_platt))) NA_real_ else ece(ifelse(pred_raw == 1, p_platt, 1 - p_platt), correct)
  e_event_none <- ece(p_raw, y)
  e_event_temp <- if (all(is.na(p_temp))) NA_real_ else ece(p_temp, y)
  e_event_platt <- if (all(is.na(p_platt))) NA_real_ else ece(p_platt, y)

  method <- choose_method_from_lofo(srow)
  e_sel <- c(none = e_none, temp = e_temp, platt = e_platt)[method]

  selected_methods[[name]] <- list(method = method, T = T, a = a, b = b)
  rows[[name]] <- data.frame(
    target = name, type = as.character(srow$type), quality = q,
    N_eval = length(y), Acc_raw = mean(correct),
    selected_method = method,
    ECE_conf_none = e_none, ECE_conf_temp = e_temp, ECE_conf_platt = e_platt,
    ECE_conf_selected = e_sel,
    ECE_event_none = e_event_none, ECE_event_temp = e_event_temp, ECE_event_platt = e_event_platt,
    ECE_event_selected = c(none = e_event_none, temp = e_event_temp, platt = e_event_platt)[method],
    delta_conf_vs_none = e_none - e_sel,
    delta_event_vs_none = e_event_none - c(none = e_event_none, temp = e_event_temp, platt = e_event_platt)[method],
    param_source = param_source,
    param_T = T, param_a = a, param_b = b,
    stringsAsFactors = FALSE
  )
}

for (name in names(softmax2_heads)) {
  h <- softmax2_heads[[name]]
  srow <- get_row(h$target)
  if (is.null(srow)) stop("Missing calibration row for target: ", h$target, call. = FALSE)

  pmat <- parse_prob_vec(eval_df[[h$val_col]])
  ymat <- as.matrix(eval_df[, h$true_cols])
  storage.mode(ymat) <- "numeric"
  ok <- apply(ymat, 1, function(r) !any(is.na(r) | r == -9999)) & rowSums(ymat) > 0
  pmat <- pmat[ok, , drop = FALSE]
  ymat <- ymat[ok, , drop = FALSE]

  pred_raw_idx <- apply(pmat, 1, which.max)
  true_idx <- apply(ymat, 1, which.max)
  correct <- as.numeric(pred_raw_idx == true_idx)

  T <- pick_param(srow, "T", param_source)
  a <- pick_param(srow, "a", param_source)
  b <- pick_param(srow, "b", param_source)
  q <- as.character(srow$quality)

  p_temp <- if (is.finite(T) && T > 0) softmax_t_rows(pmat, T) else matrix(NA_real_, nrow(pmat), ncol(pmat))
  p2_platt <- if (is.finite(a) && is.finite(b)) clip01(sigmoid(a * logit(pmat[, 2]) + b)) else rep(NA_real_, nrow(pmat))
  p_platt <- if (all(is.na(p2_platt))) matrix(NA_real_, nrow(pmat), ncol(pmat)) else cbind(1 - p2_platt, p2_platt)

  conf_none <- pmat[cbind(seq_len(nrow(pmat)), pred_raw_idx)]
  conf_temp <- if (all(is.na(p_temp))) rep(NA_real_, nrow(pmat)) else p_temp[cbind(seq_len(nrow(pmat)), pred_raw_idx)]
  conf_platt <- if (all(is.na(p_platt))) rep(NA_real_, nrow(pmat)) else p_platt[cbind(seq_len(nrow(pmat)), pred_raw_idx)]
  y_event <- as.numeric(true_idx == 2)
  p2_none <- pmat[, 2]
  p2_temp <- if (all(is.na(p_temp))) rep(NA_real_, nrow(pmat)) else p_temp[, 2]
  p2_platt <- if (all(is.na(p_platt))) rep(NA_real_, nrow(pmat)) else p_platt[, 2]

  e_none <- ece(conf_none, correct)
  e_temp <- if (all(is.na(conf_temp))) NA_real_ else ece(conf_temp, correct)
  e_platt <- if (all(is.na(conf_platt))) NA_real_ else ece(conf_platt, correct)
  e_event_none <- ece(p2_none, y_event)
  e_event_temp <- if (all(is.na(p2_temp))) NA_real_ else ece(p2_temp, y_event)
  e_event_platt <- if (all(is.na(p2_platt))) NA_real_ else ece(p2_platt, y_event)

  method <- choose_method_from_lofo(srow)
  e_sel <- c(none = e_none, temp = e_temp, platt = e_platt)[method]

  selected_methods[[name]] <- list(method = method, T = T, a = a, b = b)
  rows[[name]] <- data.frame(
    target = name, type = as.character(srow$type), quality = q,
    N_eval = nrow(pmat), Acc_raw = mean(correct),
    selected_method = method,
    ECE_conf_none = e_none, ECE_conf_temp = e_temp, ECE_conf_platt = e_platt,
    ECE_conf_selected = e_sel,
    ECE_event_none = e_event_none, ECE_event_temp = e_event_temp, ECE_event_platt = e_event_platt,
    ECE_event_selected = c(none = e_event_none, temp = e_event_temp, platt = e_event_platt)[method],
    delta_conf_vs_none = e_none - e_sel,
    delta_event_vs_none = e_event_none - c(none = e_event_none, temp = e_event_temp, platt = e_event_platt)[method],
    param_source = param_source,
    param_T = T, param_a = a, param_b = b,
    stringsAsFactors = FALSE
  )
}

for (name in names(multiclass_heads)) {
  h <- multiclass_heads[[name]]
  srow <- get_row(h$target)
  if (is.null(srow)) stop("Missing calibration row for target: ", h$target, call. = FALSE)

  pmat <- parse_prob_vec(eval_df[[h$val_col]])
  ymat <- as.matrix(eval_df[, h$true_cols])
  storage.mode(ymat) <- "numeric"
  ok <- apply(ymat, 1, function(r) !any(is.na(r) | r == -9999)) & rowSums(ymat) > 0
  pmat <- pmat[ok, , drop = FALSE]
  ymat <- ymat[ok, , drop = FALSE]

  pred_raw_idx <- apply(pmat, 1, which.max)
  true_idx <- apply(ymat, 1, which.max)
  correct <- as.numeric(pred_raw_idx == true_idx)

  T <- pick_param(srow, "T", param_source)
  q <- as.character(srow$quality)

  p_temp <- if (is.finite(T) && T > 0) softmax_t_rows(pmat, T) else matrix(NA_real_, nrow(pmat), ncol(pmat))
  conf_none <- pmat[cbind(seq_len(nrow(pmat)), pred_raw_idx)]
  conf_temp <- if (all(is.na(p_temp))) rep(NA_real_, nrow(pmat)) else p_temp[cbind(seq_len(nrow(pmat)), pred_raw_idx)]

  e_none <- ece(conf_none, correct)
  e_temp <- if (all(is.na(conf_temp))) NA_real_ else ece(conf_temp, correct)
  e_event_none <- ece(as.vector(pmat), as.vector(ymat))
  e_event_temp <- if (all(is.na(p_temp))) NA_real_ else ece(as.vector(p_temp), as.vector(ymat))

  method <- choose_method_from_lofo(srow)
  e_sel <- c(none = e_none, temp = e_temp)[method]

  selected_methods[[name]] <- list(method = method, T = T, a = NA_real_, b = NA_real_)
  rows[[name]] <- data.frame(
    target = name, type = as.character(srow$type), quality = q,
    N_eval = nrow(pmat), Acc_raw = mean(correct),
    selected_method = method,
    ECE_conf_none = e_none, ECE_conf_temp = e_temp, ECE_conf_platt = NA_real_,
    ECE_conf_selected = e_sel,
    ECE_event_none = e_event_none, ECE_event_temp = e_event_temp, ECE_event_platt = NA_real_,
    ECE_event_selected = c(none = e_event_none, temp = e_event_temp)[method],
    delta_conf_vs_none = e_none - e_sel,
    delta_event_vs_none = e_event_none - c(none = e_event_none, temp = e_event_temp)[method],
    param_source = param_source,
    param_T = T, param_a = NA_real_, param_b = NA_real_,
    stringsAsFactors = FALSE
  )
}

for (name in names(patho_heads)) {
  h <- patho_heads[[name]]
  srow <- get_row(h$target)
  if (is.null(srow)) stop("Missing calibration row for target: ", h$target, call. = FALSE)

  raw <- as.numeric(eval_df[[h$val_col]])
  y_raw <- as.numeric(eval_df[[h$true_col]])
  ok <- !is.na(y_raw) & y_raw != -9999 & is.finite(raw)
  raw <- raw[ok]
  y <- as.numeric(y_raw[ok] > 0)
  pred_raw <- as.numeric(raw > h$threshold)
  correct <- as.numeric(pred_raw == y)

  a <- pick_param(srow, "a", param_source)
  b <- pick_param(srow, "b", param_source)
  q <- as.character(srow$quality)

  # Legacy deployed confidence (uncalibrated baseline).
  conf_none <- (tanh(5 * abs(raw - h$threshold)) + 1) / 2
  p_platt <- if (is.finite(a) && is.finite(b)) clip01(sigmoid(a * raw + b)) else rep(NA_real_, length(raw))
  conf_platt <- if (all(is.na(p_platt))) rep(NA_real_, length(raw)) else ifelse(pred_raw == 1, p_platt, 1 - p_platt)

  e_none <- ece(conf_none, correct)
  e_platt <- if (all(is.na(conf_platt))) NA_real_ else ece(conf_platt, correct)
  e_event_none <- ece(sigmoid(raw - h$threshold), y)
  e_event_platt <- if (all(is.na(p_platt))) NA_real_ else ece(p_platt, y)

  method <- choose_method_from_lofo(srow)
  e_sel <- c(none = e_none, platt_linear = e_platt)[method]

  selected_methods[[name]] <- list(method = method,
                                   T = NA_real_, a = a, b = b)
  rows[[name]] <- data.frame(
    target = name, type = as.character(srow$type), quality = q,
    N_eval = length(y), Acc_raw = mean(correct),
    selected_method = method,
    ECE_conf_none = e_none, ECE_conf_temp = NA_real_, ECE_conf_platt = e_platt,
    ECE_conf_selected = e_sel,
    ECE_event_none = e_event_none, ECE_event_temp = NA_real_, ECE_event_platt = e_event_platt,
    ECE_event_selected = c(none = e_event_none, platt_linear = e_event_platt)[method],
    delta_conf_vs_none = e_none - e_sel,
    delta_event_vs_none = e_event_none - c(none = e_event_none, platt_linear = e_event_platt)[method],
    param_source = param_source,
    param_T = NA_real_, param_a = a, param_b = b,
    stringsAsFactors = FALSE
  )
}

ece_summary <- bind_rows(rows) %>%
  arrange(desc(delta_conf_vs_none), target) %>%
  mutate(
    calibration_rds = calibration_rds,
    calibration_fit_scope = cal_fit_scope,
    eval_split = "valtest",
    taxonomy_filter = taxonomy_filter,
    objective = "lofo_event_ece_method_selection",
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )

excluded_rows <- data.frame(
  target = c("biosafety", "cellsize"),
  type = c("excluded", "excluded"),
  quality = c("EXCLUDED", "EXCLUDED"),
  N_eval = as.numeric(NA),
  Acc_raw = as.numeric(NA),
  selected_method = c("none", "none"),
  ECE_conf_none = as.numeric(NA),
  ECE_conf_temp = as.numeric(NA),
  ECE_conf_platt = as.numeric(NA),
  ECE_conf_selected = as.numeric(NA),
  ECE_event_none = as.numeric(NA),
  ECE_event_temp = as.numeric(NA),
  ECE_event_platt = as.numeric(NA),
  ECE_event_selected = as.numeric(NA),
  delta_conf_vs_none = as.numeric(NA),
  delta_event_vs_none = as.numeric(NA),
  param_source = param_source,
  param_T = as.numeric(NA),
  param_a = as.numeric(NA),
  param_b = as.numeric(NA),
  calibration_rds = calibration_rds,
  calibration_fit_scope = cal_fit_scope,
  eval_split = "valtest",
  taxonomy_filter = taxonomy_filter,
  objective = "lofo_event_ece_method_selection",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)

ece_summary <- bind_rows(ece_summary, excluded_rows) %>% arrange(target)

cat("\nSelected methods (LOFO event-ECE objective from calibration_logo):\n")
print(ece_summary %>% select(target, quality, selected_method, delta_conf_vs_none, delta_event_vs_none), row.names = FALSE)
cat("Parameter source:", param_source, "\n")

cat("\nApplying selected calibration to _prob only (labels/decisions unchanged)...\n")
cal_full <- raw_full

# Keep raw _val for provenance and add calibrated counterparts.
all_heads <- c(binary_heads, softmax2_heads, multiclass_heads, patho_heads)
for (name in names(all_heads)) {
  h <- all_heads[[name]]
  val_cal_col <- sub("_val$", "_val_cal", h$val_col)
  if (!identical(val_cal_col, h$val_col)) {
    cal_full[[val_cal_col]] <- cal_full[[h$val_col]]
  }
}

for (name in names(binary_heads)) {
  h <- binary_heads[[name]]
  sel <- selected_methods[[name]]
  if (is.null(sel) || identical(sel$method, "none")) next

  p_raw <- clip01(as.numeric(cal_full[[h$val_col]]))
  pred_lbl <- to_bool(cal_full[[h$label_col]])
  pred_raw <- ifelse(is.na(pred_lbl), p_raw > 0.5, pred_lbl)

  p_cal <- p_raw
  if (identical(sel$method, "temp") && is.finite(sel$T) && sel$T > 0) {
    p_cal <- clip01(sigmoid(logit(p_raw) / sel$T))
  } else if (identical(sel$method, "platt") && is.finite(sel$a) && is.finite(sel$b)) {
    p_cal <- clip01(sigmoid(sel$a * logit(p_raw) + sel$b))
  } else {
    next
  }

  cal_full[[sub("_val$", "_val_cal", h$val_col)]] <- p_cal
  cal_full[[h$prob_col]] <- clip01(ifelse(pred_raw, p_cal, 1 - p_cal))
}

for (name in names(softmax2_heads)) {
  h <- softmax2_heads[[name]]
  sel <- selected_methods[[name]]
  if (is.null(sel) || identical(sel$method, "none")) next

  pmat <- parse_prob_vec(cal_full[[h$val_col]])
  pmat <- pmat / rowSums(pmat)
  pred_idx <- match(as.character(cal_full[[h$label_col]]), h$labels)
  if (any(is.na(pred_idx))) {
    stop("Unmapped predicted labels for target: ", name, call. = FALSE)
  }

  p_cal <- pmat
  if (identical(sel$method, "temp") && is.finite(sel$T) && sel$T > 0) {
    p_cal <- softmax_t_rows(pmat, sel$T)
  } else if (identical(sel$method, "platt") && is.finite(sel$a) && is.finite(sel$b)) {
    p2 <- clip01(sigmoid(sel$a * logit(pmat[, 2]) + sel$b))
    p_cal <- cbind(1 - p2, p2)
  } else {
    next
  }

  cal_full[[sub("_val$", "_val_cal", h$val_col)]] <- apply(p_cal, 1, function(v) paste(v, collapse = ", "))
  cal_full[[h$prob_col]] <- clip01(p_cal[cbind(seq_len(nrow(p_cal)), pred_idx)])
}

for (name in names(multiclass_heads)) {
  h <- multiclass_heads[[name]]
  sel <- selected_methods[[name]]
  if (is.null(sel) || identical(sel$method, "none")) next
  if (!(identical(sel$method, "temp") && is.finite(sel$T) && sel$T > 0)) next

  pmat <- parse_prob_vec(cal_full[[h$val_col]])
  pmat <- pmat / rowSums(pmat)
  pred_idx <- match(as.character(cal_full[[h$label_col]]), h$labels)
  if (any(is.na(pred_idx))) {
    stop("Unmapped predicted labels for target: ", name, call. = FALSE)
  }

  p_cal <- softmax_t_rows(pmat, sel$T)
  cal_full[[sub("_val$", "_val_cal", h$val_col)]] <- apply(p_cal, 1, function(v) paste(v, collapse = ", "))
  cal_full[[h$prob_col]] <- clip01(p_cal[cbind(seq_len(nrow(p_cal)), pred_idx)])
}

for (name in names(patho_heads)) {
  h <- patho_heads[[name]]
  sel <- selected_methods[[name]]
  if (is.null(sel) || !identical(sel$method, "platt_linear")) next
  if (!is.finite(sel$a) || !is.finite(sel$b)) next

  raw <- as.numeric(cal_full[[h$val_col]])
  pred_lbl <- to_bool(cal_full[[h$label_col]])
  pred_raw <- ifelse(is.na(pred_lbl), raw > h$threshold, pred_lbl)
  p_cal <- clip01(sigmoid(sel$a * raw + sel$b))
  cal_full[[sub("_val$", "_val_cal", h$val_col)]] <- p_cal
  cal_full[[h$prob_col]] <- clip01(ifelse(pred_raw, p_cal, 1 - p_cal))
}

to_integrate_subset <- function(df_in) {
  df2 <- df_in[, !grepl("_val$", names(df_in)) & !grepl("^rg", names(df_in))]
  df2 <- df2[, c(
    "ability_spore", "ability_spore_prob",
    "oxygen_growth", "oxygen_growth_prob",
    "oxygen_obligate", "oxygen_obligate_prob",
    "oxygen_microaerophile", "oxygen_microaerophile_prob",
    "oxygen_facultative", "oxygen_facultative_prob",
    "gram", "gram_prob",
    "is_motile", "is_motile_prob",
    "biosafety", "biosafety_prob",
    "file"
  )]
  df2
}

cat("Raw predictions source:", raw_full_path, "\n")

cal_integrate <- to_integrate_subset(cal_full)

write.csv(cal_full, out_full_cal, row.names = FALSE)
write.csv(cal_integrate, out_integrate_cal, row.names = FALSE)
write.csv(ece_summary, out_ece, row.names = FALSE)

cat("\nSaved:\n")
cat("  ", out_full_cal, "\n")
cat("  ", out_integrate_cal, "\n")
cat("  ", out_ece, "\n")
