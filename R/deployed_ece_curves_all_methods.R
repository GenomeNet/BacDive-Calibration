# Deployed calibration audit with both confidence and event-probability views.
# Produces one PDF with both reliability curves and one CSV with both ECE types.
# Usage: Rscript R/deployed_ece_curves_all_methods.R

library(dplyr)
library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg)) else normalizePath(getwd())
script_dir <- if (dir.exists(script_path)) script_path else dirname(script_path)
proj_dir <- if (basename(script_dir) == "R") dirname(script_dir) else script_dir

source(file.path(proj_dir, "R", "config_paths.R"))
source(file.path(proj_dir, "R", "calibration_common.R"))

cfg <- load_bacdive_config(proj_dir, required = c("bd_pred_csv", "bd_labels_csv", "bd_splits_csv"))

deployment_dir <- file.path(proj_dir, "data", "deployment")
raw_csv <- cfg$bd_pred_csv
param_source <- Sys.getenv("BACDIVE_DEPLOY_PARAM_SOURCE", unset = "lofo_median")

split_target <- Sys.getenv("BACDIVE_AUDIT_SPLIT", unset = "test")
default_cal_rds <- if (split_target == "valtest") {
  file.path(proj_dir, "output", "calibration_val_test.rds")
} else {
  file.path(proj_dir, "output", "calibration_logo.rds")
}
cal_rds <- Sys.getenv("BACDIVE_DEPLOY_CALIBRATION_RDS", unset = default_cal_rds)

keep_split <- if (split_target == "lofo") {
  c("validation", "test")
} else if (split_target == "valtest") {
  "test"
} else {
  split_target
}

clip01 <- function(p, eps = 1e-9) pmax(pmin(p, 1 - eps), eps)

as_bool <- function(x) {
  lx <- tolower(as.character(x))
  ifelse(lx %in% c("true", "1", "t", "yes", "y"), TRUE,
         ifelse(lx %in% c("false", "0", "f", "no", "n"), FALSE, NA))
}

safe_rel <- function(p, y) {
  if (length(p) == 0 || length(y) == 0) return(NULL)
  reliability_data(clip01(as.numeric(p)), as.numeric(y))
}

safe_ece <- function(p, y) {
  if (length(p) == 0 || length(y) == 0) return(NA_real_)
  ece(clip01(as.numeric(p)), as.numeric(y))
}

softmax_from_probs_with_temp <- function(p_mat, T) {
  z <- log(clip01(p_mat))
  softmax_t(z, T)
}

method_order <- c("Raw", "Temp", "Platt", "Selected")

infer_calibration_fit_scope <- function(tbl) {
  nms <- names(tbl)
  if (all(c("ECE_lofo_temp", "median_lofo_T", "global_T") %in% nms)) return("valtest_lofo")
  if (all(c("ECE_tc_temp", "tc_T") %in% nms) || any(startsWith(nms, "tc_"))) return("test_7030_tc")
  if (all(c("ECE_val_temp", "val_T") %in% nms) || any(startsWith(nms, "val_"))) return("validation_to_test")
  if (all(c("ECE_temp", "T", "a", "b") %in% nms)) return("single_split_7030")
  "unknown"
}

pick_param <- function(row, which, split_target, source = "lofo_median") {
  col <- if (split_target == "valtest") {
    paste0("val_", which)
  } else if (source == "lofo_median") {
    paste0("median_lofo_", which)
  } else {
    paste0("global_", which)
  }
  suppressWarnings(as.numeric(row[[col]][1]))
}

choose_selected_method <- function(target, row, split_target, sel_tbl) {
  get_num <- function(col) suppressWarnings(as.numeric(row[[col]][1]))
  if (split_target == "valtest") {
    type <- as.character(row$type[1])
    cands <- if (type %in% c("binary_sigmoid", "softmax_2class")) {
      c(Raw = get_num("ECE_uncal"), Temp = get_num("ECE_val_temp"), Platt = get_num("ECE_val_platt"))
    } else if (startsWith(type, "multiclass_")) {
      c(Raw = get_num("ECE_uncal"), Temp = get_num("ECE_val_temp"))
    } else if (identical(type, "linear_platt")) {
      c(Raw = get_num("ECE_uncal"), Platt = get_num("ECE_val_platt"))
    } else {
      stop("Unsupported type for valtest selected-method mapping: ", type, call. = FALSE)
    }
    cands <- cands[is.finite(cands)]
    if (length(cands) == 0) stop("No finite valtest ECE candidates for target: ", target, call. = FALSE)
    return(names(cands)[which.min(cands)][1])
  }

  r <- sel_tbl[sel_tbl$target == target, , drop = FALSE]
  if (nrow(r) == 0) stop("Missing selected_method row for target: ", target, call. = FALSE)
  sm <- as.character(r$selected_method[1])
  if (sm == "none") return("Raw")
  if (sm %in% c("platt", "platt_linear")) return("Platt")
  if (sm == "temp") return("Temp")
  stop("Unsupported selected_method token: ", sm, " for target: ", target, call. = FALSE)
}

raw_df <- read.csv(raw_csv, check.names = FALSE)
labels <- read.csv(cfg$bd_labels_csv, check.names = FALSE)
splits <- read.csv(cfg$bd_splits_csv, check.names = FALSE)
cal_obj <- readRDS(cal_rds)
cal_tbl <- cal_obj$summary
cal_fit_scope <- infer_calibration_fit_scope(cal_tbl)
ece_candidates <- sort(list.files(deployment_dir, pattern = "^bd_calibration_ece_.*\\.csv$",
                                  full.names = TRUE), decreasing = TRUE)
if (length(ece_candidates) == 0) stop("No ECE summary CSV found in ", deployment_dir, call. = FALSE)
summary_csv <- ece_candidates[1]
cat("Using ECE summary:", summary_csv, "\n")
sel_tbl <- read.csv(summary_csv, check.names = FALSE)

raw_df$genome <- basename(raw_df$file)
labels$genome <- labels$file
splits$genome <- splits$file

raw_eval <- merge(raw_df, labels, by = "genome", suffixes = c("_pred", "_true"))
raw_eval <- merge(raw_eval, splits[, c("genome", "type")], by = "genome") %>%
  filter(type %in% keep_split)

if (split_target %in% c("test", "lofo")) {
  tax_tbl <- cal_obj$genome_taxonomy
  tax_tbl$genome <- basename(as.character(tax_tbl$file))
  keep_genomes_tax <- unique(tax_tbl$genome[!is.na(tax_tbl$Phylum)])
  n_before <- nrow(raw_eval)
  raw_eval <- raw_eval %>% filter(genome %in% keep_genomes_tax)
  taxonomy_filter <- "phylum_non_na"
  cat(
    "Applied taxonomy filter (Phylum non-NA) from calibration RDS:",
    n_before, "->", nrow(raw_eval), "\n"
  )
} else {
  taxonomy_filter <- "none"
}

heads <- list(
  motility = list(kind = "binary", target = "motility", label_col = "is_motile_pred",
                  val_col = "is_motile_val", val_cal_col = "is_motile_val_cal",
                  prob_col = "is_motile_prob", prob_cal_col = "is_motile_prob",
                  true_col = "is_motile_true"),
  spore = list(kind = "binary", target = "spore", label_col = "ability_spore_pred",
               val_col = "ability_spore_val", val_cal_col = "ability_spore_val_cal",
               prob_col = "ability_spore_prob", prob_cal_col = "ability_spore_prob",
               true_col = "ability_spore_true"),
  oxygen_microaerophile = list(kind = "binary", target = "oxygen_microaerophile",
                               label_col = "oxygen_microaerophile",
                               val_col = "oxygen_microaerophile_val", val_cal_col = "oxygen_microaerophile_val_cal",
                               prob_col = "oxygen_microaerophile_prob", prob_cal_col = "oxygen_microaerophile_prob",
                               true_col = "microaerophile"),
  oxygen_growth = list(kind = "soft2", target = "oxygen_growth",
                       label_col = "oxygen_growth",
                       labels = c("aerobe", "anaerobe"),
                       val_col = "oxygen_growth_val", val_cal_col = "oxygen_growth_val_cal",
                       prob_col = "oxygen_growth_prob", prob_cal_col = "oxygen_growth_prob",
                       true_cols = c("aerobe", "anaerobe")),
  oxygen_facultative = list(kind = "soft2", target = "oxygen_facultative",
                            label_col = "oxygen_facultative",
                            labels = c("facultative.aerobe", "facultative.anaerobe"),
                            val_col = "oxygen_facultative_val", val_cal_col = "oxygen_facultative_val_cal",
                            prob_col = "oxygen_facultative_prob", prob_cal_col = "oxygen_facultative_prob",
                            true_cols = c("facultative.aerobe", "facultative.anaerobe")),
  oxygen_obligate = list(kind = "soft2", target = "oxygen_obligate",
                         label_col = "oxygen_obligate",
                         labels = c("obligate.aerobe", "obligate.anaerobe"),
                         val_col = "oxygen_obligate_val", val_cal_col = "oxygen_obligate_val_cal",
                         prob_col = "oxygen_obligate_prob", prob_cal_col = "oxygen_obligate_prob",
                         true_cols = c("obligate.aerobe", "obligate.anaerobe")),
  gram = list(kind = "multi", target = "gram",
              label_col = "gram",
              labels = c("positive", "variable", "negative"),
              val_col = "gram_val", val_cal_col = "gram_val_cal",
              prob_col = "gram_prob", prob_cal_col = "gram_prob",
              true_cols = c("is_gram_stain_positive", "is_gram_stain_variable", "is_gram_stain_negative")),
  pathogenicity_human = list(kind = "patho", target = "pathogenicity_human",
                             val_col = "pathogenicity_human_val",
                             prob_col = "pathogenicity_human_prob",
                             label_col = "pathogenicity_human",
                             true_col = "pathogenicity_human_true",
                             threshold = 1.4),
  pathogenicity_animal = list(kind = "patho", target = "pathogenicity_animal",
                              val_col = "pathogenicity_animal_val",
                              prob_col = "pathogenicity_animal_prob",
                              label_col = "pathogenicity_animal",
                              true_col = "pathogenicity_animal_true",
                              threshold = 1.4),
  pathogenicity_plant = list(kind = "patho", target = "pathogenicity_plant",
                             val_col = "pathogenicity_plant_val",
                             prob_col = "pathogenicity_plant_prob",
                             label_col = "pathogenicity_plant",
                             true_col = "pathogenicity_plant_true",
                             threshold = 1.0)
)

summary_rows <- list()
plot_rows <- list()
selected_method_map <- list()

for (name in names(heads)) {
  h <- heads[[name]]
  crow <- cal_tbl[cal_tbl$target == h$target, , drop = FALSE]
  if (nrow(crow) == 0) stop("Missing calibration row for target: ", h$target, call. = FALSE)
  T <- pick_param(crow, "T", split_target, param_source)
  a <- pick_param(crow, "a", split_target, param_source)
  b <- pick_param(crow, "b", split_target, param_source)
  sel_method <- choose_selected_method(name, crow, split_target, sel_tbl)

  conf_methods <- list()
  event_methods <- list()
  acc_ref <- NA_real_
  n_ref <- 0L

  if (h$kind == "binary") {
    y <- as.numeric(raw_eval[[h$true_col]])
    p_raw <- clip01(as.numeric(raw_eval[[h$val_col]]))
    pred <- as_bool(raw_eval[[h$label_col]])
    ok <- !is.na(y) & y != -9999 & !is.na(pred) & is.finite(p_raw)
    y <- y[ok]; p_raw <- p_raw[ok]; pred <- pred[ok]
    correct <- as.numeric(pred == y)
    n_ref <- length(y); acc_ref <- mean(correct)

    conf_methods$Raw <- ifelse(pred, p_raw, 1 - p_raw)
    event_methods$Raw <- list(p = p_raw, y = y)
    if (is.finite(T) && T > 0) {
      p_temp <- clip01(sigmoid(logit(p_raw) / T))
      conf_methods$Temp <- ifelse(pred, p_temp, 1 - p_temp)
      event_methods$Temp <- list(p = p_temp, y = y)
    }
    if (is.finite(a) && is.finite(b)) {
      p_platt <- clip01(sigmoid(a * logit(p_raw) + b))
      conf_methods$Platt <- ifelse(pred, p_platt, 1 - p_platt)
      event_methods$Platt <- list(p = p_platt, y = y)
    }
    correct_ref <- correct
  } else if (h$kind == "soft2") {
    pmat <- parse_prob_vec(raw_eval[[h$val_col]])
    ymat <- as.matrix(raw_eval[, h$true_cols]); storage.mode(ymat) <- "numeric"
    pred_idx <- match(as.character(raw_eval[[h$label_col]]), h$labels)
    ok <- apply(ymat, 1, function(r) !any(is.na(r) | r == -9999)) & rowSums(ymat) > 0 &
      !is.na(pred_idx)
    pmat <- pmat[ok, , drop = FALSE]
    ymat <- ymat[ok, , drop = FALSE]
    pred_idx <- pred_idx[ok]
    true_idx <- apply(ymat, 1, which.max)
    y_ev <- as.numeric(true_idx == 2)
    correct <- as.numeric(pred_idx == true_idx)
    n_ref <- length(y_ev); acc_ref <- mean(correct)

    conf_methods$Raw <- pmat[cbind(seq_len(nrow(pmat)), pred_idx)]
    event_methods$Raw <- list(p = pmat[, 2], y = y_ev)
    if (is.finite(T) && T > 0) {
      p_temp <- softmax_from_probs_with_temp(pmat, T)
      conf_methods$Temp <- p_temp[cbind(seq_len(nrow(p_temp)), pred_idx)]
      event_methods$Temp <- list(p = p_temp[, 2], y = y_ev)
    }
    if (is.finite(a) && is.finite(b)) {
      p2 <- clip01(sigmoid(a * logit(pmat[, 2]) + b))
      p_platt <- cbind(1 - p2, p2)
      conf_methods$Platt <- p_platt[cbind(seq_len(nrow(p_platt)), pred_idx)]
      event_methods$Platt <- list(p = p2, y = y_ev)
    }
    correct_ref <- correct
  } else if (h$kind == "multi") {
    pmat <- parse_prob_vec(raw_eval[[h$val_col]])
    ymat <- as.matrix(raw_eval[, h$true_cols]); storage.mode(ymat) <- "numeric"
    pred_idx <- match(as.character(raw_eval[[h$label_col]]), h$labels)
    ok <- apply(ymat, 1, function(r) !any(is.na(r) | r == -9999)) & rowSums(ymat) > 0 &
      !is.na(pred_idx)
    pmat <- pmat[ok, , drop = FALSE]
    ymat <- ymat[ok, , drop = FALSE]
    pred_idx <- pred_idx[ok]
    true_idx <- apply(ymat, 1, which.max)
    correct <- as.numeric(pred_idx == true_idx)
    n_ref <- nrow(pmat); acc_ref <- mean(correct)

    conf_methods$Raw <- pmat[cbind(seq_len(nrow(pmat)), pred_idx)]
    event_methods$Raw <- list(p = as.vector(pmat), y = as.vector(ymat))
    if (is.finite(T) && T > 0) {
      p_temp <- softmax_from_probs_with_temp(pmat, T)
      conf_methods$Temp <- p_temp[cbind(seq_len(nrow(p_temp)), pred_idx)]
      event_methods$Temp <- list(p = as.vector(p_temp), y = as.vector(ymat))
    }
    correct_ref <- correct
  } else if (h$kind == "patho") {
    raw_val <- as.numeric(raw_eval[[h$val_col]])
    y_raw <- as.numeric(raw_eval[[h$true_col]])
    ok <- !is.na(y_raw) & y_raw != -9999 & is.finite(raw_val)
    raw_val <- raw_val[ok]
    y <- as.numeric(y_raw[ok] > 0)
    pred_raw <- as.numeric(raw_val > h$threshold)
    correct <- as.numeric(pred_raw == y)
    n_ref <- length(y); acc_ref <- mean(correct)

    conf_none <- (tanh(5 * abs(raw_val - h$threshold)) + 1) / 2
    conf_methods$Raw <- conf_none
    event_methods$Raw <- list(p = clip01(sigmoid(raw_val - h$threshold)), y = y)
    if (is.finite(a) && is.finite(b)) {
      p_platt <- clip01(sigmoid(a * raw_val + b))
      conf_platt <- ifelse(pred_raw == 1, p_platt, 1 - p_platt)
      conf_methods$Platt <- conf_platt
      event_methods$Platt <- list(p = p_platt, y = y)
    }
    correct_ref <- correct
  } else {
    stop("Unsupported head kind: ", h$kind, " for target: ", name, call. = FALSE)
  }

  if (!(sel_method %in% names(conf_methods))) {
    stop("Selected method ", sel_method, " unavailable for target: ", name, call. = FALSE)
  }
  conf_methods$Selected <- conf_methods[[sel_method]]
  if (sel_method %in% names(event_methods)) {
    event_methods$Selected <- event_methods[[sel_method]]
  }
  selected_method_map[[name]] <- sel_method

  conf_available <- intersect(method_order, names(conf_methods))
  event_available <- intersect(method_order, names(event_methods))

  for (m in conf_available) {
    e_conf <- safe_ece(conf_methods[[m]], correct_ref)
    if (m %in% event_available) {
      e_event <- safe_ece(event_methods[[m]]$p, event_methods[[m]]$y)
      n_event <- length(event_methods[[m]]$y)
    } else {
      e_event <- NA_real_
      n_event <- NA_real_
    }

    summary_rows[[paste(name, m, sep = "_")]] <- data.frame(
      target = name,
      method = m,
      N_conf = n_ref,
      Acc = acc_ref,
      ECE_conf = e_conf,
      N_event = n_event,
      ECE_event = e_event,
      stringsAsFactors = FALSE
    )
  }

  rd_conf <- bind_rows(lapply(conf_available, function(m) {
    rr <- safe_rel(conf_methods[[m]], correct_ref)
    if (is.null(rr)) return(NULL)
    rr %>% mutate(method = m, view = "Confidence")
  }))

  if (length(event_available) > 0) {
    rd_event <- bind_rows(lapply(event_available, function(m) {
      rr <- safe_rel(event_methods[[m]]$p, event_methods[[m]]$y)
      if (is.null(rr)) return(NULL)
      rr %>% mutate(method = m, view = "Event")
    }))
  } else {
    rd_event <- NULL
  }

  if (!is.null(rd_conf)) plot_rows[[paste0(name, "_conf")]] <- rd_conf %>% mutate(target = name)
  if (!is.null(rd_event)) plot_rows[[paste0(name, "_event")]] <- rd_event %>% mutate(target = name)
}

summary_df <- bind_rows(summary_rows) %>%
  mutate(
    method = factor(method, levels = method_order),
    eval_split = split_target,
    param_source = param_source,
    calibration_fit_scope = cal_fit_scope,
    calibration_rds = basename(cal_rds),
    taxonomy_filter = taxonomy_filter
  ) %>%
  arrange(target, method)

split_tag <- if (split_target == "validation") "val" else split_target
pdf_path <- file.path(cfg$output_dir, paste0("deployed_ece_curves_all_methods_", split_tag, ".pdf"))
csv_path <- file.path(cfg$output_dir, paste0("deployed_ece_all_methods_", split_tag, ".csv"))
write.csv(summary_df, csv_path, row.names = FALSE)

pdf(pdf_path, width = 14, height = 7, onefile = TRUE)

theme_tbl <- ttheme_minimal(
  core = list(fg_params = list(fontsize = 8, hjust = 0, x = 0.02),
              bg_params = list(fill = "white")),
  colhead = list(fg_params = list(fontsize = 8, fontface = "bold", hjust = 0, x = 0.02),
                 bg_params = list(fill = "#E8E8E8"))
)

grid::grid.newpage()
grid::grid.text("Deployed Calibration Audit: Confidence and Event ECE (All Methods)",
                x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 14, fontface = "bold"))
grid::grid.text(
  paste0(
    "Split: ", split_target,
    " | Calibration fit: ", cal_fit_scope,
    " | Param source: ", param_source,
    " | Taxonomy filter: ", taxonomy_filter,
    " | RDS: ", basename(cal_rds)
  ),
  x = 0.5, y = 0.90, gp = grid::gpar(fontsize = 9, col = "grey40")
)
tg <- tableGrob(summary_df %>% mutate(across(where(is.numeric), ~round(., 4))),
                rows = NULL, theme = theme_tbl)
tg$widths <- grid::unit(rep(1 / ncol(summary_df), ncol(summary_df)), "npc")
grid::pushViewport(grid::viewport(y = 0.42, height = 0.75))
grid::grid.draw(tg)
grid::popViewport()

for (tgt in unique(summary_df$target)) {
  dfc <- bind_rows(plot_rows[paste0(tgt, "_conf")])
  dfe <- bind_rows(plot_rows[paste0(tgt, "_event")])
  if (!is.null(dfc) && nrow(dfc) > 0) dfc <- dfc %>% filter(method != "Selected")
  if (!is.null(dfe) && nrow(dfe) > 0) dfe <- dfe %>% filter(method != "Selected")

  sel_method <- selected_method_map[[tgt]]

  p_conf <- ggplot(dfc, aes(x = mean_pred, y = obs_freq, color = method, size = n)) +
    geom_point(alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    xlim(0, 1) + ylim(0, 1) +
    labs(
      title = paste0(tgt, " - Confidence (", split_target, ")"),
      subtitle = paste0(
        "Selected method: ", sel_method,
        " | Param source: ", param_source,
        " | Cal fit: ", cal_fit_scope
      ),
      x = "Predicted confidence", y = "Observed correctness"
    ) +
    theme_minimal() + theme(legend.position = "bottom")

  p_event <- ggplot(dfe, aes(x = mean_pred, y = obs_freq, color = method, size = n)) +
    geom_point(alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    xlim(0, 1) + ylim(0, 1) +
    labs(
      title = paste0(tgt, " - Event Probability (", split_target, ")"),
      subtitle = paste0(
        "Selected method: ", sel_method,
        " | Param source: ", param_source,
        " | Cal fit: ", cal_fit_scope
      ),
      x = "Predicted event probability", y = "Observed frequency"
    ) +
    theme_minimal() + theme(legend.position = "bottom")

  grid.arrange(p_conf, p_event, ncol = 2)
}

dev.off()

# Canonical unsuffixed names represent LOFO pooled val+test output.
if (split_target == "lofo") {
  file.copy(pdf_path, file.path(cfg$output_dir, "deployed_ece_curves_all_methods.pdf"), overwrite = TRUE)
  file.copy(csv_path, file.path(cfg$output_dir, "deployed_ece_all_methods.csv"), overwrite = TRUE)
}

cat("Saved:\n")
cat("  PDF:", pdf_path, "\n")
cat("  CSV:", csv_path, "\n")
