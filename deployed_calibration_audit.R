# ── Deployed calibration audit ─────────────────────────────────────────────────
# How well-calibrated are the _prob columns (deployed on BacDive website)?
# These use custom transformations from eval_non_train_self_pad.R:
#   Binary: ifelse(m > 0.5, m, 1-m)  → confidence (always >= 0.5)
#   Multiclass ≥3: 0.5 + 0.5 * log(N * max(m)) / log(N)  → log-scaled
#   2-class softmax: max(m)  → raw max
#   Pathogenicity: (tanh(5 * abs(m - threshold)) + 1) / 2  → tanh-distance
#   Biosafety: threshold-based with tanh
# Usage: conda run -n genome Rscript deployed_calibration_audit.R

library(ggplot2)
library(dplyr)
library(gridExtra)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg)) else normalizePath(getwd())
proj_dir <- if (dir.exists(script_path)) script_path else dirname(script_path)
source(file.path(proj_dir, "R", "config_paths.R"))
source(file.path(proj_dir, "R", "calibration_common.R"))
cfg <- load_bacdive_config(
  proj_dir,
  required = c("bd_pred_csv", "bd_labels_csv", "bd_splits_csv")
)

# ── Load data ─────────────────────────────────────────────────────────────────
bdpred <- read.csv(cfg$bd_pred_csv)
labels_true <- read.csv(cfg$bd_labels_csv)
splits <- read.csv(cfg$bd_splits_csv)

bdpred$genome <- basename(bdpred$file)
labels_true$genome <- labels_true$file
splits$genome <- splits$file

df <- merge(bdpred, labels_true, by = "genome", suffixes = c("_pred", "_true"))
df <- merge(df, splits[, c("genome", "type")], by = "genome")

test_df <- df %>% filter(type == "test")
val_df  <- df %>% filter(type == "validation")
cat("Test genomes:", nrow(test_df), "\n")
cat("Val genomes:", nrow(val_df), "\n\n")

# ── Audit each head ──────────────────────────────────────────────────────────
results <- list()
plot_list <- list()

# --------------------------------------------------------------------------
# 1. Binary sigmoid heads
#    _val = raw sigmoid output (positive-class prob)
#    _prob = max(p, 1-p) = deployed confidence (always >= 0.5)
#    Prediction: y_hat = (p_val > 0.5)
#    Correct: y_hat == y_true
#    Raw conf: max(p_val, 1-p_val) — same as _prob for binary
# --------------------------------------------------------------------------
binary_heads <- list(
  motility = list(val = "is_motile_val", prob = "is_motile_prob",
                  true_col = "is_motile_true"),
  oxygen_microaerophile = list(val = "oxygen_microaerophile_val",
                               prob = "oxygen_microaerophile_prob",
                               true_col = "microaerophile"),
  spore = list(val = "ability_spore_val", prob = "ability_spore_prob",
               true_col = "ability_spore_true")
)

for (name in names(binary_heads)) {
  h <- binary_heads[[name]]
  p_raw <- as.numeric(test_df[[h$val]])
  p_prob <- as.numeric(test_df[[h$prob]])
  y <- as.numeric(test_df[[h$true_col]])
  ok <- !is.na(y) & y != -9999
  p_raw <- p_raw[ok]; p_prob <- p_prob[ok]; y <- y[ok]

  y_hat <- as.numeric(p_raw > 0.5)
  correct <- as.numeric(y_hat == y)
  raw_conf <- pmax(p_raw, 1 - p_raw)
  acc <- mean(correct)

  ece_raw <- ece(raw_conf, correct)
  ece_dep <- ece(p_prob, correct)

  results[[name]] <- data.frame(
    target = name, type = "binary_sigmoid",
    N = length(y), Acc = round(acc, 3),
    ECE_raw_conf = round(ece_raw, 4),
    ECE_deployed = round(ece_dep, 4),
    stringsAsFactors = FALSE)

  rd <- bind_rows(
    reliability_data(raw_conf, correct) %>% mutate(method = "max(p, 1-p)"),
    reliability_data(p_prob, correct) %>% mutate(method = "Deployed _prob")
  )
  rd$method <- factor(rd$method, levels = unique(rd$method))
  plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                       color = method, size = n)) +
    geom_point(alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = paste0(name, " [binary]  N=", length(y), " Acc=", round(acc, 3)),
         subtitle = paste0("ECE raw: ", round(ece_raw, 4),
                           " | ECE deployed: ", round(ece_dep, 4)),
         x = "Confidence", y = "Accuracy") +
    theme_minimal() + theme(legend.position = "bottom")
}

# --------------------------------------------------------------------------
# 2. Two-class softmax heads (oxygen_growth, oxygen_facultative, oxygen_obligate)
#    _val = "p0, p1" softmax probabilities
#    _prob = max(p0, p1) = deployed confidence
#    Prediction: which.max
#    Raw conf: max(softmax)
# --------------------------------------------------------------------------
softmax2_heads <- list(
  oxygen_growth = list(val = "oxygen_growth_val", prob = "oxygen_growth_prob",
                       true_cols = c("aerobe", "anaerobe")),
  oxygen_facultative = list(val = "oxygen_facultative_val",
                            prob = "oxygen_facultative_prob",
                            true_cols = c("facultative.aerobe",
                                          "facultative.anaerobe")),
  oxygen_obligate = list(val = "oxygen_obligate_val",
                         prob = "oxygen_obligate_prob",
                         true_cols = c("obligate.aerobe", "obligate.anaerobe"))
)

for (name in names(softmax2_heads)) {
  h <- softmax2_heads[[name]]
  pmat <- parse_prob_vec(test_df[[h$val]])
  p_prob <- as.numeric(test_df[[h$prob]])
  ymat <- as.matrix(test_df[, h$true_cols])
  storage.mode(ymat) <- "numeric"
  ok <- apply(ymat, 1, function(r) !any(is.na(r) | r == -9999)) &
    rowSums(ymat) > 0
  pmat <- pmat[ok, ]; p_prob <- p_prob[ok]; ymat <- ymat[ok, ]

  y_hat <- apply(pmat, 1, which.max)
  y_true <- apply(ymat, 1, which.max)
  correct <- as.numeric(y_hat == y_true)
  raw_conf <- apply(pmat, 1, max)
  acc <- mean(correct)

  ece_raw <- ece(raw_conf, correct)
  ece_dep <- ece(p_prob, correct)

  results[[name]] <- data.frame(
    target = name, type = "softmax_2class",
    N = nrow(pmat), Acc = round(acc, 3),
    ECE_raw_conf = round(ece_raw, 4),
    ECE_deployed = round(ece_dep, 4),
    stringsAsFactors = FALSE)

  rd <- bind_rows(
    reliability_data(raw_conf, correct) %>% mutate(method = "max(softmax)"),
    reliability_data(p_prob, correct) %>% mutate(method = "Deployed _prob")
  )
  rd$method <- factor(rd$method, levels = unique(rd$method))
  plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                       color = method, size = n)) +
    geom_point(alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = paste0(name, " [2-class]  N=", nrow(pmat),
                        " Acc=", round(acc, 3)),
         subtitle = paste0("ECE raw: ", round(ece_raw, 4),
                           " | ECE deployed: ", round(ece_dep, 4)),
         x = "Confidence", y = "Accuracy") +
    theme_minimal() + theme(legend.position = "bottom")
}

# --------------------------------------------------------------------------
# 3. Multiclass softmax heads (gram, cellshape, flagellum)
#    _val = "p0, p1, ..." softmax probabilities
#    _prob = 0.5 + 0.5 * log(N * max(p)) / log(N)  (log-scaled)
#    Prediction: which.max
#    Raw conf: max(softmax)
# --------------------------------------------------------------------------
multi_heads <- list(
  gram = list(val = "gram_val", prob = "gram_prob",
              true_cols = c("is_gram_stain_positive", "is_gram_stain_variable",
                            "is_gram_stain_negative")),
  cellshape = list(val = "cellshape_val", prob = "cellshape_prob",
                   true_cols = c("is_cell_shape_rod.shaped",
                     "is_cell_shape_other", "is_cell_shape_coccus.shaped",
                     "is_cell_shape_vibrio.shaped",
                     "is_cell_shape_filament.shaped",
                     "is_cell_shape_sphere.shaped",
                     "is_cell_shape_ovoid.shaped",
                     "is_cell_shape_pleomorphic.shaped",
                     "is_cell_shape_spiral.shaped",
                     "is_cell_shape_curved.shaped",
                     "is_cell_shape_oval.shaped")),
  flagellum = list(val = "flagellum_val", prob = "flagellum_prob",
                   true_cols = c("is_flagellum_arrangement_monotrichous",
                     "is_flagellum_arrangement_monotrichous_polar",
                     "is_flagellum_arrangement_polar",
                     "is_flagellum_arrangement_peritrichous",
                     "is_flagellum_arrangement_lophotrichous",
                     "is_flagellum_arrangement_gliding"))
)

for (name in names(multi_heads)) {
  h <- multi_heads[[name]]
  pmat <- parse_prob_vec(test_df[[h$val]])
  p_prob <- as.numeric(test_df[[h$prob]])
  ymat <- as.matrix(test_df[, h$true_cols])
  storage.mode(ymat) <- "numeric"
  ok <- apply(ymat, 1, function(r) !any(is.na(r) | r == -9999)) &
    rowSums(ymat) > 0
  pmat <- pmat[ok, ]; p_prob <- p_prob[ok]; ymat <- ymat[ok, ]

  y_hat <- apply(pmat, 1, which.max)
  y_true <- apply(ymat, 1, which.max)
  correct <- as.numeric(y_hat == y_true)
  raw_conf <- apply(pmat, 1, max)
  acc <- mean(correct)

  N_classes <- ncol(pmat)
  ece_raw <- ece(raw_conf, correct)
  ece_dep <- ece(p_prob, correct)

  results[[name]] <- data.frame(
    target = name, type = paste0("multiclass_", N_classes),
    N = nrow(pmat), Acc = round(acc, 3),
    ECE_raw_conf = round(ece_raw, 4),
    ECE_deployed = round(ece_dep, 4),
    stringsAsFactors = FALSE)

  rd <- bind_rows(
    reliability_data(raw_conf, correct) %>% mutate(method = "max(softmax)"),
    reliability_data(p_prob, correct) %>%
      mutate(method = paste0("Deployed (log-scaled, K=", N_classes, ")"))
  )
  rd$method <- factor(rd$method, levels = unique(rd$method))
  plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                       color = method, size = n)) +
    geom_point(alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = paste0(name, " [", N_classes, "-class]  N=", nrow(pmat),
                        " Acc=", round(acc, 3)),
         subtitle = paste0("ECE raw: ", round(ece_raw, 4),
                           " | ECE deployed: ", round(ece_dep, 4)),
         x = "Confidence", y = "Accuracy") +
    theme_minimal() + theme(legend.position = "bottom")
}

# --------------------------------------------------------------------------
# 4. Pathogenicity heads
#    _val = raw linear output
#    _prob = (tanh(5 * abs(m - threshold)) + 1) / 2 where threshold differs
#    Prediction: m > threshold (matches deployment code)
#    We compare threshold-aware raw confidence vs deployed _prob
# --------------------------------------------------------------------------
patho_heads <- list(
  pathogenicity_human  = list(val = "pathogenicity_human_val",
                              prob = "pathogenicity_human_prob",
                              true_col = "pathogenicity_human_true",
                              threshold = 1.4),
  pathogenicity_animal = list(val = "pathogenicity_animal_val",
                              prob = "pathogenicity_animal_prob",
                              true_col = "pathogenicity_animal_true",
                              threshold = 1.4),
  pathogenicity_plant  = list(val = "pathogenicity_plant_val",
                              prob = "pathogenicity_plant_prob",
                              true_col = "pathogenicity_plant_true",
                              threshold = 1.0)
)

for (name in names(patho_heads)) {
  h <- patho_heads[[name]]
  raw <- as.numeric(test_df[[h$val]])
  p_prob <- as.numeric(test_df[[h$prob]])
  y_raw <- as.numeric(test_df[[h$true_col]])
  ok <- !is.na(y_raw) & y_raw != -9999
  raw <- raw[ok]; p_prob <- p_prob[ok]; y_raw <- y_raw[ok]
  y <- as.numeric(y_raw > 0)

  # Prediction threshold mirrors deployment logic.
  y_hat <- as.numeric(raw > h$threshold)
  correct <- as.numeric(y_hat == y)
  # Raw confidence: distance from threshold in shifted-sigmoid space.
  raw_conf <- abs(sigmoid(raw - h$threshold) - 0.5) * 2
  raw_conf <- raw_conf * 0.5 + 0.5
  acc <- mean(correct)

  ece_raw <- ece(raw_conf, correct)
  ece_dep <- ece(p_prob, correct)

  results[[name]] <- data.frame(
    target = name, type = "linear_patho",
    N = length(y), Acc = round(acc, 3),
    ECE_raw_conf = round(ece_raw, 4),
    ECE_deployed = round(ece_dep, 4),
    stringsAsFactors = FALSE)

  rd <- bind_rows(
    reliability_data(raw_conf, correct) %>%
      mutate(method = "sigmoid(|raw-thr|) conf"),
    reliability_data(p_prob, correct) %>%
      mutate(method = "Deployed tanh-prob")
  )
  rd$method <- factor(rd$method, levels = unique(rd$method))
  plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                       color = method, size = n)) +
    geom_point(alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = paste0(name, " [linear->patho]  N=", length(y),
                        " Acc=", round(acc, 3)),
         subtitle = paste0("ECE raw: ", round(ece_raw, 4),
                           " | ECE deployed: ", round(ece_dep, 4)),
         x = "Confidence", y = "Accuracy") +
    theme_minimal() + theme(legend.position = "bottom")
}

# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════
summary_df <- do.call(rbind, results) %>% arrange(desc(Acc))
cat("\n=== Deployed Calibration Audit (Test Set) ===\n")
print(summary_df, row.names = FALSE)

cat("\nMean ECE raw confidence:", round(mean(summary_df$ECE_raw_conf), 4), "\n")
cat("Mean ECE deployed _prob:", round(mean(summary_df$ECE_deployed), 4), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# Save PDF
# ══════════════════════════════════════════════════════════════════════════════
pdf_path <- file.path(cfg$reports_dir, "deployed_calibration_audit.pdf")
pdf(pdf_path, width = 14, height = 5, onefile = TRUE)

# Page 1: summary table
theme_tbl <- ttheme_minimal(
  core = list(fg_params = list(fontsize = 8, hjust = 0, x = 0.02),
              bg_params = list(fill = "white")),
  colhead = list(fg_params = list(fontsize = 8, fontface = "bold",
                                  hjust = 0, x = 0.02),
                 bg_params = list(fill = "#E8E8E8"))
)
grid::grid.newpage()
grid::grid.text("Deployed Calibration Audit (Test Set)",
                x = 0.5, y = 0.95,
                gp = grid::gpar(fontsize = 14, fontface = "bold"))
grid::grid.text(
  paste0("Comparing raw model confidence vs deployed _prob transformations | ",
         "Test N=", nrow(test_df)),
  x = 0.5, y = 0.90, gp = grid::gpar(fontsize = 9, col = "grey40"))
tg <- tableGrob(summary_df %>%
                  mutate(across(where(is.numeric), ~round(., 4))),
                rows = NULL, theme = theme_tbl)
tg$widths <- unit(rep(1 / ncol(summary_df), ncol(summary_df)), "npc")
grid::pushViewport(grid::viewport(y = 0.42, height = 0.75))
grid::grid.draw(tg)
grid::popViewport()

# Reliability diagrams
for (p in plot_list) print(p)
dev.off()

cat("\nPDF saved to:", pdf_path, "\n")
