# ── Calibration: within-test 70/30 split vs cross-population (val→test) ───────
# The val set shares 90.6% genus overlap with train (data leakage), while test
# is a proper genus-level holdout (75% unseen genera). We therefore fit
# calibration on 70% of the test set and evaluate on the remaining 30%.
# Val-based calibration is shown alongside for comparison.
# Usage: conda run -n genome Rscript R/calibration_val_test.R

library(ggplot2)
library(dplyr)
library(gridExtra)

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

ACC_THRESHOLD <- 0.6
N_FILTER <- 100
CAL_FRAC <- 0.7
SEED <- 42
ENABLE_DIRICHLET <- tolower(Sys.getenv("ENABLE_DIRICHLET_MULTICLASS", "false")) %in%
  c("1", "true", "yes")

# ── Load predictions and ground truth ────────────────────────────────────────
bdpred <- read.csv(cfg$bd_pred_csv)
labels_true <- read.csv(cfg$bd_labels_csv)
splits <- read.csv(cfg$bd_splits_csv)

bdpred$genome <- basename(bdpred$file)
labels_true$genome <- labels_true$file
splits$genome <- splits$file

df <- merge(bdpred, labels_true, by = "genome", suffixes = c("_pred", "_true"))
df <- merge(df, splits[, c("genome", "type")], by = "genome")

val_df  <- df %>% filter(type == "validation")
test_df <- df %>% filter(type == "test")

# ── Split test set 70/30 ─────────────────────────────────────────────────────
set.seed(SEED)
n_test <- nrow(test_df)
cal_idx <- sample(n_test, floor(CAL_FRAC * n_test))
test_cal  <- test_df[cal_idx, ]   # 70% — fit calibration
test_eval <- test_df[-cal_idx, ]  # 30% — evaluate

cat("Validation genomes:", nrow(val_df), "\n")
cat("Test genomes total:", n_test, "\n")
cat("  Test-cal (fit):", nrow(test_cal), "\n")
cat("  Test-eval (eval):", nrow(test_eval), "\n")
cat("Calibration fraction:", CAL_FRAC, "| Seed:", SEED, "\n")
cat("Accuracy threshold:", ACC_THRESHOLD, "\n\n")

# ── Collect results ──────────────────────────────────────────────────────────
all_rows <- list()
plot_list <- list()

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Binary sigmoid heads
# ══════════════════════════════════════════════════════════════════════════════

binary_heads <- list(
  motility = list(
    pred_val = "is_motile_val", true_col = "is_motile_true"),
  oxygen_microaerophile = list(
    pred_val = "oxygen_microaerophile_val", true_col = "microaerophile"),
  spore = list(
    pred_val = "ability_spore_val", true_col = "ability_spore_true")
)

for (name in names(binary_heads)) {
  h <- binary_heads[[name]]

  # Val (cross-population fit)
  p_v <- as.numeric(val_df[[h$pred_val]])
  y_v <- as.numeric(val_df[[h$true_col]])
  ok_v <- !is.na(y_v) & y_v != -9999
  p_v <- p_v[ok_v]; y_v <- y_v[ok_v]

  # Test-cal (within-population fit)
  p_c <- as.numeric(test_cal[[h$pred_val]])
  y_c <- as.numeric(test_cal[[h$true_col]])
  ok_c <- !is.na(y_c) & y_c != -9999
  p_c <- p_c[ok_c]; y_c <- y_c[ok_c]

  # Test-eval (held-out evaluation)
  p_e <- as.numeric(test_eval[[h$pred_val]])
  y_e <- as.numeric(test_eval[[h$true_col]])
  ok_e <- !is.na(y_e) & y_e != -9999
  p_e <- p_e[ok_e]; y_e <- y_e[ok_e]

  z_v <- logit(p_v); z_c <- logit(p_c); z_e <- logit(p_e)

  # Fit on val (cross-population)
  val_T <- optimize(nll_binary_temp, c(0.01, 20), logits = z_v, y = y_v)$minimum
  val_P <- optim(c(1, 0), nll_binary_platt, logits = z_v, y = y_v,
                 method = "BFGS")$par

  # Fit on test-cal (within-population)
  tc_T <- optimize(nll_binary_temp, c(0.01, 20), logits = z_c, y = y_c)$minimum
  tc_P <- optim(c(1, 0), nll_binary_platt, logits = z_c, y = y_c,
                method = "BFGS")$par

  # Evaluate on test-eval
  p_e_uncal    <- p_e
  p_e_val_temp <- sigmoid(z_e / val_T)
  p_e_val_pla  <- sigmoid(val_P[1] * z_e + val_P[2])
  p_e_tc_temp  <- sigmoid(z_e / tc_T)
  p_e_tc_pla   <- sigmoid(tc_P[1] * z_e + tc_P[2])

  acc <- mean((p_e > 0.5) == y_e)

  row <- data.frame(
    target = name, type = "binary_sigmoid",
    N_cal = length(y_c), N_eval = length(y_e), Acc = round(acc, 3),
    ECE_uncal = round(ece(p_e, y_e), 4),
    ECE_val_temp  = round(ece(p_e_val_temp, y_e), 4),
    ECE_val_platt = round(ece(p_e_val_pla, y_e), 4),
    ECE_val_dirichlet = NA,
    ECE_tc_temp   = round(ece(p_e_tc_temp, y_e), 4),
    ECE_tc_platt  = round(ece(p_e_tc_pla, y_e), 4),
    ECE_tc_dirichlet = NA,
    NLL_uncal = round(nll_binary_temp(1, z_e, y_e), 4),
    NLL_val_platt = round(nll_binary_platt(val_P, z_e, y_e), 4),
    NLL_tc_platt  = round(nll_binary_platt(tc_P, z_e, y_e), 4),
    NLL_val_dirichlet = NA,
    NLL_tc_dirichlet = NA,
    val_T = round(val_T, 3), val_a = round(val_P[1], 3), val_b = round(val_P[2], 3),
    tc_T = round(tc_T, 3), tc_a = round(tc_P[1], 3), tc_b = round(tc_P[2], 3),
    stringsAsFactors = FALSE)
  all_rows[[name]] <- row

  if (acc >= ACC_THRESHOLD && length(y_e) >= N_FILTER) {
    rd <- bind_rows(
      reliability_data(p_e, y_e) %>% mutate(method = "Uncalibrated"),
      reliability_data(p_e_val_pla, y_e) %>%
        mutate(method = paste0("Val-Platt (a=", round(val_P[1], 2),
                               " b=", round(val_P[2], 2), ")")),
      reliability_data(p_e_tc_pla, y_e) %>%
        mutate(method = paste0("TestCal-Platt (a=", round(tc_P[1], 2),
                               " b=", round(tc_P[2], 2), ")"))
    )
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                        color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, "  (N_eval=", length(y_e),
                          ", Acc=", round(acc, 3), ")"),
           subtitle = paste0("ECE uncal:", round(ece(p_e, y_e), 4),
                             " | val-Platt:", round(ece(p_e_val_pla, y_e), 4),
                             " | tc-Platt:", round(ece(p_e_tc_pla, y_e), 4)),
           x = "Mean Predicted Prob", y = "Observed Frequency") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: 2-class softmax heads
# ══════════════════════════════════════════════════════════════════════════════

softmax2_heads <- list(
  oxygen_growth = list(
    pred_val = "oxygen_growth_val",
    true_cols = c("aerobe", "anaerobe")),
  oxygen_facultative = list(
    pred_val = "oxygen_facultative_val",
    true_cols = c("facultative.aerobe", "facultative.anaerobe")),
  oxygen_obligate = list(
    pred_val = "oxygen_obligate_val",
    true_cols = c("obligate.aerobe", "obligate.anaerobe"))
)

for (name in names(softmax2_heads)) {
  h <- softmax2_heads[[name]]

  extract_s2 <- function(d) {
    pmat <- parse_prob_vec(d[[h$pred_val]])
    y1 <- as.numeric(d[[h$true_cols[2]]])
    ok <- !is.na(y1) & y1 != -9999 &
      !is.na(as.numeric(d[[h$true_cols[1]]])) &
      (as.numeric(d[[h$true_cols[1]]]) + y1) > 0
    list(p = pmat[ok, 2], y = y1[ok])
  }

  dv <- extract_s2(val_df)
  dc <- extract_s2(test_cal)
  de <- extract_s2(test_eval)

  z_v <- logit(dv$p); z_c <- logit(dc$p); z_e <- logit(de$p)

  # Fit on val
  val_T <- optimize(nll_binary_temp, c(0.01, 20), logits = z_v, y = dv$y)$minimum
  val_P <- optim(c(1, 0), nll_binary_platt, logits = z_v, y = dv$y,
                 method = "BFGS")$par
  # Fit on test-cal
  tc_T <- optimize(nll_binary_temp, c(0.01, 20), logits = z_c, y = dc$y)$minimum
  tc_P <- optim(c(1, 0), nll_binary_platt, logits = z_c, y = dc$y,
                method = "BFGS")$par

  p_e_uncal    <- de$p
  p_e_val_pla  <- sigmoid(val_P[1] * z_e + val_P[2])
  p_e_tc_pla   <- sigmoid(tc_P[1] * z_e + tc_P[2])

  acc <- mean((de$p > 0.5) == de$y)

  row <- data.frame(
    target = name, type = "softmax_2class",
    N_cal = length(dc$y), N_eval = length(de$y), Acc = round(acc, 3),
    ECE_uncal = round(ece(de$p, de$y), 4),
    ECE_val_temp  = round(ece(sigmoid(z_e / val_T), de$y), 4),
    ECE_val_platt = round(ece(p_e_val_pla, de$y), 4),
    ECE_val_dirichlet = NA,
    ECE_tc_temp   = round(ece(sigmoid(z_e / tc_T), de$y), 4),
    ECE_tc_platt  = round(ece(p_e_tc_pla, de$y), 4),
    ECE_tc_dirichlet = NA,
    NLL_uncal = round(nll_binary_temp(1, z_e, de$y), 4),
    NLL_val_platt = round(nll_binary_platt(val_P, z_e, de$y), 4),
    NLL_tc_platt  = round(nll_binary_platt(tc_P, z_e, de$y), 4),
    NLL_val_dirichlet = NA,
    NLL_tc_dirichlet = NA,
    val_T = round(val_T, 3), val_a = round(val_P[1], 3), val_b = round(val_P[2], 3),
    tc_T = round(tc_T, 3), tc_a = round(tc_P[1], 3), tc_b = round(tc_P[2], 3),
    stringsAsFactors = FALSE)
  all_rows[[name]] <- row

  if (acc >= ACC_THRESHOLD && length(de$y) >= N_FILTER) {
    rd <- bind_rows(
      reliability_data(de$p, de$y) %>% mutate(method = "Uncalibrated"),
      reliability_data(p_e_val_pla, de$y) %>%
        mutate(method = paste0("Val-Platt (a=", round(val_P[1], 2),
                               " b=", round(val_P[2], 2), ")")),
      reliability_data(p_e_tc_pla, de$y) %>%
        mutate(method = paste0("TestCal-Platt (a=", round(tc_P[1], 2),
                               " b=", round(tc_P[2], 2), ")"))
    )
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                        color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, "  (N_eval=", length(de$y),
                          ", Acc=", round(acc, 3), ")"),
           subtitle = paste0("ECE uncal:", round(ece(de$p, de$y), 4),
                             " | val-Platt:", round(ece(p_e_val_pla, de$y), 4),
                             " | tc-Platt:", round(ece(p_e_tc_pla, de$y), 4)),
           x = "Mean Predicted Prob", y = "Observed Frequency") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Multiclass softmax heads (gram, cellshape, flagellum)
# ══════════════════════════════════════════════════════════════════════════════

multi_heads <- list(
  gram = list(
    pred_val = "gram_val",
    true_cols = c("is_gram_stain_positive", "is_gram_stain_variable",
                  "is_gram_stain_negative")),
  cellshape = list(
    pred_val = "cellshape_val",
    true_cols = c("is_cell_shape_rod.shaped", "is_cell_shape_other",
      "is_cell_shape_coccus.shaped", "is_cell_shape_vibrio.shaped",
      "is_cell_shape_filament.shaped", "is_cell_shape_sphere.shaped",
      "is_cell_shape_ovoid.shaped", "is_cell_shape_pleomorphic.shaped",
      "is_cell_shape_spiral.shaped", "is_cell_shape_curved.shaped",
      "is_cell_shape_oval.shaped")),
  flagellum = list(
    pred_val = "flagellum_val",
    true_cols = c("is_flagellum_arrangement_monotrichous",
      "is_flagellum_arrangement_monotrichous_polar",
      "is_flagellum_arrangement_polar",
      "is_flagellum_arrangement_peritrichous",
      "is_flagellum_arrangement_lophotrichous",
      "is_flagellum_arrangement_gliding"))
)

for (name in names(multi_heads)) {
  h <- multi_heads[[name]]

  extract_mc <- function(d) {
    pmat <- parse_prob_vec(d[[h$pred_val]])
    ymat <- as.matrix(d[, h$true_cols])
    storage.mode(ymat) <- "numeric"
    ok <- apply(ymat, 1, function(r) !any(is.na(r) | r == -9999)) &
      rowSums(ymat) > 0
    list(pmat = pmat[ok, ], ymat = ymat[ok, ])
  }

  dv <- extract_mc(val_df)
  dc <- extract_mc(test_cal)
  de <- extract_mc(test_eval)

  eps <- 1e-7
  zv <- log(pmax(dv$pmat, eps))
  zc <- log(pmax(dc$pmat, eps))
  ze <- log(pmax(de$pmat, eps))

  # Fit on val
  val_T <- optimize(nll_multi_temp, c(0.01, 20),
                    z_mat = zv, y_mat = dv$ymat)$minimum
  if (isTRUE(ENABLE_DIRICHLET)) {
    val_D <- fit_dirichlet_calibration(dv$pmat, dv$ymat)
  }
  # Fit on test-cal
  tc_T <- optimize(nll_multi_temp, c(0.01, 20),
                   z_mat = zc, y_mat = dc$ymat)$minimum
  if (isTRUE(ENABLE_DIRICHLET)) {
    tc_D <- fit_dirichlet_calibration(dc$pmat, dc$ymat)
  }

  pe_uncal    <- de$pmat
  pe_val_temp <- softmax_t(ze, val_T)
  if (isTRUE(ENABLE_DIRICHLET)) {
    pe_val_dir <- predict_dirichlet_calibration(val_D, de$pmat)
  } else {
    pe_val_dir <- de$pmat
  }
  pe_tc_temp  <- softmax_t(ze, tc_T)
  if (isTRUE(ENABLE_DIRICHLET)) {
    pe_tc_dir <- predict_dirichlet_calibration(tc_D, de$pmat)
  } else {
    pe_tc_dir <- de$pmat
  }

  acc <- mean(apply(de$pmat, 1, which.max) == apply(de$ymat, 1, which.max))
  ece_uncal    <- ece_confidence(pe_uncal, de$ymat)
  ece_val_temp <- ece_confidence(pe_val_temp, de$ymat)
  ece_val_dir <- if (isTRUE(ENABLE_DIRICHLET)) ece_confidence(pe_val_dir, de$ymat) else NA_real_
  ece_tc_temp  <- ece_confidence(pe_tc_temp, de$ymat)
  ece_tc_dir <- if (isTRUE(ENABLE_DIRICHLET)) ece_confidence(pe_tc_dir, de$ymat) else NA_real_
  nll_uncal    <- nll_multi_temp(1, ze, de$ymat)
  nll_val_temp <- nll_multi_temp(val_T, ze, de$ymat)
  nll_tc_temp  <- nll_multi_temp(tc_T, ze, de$ymat)
  nll_val_dir <- if (isTRUE(ENABLE_DIRICHLET)) nll_multi_temp(1, log(pmax(pe_val_dir, 1e-9)), de$ymat) else NA_real_
  nll_tc_dir <- if (isTRUE(ENABLE_DIRICHLET)) nll_multi_temp(1, log(pmax(pe_tc_dir, 1e-9)), de$ymat) else NA_real_

  row <- data.frame(
    target = name, type = paste0("multiclass_", ncol(de$pmat)),
    N_cal = nrow(dc$pmat), N_eval = nrow(de$pmat), Acc = round(acc, 3),
    ECE_uncal = round(ece_uncal, 4),
    ECE_val_temp = round(ece_val_temp, 4),
    ECE_val_platt = NA,
    ECE_val_dirichlet = ifelse(is.na(ece_val_dir), NA, round(ece_val_dir, 4)),
    ECE_tc_temp = round(ece_tc_temp, 4),
    ECE_tc_platt = NA,
    ECE_tc_dirichlet = ifelse(is.na(ece_tc_dir), NA, round(ece_tc_dir, 4)),
    NLL_uncal = round(nll_uncal, 4),
    NLL_val_platt = NA,
    NLL_tc_platt = NA,
    NLL_val_dirichlet = ifelse(is.na(nll_val_dir), NA, round(nll_val_dir, 4)),
    NLL_tc_dirichlet = ifelse(is.na(nll_tc_dir), NA, round(nll_tc_dir, 4)),
    val_T = round(val_T, 3), val_a = NA, val_b = NA,
    tc_T = round(tc_T, 3), tc_a = NA, tc_b = NA,
    stringsAsFactors = FALSE)
  all_rows[[name]] <- row

  if (acc >= ACC_THRESHOLD && nrow(de$pmat) >= N_FILTER) {
    conf_uncal <- apply(pe_uncal, 1, max)
    conf_val   <- apply(pe_val_temp, 1, max)
    conf_val_dir <- apply(pe_val_dir, 1, max)
    conf_tc    <- apply(pe_tc_temp, 1, max)
    conf_tc_dir <- apply(pe_tc_dir, 1, max)
    correct <- (apply(de$pmat, 1, which.max) == apply(de$ymat, 1, which.max)) * 1

    rd <- bind_rows(
      reliability_data(conf_uncal, correct) %>% mutate(method = "Uncalibrated"),
      reliability_data(conf_val, correct) %>%
        mutate(method = paste0("Val-Temp (T=", round(val_T, 2), ")")),
      reliability_data(conf_tc, correct) %>%
        mutate(method = paste0("TestCal-Temp (T=", round(tc_T, 2), ")"))
    )
    if (isTRUE(ENABLE_DIRICHLET)) {
      rd <- bind_rows(
        rd,
        reliability_data(conf_val_dir, correct) %>% mutate(method = "Val-Dirichlet"),
        reliability_data(conf_tc_dir, correct) %>% mutate(method = "TestCal-Dirichlet")
      )
    }
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                        color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, " [", ncol(de$pmat), "-class]  (N_eval=",
                          nrow(de$pmat), ", Acc=", round(acc, 3), ")"),
           subtitle = if (isTRUE(ENABLE_DIRICHLET)) {
             paste0("ECE uncal:", round(ece_uncal, 4),
                    " | val-T:", round(ece_val_temp, 4),
                    " | val-Dir:", round(ece_val_dir, 4),
                    " | tc-T:", round(ece_tc_temp, 4),
                    " | tc-Dir:", round(ece_tc_dir, 4))
           } else {
             paste0("ECE uncal:", round(ece_uncal, 4),
                    " | val-T:", round(ece_val_temp, 4),
                    " | tc-T:", round(ece_tc_temp, 4),
                    " | Dir: disabled")
           },
           x = "Confidence (max prob)", y = "Accuracy") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Pathogenicity (linear -> Platt)
# ══════════════════════════════════════════════════════════════════════════════

patho_heads <- list(
  pathogenicity_human  = list(
    pred_val = "pathogenicity_human_val", true_col = "pathogenicity_human_true"),
  pathogenicity_animal = list(
    pred_val = "pathogenicity_animal_val", true_col = "pathogenicity_animal_true"),
  pathogenicity_plant  = list(
    pred_val = "pathogenicity_plant_val", true_col = "pathogenicity_plant_true")
)

for (name in names(patho_heads)) {
  h <- patho_heads[[name]]

  extract_patho <- function(d) {
    raw <- as.numeric(d[[h$pred_val]])
    y_raw <- as.numeric(d[[h$true_col]])
    ok <- !is.na(y_raw) & y_raw != -9999
    list(raw = raw[ok], y = as.numeric(y_raw[ok] > 0))
  }

  dv <- extract_patho(val_df)
  dc <- extract_patho(test_cal)
  de <- extract_patho(test_eval)

  # Fit on val
  val_P <- optim(c(1, 0), nll_binary_platt, logits = dv$raw, y = dv$y,
                 method = "BFGS")$par
  # Fit on test-cal
  tc_P <- optim(c(1, 0), nll_binary_platt, logits = dc$raw, y = dc$y,
                method = "BFGS")$par

  p_e_val  <- sigmoid(val_P[1] * de$raw + val_P[2])
  p_e_tc   <- sigmoid(tc_P[1] * de$raw + tc_P[2])
  acc_val  <- mean((p_e_val > 0.5) == de$y)
  acc_tc   <- mean((p_e_tc > 0.5) == de$y)

  row <- data.frame(
    target = name, type = "linear_platt",
    N_cal = length(dc$y), N_eval = length(de$y), Acc = round(acc_tc, 3),
    ECE_uncal = NA,
    ECE_val_temp = NA,
    ECE_val_platt = round(ece(p_e_val, de$y), 4),
    ECE_val_dirichlet = NA,
    ECE_tc_temp = NA,
    ECE_tc_platt = round(ece(p_e_tc, de$y), 4),
    ECE_tc_dirichlet = NA,
    NLL_uncal = NA,
    NLL_val_platt = round(nll_binary_platt(val_P, de$raw, de$y), 4),
    NLL_tc_platt  = round(nll_binary_platt(tc_P, de$raw, de$y), 4),
    NLL_val_dirichlet = NA,
    NLL_tc_dirichlet = NA,
    val_T = NA, val_a = round(val_P[1], 3), val_b = round(val_P[2], 3),
    tc_T = NA, tc_a = round(tc_P[1], 3), tc_b = round(tc_P[2], 3),
    stringsAsFactors = FALSE)
  all_rows[[name]] <- row

  if (acc_tc >= ACC_THRESHOLD && length(de$y) >= N_FILTER) {
    rd <- bind_rows(
      reliability_data(p_e_val, de$y) %>%
        mutate(method = paste0("Val-Platt (a=", round(val_P[1], 2),
                               " b=", round(val_P[2], 2), ")")),
      reliability_data(p_e_tc, de$y) %>%
        mutate(method = paste0("TestCal-Platt (a=", round(tc_P[1], 2),
                               " b=", round(tc_P[2], 2), ")"))
    )
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                        color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, " [linear->Platt]  (N_eval=", length(de$y),
                          ", Acc=", round(acc_tc, 3), ")"),
           subtitle = paste0("Val: a=", round(val_P[1], 2), " b=", round(val_P[2], 2),
                             " ECE=", round(ece(p_e_val, de$y), 4),
                             " | TC: a=", round(tc_P[1], 2), " b=", round(tc_P[2], 2),
                             " ECE=", round(ece(p_e_tc, de$y), 4)),
           x = "Platt Predicted Prob", y = "Observed Frequency") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# BUILD SUMMARY TABLE
# ══════════════════════════════════════════════════════════════════════════════

summary_df <- do.call(rbind, all_rows)
summary_df <- summary_df %>%
  arrange(desc(Acc)) %>%
  mutate(quality = ifelse(Acc >= ACC_THRESHOLD, "PASS", "FAIL"))

cat("\n")
print(summary_df, row.names = FALSE)

# ── Compact comparison table ─────────────────────────────────────────────────
cat("\n\n=== Val-based vs Within-test Calibration (ECE on test-eval) ===\n")
compact <- summary_df %>%
  transmute(
    target, type, N_eval, Acc,
    ECE_uncal,
    ECE_val = coalesce(ECE_val_platt, ECE_val_dirichlet, ECE_val_temp),
    ECE_tc  = coalesce(ECE_tc_platt, ECE_tc_dirichlet, ECE_tc_temp),
    val_params = ifelse(is.na(val_T),
                        paste0("a=", val_a, " b=", val_b),
                        paste0("T=", val_T)),
    tc_params  = ifelse(is.na(tc_T),
                        paste0("a=", tc_a, " b=", tc_b),
                        paste0("T=", tc_T))
  )
print(compact, row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# PARAMETER COMPARISON PLOTS
# ══════════════════════════════════════════════════════════════════════════════

# Filter to N_eval >= N_FILTER
sf <- summary_df %>% filter(N_eval >= N_FILTER)

# --- Temperature plot: val_T vs tc_T for heads that use temperature scaling ---
temp_df <- sf %>%
  filter(!is.na(val_T)) %>%
  select(target, N_eval, Acc, val_T, tc_T)

p_temp <- NULL
if (nrow(temp_df) > 0) {
  p_temp <- ggplot(temp_df, aes(x = val_T, y = tc_T)) +
    geom_point(aes(size = N_eval, color = Acc), alpha = 0.8) +
    geom_text(aes(label = target), vjust = -1.2, size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "blue", alpha = 0.5) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "blue", alpha = 0.5) +
    scale_color_gradient(low = "orange", high = "darkgreen") +
    labs(title = "Temperature parameter: Val-based vs Within-test",
         subtitle = "T<1 = sharpen (model underconfident), T>1 = soften (model overconfident)",
         x = "Val-fitted T", y = "Test-cal-fitted T") +
    theme_minimal() +
    coord_fixed(xlim = c(min(c(temp_df$val_T, temp_df$tc_T)) - 0.1,
                         max(c(temp_df$val_T, temp_df$tc_T)) + 0.1),
                ylim = c(min(c(temp_df$val_T, temp_df$tc_T)) - 0.1,
                         max(c(temp_df$val_T, temp_df$tc_T)) + 0.1))
}

# --- Platt a,b plot: val vs tc for heads that use Platt scaling ---
platt_df <- sf %>%
  filter(!is.na(val_a)) %>%
  select(target, N_eval, Acc, val_a, val_b, tc_a, tc_b)

p_platt_a <- NULL
p_platt_b <- NULL
if (nrow(platt_df) > 0) {
  p_platt_a <- ggplot(platt_df, aes(x = val_a, y = tc_a)) +
    geom_point(aes(size = N_eval, color = Acc), alpha = 0.8) +
    geom_text(aes(label = target), vjust = -1.2, size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    scale_color_gradient(low = "orange", high = "darkgreen") +
    labs(title = "Platt slope (a): Val-based vs Within-test",
         subtitle = "a=1 means identity scaling; sign flip = different population",
         x = "Val-fitted a", y = "Test-cal-fitted a") +
    theme_minimal()

  p_platt_b <- ggplot(platt_df, aes(x = val_b, y = tc_b)) +
    geom_point(aes(size = N_eval, color = Acc), alpha = 0.8) +
    geom_text(aes(label = target), vjust = -1.2, size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    scale_color_gradient(low = "orange", high = "darkgreen") +
    labs(title = "Platt intercept (b): Val-based vs Within-test",
         subtitle = "b=0 means no bias shift; large difference = population mismatch",
         x = "Val-fitted b", y = "Test-cal-fitted b") +
    theme_minimal()
}

# --- ECE comparison: uncal vs val vs tc ---
ece_df <- sf %>%
  transmute(
    target, N_eval,
    ECE_uncal,
    ECE_val = coalesce(ECE_val_platt, ECE_val_dirichlet, ECE_val_temp),
    ECE_tc  = coalesce(ECE_tc_platt, ECE_tc_dirichlet, ECE_tc_temp)
  ) %>%
  tidyr::pivot_longer(cols = starts_with("ECE_"),
                      names_to = "method", values_to = "ECE") %>%
  filter(!is.na(ECE)) %>%
  mutate(method = factor(method,
                         levels = c("ECE_uncal", "ECE_val", "ECE_tc"),
                         labels = c("Uncalibrated", "Val-based", "Test-cal (70%)")))

p_ece <- ggplot(ece_df, aes(x = reorder(target, -ECE), y = ECE, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(aes(label = round(ECE, 3)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("Uncalibrated" = "#E8808080",
                                "Val-based" = "#66C2A5",
                                "Test-cal (70%)" = "#8DA0CB")) +
  labs(title = paste0("ECE comparison on test-eval (N_eval >= ", N_FILTER, ")"),
       subtitle = "Lower is better. Val-based calibration can hurt when populations differ.",
       x = NULL, y = "Expected Calibration Error") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom")

# ══════════════════════════════════════════════════════════════════════════════
# SAVE PDF
# ══════════════════════════════════════════════════════════════════════════════

pdf_path <- file.path(cfg$output_dir, "calibration_val_test.pdf")
pdf(pdf_path, width = 14, height = 9, onefile = TRUE)

# Page 1: compact comparison table
grid::grid.newpage()
grid::grid.text(
  "Calibration: val-based vs within-test (70/30) on test-eval",
  x = 0.5, y = 0.96,
  gp = grid::gpar(fontsize = 14, fontface = "bold"))
grid::grid.text(
  paste0("Model: bacdive_attn_lstm_sum_28 Ep.23 | ",
         "Source: bd_pred_new_full.csv | ",
         "Test-cal: ", nrow(test_cal), " | Test-eval: ", nrow(test_eval),
         " | Seed: ", SEED, " | N_filter: ", N_FILTER),
  x = 0.5, y = 0.92, gp = grid::gpar(fontsize = 9, col = "grey40"))
grid::grid.text(
  paste0("NOTE: Val has 90.6% genus overlap with train (data leakage). ",
         "Test is a genus-level holdout (75% unseen genera)."),
  x = 0.5, y = 0.88, gp = grid::gpar(fontsize = 9, col = "red3"))

theme_tbl <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize = 7, hjust = 0, x = 0.02),
    bg_params = list(fill = ifelse(compact$Acc < ACC_THRESHOLD,
                                   "#FFE0E0", "white"))),
  colhead = list(
    fg_params = list(fontsize = 7, fontface = "bold", hjust = 0, x = 0.02),
    bg_params = list(fill = "#E8E8E8"))
)

tbl_display <- compact %>%
  mutate(across(where(is.numeric), ~ifelse(is.na(.), "-", as.character(.))))
tg <- tableGrob(tbl_display, rows = NULL, theme = theme_tbl)
tg$widths <- unit(rep(1 / ncol(tbl_display), ncol(tbl_display)), "npc")
grid::pushViewport(grid::viewport(y = 0.42, height = 0.75))
grid::grid.draw(tg)
grid::popViewport()

# Page 2: ECE bar comparison
print(p_ece)

# Page 3: Temperature parameter comparison
if (!is.null(p_temp)) print(p_temp)

# Page 4-5: Platt a,b parameter comparison
if (!is.null(p_platt_a)) print(p_platt_a)
if (!is.null(p_platt_b)) print(p_platt_b)

# Reliability diagrams (N_eval >= N_FILTER only)
for (p in plot_list) print(p)

dev.off()
cat("\nReport saved to:", pdf_path, "\n")
cat("Passing heads:", paste(names(plot_list), collapse = ", "), "\n")

# ── Save results ─────────────────────────────────────────────────────────────
all_results <- list(
  summary = summary_df,
  compact = compact,
  n_val = nrow(val_df),
  n_test_cal = nrow(test_cal),
  n_test_eval = nrow(test_eval),
  cal_frac = CAL_FRAC,
  seed = SEED,
  acc_threshold = ACC_THRESHOLD,
  source = "bd_pred_new_full.csv"
)
rds_path <- file.path(cfg$output_dir, "calibration_val_test.rds")
saveRDS(all_results, rds_path)
cat("Results saved to:", rds_path, "\n")
