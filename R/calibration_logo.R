# ── Leave-One-Family-Out (LOFO) Calibration ──────────────────────────────────
# Pool val + test, fit calibration via LOFO CV.
# For each family, hold it out, fit Platt/temperature on the rest,
# predict on the held-out family. Aggregate held-out predictions for ECE.
#
# Prereq: run build_genome_family_map.R to create genome_family_map.csv
# Usage: conda run -n genome Rscript R/calibration_logo.R

library(ggplot2)
library(dplyr)
library(gridExtra)
library(data.table)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg)) else normalizePath(getwd())
script_dir <- if (dir.exists(script_path)) script_path else dirname(script_path)
proj_dir <- if (basename(script_dir) == "R") dirname(script_dir) else script_dir
source(file.path(proj_dir, "R", "config_paths.R"))
source(file.path(proj_dir, "R", "calibration_common.R"))
cfg <- load_bacdive_config(
  proj_dir,
  required = c("bd_pred_csv", "bd_labels_csv", "bd_splits_csv", "genome_family_map_csv")
)

ACC_THRESHOLD <- 0.6
N_FILTER <- 100
ENABLE_DIRICHLET <- tolower(Sys.getenv("ENABLE_DIRICHLET_MULTICLASS", "false")) %in%
  c("1", "true", "yes")

# ── Load predictions and ground truth (fread for speed) ──────────────────────
cat("Loading data...\n"); t0 <- proc.time()
bdpred <- fread(cfg$bd_pred_csv, data.table = FALSE)
labels_true <- fread(cfg$bd_labels_csv, data.table = FALSE)
splits <- fread(cfg$bd_splits_csv, data.table = FALSE)
cat("CSVs loaded in", round((proc.time() - t0)[3], 1), "s\n")

bdpred$genome <- basename(bdpred$file)
labels_true$genome <- labels_true$file
splits$genome <- splits$file

df <- merge(bdpred, labels_true, by = "genome", suffixes = c("_pred", "_true"))
df <- merge(df, splits[, c("genome", "type")], by = "genome")

pool_df <- df %>% filter(type %in% c("validation", "test"))
cat("Pooled val+test genomes:", nrow(pool_df), "\n")

# ── Load precomputed genome → family mapping ────────────────────────────────
fam_map <- fread(cfg$genome_family_map_csv, data.table = FALSE)
pool_df <- merge(pool_df, fam_map[, c("file", "genus", "Family", "Order", "Phylum")],
                 by.x = "genome", by.y = "file")
n_na <- sum(is.na(pool_df$Phylum))
cat("Genomes with Phylum:", nrow(pool_df) - n_na, "\n")
cat("Genomes missing Phylum:", n_na, "(dropping)\n")
pool_df <- pool_df[!is.na(pool_df$Phylum), ]

group_counts <- table(pool_df$Phylum)
cat("Unique phyla:", length(group_counts), "\n")
cat("Phyla with >=5 samples:", sum(group_counts >= 5), "\n")
cat("Phyla with 1 sample:", sum(group_counts == 1), "\n\n")

# Safe optimizer wrappers (catch errors from degenerate folds)
safe_opt_temp <- function(logits, y) {
  tryCatch({
    optimize(nll_binary_temp, c(0.01, 20), logits = logits, y = y)$minimum
  }, error = function(e) 1.0)
}
safe_opt_platt <- function(logits, y) {
  tryCatch({
    optim(c(1, 0), nll_binary_platt, logits = logits, y = y,
          method = "BFGS")$par
  }, error = function(e) c(1, 0))
}
safe_opt_multi_temp <- function(z_mat, y_mat) {
  tryCatch({
    optimize(nll_multi_temp, c(0.01, 20), z_mat = z_mat, y_mat = y_mat)$minimum
  }, error = function(e) 1.0)
}

# ── Collect results ──────────────────────────────────────────────────────────
all_rows <- list()
plot_list <- list()

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Binary sigmoid heads — LOFO CV
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
  cat("LOFO CV:", name, "...\n")

  p_all <- as.numeric(pool_df[[h$pred_val]])
  y_all <- as.numeric(pool_df[[h$true_col]])
  fam_all <- pool_df$Phylum
  ok <- !is.na(y_all) & y_all != -9999
  p_all <- p_all[ok]; y_all <- y_all[ok]; fam_all <- fam_all[ok]
  z_all <- logit(p_all)

  families <- unique(fam_all)
  # LOFO: for each family, hold it out
  p_lofo_uncal <- numeric(length(y_all))
  p_lofo_platt <- numeric(length(y_all))
  p_lofo_temp  <- numeric(length(y_all))
  lofo_a <- numeric(length(families)); lofo_b <- numeric(length(families))
  lofo_T <- numeric(length(families))

  for (fi in seq_along(families)) {
    fam <- families[fi]
    held <- fam_all == fam
    train_idx <- !held

    if (length(unique(y_all[train_idx])) < 2) {
      p_lofo_uncal[held] <- p_all[held]
      p_lofo_platt[held] <- p_all[held]
      p_lofo_temp[held]  <- p_all[held]
      lofo_a[fi] <- 1; lofo_b[fi] <- 0; lofo_T[fi] <- 1
      next
    }

    T_fit <- safe_opt_temp(z_all[train_idx], y_all[train_idx])
    P_fit <- safe_opt_platt(z_all[train_idx], y_all[train_idx])

    p_lofo_uncal[held] <- p_all[held]
    p_lofo_temp[held]  <- sigmoid(z_all[held] / T_fit)
    p_lofo_platt[held] <- sigmoid(P_fit[1] * z_all[held] + P_fit[2])
    lofo_a[fi] <- P_fit[1]; lofo_b[fi] <- P_fit[2]; lofo_T[fi] <- T_fit
  }

  acc <- mean((p_lofo_uncal > 0.5) == y_all)
  ece_uncal <- ece(p_lofo_uncal, y_all)
  ece_temp  <- ece(p_lofo_temp, y_all)
  ece_platt <- ece(p_lofo_platt, y_all)

  global_T <- safe_opt_temp(z_all, y_all)
  global_P <- safe_opt_platt(z_all, y_all)

  row <- data.frame(
    target = name, type = "binary_sigmoid",
    N_total = length(y_all), N_families = length(families),
    Acc = round(acc, 3),
    ECE_uncal = round(ece_uncal, 4),
    ECE_lofo_temp = round(ece_temp, 4),
    ECE_lofo_platt = round(ece_platt, 4),
    ECE_lofo_dirichlet = NA,
    global_T = round(global_T, 3),
    global_a = round(global_P[1], 3),
    global_b = round(global_P[2], 3),
    median_lofo_T = round(median(lofo_T), 3),
    median_lofo_a = round(median(lofo_a), 3),
    median_lofo_b = round(median(lofo_b), 3),
    stringsAsFactors = FALSE)
  all_rows[[name]] <- row

  if (acc >= ACC_THRESHOLD && length(y_all) >= N_FILTER) {
    rd <- bind_rows(
      reliability_data(p_lofo_uncal, y_all) %>% mutate(method = "Uncalibrated"),
      reliability_data(p_lofo_temp, y_all) %>%
        mutate(method = paste0("LOFO-Temp (med T=", round(median(lofo_T), 2), ")")),
      reliability_data(p_lofo_platt, y_all) %>%
        mutate(method = paste0("LOFO-Platt (med a=", round(median(lofo_a), 2),
                               " b=", round(median(lofo_b), 2), ")"))
    )
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                        color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, " [LOFO CV]  (N=", length(y_all),
                          ", ", length(families), " families, Acc=", round(acc, 3), ")"),
           subtitle = paste0("ECE uncal:", round(ece_uncal, 4),
                             " | LOFO-Temp:", round(ece_temp, 4),
                             " | LOFO-Platt:", round(ece_platt, 4)),
           x = "Mean Predicted Prob", y = "Observed Frequency") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: 2-class softmax heads — LOFO CV
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
  cat("LOFO CV:", name, "...\n")

  extract_s2 <- function(d) {
    pmat <- parse_prob_vec(d[[h$pred_val]])
    y1 <- as.numeric(d[[h$true_cols[2]]])
    ok <- !is.na(y1) & y1 != -9999 &
      !is.na(as.numeric(d[[h$true_cols[1]]])) &
      (as.numeric(d[[h$true_cols[1]]]) + y1) > 0
    list(p = pmat[ok, 2], y = y1[ok], family = d$Phylum[ok])
  }

  dd <- extract_s2(pool_df)
  z_all <- logit(dd$p)

  families <- unique(dd$family)
  p_lofo_uncal <- numeric(length(dd$y))
  p_lofo_platt <- numeric(length(dd$y))
  p_lofo_temp  <- numeric(length(dd$y))
  lofo_a <- numeric(length(families)); lofo_b <- numeric(length(families))
  lofo_T <- numeric(length(families))

  for (fi in seq_along(families)) {
    fam <- families[fi]
    held <- dd$family == fam
    train_idx <- !held

    if (length(unique(dd$y[train_idx])) < 2) {
      p_lofo_uncal[held] <- dd$p[held]
      p_lofo_platt[held] <- dd$p[held]
      p_lofo_temp[held]  <- dd$p[held]
      lofo_a[fi] <- 1; lofo_b[fi] <- 0; lofo_T[fi] <- 1
      next
    }

    T_fit <- safe_opt_temp(z_all[train_idx], dd$y[train_idx])
    P_fit <- safe_opt_platt(z_all[train_idx], dd$y[train_idx])

    p_lofo_uncal[held] <- dd$p[held]
    p_lofo_temp[held]  <- sigmoid(z_all[held] / T_fit)
    p_lofo_platt[held] <- sigmoid(P_fit[1] * z_all[held] + P_fit[2])
    lofo_a[fi] <- P_fit[1]; lofo_b[fi] <- P_fit[2]; lofo_T[fi] <- T_fit
  }

  acc <- mean((p_lofo_uncal > 0.5) == dd$y)
  ece_uncal <- ece(p_lofo_uncal, dd$y)
  ece_temp  <- ece(p_lofo_temp, dd$y)
  ece_platt <- ece(p_lofo_platt, dd$y)

  global_T <- safe_opt_temp(z_all, dd$y)
  global_P <- safe_opt_platt(z_all, dd$y)

  row <- data.frame(
    target = name, type = "softmax_2class",
    N_total = length(dd$y), N_families = length(families),
    Acc = round(acc, 3),
    ECE_uncal = round(ece_uncal, 4),
    ECE_lofo_temp = round(ece_temp, 4),
    ECE_lofo_platt = round(ece_platt, 4),
    ECE_lofo_dirichlet = NA,
    global_T = round(global_T, 3),
    global_a = round(global_P[1], 3),
    global_b = round(global_P[2], 3),
    median_lofo_T = round(median(lofo_T), 3),
    median_lofo_a = round(median(lofo_a), 3),
    median_lofo_b = round(median(lofo_b), 3),
    stringsAsFactors = FALSE)
  all_rows[[name]] <- row

  if (acc >= ACC_THRESHOLD && length(dd$y) >= N_FILTER) {
    rd <- bind_rows(
      reliability_data(p_lofo_uncal, dd$y) %>% mutate(method = "Uncalibrated"),
      reliability_data(p_lofo_temp, dd$y) %>%
        mutate(method = paste0("LOFO-Temp (med T=", round(median(lofo_T), 2), ")")),
      reliability_data(p_lofo_platt, dd$y) %>%
        mutate(method = paste0("LOFO-Platt (med a=", round(median(lofo_a), 2),
                               " b=", round(median(lofo_b), 2), ")"))
    )
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                        color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, " [LOFO CV]  (N=", length(dd$y),
                          ", ", length(families), " families, Acc=", round(acc, 3), ")"),
           subtitle = paste0("ECE uncal:", round(ece_uncal, 4),
                             " | LOFO-Temp:", round(ece_temp, 4),
                             " | LOFO-Platt:", round(ece_platt, 4)),
           x = "Mean Predicted Prob", y = "Observed Frequency") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Multiclass softmax heads (gram, cellshape, flagellum) — LOFO CV
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
  cat("LOFO CV:", name, "...\n")

  extract_mc <- function(d) {
    pmat <- parse_prob_vec(d[[h$pred_val]])
    ymat <- as.matrix(d[, h$true_cols])
    storage.mode(ymat) <- "numeric"
    ok <- apply(ymat, 1, function(r) !any(is.na(r) | r == -9999)) &
      rowSums(ymat) > 0
    list(pmat = pmat[ok, ], ymat = ymat[ok, ], family = d$Phylum[ok])
  }

  dd <- extract_mc(pool_df)
  eps <- 1e-7
  z_all <- log(pmax(dd$pmat, eps))

  families <- unique(dd$family)
  K <- ncol(dd$pmat)

  pmat_lofo_uncal <- dd$pmat
  pmat_lofo_temp  <- matrix(0, nrow(dd$pmat), K)
  pmat_lofo_dir   <- matrix(0, nrow(dd$pmat), K)
  lofo_T <- numeric(length(families))

  for (fi in seq_along(families)) {
    fam <- families[fi]
    held <- dd$family == fam
    train_idx <- !held

    if (length(unique(apply(dd$ymat[train_idx, , drop = FALSE], 1, which.max))) < 2) {
      pmat_lofo_temp[held, ] <- dd$pmat[held, ]
      pmat_lofo_dir[held, ] <- dd$pmat[held, ]
      lofo_T[fi] <- 1
      next
    }

    T_fit <- safe_opt_multi_temp(z_all[train_idx, , drop = FALSE],
                                  dd$ymat[train_idx, , drop = FALSE])
    if (isTRUE(ENABLE_DIRICHLET)) {
      D_fit <- fit_dirichlet_calibration(
        dd$pmat[train_idx, , drop = FALSE],
        dd$ymat[train_idx, , drop = FALSE]
      )
    }
    pmat_lofo_temp[held, ] <- softmax_t(z_all[held, , drop = FALSE], T_fit)
    if (isTRUE(ENABLE_DIRICHLET)) {
      pmat_lofo_dir[held, ] <- predict_dirichlet_calibration(D_fit, dd$pmat[held, , drop = FALSE])
    } else {
      pmat_lofo_dir[held, ] <- dd$pmat[held, ]
    }
    lofo_T[fi] <- T_fit
  }

  acc <- mean(apply(dd$pmat, 1, which.max) == apply(dd$ymat, 1, which.max))
  ece_uncal <- ece_confidence(pmat_lofo_uncal, dd$ymat)
  ece_temp  <- ece_confidence(pmat_lofo_temp, dd$ymat)
  ece_dir <- if (isTRUE(ENABLE_DIRICHLET)) ece_confidence(pmat_lofo_dir, dd$ymat) else NA_real_

  global_T <- safe_opt_multi_temp(z_all, dd$ymat)

  row <- data.frame(
    target = name, type = paste0("multiclass_", K),
    N_total = nrow(dd$pmat), N_families = length(families),
    Acc = round(acc, 3),
    ECE_uncal = round(ece_uncal, 4),
    ECE_lofo_temp = round(ece_temp, 4),
    ECE_lofo_platt = NA,
    ECE_lofo_dirichlet = ifelse(is.na(ece_dir), NA, round(ece_dir, 4)),
    global_T = round(global_T, 3),
    global_a = NA, global_b = NA,
    median_lofo_T = round(median(lofo_T), 3),
    median_lofo_a = NA, median_lofo_b = NA,
    stringsAsFactors = FALSE)
  all_rows[[name]] <- row

  if (acc >= ACC_THRESHOLD && nrow(dd$pmat) >= N_FILTER) {
    conf_uncal <- apply(pmat_lofo_uncal, 1, max)
    conf_temp  <- apply(pmat_lofo_temp, 1, max)
    conf_dir   <- apply(pmat_lofo_dir, 1, max)
    correct <- (apply(dd$pmat, 1, which.max) == apply(dd$ymat, 1, which.max)) * 1

    rd <- bind_rows(
      reliability_data(conf_uncal, correct) %>% mutate(method = "Uncalibrated"),
      reliability_data(conf_temp, correct) %>%
        mutate(method = paste0("LOFO-Temp (med T=", round(median(lofo_T), 2), ")"))
    )
    if (isTRUE(ENABLE_DIRICHLET)) {
      rd <- bind_rows(
        rd,
        reliability_data(conf_dir, correct) %>% mutate(method = "LOFO-Dirichlet")
      )
    }
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                        color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, " [", K, "-class, LOFO CV]  (N=",
                          nrow(dd$pmat), ", ", length(families),
                          " families, Acc=", round(acc, 3), ")"),
           subtitle = if (isTRUE(ENABLE_DIRICHLET)) {
             paste0("ECE uncal:", round(ece_uncal, 4),
                    " | LOFO-Temp:", round(ece_temp, 4),
                    " | LOFO-Dir:", round(ece_dir, 4))
           } else {
             paste0("ECE uncal:", round(ece_uncal, 4),
                    " | LOFO-Temp:", round(ece_temp, 4),
                    " | LOFO-Dir: disabled")
           },
           x = "Confidence (max prob)", y = "Accuracy") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Pathogenicity (linear -> Platt) — LOFO CV
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
  cat("LOFO CV:", name, "...\n")

  raw_all <- as.numeric(pool_df[[h$pred_val]])
  y_raw   <- as.numeric(pool_df[[h$true_col]])
  fam_all <- pool_df$Phylum
  ok <- !is.na(y_raw) & y_raw != -9999
  raw_all <- raw_all[ok]; y_all <- as.numeric(y_raw[ok] > 0); fam_all <- fam_all[ok]

  families <- unique(fam_all)
  p_lofo_platt <- numeric(length(y_all))
  lofo_a <- numeric(length(families)); lofo_b <- numeric(length(families))

  for (fi in seq_along(families)) {
    fam <- families[fi]
    held <- fam_all == fam
    train_idx <- !held

    if (length(unique(y_all[train_idx])) < 2) {
      p_lofo_platt[held] <- sigmoid(raw_all[held])
      lofo_a[fi] <- 1; lofo_b[fi] <- 0
      next
    }

    P_fit <- safe_opt_platt(raw_all[train_idx], y_all[train_idx])
    p_lofo_platt[held] <- sigmoid(P_fit[1] * raw_all[held] + P_fit[2])
    lofo_a[fi] <- P_fit[1]; lofo_b[fi] <- P_fit[2]
  }

  acc <- mean((p_lofo_platt > 0.5) == y_all)
  ece_platt <- ece(p_lofo_platt, y_all)

  global_P <- safe_opt_platt(raw_all, y_all)

  row <- data.frame(
    target = name, type = "linear_platt",
    N_total = length(y_all), N_families = length(families),
    Acc = round(acc, 3),
    ECE_uncal = NA,
    ECE_lofo_temp = NA,
    ECE_lofo_platt = round(ece_platt, 4),
    ECE_lofo_dirichlet = NA,
    global_T = NA,
    global_a = round(global_P[1], 3),
    global_b = round(global_P[2], 3),
    median_lofo_T = NA,
    median_lofo_a = round(median(lofo_a), 3),
    median_lofo_b = round(median(lofo_b), 3),
    stringsAsFactors = FALSE)
  all_rows[[name]] <- row

  if (acc >= ACC_THRESHOLD && length(y_all) >= N_FILTER) {
    rd <- reliability_data(p_lofo_platt, y_all) %>%
      mutate(method = paste0("LOFO-Platt (med a=", round(median(lofo_a), 2),
                             " b=", round(median(lofo_b), 2), ")"))
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                        color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, " [LOFO CV, linear->Platt]  (N=", length(y_all),
                          ", ", length(families), " families, Acc=", round(acc, 3), ")"),
           subtitle = paste0("LOFO-Platt ECE:", round(ece_platt, 4)),
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
cat("\n\n=== LOFO CV Calibration Results ===\n")
compact <- summary_df %>%
  transmute(
    target, type, N_total, N_families, Acc,
    ECE_uncal,
    ECE_lofo = coalesce(ECE_lofo_platt, ECE_lofo_dirichlet, ECE_lofo_temp),
    global_params = case_when(
      !is.na(global_T) & !is.na(global_a) ~
        paste0("T=", global_T, " a=", global_a, " b=", global_b),
      !is.na(global_T) ~ paste0("T=", global_T),
      TRUE ~ paste0("a=", global_a, " b=", global_b)
    ),
    median_lofo_params = case_when(
      !is.na(median_lofo_T) & !is.na(median_lofo_a) ~
        paste0("T=", median_lofo_T, " a=", median_lofo_a, " b=", median_lofo_b),
      !is.na(median_lofo_T) ~ paste0("T=", median_lofo_T),
      TRUE ~ paste0("a=", median_lofo_a, " b=", median_lofo_b)
    )
  )
print(compact, row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# PARAMETER DISTRIBUTION PLOTS
# ══════════════════════════════════════════════════════════════════════════════

# Filter to passing heads
sf <- summary_df %>% filter(N_total >= N_FILTER)

# --- ECE comparison: uncal vs LOGO ---
ece_df <- sf %>%
  transmute(
    target, N_total,
    ECE_uncal,
    ECE_lofo = coalesce(ECE_lofo_platt, ECE_lofo_dirichlet, ECE_lofo_temp)
  ) %>%
  tidyr::pivot_longer(cols = starts_with("ECE_"),
                      names_to = "method", values_to = "ECE") %>%
  filter(!is.na(ECE)) %>%
  mutate(method = factor(method,
                         levels = c("ECE_uncal", "ECE_lofo"),
                         labels = c("Uncalibrated", "LOFO CV")))

p_ece <- ggplot(ece_df, aes(x = reorder(target, -ECE), y = ECE, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(aes(label = round(ECE, 3)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("Uncalibrated" = "#E8808080",
                                "LOFO CV" = "#FC8D62")) +
  labs(title = "ECE comparison: Uncalibrated vs LOFO CV",
       subtitle = "Lower is better. LOFO CV = leave-one-family-out cross-validated calibration.",
       x = NULL, y = "Expected Calibration Error") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom")

# ══════════════════════════════════════════════════════════════════════════════
# SAVE PDF
# ══════════════════════════════════════════════════════════════════════════════

pdf_path <- file.path(cfg$output_dir, "calibration_logo.pdf")
pdf(pdf_path, width = 14, height = 9, onefile = TRUE)

# Page 1: compact comparison table
grid::grid.newpage()
grid::grid.text(
  "Calibration: Leave-One-Family-Out Cross-Validation (LOFO CV)",
  x = 0.5, y = 0.96,
  gp = grid::gpar(fontsize = 14, fontface = "bold"))
grid::grid.text(
  paste0("Model: bacdive_attn_lstm_sum_28 Ep.23 | ",
         "Source: bd_pred_new_full.csv | ",
         "Pool: val+test = ", nrow(pool_df), " genomes | ",
         "Phyla: ", length(unique(pool_df$Phylum))),
  x = 0.5, y = 0.92, gp = grid::gpar(fontsize = 9, col = "grey40"))
grid::grid.text(
  paste0("For each family: hold out, fit calibration on rest, predict held-out. ",
         "Aggregated held-out predictions -> ECE. Genus from FASTA headers, Family from microbe.cards S1."),
  x = 0.5, y = 0.88, gp = grid::gpar(fontsize = 9, col = "blue3"))

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

# Reliability diagrams
for (p in plot_list) print(p)

dev.off()
cat("\nReport saved to:", pdf_path, "\n")
cat("Passing heads:", paste(names(plot_list), collapse = ", "), "\n")

# ── Save results ─────────────────────────────────────────────────────────────
all_results <- list(
  summary = summary_df,
  compact = compact,
  passing_heads = names(plot_list),
  n_pool = nrow(pool_df),
  n_families = length(unique(pool_df$Phylum)),  # now Order-level
  acc_threshold = ACC_THRESHOLD,
  source = "bd_pred_new_full.csv",
  genome_taxonomy = fam_map
)
rds_path <- file.path(cfg$output_dir, "calibration_logo.rds")
saveRDS(all_results, rds_path)
cat("Results saved to:", rds_path, "\n")
