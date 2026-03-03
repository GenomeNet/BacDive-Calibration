# ── IPW-Platt Scaling: Selection-Bias-Corrected Calibration ─────────────────
# Adds inverse probability weighted Platt scaling to the BacDive genome
# phenotype prediction calibration pipeline. Compares plain vs IPW-weighted
# Platt parameters and runs diagnostic checks.
#
# Prereqs:
#   - local selection artifacts: data/selection/selection_probs.rds +
#     data/selection/bugphyzz_with_proxies.rds
#   - bacdive: genome_family_map.csv, bd_pred_new_full.csv, labels, splits
#
# Usage: Rscript R/calibration_ipw_platt.R

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
  required = c("bd_pred_csv", "bd_labels_csv", "bd_splits_csv", "genome_family_map_csv"),
  require_selection = TRUE
)

cat("=== IPW-Platt Scaling for BacDive Calibration ===\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Load BacDive data (same as calibration_logo.R)
# ══════════════════════════════════════════════════════════════════════════════

ACC_THRESHOLD <- 0.6
N_FILTER <- 100

cat("Loading BacDive data...\n"); t0 <- proc.time()
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

# Load genome -> family/phylum mapping
fam_map <- fread(cfg$genome_family_map_csv, data.table = FALSE)
pool_df <- merge(pool_df, fam_map[, c("file", "genus", "Family", "Order", "Phylum")],
                 by.x = "genome", by.y = "file")
pool_df <- pool_df[!is.na(pool_df$Phylum), ]
cat("Genomes with Phylum:", nrow(pool_df), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Load phenomnar selection probabilities & build genus-level medians
# ══════════════════════════════════════════════════════════════════════════════

cat("\nLoading phenomnar selection probabilities...\n")
sel_probs <- readRDS(cfg$selection_probs_rds)
bp <- readRDS(cfg$bugphyzz_proxies_rds)

# sel_probs is aligned with bugphyzz rows where !is.na(year_valid_pub)
idx_valid <- !is.na(bp$year_valid_pub)
bp_valid <- bp[idx_valid, ]
stopifnot(nrow(bp_valid) == nrow(sel_probs))

# Extract genus name from binomial species name
bp_valid$genus_name <- sub(" .*", "", bp_valid$Species)

# Also get proxy columns for covariate balance diagnostics
proxy_cols_for_smd <- c("log_pub_count", "log_assembly_count", "log_wiki_size")
bp_valid_proxies <- bp_valid[, c("genus_name", proxy_cols_for_smd), drop = FALSE]

# Build genus-level median selection probabilities
sel_prob_cols <- colnames(sel_probs)
genus_sel <- cbind(genus_name = bp_valid$genus_name, as.data.frame(sel_probs))
genus_median <- genus_sel %>%
  group_by(genus_name) %>%
  summarise(across(all_of(sel_prob_cols), \(x) median(x, na.rm = TRUE)), .groups = "drop")
cat("Genus-level sel_prob table:", nrow(genus_median), "genera\n")

# Also build genus-level median proxies (for SMD diagnostic)
genus_proxy <- cbind(genus_name = bp_valid$genus_name,
                     bp_valid[, proxy_cols_for_smd, drop = FALSE]) %>%
  as.data.frame() %>%
  group_by(genus_name) %>%
  summarise(across(all_of(proxy_cols_for_smd), function(x) median(as.numeric(x), na.rm = TRUE)),
            .groups = "drop")

# Join genus sel_probs to pool_df
pool_df <- merge(pool_df, genus_median, by.x = "genus", by.y = "genus_name", all.x = TRUE)
pool_df <- merge(pool_df, genus_proxy, by.x = "genus", by.y = "genus_name", all.x = TRUE)
n_matched <- sum(!is.na(pool_df$sel_prob_motility))
cat("Genomes with sel_prob match:", n_matched, "/", nrow(pool_df),
    "(", round(n_matched/nrow(pool_df)*100, 1), "%)\n")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Head definitions and phenotype → sel_prob mapping
# ══════════════════════════════════════════════════════════════════════════════

# Mapping: BacDive head -> phenomnar sel_prob column
head_sel_map <- list(
  motility              = "sel_prob_motility",
  oxygen_microaerophile = "sel_prob_aerophilicity",
  spore                 = "sel_prob_spore_formation",
  oxygen_growth         = "sel_prob_aerophilicity",
  oxygen_facultative    = "sel_prob_aerophilicity",
  oxygen_obligate       = "sel_prob_aerophilicity",
  pathogenicity_human   = "sel_prob_animal_pathogenicity",
  pathogenicity_animal  = "sel_prob_animal_pathogenicity",
  pathogenicity_plant   = "sel_prob_plant_pathogenicity"
)

binary_heads <- list(
  motility = list(pred_val = "is_motile_val", true_col = "is_motile_true"),
  oxygen_microaerophile = list(pred_val = "oxygen_microaerophile_val", true_col = "microaerophile"),
  spore = list(pred_val = "ability_spore_val", true_col = "ability_spore_true")
)

softmax2_heads <- list(
  oxygen_growth = list(pred_val = "oxygen_growth_val", true_cols = c("aerobe", "anaerobe")),
  oxygen_facultative = list(pred_val = "oxygen_facultative_val",
                            true_cols = c("facultative.aerobe", "facultative.anaerobe")),
  oxygen_obligate = list(pred_val = "oxygen_obligate_val",
                         true_cols = c("obligate.aerobe", "obligate.anaerobe"))
)

patho_heads <- list(
  pathogenicity_human  = list(pred_val = "pathogenicity_human_val",
                              true_col = "pathogenicity_human_true"),
  pathogenicity_animal = list(pred_val = "pathogenicity_animal_val",
                              true_col = "pathogenicity_animal_true"),
  pathogenicity_plant  = list(pred_val = "pathogenicity_plant_val",
                              true_col = "pathogenicity_plant_true")
)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Helper functions
# ══════════════════════════════════════════════════════════════════════════════

nll_binary_platt_ipw <- function(params, logits, y, w) {
  a <- params[1]; b <- params[2]
  p <- sigmoid(a * logits + b)
  eps <- 1e-7; p <- pmax(pmin(p, 1 - eps), eps)
  nll_i <- -(y * log(p) + (1 - y) * log(1 - p))
  sum(w * nll_i) / sum(w)
}

safe_opt_platt <- function(logits, y) {
  tryCatch({
    optim(c(1, 0), nll_binary_platt, logits = logits, y = y,
          method = "BFGS")$par
  }, error = function(e) c(1, 0))
}

safe_opt_platt_ipw <- function(logits, y, w) {
  tryCatch({
    optim(c(1, 0), nll_binary_platt_ipw, logits = logits, y = y, w = w,
          method = "BFGS")$par
  }, error = function(e) c(1, 0))
}

# Compute IPW weights for a given sel_prob column
compute_ipw_weights <- function(pi_hat, trim_quantiles = c(0.01, 0.99)) {
  # Unmapped genomes (NA) get w = 1
  w <- rep(1.0, length(pi_hat))
  has_pi <- !is.na(pi_hat) & pi_hat > 0
  if (sum(has_pi) == 0) return(w)

  # Stabilized Hajek weights: w = P(R=1) / pi_hat
  prev_R <- mean(has_pi)  # approximation: fraction with sel_prob
  w_raw <- rep(1.0, length(pi_hat))
  w_raw[has_pi] <- prev_R / pi_hat[has_pi]

  # Percentile-based trimming
  bounds <- quantile(w_raw[has_pi], trim_quantiles, na.rm = TRUE)
  w_raw[has_pi] <- pmax(pmin(w_raw[has_pi], bounds[2]), bounds[1])

  # Normalize so mean(w) = 1
  w_raw / mean(w_raw)
}

# Compute IPW weights with gamma sensitivity parameter
compute_ipw_weights_gamma <- function(pi_hat, gamma, trim_quantiles = c(0.01, 0.99)) {
  w <- rep(1.0, length(pi_hat))
  has_pi <- !is.na(pi_hat) & pi_hat > 0
  if (sum(has_pi) == 0 || gamma == 0) return(w)

  prev_R <- mean(has_pi)
  w_raw <- rep(1.0, length(pi_hat))
  w_raw[has_pi] <- prev_R / (pi_hat[has_pi]^gamma)

  bounds <- quantile(w_raw[has_pi], trim_quantiles, na.rm = TRUE)
  w_raw[has_pi] <- pmax(pmin(w_raw[has_pi], bounds[2]), bounds[1])
  w_raw / mean(w_raw)
}

# SMD: standardized mean difference
smd <- function(x, group) {
  m1 <- mean(x[group == 1], na.rm = TRUE)
  m0 <- mean(x[group == 0], na.rm = TRUE)
  v1 <- var(x[group == 1], na.rm = TRUE)
  v0 <- var(x[group == 0], na.rm = TRUE)
  abs(m1 - m0) / sqrt((v1 + v0) / 2)
}

# Weighted SMD
weighted_smd <- function(x, group, w) {
  w1 <- w[group == 1]; x1 <- x[group == 1]
  w0 <- w[group == 0]; x0 <- x[group == 0]
  wm1 <- weighted.mean(x1, w1, na.rm = TRUE)
  wm0 <- weighted.mean(x0, w0, na.rm = TRUE)
  # Weighted variance
  wv <- function(xx, ww) {
    ok <- !is.na(xx)
    xx <- xx[ok]; ww <- ww[ok]
    ww <- ww / sum(ww)
    sum(ww * (xx - weighted.mean(xx, ww))^2) / (1 - sum(ww^2))
  }
  wv1 <- wv(x1, w1); wv0 <- wv(x0, w0)
  abs(wm1 - wm0) / sqrt((wv1 + wv0) / 2)
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Extract data + weights per head (unified function)
# ══════════════════════════════════════════════════════════════════════════════

# Returns list(z, y, w, family, pi_hat, proxies) for a binary/linear head
extract_binary_head <- function(pool_df, pred_col, true_col, sel_col, is_linear = FALSE) {
  p_all <- as.numeric(pool_df[[pred_col]])
  y_raw <- as.numeric(pool_df[[true_col]])
  fam   <- pool_df$Phylum
  pi    <- pool_df[[sel_col]]

  if (is_linear) {
    ok <- !is.na(y_raw) & y_raw != -9999
    z <- p_all[ok]  # raw linear output
    y <- as.numeric(y_raw[ok] > 0)
  } else {
    ok <- !is.na(y_raw) & y_raw != -9999
    z <- logit(p_all[ok])
    y <- y_raw[ok]
  }

  w <- compute_ipw_weights(pi[ok])
  proxies <- pool_df[ok, proxy_cols_for_smd, drop = FALSE]
  list(z = z, y = y, w = w, family = fam[ok], pi_hat = pi[ok], proxies = proxies)
}

extract_softmax2_head <- function(pool_df, pred_col, true_cols, sel_col) {
  pmat <- parse_prob_vec(pool_df[[pred_col]])
  y1 <- as.numeric(pool_df[[true_cols[2]]])
  ok <- !is.na(y1) & y1 != -9999 &
    !is.na(as.numeric(pool_df[[true_cols[1]]])) &
    (as.numeric(pool_df[[true_cols[1]]]) + y1) > 0

  z <- logit(pmat[ok, 2])
  y <- y1[ok]
  pi <- pool_df[[sel_col]]
  w <- compute_ipw_weights(pi[ok])
  proxies <- pool_df[ok, proxy_cols_for_smd, drop = FALSE]
  list(z = z, y = y, w = w, family = pool_df$Phylum[ok], pi_hat = pi[ok], proxies = proxies)
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: DIAGNOSTIC D1 — Weight sanity
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== D1: Weight Diagnostics ===\n")
weight_diag_rows <- list()
weight_hist_data <- list()

all_heads <- c(names(binary_heads), names(softmax2_heads), names(patho_heads))

for (name in all_heads) {
  sel_col <- head_sel_map[[name]]
  pi <- pool_df[[sel_col]]
  has_pi <- !is.na(pi) & pi > 0
  if (sum(has_pi) < 10) next

  prev_R <- mean(has_pi)
  w_raw <- rep(1.0, length(pi))
  w_raw[has_pi] <- prev_R / pi[has_pi]

  bounds <- quantile(w_raw[has_pi], c(0.01, 0.99), na.rm = TRUE)
  w_trim <- w_raw
  w_trim[has_pi] <- pmax(pmin(w_raw[has_pi], bounds[2]), bounds[1])
  w_trim <- w_trim / mean(w_trim)

  n_eff <- sum(w_trim)^2 / sum(w_trim^2)

  row <- data.frame(
    head = name, n = length(pi), n_mapped = sum(has_pi),
    w_min = round(min(w_trim), 3), w_p05 = round(quantile(w_trim, 0.05), 3),
    w_p25 = round(quantile(w_trim, 0.25), 3), w_median = round(median(w_trim), 3),
    w_p75 = round(quantile(w_trim, 0.75), 3), w_p95 = round(quantile(w_trim, 0.95), 3),
    w_max = round(max(w_trim), 3),
    n_eff = round(n_eff, 0), n_eff_ratio = round(n_eff / length(pi), 3),
    stringsAsFactors = FALSE
  )
  weight_diag_rows[[name]] <- row
  weight_hist_data[[name]] <- data.frame(head = name, w = w_trim[has_pi])

  cat(sprintf("  %-25s  n=%d  mapped=%d  w: [%.2f, %.2f]  n_eff/n=%.3f\n",
              name, length(pi), sum(has_pi), min(w_trim), max(w_trim), n_eff/length(pi)))
}

weight_diag_df <- do.call(rbind, weight_diag_rows)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: DIAGNOSTIC D3 — Covariate balance (SMD)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== D3: Covariate Balance (SMD) ===\n")
smd_rows <- list()

for (name in all_heads) {
  sel_col <- head_sel_map[[name]]

  # Extract head-specific data
  if (name %in% names(binary_heads)) {
    h <- binary_heads[[name]]
    dd <- extract_binary_head(pool_df, h$pred_val, h$true_col, sel_col)
  } else if (name %in% names(softmax2_heads)) {
    h <- softmax2_heads[[name]]
    dd <- extract_softmax2_head(pool_df, h$pred_val, h$true_cols, sel_col)
  } else {
    h <- patho_heads[[name]]
    dd <- extract_binary_head(pool_df, h$pred_val, h$true_col, sel_col, is_linear = TRUE)
  }

  # Use y as group indicator (annotated = 1)
  # In BacDive context, all samples have labels, so group by propensity tertile instead
  # Actually: SMD between high-pi and low-pi groups (above/below median pi_hat)
  has_pi <- !is.na(dd$pi_hat)
  if (sum(has_pi) < 20) next

  pi_med <- median(dd$pi_hat[has_pi], na.rm = TRUE)
  group <- as.integer(dd$pi_hat >= pi_med)  # 1 = high sel_prob, 0 = low

  for (pcol in proxy_cols_for_smd) {
    x <- as.numeric(dd$proxies[[pcol]])
    smd_unw <- smd(x[has_pi], group[has_pi])
    smd_w   <- weighted_smd(x[has_pi], group[has_pi], dd$w[has_pi])

    smd_rows[[paste(name, pcol)]] <- data.frame(
      head = name, covariate = pcol,
      smd_unweighted = round(smd_unw, 4),
      smd_weighted = round(smd_w, 4),
      stringsAsFactors = FALSE
    )
  }
}

smd_df <- do.call(rbind, smd_rows)
cat(sprintf("  Mean unweighted SMD: %.4f\n  Mean weighted SMD:   %.4f\n",
            mean(smd_df$smd_unweighted, na.rm = TRUE),
            mean(smd_df$smd_weighted, na.rm = TRUE)))

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8: LOFO CV — Plain vs IPW Platt (binary + softmax2 + pathogenicity)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== LOFO CV: Plain vs IPW Platt ===\n")
all_results <- list()
plot_list <- list()

# --- Binary sigmoid heads ---
for (name in names(binary_heads)) {
  h <- binary_heads[[name]]
  sel_col <- head_sel_map[[name]]
  cat("LOFO CV:", name, "...\n")

  dd <- extract_binary_head(pool_df, h$pred_val, h$true_col, sel_col)

  families <- unique(dd$family)
  n <- length(dd$y)
  p_plain <- numeric(n); p_ipw <- numeric(n)
  lofo_a <- numeric(length(families)); lofo_b <- numeric(length(families))
  lofo_a_ipw <- numeric(length(families)); lofo_b_ipw <- numeric(length(families))

  for (fi in seq_along(families)) {
    fam <- families[fi]
    held <- dd$family == fam
    tr <- !held

    if (length(unique(dd$y[tr])) < 2) {
      p_plain[held] <- sigmoid(dd$z[held])
      p_ipw[held]   <- sigmoid(dd$z[held])
      lofo_a[fi] <- 1; lofo_b[fi] <- 0
      lofo_a_ipw[fi] <- 1; lofo_b_ipw[fi] <- 0
      next
    }

    P_fit <- safe_opt_platt(dd$z[tr], dd$y[tr])
    P_ipw <- safe_opt_platt_ipw(dd$z[tr], dd$y[tr], dd$w[tr])

    p_plain[held] <- sigmoid(P_fit[1] * dd$z[held] + P_fit[2])
    p_ipw[held]   <- sigmoid(P_ipw[1] * dd$z[held] + P_ipw[2])
    lofo_a[fi] <- P_fit[1]; lofo_b[fi] <- P_fit[2]
    lofo_a_ipw[fi] <- P_ipw[1]; lofo_b_ipw[fi] <- P_ipw[2]
  }

  acc <- mean((sigmoid(dd$z) > 0.5) == dd$y)
  ece_plain <- ece(p_plain, dd$y)
  ece_ipw   <- ece(p_ipw, dd$y)
  brier_plain <- mean((p_plain - dd$y)^2)
  brier_ipw   <- mean((p_ipw - dd$y)^2)

  global_P <- safe_opt_platt(dd$z, dd$y)
  global_I <- safe_opt_platt_ipw(dd$z, dd$y, dd$w)

  all_results[[name]] <- list(
    type = "binary_sigmoid", n = n, n_families = length(families), acc = acc,
    global_a = global_P[1], global_b = global_P[2],
    global_a_ipw = global_I[1], global_b_ipw = global_I[2],
    ece_plain = ece_plain, ece_ipw = ece_ipw,
    brier_plain = brier_plain, brier_ipw = brier_ipw,
    median_a = median(lofo_a), median_b = median(lofo_b),
    median_a_ipw = median(lofo_a_ipw), median_b_ipw = median(lofo_b_ipw),
    p_plain = p_plain, p_ipw = p_ipw, y = dd$y, pi_hat = dd$pi_hat
  )

  cat(sprintf("  N=%d  Acc=%.3f  ECE: plain=%.4f  IPW=%.4f  delta=%.4f\n",
              n, acc, ece_plain, ece_ipw, ece_plain - ece_ipw))

  if (acc >= ACC_THRESHOLD && n >= N_FILTER) {
    rd <- bind_rows(
      reliability_data(p_plain, dd$y) %>%
        mutate(method = paste0("Plain (a=", round(median(lofo_a), 2),
                               " b=", round(median(lofo_b), 2), ")")),
      reliability_data(p_ipw, dd$y) %>%
        mutate(method = paste0("IPW (a=", round(median(lofo_a_ipw), 2),
                               " b=", round(median(lofo_b_ipw), 2), ")"))
    )
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                        color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, " [LOFO]  ECE: plain=", round(ece_plain, 4),
                          " IPW=", round(ece_ipw, 4)),
           x = "Mean Predicted Prob", y = "Observed Frequency") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# --- Softmax-2 heads ---
for (name in names(softmax2_heads)) {
  h <- softmax2_heads[[name]]
  sel_col <- head_sel_map[[name]]
  cat("LOFO CV:", name, "...\n")

  dd <- extract_softmax2_head(pool_df, h$pred_val, h$true_cols, sel_col)

  families <- unique(dd$family)
  n <- length(dd$y)
  p_plain <- numeric(n); p_ipw <- numeric(n)
  lofo_a <- numeric(length(families)); lofo_b <- numeric(length(families))
  lofo_a_ipw <- numeric(length(families)); lofo_b_ipw <- numeric(length(families))

  for (fi in seq_along(families)) {
    fam <- families[fi]
    held <- dd$family == fam
    tr <- !held

    if (length(unique(dd$y[tr])) < 2) {
      p_plain[held] <- sigmoid(dd$z[held])
      p_ipw[held]   <- sigmoid(dd$z[held])
      lofo_a[fi] <- 1; lofo_b[fi] <- 0
      lofo_a_ipw[fi] <- 1; lofo_b_ipw[fi] <- 0
      next
    }

    P_fit <- safe_opt_platt(dd$z[tr], dd$y[tr])
    P_ipw <- safe_opt_platt_ipw(dd$z[tr], dd$y[tr], dd$w[tr])

    p_plain[held] <- sigmoid(P_fit[1] * dd$z[held] + P_fit[2])
    p_ipw[held]   <- sigmoid(P_ipw[1] * dd$z[held] + P_ipw[2])
    lofo_a[fi] <- P_fit[1]; lofo_b[fi] <- P_fit[2]
    lofo_a_ipw[fi] <- P_ipw[1]; lofo_b_ipw[fi] <- P_ipw[2]
  }

  acc <- mean((sigmoid(dd$z) > 0.5) == dd$y)
  ece_plain <- ece(p_plain, dd$y)
  ece_ipw   <- ece(p_ipw, dd$y)
  brier_plain <- mean((p_plain - dd$y)^2)
  brier_ipw   <- mean((p_ipw - dd$y)^2)

  global_P <- safe_opt_platt(dd$z, dd$y)
  global_I <- safe_opt_platt_ipw(dd$z, dd$y, dd$w)

  all_results[[name]] <- list(
    type = "softmax_2class", n = n, n_families = length(families), acc = acc,
    global_a = global_P[1], global_b = global_P[2],
    global_a_ipw = global_I[1], global_b_ipw = global_I[2],
    ece_plain = ece_plain, ece_ipw = ece_ipw,
    brier_plain = brier_plain, brier_ipw = brier_ipw,
    median_a = median(lofo_a), median_b = median(lofo_b),
    median_a_ipw = median(lofo_a_ipw), median_b_ipw = median(lofo_b_ipw),
    p_plain = p_plain, p_ipw = p_ipw, y = dd$y, pi_hat = dd$pi_hat
  )

  cat(sprintf("  N=%d  Acc=%.3f  ECE: plain=%.4f  IPW=%.4f  delta=%.4f\n",
              n, acc, ece_plain, ece_ipw, ece_plain - ece_ipw))

  if (acc >= ACC_THRESHOLD && n >= N_FILTER) {
    rd <- bind_rows(
      reliability_data(p_plain, dd$y) %>%
        mutate(method = paste0("Plain (a=", round(median(lofo_a), 2),
                               " b=", round(median(lofo_b), 2), ")")),
      reliability_data(p_ipw, dd$y) %>%
        mutate(method = paste0("IPW (a=", round(median(lofo_a_ipw), 2),
                               " b=", round(median(lofo_b_ipw), 2), ")"))
    )
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                        color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, " [LOFO]  ECE: plain=", round(ece_plain, 4),
                          " IPW=", round(ece_ipw, 4)),
           x = "Mean Predicted Prob", y = "Observed Frequency") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# --- Pathogenicity (linear -> Platt) heads ---
for (name in names(patho_heads)) {
  h <- patho_heads[[name]]
  sel_col <- head_sel_map[[name]]
  cat("LOFO CV:", name, "...\n")

  dd <- extract_binary_head(pool_df, h$pred_val, h$true_col, sel_col, is_linear = TRUE)

  families <- unique(dd$family)
  n <- length(dd$y)
  p_plain <- numeric(n); p_ipw <- numeric(n)
  lofo_a <- numeric(length(families)); lofo_b <- numeric(length(families))
  lofo_a_ipw <- numeric(length(families)); lofo_b_ipw <- numeric(length(families))

  for (fi in seq_along(families)) {
    fam <- families[fi]
    held <- dd$family == fam
    tr <- !held

    if (length(unique(dd$y[tr])) < 2) {
      p_plain[held] <- sigmoid(dd$z[held])
      p_ipw[held]   <- sigmoid(dd$z[held])
      lofo_a[fi] <- 1; lofo_b[fi] <- 0
      lofo_a_ipw[fi] <- 1; lofo_b_ipw[fi] <- 0
      next
    }

    P_fit <- safe_opt_platt(dd$z[tr], dd$y[tr])
    P_ipw <- safe_opt_platt_ipw(dd$z[tr], dd$y[tr], dd$w[tr])

    p_plain[held] <- sigmoid(P_fit[1] * dd$z[held] + P_fit[2])
    p_ipw[held]   <- sigmoid(P_ipw[1] * dd$z[held] + P_ipw[2])
    lofo_a[fi] <- P_fit[1]; lofo_b[fi] <- P_fit[2]
    lofo_a_ipw[fi] <- P_ipw[1]; lofo_b_ipw[fi] <- P_ipw[2]
  }

  acc <- mean((p_plain > 0.5) == dd$y)
  ece_plain <- ece(p_plain, dd$y)
  ece_ipw   <- ece(p_ipw, dd$y)
  brier_plain <- mean((p_plain - dd$y)^2)
  brier_ipw   <- mean((p_ipw - dd$y)^2)

  global_P <- safe_opt_platt(dd$z, dd$y)
  global_I <- safe_opt_platt_ipw(dd$z, dd$y, dd$w)

  all_results[[name]] <- list(
    type = "linear_platt", n = n, n_families = length(families), acc = acc,
    global_a = global_P[1], global_b = global_P[2],
    global_a_ipw = global_I[1], global_b_ipw = global_I[2],
    ece_plain = ece_plain, ece_ipw = ece_ipw,
    brier_plain = brier_plain, brier_ipw = brier_ipw,
    median_a = median(lofo_a), median_b = median(lofo_b),
    median_a_ipw = median(lofo_a_ipw), median_b_ipw = median(lofo_b_ipw),
    p_plain = p_plain, p_ipw = p_ipw, y = dd$y, pi_hat = dd$pi_hat
  )

  cat(sprintf("  N=%d  Acc=%.3f  ECE: plain=%.4f  IPW=%.4f  delta=%.4f\n",
              n, acc, ece_plain, ece_ipw, ece_plain - ece_ipw))

  if (acc >= ACC_THRESHOLD && n >= N_FILTER) {
    rd <- bind_rows(
      reliability_data(p_plain, dd$y) %>%
        mutate(method = paste0("Plain (a=", round(median(lofo_a), 2),
                               " b=", round(median(lofo_b), 2), ")")),
      reliability_data(p_ipw, dd$y) %>%
        mutate(method = paste0("IPW (a=", round(median(lofo_a_ipw), 2),
                               " b=", round(median(lofo_b_ipw), 2), ")"))
    )
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq,
                                        color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, " [LOFO]  ECE: plain=", round(ece_plain, 4),
                          " IPW=", round(ece_ipw, 4)),
           x = "Platt Predicted Prob", y = "Observed Frequency") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8b: TRIMMING SENSITIVITY — Compare P90, P95, P99 upper bounds
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== Trimming Sensitivity (P90 / P95 / P99 upper bound) ===\n")
trim_levels <- c(0.90, 0.95, 0.99)
trim_rows <- list()

# Helper: run LOFO CV for one head with given trim quantiles, return summary
run_lofo_ipw <- function(dd, trim_q) {
  # Recompute weights with the given trim quantiles
  w <- compute_ipw_weights(dd$pi_hat, trim_quantiles = c(1 - trim_q, trim_q))
  families <- unique(dd$family)
  n <- length(dd$y)
  p_ipw <- numeric(n)
  lofo_a_ipw <- numeric(length(families))
  lofo_b_ipw <- numeric(length(families))

  for (fi in seq_along(families)) {
    fam <- families[fi]
    held <- dd$family == fam
    tr <- !held
    if (length(unique(dd$y[tr])) < 2) {
      p_ipw[held] <- sigmoid(dd$z[held])
      lofo_a_ipw[fi] <- 1; lofo_b_ipw[fi] <- 0
      next
    }
    P_ipw <- safe_opt_platt_ipw(dd$z[tr], dd$y[tr], w[tr])
    p_ipw[held] <- sigmoid(P_ipw[1] * dd$z[held] + P_ipw[2])
    lofo_a_ipw[fi] <- P_ipw[1]; lofo_b_ipw[fi] <- P_ipw[2]
  }

  n_eff <- sum(w)^2 / sum(w^2)
  global_I <- safe_opt_platt_ipw(dd$z, dd$y, w)

  list(
    ece = ece(p_ipw, dd$y),
    brier = mean((p_ipw - dd$y)^2),
    n_eff = n_eff, n_eff_ratio = n_eff / n,
    w_max = max(w), w_min = min(w),
    global_a = global_I[1], global_b = global_I[2],
    median_a = median(lofo_a_ipw), median_b = median(lofo_b_ipw)
  )
}

for (name in all_heads) {
  sel_col <- head_sel_map[[name]]

  # Extract head data (without baking in any trim level)
  if (name %in% names(binary_heads)) {
    h <- binary_heads[[name]]
    p_all <- as.numeric(pool_df[[h$pred_val]])
    y_raw <- as.numeric(pool_df[[h$true_col]])
    ok <- !is.na(y_raw) & y_raw != -9999
    dd <- list(z = logit(p_all[ok]), y = y_raw[ok],
               family = pool_df$Phylum[ok], pi_hat = pool_df[[sel_col]][ok])
  } else if (name %in% names(softmax2_heads)) {
    h <- softmax2_heads[[name]]
    pmat <- parse_prob_vec(pool_df[[h$pred_val]])
    y1 <- as.numeric(pool_df[[h$true_cols[2]]])
    ok <- !is.na(y1) & y1 != -9999 &
      !is.na(as.numeric(pool_df[[h$true_cols[1]]])) &
      (as.numeric(pool_df[[h$true_cols[1]]]) + y1) > 0
    dd <- list(z = logit(pmat[ok, 2]), y = y1[ok],
               family = pool_df$Phylum[ok], pi_hat = pool_df[[sel_col]][ok])
  } else {
    h <- patho_heads[[name]]
    p_all <- as.numeric(pool_df[[h$pred_val]])
    y_raw <- as.numeric(pool_df[[h$true_col]])
    ok <- !is.na(y_raw) & y_raw != -9999
    dd <- list(z = p_all[ok], y = as.numeric(y_raw[ok] > 0),
               family = pool_df$Phylum[ok], pi_hat = pool_df[[sel_col]][ok])
  }

  # Plain Platt ECE (from existing results)
  ece_plain <- all_results[[name]]$ece_plain
  brier_plain <- all_results[[name]]$brier_plain

  for (tl in trim_levels) {
    res <- run_lofo_ipw(dd, tl)
    tag <- paste0("P", tl * 100)
    trim_rows[[paste(name, tag)]] <- data.frame(
      head = name, trim = tag, trim_upper = tl,
      n = length(dd$y),
      w_min = round(res$w_min, 3), w_max = round(res$w_max, 3),
      n_eff = round(res$n_eff, 0), n_eff_ratio = round(res$n_eff_ratio, 3),
      ECE_plain = round(ece_plain, 4), ECE_ipw = round(res$ece, 4),
      delta_ECE = round(ece_plain - res$ece, 4),
      Brier_plain = round(brier_plain, 4), Brier_ipw = round(res$brier, 4),
      a_ipw = round(res$global_a, 4), b_ipw = round(res$global_b, 4),
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf("  %-25s  ECE plain=%.4f | P90=%.4f  P95=%.4f  P99=%.4f | n_eff/n: P90=%.3f  P95=%.3f  P99=%.3f\n",
              name, ece_plain,
              trim_rows[[paste(name, "P90")]]$ECE_ipw,
              trim_rows[[paste(name, "P95")]]$ECE_ipw,
              trim_rows[[paste(name, "P99")]]$ECE_ipw,
              trim_rows[[paste(name, "P90")]]$n_eff_ratio,
              trim_rows[[paste(name, "P95")]]$n_eff_ratio,
              trim_rows[[paste(name, "P99")]]$n_eff_ratio))
}

trim_df <- do.call(rbind, trim_rows)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9: DIAGNOSTIC D4 — Stratified ECE by propensity tertile
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== D4: Stratified ECE by Propensity Tertile ===\n")
strat_ece_rows <- list()

for (name in names(all_results)) {
  r <- all_results[[name]]
  pi <- r$pi_hat
  has_pi <- !is.na(pi)
  if (sum(has_pi) < 30) next

  # Tertiles of pi_hat
  q33 <- quantile(pi[has_pi], c(1/3, 2/3), na.rm = TRUE)
  tertile <- rep(NA_character_, length(pi))
  tertile[has_pi & pi <= q33[1]] <- "Low"
  tertile[has_pi & pi > q33[1] & pi <= q33[2]] <- "Mid"
  tertile[has_pi & pi > q33[2]] <- "High"

  for (tert in c("Low", "Mid", "High")) {
    idx <- which(tertile == tert)
    if (length(idx) < 10) next
    ece_p <- ece(r$p_plain[idx], r$y[idx])
    ece_i <- ece(r$p_ipw[idx], r$y[idx])
    strat_ece_rows[[paste(name, tert)]] <- data.frame(
      head = name, tertile = tert, n = length(idx),
      pi_range = paste0("[", round(min(pi[idx], na.rm=TRUE), 3), ", ",
                        round(max(pi[idx], na.rm=TRUE), 3), "]"),
      ece_plain = round(ece_p, 4), ece_ipw = round(ece_i, 4),
      delta_ece = round(ece_p - ece_i, 4),
      stringsAsFactors = FALSE
    )
  }
}

strat_ece_df <- do.call(rbind, strat_ece_rows)
cat("  Stratified ECE (Low tertile = understudied species):\n")
print(strat_ece_df %>% filter(tertile == "Low") %>%
        select(head, n, pi_range, ece_plain, ece_ipw, delta_ece), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 10: DIAGNOSTIC D5 — Tipping point analysis (gamma sweep)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== D5: Tipping Point Analysis (gamma sweep) ===\n")
gammas <- c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5)
tipping_rows <- list()

for (name in names(all_results)) {
  r <- all_results[[name]]
  sel_col <- head_sel_map[[name]]

  # Re-extract data for global fit
  if (name %in% names(binary_heads)) {
    h <- binary_heads[[name]]
    dd <- extract_binary_head(pool_df, h$pred_val, h$true_col, sel_col)
  } else if (name %in% names(softmax2_heads)) {
    h <- softmax2_heads[[name]]
    dd <- extract_softmax2_head(pool_df, h$pred_val, h$true_cols, sel_col)
  } else {
    h <- patho_heads[[name]]
    dd <- extract_binary_head(pool_df, h$pred_val, h$true_col, sel_col, is_linear = TRUE)
  }

  for (g in gammas) {
    w_g <- compute_ipw_weights_gamma(dd$pi_hat, g)
    if (g == 0) {
      P <- safe_opt_platt(dd$z, dd$y)
    } else {
      P <- safe_opt_platt_ipw(dd$z, dd$y, w_g)
    }
    p_cal <- sigmoid(P[1] * dd$z + P[2])
    ece_g <- ece(p_cal, dd$y)

    tipping_rows[[paste(name, g)]] <- data.frame(
      head = name, gamma = g,
      a = round(P[1], 4), b = round(P[2], 4),
      ece = round(ece_g, 4),
      stringsAsFactors = FALSE
    )
  }
  cat(sprintf("  %-25s  a range: [%.3f, %.3f]  b range: [%.3f, %.3f]\n",
              name,
              min(sapply(gammas, function(g) tipping_rows[[paste(name, g)]]$a)),
              max(sapply(gammas, function(g) tipping_rows[[paste(name, g)]]$a)),
              min(sapply(gammas, function(g) tipping_rows[[paste(name, g)]]$b)),
              max(sapply(gammas, function(g) tipping_rows[[paste(name, g)]]$b))))
}

tipping_df <- do.call(rbind, tipping_rows)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 11: DIAGNOSTIC D6 — Parameter comparison table
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== D6: Parameter Comparison (Plain vs IPW) ===\n")
comparison_rows <- list()

for (name in names(all_results)) {
  r <- all_results[[name]]
  comparison_rows[[name]] <- data.frame(
    head = name, type = r$type, N = r$n,
    a_plain = round(r$global_a, 4), b_plain = round(r$global_b, 4),
    a_ipw = round(r$global_a_ipw, 4), b_ipw = round(r$global_b_ipw, 4),
    delta_a = round(r$global_a_ipw - r$global_a, 4),
    delta_b = round(r$global_b_ipw - r$global_b, 4),
    ECE_plain = round(r$ece_plain, 4), ECE_ipw = round(r$ece_ipw, 4),
    delta_ECE = round(r$ece_plain - r$ece_ipw, 4),
    Brier_plain = round(r$brier_plain, 4), Brier_ipw = round(r$brier_ipw, 4),
    delta_Brier = round(r$brier_plain - r$brier_ipw, 4),
    stringsAsFactors = FALSE
  )
}

comparison_df <- do.call(rbind, comparison_rows)
cat("\n")
print(comparison_df, row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 12: Save outputs
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== Saving outputs ===\n")

csv_comparison <- file.path(cfg$output_dir, "ipw_platt_comparison.csv")
csv_weight_diag <- file.path(cfg$output_dir, "ipw_weight_diagnostics.csv")
csv_smd <- file.path(cfg$output_dir, "ipw_smd_balance.csv")
csv_tipping <- file.path(cfg$output_dir, "ipw_tipping_point.csv")
csv_strat <- file.path(cfg$output_dir, "ipw_stratified_ece.csv")
csv_trim <- file.path(cfg$output_dir, "ipw_trimming_sensitivity.csv")
rds_out <- file.path(cfg$output_dir, "ipw_platt_results.rds")

pdf_weights <- file.path(cfg$output_dir, "fig_ipw_weights.pdf")
pdf_calibration <- file.path(cfg$output_dir, "fig_ipw_calibration.pdf")
pdf_trimming <- file.path(cfg$output_dir, "fig_ipw_trimming.pdf")
pdf_tipping <- file.path(cfg$output_dir, "fig_ipw_tipping.pdf")

# CSV tables
write.csv(comparison_df, csv_comparison, row.names = FALSE)
write.csv(weight_diag_df, csv_weight_diag, row.names = FALSE)
write.csv(smd_df, csv_smd, row.names = FALSE)
write.csv(tipping_df, csv_tipping, row.names = FALSE)
write.csv(strat_ece_df, csv_strat, row.names = FALSE)
write.csv(trim_df, csv_trim, row.names = FALSE)

cat("  CSVs saved to:", cfg$output_dir, "\n")

# PDF: Weight distributions
pdf(pdf_weights, width = 12, height = 8)
whist <- do.call(rbind, weight_hist_data)
print(
  ggplot(whist, aes(x = w)) +
    geom_histogram(bins = 50, fill = "#4292C6", alpha = 0.7) +
    facet_wrap(~head, scales = "free_y") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    labs(title = "IPW Weight Distributions (trimmed, stabilized)",
         x = "Weight", y = "Count") +
    theme_minimal()
)
dev.off()
cat("  fig_ipw_weights.pdf saved\n")

# PDF: Reliability diagrams
pdf(pdf_calibration, width = 14, height = 9, onefile = TRUE)
for (p in plot_list) print(p)
dev.off()
cat("  fig_ipw_calibration.pdf saved\n")

# PDF: Trimming sensitivity
pdf(pdf_trimming, width = 14, height = 10)

# Panel 1: ECE by trim level
trim_ece_long <- trim_df %>%
  select(head, trim, ECE_plain, ECE_ipw) %>%
  tidyr::pivot_longer(cols = c(ECE_plain, ECE_ipw),
                      names_to = "method", values_to = "ECE") %>%
  mutate(method = ifelse(method == "ECE_plain", "Plain", paste0("IPW-", trim)),
         method = factor(method))
# Keep only IPW rows + one plain per head
trim_plot_df <- bind_rows(
  trim_df %>% filter(trim == "P99") %>%
    transmute(head, trim = "Plain", ECE = ECE_plain, n_eff_ratio = 1.0, w_max = 1.0),
  trim_df %>%
    transmute(head, trim, ECE = ECE_ipw, n_eff_ratio, w_max)
)
trim_plot_df$trim <- factor(trim_plot_df$trim, levels = c("Plain", "P90", "P95", "P99"))

p_trim_ece <- ggplot(trim_plot_df, aes(x = head, y = ECE, fill = trim)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(aes(label = round(ECE, 4)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 2.2) +
  scale_fill_manual(values = c("Plain" = "#999999", "P90" = "#E6550D",
                                "P95" = "#FD8D3C", "P99" = "#FDAE6B")) +
  labs(title = "Trimming Sensitivity: ECE by upper percentile bound",
       subtitle = "P90 = most aggressive trimming (higher n_eff, weaker correction)",
       x = NULL, y = "ECE (LOFO CV)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
print(p_trim_ece)

# Panel 2: n_eff/n by trim level
p_trim_neff <- ggplot(trim_df, aes(x = head, y = n_eff_ratio, fill = trim)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(aes(label = round(n_eff_ratio, 3)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 2.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("P90" = "#E6550D", "P95" = "#FD8D3C", "P99" = "#FDAE6B")) +
  labs(title = "Trimming Sensitivity: Effective sample size ratio (n_eff / n)",
       subtitle = "Red dashed = 0.5 threshold. Higher is more stable.",
       x = NULL, y = "n_eff / n") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
print(p_trim_neff)

dev.off()
cat("  fig_ipw_trimming.pdf saved\n")

# PDF: Tipping point
pdf(pdf_tipping, width = 14, height = 10)
tp_long <- tipping_df %>%
  tidyr::pivot_longer(cols = c(a, b, ece), names_to = "param", values_to = "value")
print(
  ggplot(tp_long, aes(x = gamma, y = value, color = head)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    facet_wrap(~param, scales = "free_y", ncol = 1) +
    labs(title = "Tipping Point Analysis: Platt parameters vs IPW strength (gamma)",
         subtitle = "gamma=0: plain Platt, gamma=1: standard IPW",
         x = "gamma", y = "Parameter value") +
    theme_minimal() + theme(legend.position = "bottom")
)
dev.off()
cat("  fig_ipw_tipping.pdf saved\n")

# RDS: Full results
saveRDS(list(
  comparison = comparison_df,
  weight_diagnostics = weight_diag_df,
  smd_balance = smd_df,
  tipping_point = tipping_df,
  trim_sensitivity = trim_df,
  stratified_ece = strat_ece_df,
  per_head_results = all_results,
  n_pool = nrow(pool_df)
), rds_out)
cat("  ipw_platt_results.rds saved\n")

cat("\n=== IPW-Platt Scaling Complete ===\n")
