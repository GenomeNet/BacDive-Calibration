#!/usr/bin/env Rscript
# Phase 3: MAR Selection Model (Extended Proxies)
# Fits logistic regression models predicting annotation status.
# M1 = taxonomy only; M2_old = 3 original proxies; M2 = extended proxies.
# Cross-phenotype indicators (old M3) removed — conceptually invalid.
# Includes collinearity diagnostics, old-vs-new comparison, propensity plots.

library(dplyr)
library(ggplot2)
library(pROC)
library(car)  # for vif()

# ---- Paths & config ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg)) else getwd()
proj_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))
source(file.path(proj_dir, "R", "selection", "00_config.R"))

paths <- resolve_paths()
data_dir <- paths$data_dir; fig_dir <- paths$fig_dir; out_dir <- paths$out_dir

df <- readRDS(file.path(data_dir, "bugphyzz_with_proxies.rds"))
cat("Loaded:", nrow(df), "rows\n")

# ---- Determine which phenotypes to analyse ----
target_phenotypes <- parse_phenotype_args(df, pool = "all")
cat("Target phenotypes:", paste(target_phenotypes, collapse = ", "), "\n\n")

# ---- Drop species without valid description year (not imputable) ----
n_before <- nrow(df)
df <- df[!is.na(df$year_valid_pub), ]
cat(sprintf("Dropped %d species without year_valid_pub (%d → %d)\n\n",
            n_before - nrow(df), n_before, nrow(df)))

# ---- Prepare shared predictors ----
df$Phylum_coll <- collapse_rare(df$Phylum, min_count = 50)
cat("Phyla after collapsing rare (<50 spp):", nlevels(df$Phylum_coll), "levels\n")

# Check proxy availability — original and extended separately
proxy_cols_old <- c()
for (pc in PROXY_COLS_ORIGINAL) {
  if (pc %in% names(df) && !all(is.na(df[[pc]]))) proxy_cols_old <- c(proxy_cols_old, pc)
}
proxy_cols_new <- c()
for (pc in PROXY_COLS_EXTENDED) {
  if (pc %in% names(df) && !all(is.na(df[[pc]]))) proxy_cols_new <- c(proxy_cols_new, pc)
}
proxy_cols_all <- c(proxy_cols_old, proxy_cols_new)

# Imputed proxy sets (contig median / P95 + corrected checkm + indicators)
# These have no NAs (except log_species_age which was already dropped), so use all rows.
proxy_cols_imp_base <- c(proxy_cols_old, "log_species_age", "n_culture_collections",
                         "has_complete_genome", "has_assembly", "has_checkm_eval",
                         "checkm_imp")
# For contig: log1p of the imputed value
df$log_contig_imp_median <- log1p(df$contig_imp_median)
df$log_contig_imp_p95    <- log1p(df$contig_imp_p95)
proxy_cols_imp_median <- c(proxy_cols_imp_base, "log_contig_imp_median")
proxy_cols_imp_p95    <- c(proxy_cols_imp_base, "log_contig_imp_p95")

if (length(proxy_cols_old) == 0) stop("No original proxies available. Run Phase 2 first.")
cat("Original proxies:", paste(proxy_cols_old, collapse = ", "), "\n")
cat("Extended proxies:", paste(proxy_cols_new, collapse = ", "), "\n")
cat("All proxies (complete-case):", paste(proxy_cols_all, collapse = ", "), "\n")
cat("Imputed proxies (median):", paste(proxy_cols_imp_median, collapse = ", "), "\n")
cat("Imputed proxies (P95):   ", paste(proxy_cols_imp_p95, collapse = ", "), "\n\n")

# ---- Open summary sink ----
sink_file <- file.path(out_dir, "selection_model_summary.txt")
sink(sink_file, split = TRUE)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SELECTION MODEL RESULTS (Extended Proxies)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ========================================================================
# Collinearity Diagnostics (run once on all proxies)
# ========================================================================
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("PROXY COLLINEARITY DIAGNOSTICS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Correlation matrix of all continuous proxies
proxy_data <- df[, proxy_cols_all, drop = FALSE]
proxy_complete <- proxy_data[complete.cases(proxy_data), ]
cat("Rows with all proxies complete:", nrow(proxy_complete), "/", nrow(df), "\n\n")

cat("--- Correlation matrix (Pearson) ---\n")
cor_mat <- cor(proxy_complete, use = "complete.obs")
print(round(cor_mat, 3))

# Flag high correlations
cat("\nProxy pairs with |r| > 0.5:\n")
for (i in 1:(ncol(cor_mat) - 1)) {
  for (j in (i + 1):ncol(cor_mat)) {
    r <- cor_mat[i, j]
    if (abs(r) > 0.5) {
      cat(sprintf("  %s <-> %s: r = %.3f\n",
                  colnames(cor_mat)[i], colnames(cor_mat)[j], r))
    }
  }
}

# PCA on proxy space — how many dimensions do the proxies capture?
cat("\n--- PCA on proxy space ---\n")
pca <- prcomp(proxy_complete, scale. = TRUE)
var_explained <- summary(pca)$importance[2, ]  # proportion of variance
cum_var <- summary(pca)$importance[3, ]  # cumulative
cat("Variance explained by each PC:\n")
for (k in seq_along(var_explained)) {
  cat(sprintf("  PC%d: %.1f%% (cumulative: %.1f%%)\n",
              k, 100 * var_explained[k], 100 * cum_var[k]))
}

# Save PCA biplot
pdf(file.path(fig_dir, "fig03b_proxy_pca.pdf"), width = 8, height = 6)
biplot(pca, main = "PCA of Selection Proxies", cex = 0.5)
dev.off()
cat("\nSaved: fig03b_proxy_pca.pdf\n")

# VIF check — fit a test model using first phenotype
test_pheno <- target_phenotypes[1]
R_test <- paste0("R_", gsub(" ", "_", test_pheno))
df[[R_test]] <- as.integer(!is.na(df[[test_pheno]]))

test_cols <- c(R_test, "Phylum_coll", proxy_cols_all)
df_test <- df[, test_cols, drop = FALSE]
df_test <- df_test[complete.cases(df_test), ]

fml_test <- as.formula(paste0(R_test, " ~ Phylum_coll + ",
                               paste(proxy_cols_all, collapse = " + ")))
m_test <- glm(fml_test, family = binomial, data = df_test)
vif_vals <- tryCatch(vif(m_test), error = function(e) NULL)

# VIF returns GVIF for factors; extract GVIF^(1/(2*Df)) for comparability
if (!is.null(vif_vals)) {
  cat("\n--- VIF check (", test_pheno, ") ---\n")
  if (is.matrix(vif_vals)) {
    # GVIF output for models with factors
    cat("  Proxy VIFs (GVIF^(1/2Df)):\n")
    for (pc in proxy_cols_all) {
      if (pc %in% rownames(vif_vals)) {
        gvif_adj <- vif_vals[pc, "GVIF^(1/(2*Df))"]
        flag <- if (gvif_adj > 5) " *** HIGH" else ""
        cat(sprintf("    %-30s %.2f%s\n", pc, gvif_adj, flag))
      }
    }
  } else {
    # Simple VIF vector
    for (pc in proxy_cols_all) {
      if (pc %in% names(vif_vals)) {
        flag <- if (vif_vals[pc] > 5) " *** HIGH" else ""
        cat(sprintf("    %-30s %.2f%s\n", pc, vif_vals[pc], flag))
      }
    }
  }

  # Identify proxies to drop due to high VIF
  drop_proxies <- c()
  if (is.matrix(vif_vals)) {
    for (pc in proxy_cols_all) {
      if (pc %in% rownames(vif_vals) && vif_vals[pc, "GVIF^(1/(2*Df))"] > 5) {
        drop_proxies <- c(drop_proxies, pc)
      }
    }
  } else {
    for (pc in proxy_cols_all) {
      if (pc %in% names(vif_vals) && vif_vals[pc] > 5) {
        drop_proxies <- c(drop_proxies, pc)
      }
    }
  }

  if (length(drop_proxies) > 0) {
    cat("\n  Dropping due to VIF > 5:", paste(drop_proxies, collapse = ", "), "\n")
    proxy_cols_all <- setdiff(proxy_cols_all, drop_proxies)
    proxy_cols_new <- setdiff(proxy_cols_new, drop_proxies)
    cat("  Retained proxies:", paste(proxy_cols_all, collapse = ", "), "\n")
  } else {
    cat("\n  All proxies have VIF <= 5. No drops needed.\n")
  }
}

cat("\n")

# ---- Per-phenotype model fitting ----
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("PER-PHENOTYPE MODEL FITTING\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

all_roc_data <- list()
model_summary <- data.frame(stringsAsFactors = FALSE)

# Helper: 5-fold CV-AUC for a model specification
cv_auc_5fold <- function(fml, data, R_col, seed = 42) {
  set.seed(seed)
  folds <- sample(rep(1:5, length.out = nrow(data)))
  cv_aucs <- numeric(5)
  cv_preds <- rep(NA_real_, nrow(data))
  for (k in 1:5) {
    mk <- tryCatch(glm(fml, family = binomial, data = data[folds != k, ]),
                    error = function(e) NULL)
    if (!is.null(mk)) {
      pk <- predict(mk, newdata = data[folds == k, ], type = "response")
      cv_preds[folds == k] <- pk
      cv_aucs[k] <- as.numeric(auc(roc(data[[R_col]][folds == k], pk, quiet = TRUE)))
    } else cv_aucs[k] <- NA
  }
  list(mean_auc = mean(cv_aucs, na.rm = TRUE), cv_preds = cv_preds)
}

for (pheno in target_phenotypes) {
  cat(paste(rep("-", 60), collapse = ""), "\n")
  cat("Phenotype:", pheno, "\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")

  # Annotation indicator
  R_col <- paste0("R_", gsub(" ", "_", pheno))
  df[[R_col]] <- as.integer(!is.na(df[[pheno]]))

  # ---- Complete-case data (original flow) ----
  mod_cols <- c(R_col, "Phylum_coll", proxy_cols_all)
  df_cc <- df[complete.cases(df[, mod_cols, drop = FALSE]), ]

  # ---- Imputed data (all rows, no contig/checkm NAs) ----
  mod_cols_med <- c(R_col, "Phylum_coll", proxy_cols_imp_median)
  df_imp <- df[complete.cases(df[, mod_cols_med, drop = FALSE]), ]

  cat(sprintf("  Complete-case rows: %d  |  Imputed rows: %d (+%d)\n",
              nrow(df_cc), nrow(df_imp), nrow(df_imp) - nrow(df_cc)))

  # Formulas
  fml_base <- as.formula(paste0(R_col, " ~ Phylum_coll"))
  fml_old  <- as.formula(paste0(R_col, " ~ Phylum_coll + ",
                                 paste(proxy_cols_old, collapse = " + ")))
  fml_cc   <- as.formula(paste0(R_col, " ~ Phylum_coll + ",
                                 paste(proxy_cols_all, collapse = " + ")))
  fml_med  <- as.formula(paste0(R_col, " ~ Phylum_coll + ",
                                 paste(proxy_cols_imp_median, collapse = " + ")))
  fml_p95  <- as.formula(paste0(R_col, " ~ Phylum_coll + ",
                                 paste(proxy_cols_imp_p95, collapse = " + ")))

  # Fit models
  m1     <- glm(fml_base, family = binomial, data = df_cc)
  m2_old <- glm(fml_old,  family = binomial, data = df_cc)
  m2_cc  <- glm(fml_cc,   family = binomial, data = df_cc)
  m2_med <- glm(fml_med,  family = binomial, data = df_imp)
  m2_p95 <- glm(fml_p95,  family = binomial, data = df_imp)

  # McFadden R²
  ll0_cc  <- as.numeric(logLik(glm(as.formula(paste0(R_col, " ~ 1")),
                                    family = binomial, data = df_cc)))
  ll0_imp <- as.numeric(logLik(glm(as.formula(paste0(R_col, " ~ 1")),
                                    family = binomial, data = df_imp)))
  r2_m1     <- 1 - as.numeric(logLik(m1)) / ll0_cc
  r2_old    <- 1 - as.numeric(logLik(m2_old)) / ll0_cc
  r2_cc     <- 1 - as.numeric(logLik(m2_cc)) / ll0_cc
  r2_med    <- 1 - as.numeric(logLik(m2_med)) / ll0_imp
  r2_p95    <- 1 - as.numeric(logLik(m2_p95)) / ll0_imp

  cat(sprintf("  M1 (taxonomy)        AIC: %8.1f  R2: %.4f  (N=%d)\n", AIC(m1), r2_m1, nrow(df_cc)))
  cat(sprintf("  M2_old (3 proxy)     AIC: %8.1f  R2: %.4f  (N=%d)\n", AIC(m2_old), r2_old, nrow(df_cc)))
  cat(sprintf("  M2_cc (complete)     AIC: %8.1f  R2: %.4f  (N=%d)\n", AIC(m2_cc), r2_cc, nrow(df_cc)))
  cat(sprintf("  M2_imp_med (median)  AIC: %8.1f  R2: %.4f  (N=%d)\n", AIC(m2_med), r2_med, nrow(df_imp)))
  cat(sprintf("  M2_imp_p95 (P95)     AIC: %8.1f  R2: %.4f  (N=%d)\n", AIC(m2_p95), r2_p95, nrow(df_imp)))

  # Proxy ORs for M2_cc (reference model)
  m2_coef <- summary(m2_cc)$coefficients
  m2_ci <- confint.default(m2_cc)
  proxy_ors <- setNames(rep(NA_real_, length(proxy_cols_all)), proxy_cols_all)
  proxy_ps  <- setNames(rep(NA_real_, length(proxy_cols_all)), proxy_cols_all)
  for (pcol in proxy_cols_all) {
    if (pcol %in% rownames(m2_coef)) {
      row <- m2_coef[pcol, ]
      ci <- m2_ci[pcol, ]
      cat(sprintf("  %s: OR = %.3f [%.3f, %.3f], p = %s\n",
                  pcol, exp(row["Estimate"]),
                  exp(ci[1]), exp(ci[2]), signif(row["Pr(>|z|)"], 3)))
      proxy_ors[pcol] <- exp(row["Estimate"])
      proxy_ps[pcol]  <- row["Pr(>|z|)"]
    }
  }

  # 5-fold CV-AUC for all model variants
  cv_old <- cv_auc_5fold(fml_old, df_cc, R_col)
  cv_cc  <- cv_auc_5fold(fml_cc, df_cc, R_col)
  cv_med <- cv_auc_5fold(fml_med, df_imp, R_col)
  cv_p95 <- cv_auc_5fold(fml_p95, df_imp, R_col)

  cat(sprintf("  CV-AUC:  M2_old=%.3f  M2_cc=%.3f  M2_med=%.3f  M2_p95=%.3f\n",
              cv_old$mean_auc, cv_cc$mean_auc, cv_med$mean_auc, cv_p95$mean_auc))

  # ROC objects for plotting
  roc_m2_old <- roc(df_cc[[R_col]], predict(m2_old, type = "response"), quiet = TRUE)
  roc_m2_cc  <- roc(df_cc[[R_col]], predict(m2_cc, type = "response"), quiet = TRUE)

  # Store — use imputed median model for selection probabilities (best coverage)
  mod_rows_imp <- which(complete.cases(df[, mod_cols_med, drop = FALSE]))
  all_roc_data[[pheno]] <- list(
    roc_m2_old = roc_m2_old, roc_m2 = roc_m2_cc,
    cv_preds = cv_cc$cv_preds,
    df_cc = df_cc, df_imp = df_imp, R_col = R_col,
    m2_cc_model = m2_cc, m2_imp_model = m2_med,
    mod_rows_cc = which(complete.cases(df[, mod_cols, drop = FALSE])),
    mod_rows_imp = mod_rows_imp
  )

  row_data <- data.frame(
    phenotype = pheno, n_cc = nrow(df_cc), n_imp = nrow(df_imp),
    m1_r2 = round(r2_m1, 4), m2_old_r2 = round(r2_old, 4),
    m2_cc_r2 = round(r2_cc, 4), m2_med_r2 = round(r2_med, 4), m2_p95_r2 = round(r2_p95, 4),
    m2_old_auc = round(cv_old$mean_auc, 3),
    m2_cc_auc  = round(cv_cc$mean_auc, 3),
    m2_med_auc = round(cv_med$mean_auc, 3),
    m2_p95_auc = round(cv_p95$mean_auc, 3),
    stringsAsFactors = FALSE
  )
  for (pcol in proxy_cols_all) {
    short <- sub("^log_", "", pcol)
    row_data[[paste0(short, "_or")]] <- round(proxy_ors[pcol], 3)
    row_data[[paste0(short, "_p")]]  <- signif(proxy_ps[pcol], 3)
  }
  model_summary <- rbind(model_summary, row_data)
  cat("\n")
}

# ---- Extract selection probabilities ----
# Use imputed median model for best coverage (covers all non-missing-year species)
cat("\n=== Building selection probability matrix ===\n")
sel_probs <- data.frame(row_id = seq_len(nrow(df)))
for (pheno in names(all_roc_data)) {
  rd <- all_roc_data[[pheno]]
  col_name <- paste0("sel_prob_", safe_pheno_name(pheno))
  sel_probs[[col_name]] <- NA_real_
  sel_probs[[col_name]][rd$mod_rows_imp] <- predict(rd$m2_imp_model,
    newdata = df[rd$mod_rows_imp, , drop = FALSE], type = "response")
}
sel_probs$row_id <- NULL
saveRDS(sel_probs, file.path(data_dir, "selection_probs.rds"))
cat("Saved: data/selection_probs.rds (", nrow(sel_probs), " rows, ",
    ncol(sel_probs), " phenotype columns)\n", sep = "")
cat("  Coverage:", sum(!is.na(sel_probs[[1]])), "/", nrow(sel_probs), "non-NA\n\n")

# ========================================================================
# Figure 4: ROC curves — old M2 vs new M2 (panel: one per phenotype)
# ========================================================================
n_pheno <- length(all_roc_data)
n_col <- min(n_pheno, 3)
n_row <- ceiling(n_pheno / n_col)

pdf(file.path(fig_dir, "fig04_selection_model_roc.pdf"),
    width = 5 * n_col, height = 5 * n_row)
par(mfrow = c(n_row, n_col), mar = c(4, 4, 3, 1))

for (pheno in names(all_roc_data)) {
  rd <- all_roc_data[[pheno]]
  ms <- model_summary[model_summary$phenotype == pheno, ]
  plot(rd$roc_m2, main = SHORT_NAMES[pheno],
       col = "steelblue", lwd = 2, print.auc = FALSE)
  plot(rd$roc_m2_old, add = TRUE, col = "grey60", lwd = 2, lty = 2, print.auc = FALSE)
  legend("bottomright", legend = c(
    sprintf("M2 extended CV-AUC=%.3f", ms$m2_cc_auc),
    sprintf("M2 original CV-AUC=%.3f", ms$m2_old_auc)
  ), col = c("steelblue", "grey60"), lwd = 2, lty = c(1, 2), cex = 0.7)
}

dev.off()
cat("Saved: fig04_selection_model_roc.pdf\n")

# ========================================================================
# Figure 5: Out-of-fold calibration plot (M2)
# ========================================================================
cal_all <- list()
for (pheno in names(all_roc_data)) {
  rd <- all_roc_data[[pheno]]
  pp <- rd$df_cc
  pp$pred_prob <- rd$cv_preds
  pp$actual <- pp[[rd$R_col]]
  pp <- pp[!is.na(pp$pred_prob), ]

  breaks <- unique(quantile(pp$pred_prob, probs = seq(0, 1, by = 0.1)))
  if (length(breaks) < 3) breaks <- range(pp$pred_prob)
  pp$bin <- cut(pp$pred_prob, breaks = breaks, include.lowest = TRUE)
  cal <- pp %>%
    group_by(bin) %>%
    summarise(n = n(),
              mean_pred = mean(pred_prob),
              actual_rate = mean(actual),
              .groups = "drop") %>%
    filter(n >= 20)
  cal$phenotype <- pheno
  cal_all[[pheno]] <- cal
}

cal_combined <- do.call(rbind, cal_all)

p_cal <- ggplot(cal_combined, aes(x = mean_pred, y = actual_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = n), alpha = 0.7, color = "steelblue") +
  geom_line(alpha = 0.4, color = "steelblue") +
  scale_size_continuous(range = c(1.5, 6)) +
  facet_wrap(~ phenotype, scales = "free", ncol = n_col) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean Predicted Probability (out-of-fold, quantile bin)",
       y = "Observed Annotation Rate",
       title = "Selection Model Calibration (M2 Extended, 5-fold CV)",
       subtitle = "Diagonal = perfect calibration",
       size = "N species") +
  theme_minimal(base_size = 11)

ggsave(file.path(fig_dir, "fig05_calibration.pdf"),
       p_cal, width = 5 * n_col, height = 5 * n_row)
cat("Saved: fig05_calibration.pdf\n")

# ========================================================================
# Figure 5b: Propensity score distributions + KS test + overlap
# Compute for both complete-case (M2_cc) and imputed (M2_imp) models.
# ========================================================================
cat("\n=== Propensity score distributions ===\n")

# Helper: compute KS + overlap for a set of probabilities and labels
compute_prop_stats <- function(probs, ann, pheno, model_label) {
  p_ann   <- probs[ann == 1]
  p_unann <- probs[ann == 0]
  ks <- ks.test(p_ann, p_unann)
  grid <- seq(0, 1, length.out = 512)
  if (length(p_ann) >= 2 && length(p_unann) >= 2) {
    d_ann   <- approx(density(p_ann,   from = 0, to = 1, n = 512)$x,
                       density(p_ann,   from = 0, to = 1, n = 512)$y,
                       xout = grid)$y
    d_unann <- approx(density(p_unann, from = 0, to = 1, n = 512)$x,
                       density(p_unann, from = 0, to = 1, n = 512)$y,
                       xout = grid)$y
    d_ann[is.na(d_ann)] <- 0; d_unann[is.na(d_unann)] <- 0
    overlap <- sum(pmin(d_ann, d_unann)) * (grid[2] - grid[1])
  } else overlap <- NA_real_
  data.frame(phenotype = pheno, model = model_label, N = length(probs),
             n_ann = length(p_ann), n_unann = length(p_unann),
             mean_ann = round(mean(p_ann), 3), mean_unann = round(mean(p_unann), 3),
             KS_D = round(ks$statistic, 3), KS_p = signif(ks$p.value, 3),
             overlap = round(overlap, 3), stringsAsFactors = FALSE)
}

prop_all <- list()
prop_stats <- list()

for (pheno in names(all_roc_data)) {
  rd <- all_roc_data[[pheno]]

  # Complete-case model propensity
  probs_cc <- predict(rd$m2_cc_model, type = "response")
  ann_cc   <- rd$df_cc[[rd$R_col]]
  prop_stats[[paste0(pheno, "_cc")]] <- compute_prop_stats(probs_cc, ann_cc, pheno, "M2_cc")

  # Imputed model propensity (on full imputed data)
  probs_imp <- predict(rd$m2_imp_model, type = "response")
  ann_imp   <- rd$df_imp[[rd$R_col]]
  prop_stats[[paste0(pheno, "_imp")]] <- compute_prop_stats(probs_imp, ann_imp, pheno, "M2_imp_med")

  # Density plot data — use imputed model (better coverage)
  prop_all[[pheno]] <- data.frame(
    sel_prob = probs_imp,
    annotated = factor(ann_imp, levels = c(0, 1),
                       labels = c("Unannotated", "Annotated")),
    phenotype = pheno
  )
}

prop_stats_df <- do.call(rbind, prop_stats)
rownames(prop_stats_df) <- NULL

# Pivot to wide: one row per phenotype, cc vs imp side by side
prop_cc  <- prop_stats_df[prop_stats_df$model == "M2_cc", ]
prop_imp <- prop_stats_df[prop_stats_df$model == "M2_imp_med", ]
prop_wide <- merge(
  prop_cc[, c("phenotype", "N", "KS_D", "overlap")],
  prop_imp[, c("phenotype", "N", "KS_D", "overlap")],
  by = "phenotype", suffixes = c("_cc", "_imp")
)

# Merge imputed KS/overlap into model_summary
model_summary <- merge(model_summary,
  prop_imp[, c("phenotype", "KS_D", "KS_p", "overlap")],
  by = "phenotype", all.x = TRUE)

cat("\n--- Propensity separation: complete-case vs imputed ---\n")
print(prop_wide, row.names = FALSE)

cat("\n--- Full propensity stats ---\n")
print(prop_stats_df, row.names = FALSE)

cat("\n\n=== Model comparison summary ===\n")
# Compact view: AUC columns only
auc_view <- model_summary[, c("phenotype", "n_cc", "n_imp",
                               "m2_old_auc", "m2_cc_auc", "m2_med_auc", "m2_p95_auc",
                               "KS_D", "overlap")]
print(auc_view, row.names = FALSE)

sink()
cat("\nSaved:", sink_file, "\n")

# Save CSVs
write.csv(model_summary, file.path(out_dir, "selection_model_table.csv"), row.names = FALSE)
write.csv(prop_stats_df, file.path(out_dir, "propensity_stats.csv"), row.names = FALSE)
cat("Saved: output/selection_model_table.csv\n")
cat("Saved: output/propensity_stats.csv\n")

# Density plot — imputed model
prop_combined <- do.call(rbind, prop_all)

p_prop <- ggplot(prop_combined, aes(x = sel_prob, fill = annotated)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Unannotated" = "grey60", "Annotated" = "steelblue")) +
  facet_wrap(~ phenotype, scales = "free_y", ncol = n_col) +
  labs(x = "Estimated Selection Probability",
       y = "Density",
       title = "Propensity Score Distributions (M2 Imputed, Median)",
       fill = "Status") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fig05b_propensity_distributions.pdf"),
       p_prop, width = 5 * n_col, height = 5 * n_row)
cat("Saved: fig05b_propensity_distributions.pdf\n")

cat("\nPhase 3 complete.\n")
