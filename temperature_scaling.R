library(ggplot2)
library(dplyr)
library(gridExtra)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg)) else normalizePath(getwd())
proj_dir <- if (dir.exists(script_path)) script_path else dirname(script_path)
source(file.path(proj_dir, "R", "config_paths.R"))
source(file.path(proj_dir, "R", "calibration_common.R"))
cfg <- load_bacdive_config(proj_dir, required = c("bacdive_complete_rds"))

ACC_THRESHOLD <- 0.6  # heads below this get flagged and excluded from plots

# ── Load pre-computed validation predictions ──────────────────────────────────
l_all <- readRDS(cfg$bacdive_complete_rds)
l <- l_all$l
pred_df <- l$pred_df; true_df <- l$true_df; sw_df <- l$sw_df

cat("Model:", l_all$run_name, "Epoch:", l_all$ep, "\n")
cat("Total validation samples:", nrow(pred_df), "\n")

# ── Cal/Test split (70/30) ────────────────────────────────────────────────────
set.seed(42)
n <- nrow(pred_df)
cal_idx <- sample(n, size = floor(0.7 * n))
test_idx <- setdiff(1:n, cal_idx)
cat("Cal split:", length(cal_idx), " Test split:", length(test_idx), "\n")
cat("Accuracy threshold for plots:", ACC_THRESHOLD, "\n\n")

# ── Helper functions ──────────────────────────────────────────────────────────
brier_score <- function(p, y) mean((p - y)^2)

# ── Target definitions ────────────────────────────────────────────────────────
target_split <- list(
  cellshape = c("is_cell_shape_rod.shaped","is_cell_shape_other","is_cell_shape_coccus.shaped",
    "is_cell_shape_vibrio.shaped","is_cell_shape_filament.shaped","is_cell_shape_sphere.shaped",
    "is_cell_shape_ovoid.shaped","is_cell_shape_pleomorphic.shaped","is_cell_shape_spiral.shaped",
    "is_cell_shape_curved.shaped","is_cell_shape_oval.shaped"),
  flagellum = c("is_flagellum_arrangement_monotrichous","is_flagellum_arrangement_monotrichous_polar",
    "is_flagellum_arrangement_polar","is_flagellum_arrangement_peritrichous",
    "is_flagellum_arrangement_lophotrichous","is_flagellum_arrangement_gliding"),
  gram = c("is_gram_stain_positive","is_gram_stain_variable","is_gram_stain_negative"),
  oxygen_growth = c("aerobe","anaerobe"),
  oxygen_facultative = c("facultative.aerobe","facultative.anaerobe"),
  oxygen_obligate = c("obligate.aerobe","obligate.anaerobe")
)

# ── Collect all results into one list ─────────────────────────────────────────
all_rows <- list()     # for summary table
plot_list <- list()    # reliability diagrams (only for passing heads)

# ── Helper: fit binary head (sigmoid or 2-class softmax or linear) ────────────
fit_binary <- function(name, p_all, y_all, sw_all, type_label, is_linear = FALSE) {
  p_c <- p_all[cal_idx]; y_c <- y_all[cal_idx]; sw_c <- sw_all[cal_idx]
  ok_c <- sw_c != 0; p_c <- p_c[ok_c]; y_c <- y_c[ok_c]
  p_t <- p_all[test_idx]; y_t <- y_all[test_idx]; sw_t <- sw_all[test_idx]
  ok_t <- sw_t != 0; p_t <- p_t[ok_t]; y_t <- y_t[ok_t]

  if (is_linear) {
    # Platt on raw linear output
    opt_P <- optim(c(1, 0), nll_binary_platt, logits = p_c, y = y_c, method = "BFGS")$par
    p_t_platt <- sigmoid(opt_P[1] * p_t + opt_P[2])
    acc <- mean((p_t_platt > 0.5) == y_t)
    row <- data.frame(
      target = name, type = type_label, N_test = length(y_t), Acc = round(acc, 3),
      ECE_uncal = NA, ECE_temp = NA, ECE_platt = round(ece(p_t_platt, y_t), 4),
      NLL_uncal = NA, NLL_temp = NA, NLL_platt = round(nll_binary_platt(opt_P, p_t, y_t), 4),
      T = NA, a = round(opt_P[1], 3), b = round(opt_P[2], 3),
      stringsAsFactors = FALSE)
    # Plot data
    pdata <- list(p_t = NULL, p_t_platt = p_t_platt, y_t = y_t, opt_T = NA, opt_P = opt_P)
  } else {
    z_c <- logit(p_c); z_t <- logit(p_t)
    opt_T <- optimize(nll_binary_temp, c(0.01, 20), logits = z_c, y = y_c)$minimum
    opt_P <- optim(c(1, 0), nll_binary_platt, logits = z_c, y = y_c, method = "BFGS")$par
    p_t_temp  <- sigmoid(z_t / opt_T)
    p_t_platt <- sigmoid(opt_P[1] * z_t + opt_P[2])
    acc <- mean((p_t > 0.5) == y_t)
    row <- data.frame(
      target = name, type = type_label, N_test = length(y_t), Acc = round(acc, 3),
      ECE_uncal = round(ece(p_t, y_t), 4), ECE_temp = round(ece(p_t_temp, y_t), 4),
      ECE_platt = round(ece(p_t_platt, y_t), 4),
      NLL_uncal = round(nll_binary_temp(1, z_t, y_t), 4),
      NLL_temp = round(nll_binary_temp(opt_T, z_t, y_t), 4),
      NLL_platt = round(nll_binary_platt(opt_P, z_t, y_t), 4),
      T = round(opt_T, 3), a = round(opt_P[1], 3), b = round(opt_P[2], 3),
      stringsAsFactors = FALSE)
    pdata <- list(p_t = p_t, p_t_temp = p_t_temp, p_t_platt = p_t_platt,
                  y_t = y_t, opt_T = opt_T, opt_P = opt_P)
  }
  list(row = row, pdata = pdata)
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Binary sigmoid heads
# ══════════════════════════════════════════════════════════════════════════════

binary_heads <- list(
  motility = list(pred = "is_motile_pred", true = "is_motile", sw = "motility"),
  oxygen_microaerophile = list(pred = "microaerophile_pred", true = "microaerophile", sw = "oxygen_microaerophile"),
  spore = list(pred = "ability_spore_pred", true = "ability_spore", sw = "spore")
)

for (name in names(binary_heads)) {
  h <- binary_heads[[name]]
  res <- fit_binary(name, pred_df[[h$pred]], true_df[[h$true]], sw_df[[h$sw]], "binary_sigmoid")
  all_rows[[name]] <- res$row
  pd <- res$pdata
  if (res$row$Acc >= ACC_THRESHOLD) {
    rd <- bind_rows(
      reliability_data(pd$p_t, pd$y_t) %>% mutate(method = "Uncalibrated"),
      reliability_data(pd$p_t_temp, pd$y_t) %>% mutate(method = paste0("Temp (T=", round(pd$opt_T, 2), ")")),
      reliability_data(pd$p_t_platt, pd$y_t) %>% mutate(method = paste0("Platt (a=", round(pd$opt_P[1], 2), " b=", round(pd$opt_P[2], 2), ")"))
    )
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq, color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, "  (N=", res$row$N_test, ", Acc=", res$row$Acc, ")"),
           subtitle = paste0("ECE: ", res$row$ECE_uncal, " -> Temp:", res$row$ECE_temp, " / Platt:", res$row$ECE_platt),
           x = "Mean Predicted Prob", y = "Observed Frequency") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: 2-class softmax heads
# ══════════════════════════════════════════════════════════════════════════════

softmax2_heads <- list(
  oxygen_growth = list(cols = c("aerobe", "anaerobe"), sw = "oxygen_growth"),
  oxygen_facultative = list(cols = c("facultative.aerobe", "facultative.anaerobe"), sw = "oxygen_facultative"),
  oxygen_obligate = list(cols = c("obligate.aerobe", "obligate.anaerobe"), sw = "oxygen_obligate")
)

for (name in names(softmax2_heads)) {
  h <- softmax2_heads[[name]]
  cols_t <- h$cols; cols_p <- paste0(cols_t, "_pred")
  res <- fit_binary(name, pred_df[[cols_p[2]]], true_df[[cols_t[2]]], sw_df[[h$sw]], "softmax_2class")
  all_rows[[name]] <- res$row
  pd <- res$pdata
  if (res$row$Acc >= ACC_THRESHOLD) {
    rd <- bind_rows(
      reliability_data(pd$p_t, pd$y_t) %>% mutate(method = "Uncalibrated"),
      reliability_data(pd$p_t_temp, pd$y_t) %>% mutate(method = paste0("Temp (T=", round(pd$opt_T, 2), ")")),
      reliability_data(pd$p_t_platt, pd$y_t) %>% mutate(method = paste0("Platt (a=", round(pd$opt_P[1], 2), " b=", round(pd$opt_P[2], 2), ")"))
    )
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq, color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, "  (N=", res$row$N_test, ", Acc=", res$row$Acc, ")"),
           subtitle = paste0("ECE: ", res$row$ECE_uncal, " -> Temp:", res$row$ECE_temp, " / Platt:", res$row$ECE_platt),
           x = "Mean Predicted Prob", y = "Observed Frequency") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Multiclass softmax heads
# ══════════════════════════════════════════════════════════════════════════════

multi_heads <- list(
  gram = list(cols = target_split$gram, sw = "gram"),
  cellshape = list(cols = target_split$cellshape, sw = "cellshape"),
  flagellum = list(cols = target_split$flagellum, sw = "flagellum")
)

for (name in names(multi_heads)) {
  h <- multi_heads[[name]]
  cols_t <- h$cols; cols_p <- paste0(cols_t, "_pred")
  sw_all <- sw_df[[h$sw]]

  p_mat_all <- as.matrix(pred_df[, cols_p])
  y_mat_all <- as.matrix(true_df[, cols_t])

  p_c <- p_mat_all[cal_idx, ]; y_c <- y_mat_all[cal_idx, ]; sw_c <- sw_all[cal_idx]
  ok_c <- sw_c != 0; p_c <- p_c[ok_c, ]; y_c <- y_c[ok_c, ]
  p_te <- p_mat_all[test_idx, ]; y_te <- y_mat_all[test_idx, ]; sw_te <- sw_all[test_idx]
  ok_te <- sw_te != 0; p_te <- p_te[ok_te, ]; y_te <- y_te[ok_te, ]

  eps <- 1e-7
  z_c <- log(pmax(p_c, eps)); z_te <- log(pmax(p_te, eps))
  opt_T <- optimize(nll_multi_temp, c(0.01, 20), z_mat = z_c, y_mat = y_c)$minimum
  p_te_temp <- softmax_t(z_te, opt_T)

  acc <- mean(apply(p_te, 1, which.max) == apply(y_te, 1, which.max))
  ece_conf_before <- ece_confidence(p_te, y_te)
  ece_conf_after  <- ece_confidence(p_te_temp, y_te)
  nll_before <- nll_multi_temp(1, z_te, y_te)
  nll_after  <- nll_multi_temp(opt_T, z_te, y_te)

  row <- data.frame(
    target = name, type = paste0("multiclass_", ncol(p_te)), N_test = nrow(p_te), Acc = round(acc, 3),
    ECE_uncal = round(ece_conf_before, 4), ECE_temp = round(ece_conf_after, 4), ECE_platt = NA,
    NLL_uncal = round(nll_before, 4), NLL_temp = round(nll_after, 4), NLL_platt = NA,
    T = round(opt_T, 3), a = NA, b = NA,
    stringsAsFactors = FALSE)
  all_rows[[name]] <- row

  if (acc >= ACC_THRESHOLD) {
    conf_before <- apply(p_te, 1, max)
    correct_before <- (apply(p_te, 1, which.max) == apply(y_te, 1, which.max)) * 1
    conf_after <- apply(p_te_temp, 1, max)
    correct_after <- (apply(p_te_temp, 1, which.max) == apply(y_te, 1, which.max)) * 1
    rd <- bind_rows(
      reliability_data(conf_before, correct_before) %>% mutate(method = "Uncalibrated"),
      reliability_data(conf_after, correct_after) %>% mutate(method = paste0("Temp (T=", round(opt_T, 2), ")"))
    )
    rd$method <- factor(rd$method, levels = unique(rd$method))
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq, color = method, size = n)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, " [", ncol(p_te), "-class]  (N=", nrow(p_te), ", Acc=", round(acc, 3), ")"),
           subtitle = paste0("Conf ECE: ", round(ece_conf_before, 3), " -> ", round(ece_conf_after, 3),
                             "  |  NLL: ", round(nll_before, 3), " -> ", round(nll_after, 3)),
           x = "Confidence (max prob)", y = "Accuracy") +
      theme_minimal() + theme(legend.position = "bottom")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Pathogenicity heads (linear -> Platt, binarized y > 0)
# ══════════════════════════════════════════════════════════════════════════════

patho_heads <- list(
  pathogenicity_human  = list(pred = "pathogenicity_human_pred",  true = "pathogenicity_human",  sw = "pathogenicity_human"),
  pathogenicity_animal = list(pred = "pathogenicity_animal_pred", true = "pathogenicity_animal", sw = "pathogenicity_animal"),
  pathogenicity_plant  = list(pred = "pathogenicity_plant_pred",  true = "pathogenicity_plant",  sw = "pathogenicity_plant")
)

for (name in names(patho_heads)) {
  h <- patho_heads[[name]]
  y_bin <- as.numeric(true_df[[h$true]] > 0)
  res <- fit_binary(name, pred_df[[h$pred]], y_bin, sw_df[[h$sw]], "linear_platt", is_linear = TRUE)
  all_rows[[name]] <- res$row
  pd <- res$pdata
  if (res$row$Acc >= ACC_THRESHOLD) {
    rd <- reliability_data(pd$p_t_platt, pd$y_t)
    plot_list[[name]] <- ggplot(rd, aes(x = mean_pred, y = obs_freq, size = n)) +
      geom_point(color = "#7CAE00", alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(name, " [linear->Platt]  (N=", res$row$N_test, ", Acc=", res$row$Acc, ")"),
           subtitle = paste0("ECE: ", res$row$ECE_platt, "  NLL: ", res$row$NLL_platt),
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

# ══════════════════════════════════════════════════════════════════════════════
# SAVE PDF: page 1 = table, then reliability diagrams for passing heads
# ══════════════════════════════════════════════════════════════════════════════

# Format table for display
tbl <- summary_df %>%
  mutate(across(where(is.numeric), ~ifelse(is.na(.), "-", as.character(.)))) %>%
  mutate(quality = ifelse(quality == "FAIL", "FAIL", ""))

theme_tbl <- ttheme_minimal(
  core = list(fg_params = list(fontsize = 8, hjust = 0, x = 0.02),
              bg_params = list(fill = ifelse(summary_df$quality == "FAIL", "#FFE0E0", "white"))),
  colhead = list(fg_params = list(fontsize = 8, fontface = "bold", hjust = 0, x = 0.02),
                 bg_params = list(fill = "#E8E8E8"))
)

pdf_path <- file.path(cfg$reports_dir, "calibration_results.pdf")
pdf(pdf_path, width = 14, height = max(4, 0.4 * nrow(tbl) + 1.5))

# Page 1: summary table
grid::grid.newpage()
grid::grid.text("Calibration Summary - bacdive_attn_lstm_sum_28 Ep.23",
                x = 0.5, y = 0.97, gp = grid::gpar(fontsize = 14, fontface = "bold"))
grid::grid.text(paste0("Cal/Test split: 70/30 (seed=42). Acc threshold for plots: ", ACC_THRESHOLD,
                        ". Red rows excluded from reliability diagrams."),
                x = 0.5, y = 0.93, gp = grid::gpar(fontsize = 9, col = "grey40"))
tg <- tableGrob(tbl, rows = NULL, theme = theme_tbl)
tg$widths <- unit(rep(1/ncol(tbl), ncol(tbl)), "npc")
grid::pushViewport(grid::viewport(y = 0.45, height = 0.8))
grid::grid.draw(tg)
grid::popViewport()

dev.off()

# Reopen and append plots
pdf(pdf_path, width = 9, height = 6, onefile = TRUE)

# Page 1 again (pdf onefile needs it in the same call)
grid::grid.newpage()
grid::grid.text("Calibration Summary - bacdive_attn_lstm_sum_28 Ep.23",
                x = 0.5, y = 0.97, gp = grid::gpar(fontsize = 14, fontface = "bold"))
grid::grid.text(paste0("Cal/Test split: 70/30 (seed=42). Acc threshold for plots: ", ACC_THRESHOLD,
                        ". Red rows excluded from reliability diagrams."),
                x = 0.5, y = 0.93, gp = grid::gpar(fontsize = 9, col = "grey40"))
grid::pushViewport(grid::viewport(y = 0.45, height = 0.8))
grid::grid.draw(tg)
grid::popViewport()

# Reliability diagrams for passing heads
for (p in plot_list) print(p)

dev.off()
cat("\nPlots saved to:", pdf_path, "\n")
cat("Passing heads plotted:", paste(names(plot_list), collapse = ", "), "\n")

# ── Save results ──────────────────────────────────────────────────────────────
all_results <- list(summary = summary_df, all_rows = all_rows,
                    seed = 42, cal_idx = cal_idx, test_idx = test_idx,
                    acc_threshold = ACC_THRESHOLD)
rds_path <- file.path(cfg$reports_dir, "calibration_results.rds")
saveRDS(all_results, rds_path)
cat("Results saved to:", rds_path, "\n")
