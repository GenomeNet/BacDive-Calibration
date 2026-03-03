sigmoid <- function(z) 1 / (1 + exp(-z))

logit <- function(p, eps = 1e-7) {
  log(pmax(p, eps) / pmax(1 - p, eps))
}

softmax_t <- function(z_mat, T = 1) {
  z_scaled <- z_mat / T
  e <- exp(z_scaled - apply(z_scaled, 1, max))
  e / rowSums(e)
}

log_softmax <- function(z_mat, T = 1) {
  z_scaled <- z_mat / T
  log_sum_exp <- log(rowSums(exp(z_scaled - apply(z_scaled, 1, max)))) +
    apply(z_scaled, 1, max)
  z_scaled - log_sum_exp
}

nll_binary_temp <- function(T, logits, y) {
  p <- sigmoid(logits / T)
  eps <- 1e-7
  p <- pmax(pmin(p, 1 - eps), eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

nll_binary_platt <- function(params, logits, y) {
  a <- params[1]
  b <- params[2]
  p <- sigmoid(a * logits + b)
  eps <- 1e-7
  p <- pmax(pmin(p, 1 - eps), eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

nll_multi_temp <- function(T, z_mat, y_mat) {
  lp <- log_softmax(z_mat, T)
  -mean(rowSums(y_mat * lp))
}

ece <- function(p, y, n_bins = 10) {
  bins <- cut(p, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
  bin_stats <- data.frame(p = p, y = y, bin = bins) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_p = mean(p),
      mean_y = mean(y),
      .groups = "drop"
    )
  sum(bin_stats$n / length(p) * abs(bin_stats$mean_p - bin_stats$mean_y))
}

ece_confidence <- function(p_mat, y_mat, n_bins = 10) {
  conf <- apply(p_mat, 1, max)
  correct <- (apply(p_mat, 1, which.max) == apply(y_mat, 1, which.max)) * 1
  ece(conf, correct, n_bins)
}

reliability_data <- function(p, y, n_bins = 10) {
  bins <- cut(p, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
  data.frame(p = p, y = y, bin = bins) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_pred = mean(p),
      obs_freq = mean(y),
      .groups = "drop"
    )
}

parse_prob_vec <- function(val_col) {
  do.call(rbind, lapply(strsplit(as.character(val_col), ",\\s*"), as.numeric))
}
