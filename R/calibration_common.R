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

softmax_rows <- function(z_mat) {
  z_shift <- z_mat - apply(z_mat, 1, max)
  ez <- exp(z_shift)
  ez / rowSums(ez)
}

fit_dirichlet_calibration <- function(p_mat, y_mat, lambda = 1e-4, maxit = 300) {
  eps <- 1e-9
  p <- pmax(pmin(p_mat, 1 - eps), eps)
  p <- p / rowSums(p)
  x <- log(p)
  k <- ncol(p)
  y_idx <- apply(y_mat, 1, which.max)

  if (nrow(p) < k || length(unique(y_idx)) < 2) {
    return(list(W = diag(k), b = rep(0, k), converged = FALSE, fallback = TRUE))
  }

  unpack_theta <- function(theta) {
    W <- matrix(theta[seq_len(k * k)], nrow = k, ncol = k, byrow = TRUE)
    b <- theta[(k * k + 1):(k * k + k)]
    list(W = W, b = b)
  }

  objective <- function(theta) {
    u <- unpack_theta(theta)
    z <- x %*% t(u$W) + matrix(u$b, nrow = nrow(x), ncol = k, byrow = TRUE)
    p_hat <- softmax_rows(z)
    nll <- -mean(log(p_hat[cbind(seq_len(nrow(p_hat)), y_idx)]))
    nll + lambda * mean(u$W^2)
  }

  theta0 <- c(as.vector(t(diag(k))), rep(0, k))
  opt <- tryCatch(
    optim(theta0, objective, method = "BFGS", control = list(maxit = maxit, reltol = 1e-9)),
    error = function(e) NULL
  )

  if (is.null(opt) || is.null(opt$par) || !is.finite(opt$value)) {
    return(list(W = diag(k), b = rep(0, k), converged = FALSE, fallback = TRUE))
  }

  u <- unpack_theta(opt$par)
  list(
    W = u$W,
    b = u$b,
    converged = isTRUE(opt$convergence == 0),
    fallback = FALSE
  )
}

predict_dirichlet_calibration <- function(fit, p_mat) {
  eps <- 1e-9
  p <- pmax(pmin(p_mat, 1 - eps), eps)
  p <- p / rowSums(p)
  x <- log(p)
  z <- x %*% t(fit$W) + matrix(fit$b, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  softmax_rows(z)
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
