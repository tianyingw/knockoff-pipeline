# Internal matrix helpers -----------------------------------------------------
#
# The upstream method files historically each defined their own sparse.cor()
# implementation.  Those implementations immediately converted a sparse
# n-by-p genotype matrix to a dense matrix, which defeats the purpose of using
# Matrix objects and can add several gigabytes per worker.  Keep the numerical
# formula in one place and materialise only the p-by-p result.

.kp_col_means <- function(x) {
  if (inherits(x, "Matrix")) Matrix::colMeans(x) else base::colMeans(x)
}

.kp_crossprod <- function(x, y = NULL) {
  if (inherits(x, "Matrix") || (!is.null(y) && inherits(y, "Matrix"))) {
    if (is.null(y)) Matrix::crossprod(x) else Matrix::crossprod(x, y)
  } else {
    if (is.null(y)) base::crossprod(x) else base::crossprod(x, y)
  }
}

.kp_col_sums <- function(x) {
  if (inherits(x, "Matrix")) Matrix::colSums(x) else base::colSums(x)
}

.kp_continuous_score_p <- function(x, result.prelim) {
  residual <- as.numeric(result.prelim$res)
  if (length(residual) != nrow(x))
    residual <- as.numeric(result.prelim$Y - result.prelim$nullglm$fitted.values)
  v0 <- as.numeric(result.prelim$v[1L])
  if (length(v0) != 1L || !is.finite(v0))
    v0 <- as.numeric(stats::var(residual))
  score <- as.numeric(.kp_crossprod(x, residual))
  raw_variance <- v0 * .kp_col_sums(x * x)
  predictor_cross <- v0 * as.matrix(.kp_crossprod(x, result.prelim$X0))
  projected <- rowSums(
    (predictor_cross %*% result.prelim$inv.X0) * predictor_cross
  )
  stats::pchisq(
    score^2 / (raw_variance - projected), df = 1, lower.tail = FALSE
  )
}

# Canonical genotype imputation used by all analysis modes.  The historical
# method files each carried an identical copy; keeping one implementation
# avoids load-order-dependent overrides.  Extract each affected column once,
# rather than repeatedly slicing the full matrix inside the missing-value loop.
Impute <- function(Z, impute.method) {
  methods <- c("random", "fixed", "bestguess")
  if (!is.character(impute.method) || length(impute.method) != 1L ||
      is.na(impute.method) || !impute.method %in% methods) {
    stop(
      "Error: Imputation method should be \"fixed\", \"random\" or \"bestguess\" "
    )
  }
  if (!anyNA(Z)) return(Z)
  Z <- as.matrix(Z)

  for (i in seq_len(ncol(Z))) {
    values <- Z[, i]
    missing <- which(is.na(values))
    if (length(missing) == 0L) next
    maf1 <- mean(values, na.rm = TRUE) / 2
    values[missing] <- switch(
      impute.method,
      random = stats::rbinom(length(missing), 2, maf1),
      fixed = 2 * maf1,
      bestguess = round(2 * maf1)
    )
    Z[, i] <- values
  }
  Z
}

.kp_sparse_cov_cor <- function(x, need_cov = TRUE, need_cor = TRUE) {
  if (!isTRUE(need_cov) && !isTRUE(need_cor)) return(list())
  n <- nrow(x)
  if (n < 2L) stop("At least two rows are required to compute covariance.")

  means <- .kp_col_means(x)
  covmat <- (as.matrix(.kp_crossprod(x)) - n * tcrossprod(means)) / (n - 1)
  out <- list()
  if (isTRUE(need_cov)) out$cov <- covmat
  if (isTRUE(need_cor)) {
    sdvec <- sqrt(diag(covmat))
    out$cor <- covmat / tcrossprod(sdvec)
  }
  out
}

.kp_sparse_cross_cov <- function(x, y) {
  if (nrow(x) != nrow(y))
    stop("Cross-covariance matrices must have the same number of rows.")
  n <- nrow(x)
  if (n < 2L) stop("At least two rows are required to compute covariance.")

  means_x <- .kp_col_means(x)
  means_y <- .kp_col_means(y)
  (as.matrix(.kp_crossprod(x, y)) - n * tcrossprod(means_x, means_y)) /
    (n - 1)
}

# Keep the historical internal names because the imported upstream method
# implementations call them.  Every definition now delegates to the same
# sparse-safe implementation, so file collation order no longer changes the
# memory behaviour.
sparse.cor <- function(x) {
  .kp_sparse_cov_cor(x, need_cov = TRUE, need_cor = TRUE)
}

sparse.cov.cross <- function(x, y) {
  list(cov = .kp_sparse_cross_cov(x, y))
}
