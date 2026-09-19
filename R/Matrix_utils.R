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

.kp_weighted_col_sums_sq <- function(x, weights) {
  if (length(weights) != nrow(x))
    stop("The weight vector must contain one value per matrix row.")
  .kp_col_sums(x * (weights * x))
}

.kp_continuous_score_p <- function(x, result.prelim) {
  residual <- as.numeric(result.prelim$Y - result.prelim$nullglm$fitted.values)
  v <- rep(as.numeric(stats::var(residual)), nrow(x))
  score <- as.numeric(.kp_crossprod(x, residual))
  raw_variance <- .kp_weighted_col_sums_sq(x, v)
  predictor_cross <- as.matrix(
    .kp_crossprod(x, v * result.prelim$X0)
  )
  projected <- rowSums(
    (predictor_cross %*% result.prelim$inv.X0) * predictor_cross
  )
  stats::pchisq(
    score^2 / (raw_variance - projected), df = 1, lower.tail = FALSE
  )
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
