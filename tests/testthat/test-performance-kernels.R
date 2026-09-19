test_that("sparse covariance helpers match dense reference formulas", {
  set.seed(17)
  dense <- matrix(rnorm(300), nrow = 50)
  dense[abs(dense) < 0.8] <- 0
  sparse <- Matrix::Matrix(dense, sparse = TRUE)

  observed <- KnockoffPipeline:::.kp_sparse_cov_cor(sparse)
  expect_equal(observed$cov, stats::cov(dense), tolerance = 1e-12)
  expect_equal(observed$cor, stats::cor(dense), tolerance = 1e-12)

  y_dense <- dense[, 1:2, drop = FALSE]
  y_sparse <- Matrix::Matrix(y_dense, sparse = TRUE)
  expect_equal(
    KnockoffPipeline:::.kp_sparse_cross_cov(sparse, y_sparse),
    stats::cov(dense, y_dense), tolerance = 1e-12
  )
})


test_that("imputation preserves an already complete sparse matrix", {
  x <- Matrix::Matrix(
    matrix(c(0, 1, 2, 0, 0, 1, 0, 2), nrow = 4L), sparse = TRUE
  )
  observed <- KnockoffPipeline:::Impute(x, "fixed")
  expect_s4_class(observed, "Matrix")
  expect_identical(observed, x)
})


test_that("sparse missing-value imputation matches column means", {
  x <- Matrix::Matrix(
    matrix(c(0, NA, 2, 0, 1, 2, NA, 1), nrow = 4L), sparse = TRUE
  )
  observed <- KnockoffPipeline:::Impute(x, "fixed")
  expected <- matrix(c(0, 2 / 3, 2, 0, 1, 2, 4 / 3, 1), nrow = 4L)
  expect_equal(observed, expected)
})


test_that("Single knockoff construction reuses precomputed correlations", {
  set.seed(101)
  x <- matrix(rbinom(120L * 4L, 2L, 0.25), nrow = 120L)
  x <- Matrix::Matrix(x, sparse = TRUE)
  pos <- c(100, 1000, 2500, 5000)
  cor_x <- KnockoffPipeline:::.kp_sparse_cov_cor(
    x, need_cov = FALSE, need_cor = TRUE
  )$cor
  clusters <- stats::cutree(
    stats::hclust(stats::as.dist(1 - abs(cor_x)), method = "single"),
    h = 0.25
  )
  expect_identical(unname(clusters), seq_len(ncol(x)))

  set.seed(73)
  expected <- KnockoffPipeline:::create.KS(
    x, pos, M = 2L, n.AL = 60L, thres.ultrarare = 0,
    method = "uniform", bigmemory = FALSE
  )
  set.seed(73)
  observed <- KnockoffPipeline:::create.KS(
    x, pos, M = 2L, n.AL = 60L, thres.ultrarare = 0,
    method = "uniform", bigmemory = FALSE,
    cor.X.precomputed = cor_x, preclustered = TRUE
  )
  expect_equal(observed, expected, tolerance = 0)
})


test_that("continuous score helper preserves the original formula", {
  set.seed(29)
  n <- 80L
  x <- matrix(rbinom(n * 9L, 2L, 0.25), nrow = n)
  covar <- matrix(rnorm(n), ncol = 1L)
  y <- 0.4 * covar[, 1L] + rnorm(n)
  null <- KnockoffPipeline:::Fit_null_model(y, covar, out_type = "C")

  residual <- null$Y - null$nullglm$fitted.values
  v <- rep(as.numeric(stats::var(residual)), n)
  score_sq <- (t(x) %*% residual)^2
  raw_variance <- colSums(v * x^2)
  left <- t(x) %*% (v * null$X0) %*% null$inv.X0
  right <- t(t(null$X0) %*% as.matrix(v * x))
  expected <- stats::pchisq(
    as.numeric(score_sq / (raw_variance - rowSums(left * right))),
    df = 1, lower.tail = FALSE
  )

  expect_equal(
    KnockoffPipeline:::.kp_continuous_score_p(
      Matrix::Matrix(x, sparse = TRUE), null
    ),
    expected, tolerance = 1e-12
  )
})


test_that("binary single-variant scores use the fitted null-model design", {
  set.seed(30)
  n <- 100L
  x <- matrix(rbinom(n * 4L, 2L, 0.3), nrow = n)
  y <- rbinom(n, 1L, 0.4)
  null <- KnockoffPipeline:::Fit_null_model(y, out_type = "D")

  observed <- KnockoffPipeline:::Get.p(
    Matrix::Matrix(x, sparse = TRUE), null
  )

  expect_equal(dim(observed), c(ncol(x), 1L))
  expect_true(all(is.finite(observed)))
  expect_true(all(observed >= 0 & observed <= 1))
})


test_that("chunked single-variant scores preserve column order", {
  set.seed(31)
  n <- 60L
  x <- matrix(rbinom(n * 11L, 2L, 0.3), nrow = n)
  y <- rnorm(n)
  null <- KnockoffPipeline:::Fit_null_model(y, out_type = "C")
  sparse <- Matrix::Matrix(x, sparse = TRUE)

  expected <- KnockoffPipeline:::Get.p(sparse, null)
  observed <- KnockoffPipeline:::.get_p_in_chunks(
    sparse, null, chunk_size = 3L
  )
  expect_equal(observed, expected, tolerance = 1e-12)

  index <- c(11L, 2L, 7L, 1L)
  expect_equal(
    KnockoffPipeline:::.get_p_in_chunks(
      sparse, null, column_index = index, chunk_size = 2L
    ),
    KnockoffPipeline:::Get.p(sparse[, index, drop = FALSE], null),
    tolerance = 1e-12
  )
})


test_that("cached SKAT eigenvalues preserve the direct calculation", {
  score <- matrix(c(0.2, -0.1, 0.4, 0.3), ncol = 1L)
  K <- crossprod(matrix(c(
    1.0, 0.2, 0.1, 0.0,
    0.2, 1.1, 0.3, 0.1,
    0.1, 0.3, 0.9, 0.2,
    0.0, 0.1, 0.2, 1.2
  ), nrow = 4L))
  windows <- cbind(
    first = c(1, 1, 0, 0),
    second = c(0, 1, 1, 1)
  )
  weight <- c(0.7, 1.0, 1.2, 0.8)
  q_stat <- as.vector(t(score^2) %*% (weight * windows)^2)
  weighted_k <- weight * t(weight * K)
  expected <- vapply(seq_along(q_stat), function(i) {
    member <- windows[, i] != 0
    lambda <- eigen(
      weighted_k[member, member, drop = FALSE],
      symmetric = TRUE, only.values = TRUE
    )$values
    p <- suppressWarnings(
      CompQuadForm::davies(q_stat[i], lambda, acc = 1e-6)$Qq
    )
    if (p > 1 || p <= 0) {
      p <- KnockoffPipeline:::Get_Liu_PVal.MOD.Lambda(q_stat[i], lambda)
    }
    p
  }, numeric(1))

  prepared <- KnockoffPipeline:::.prepare_skat_ks(K, windows, weight)
  observed <- suppressWarnings(
    KnockoffPipeline:::.get_p_skat_ks_prepared(
      score, windows, weight, prepared
    )
  )

  expect_equal(as.numeric(observed), expected, tolerance = 1e-12)
})


test_that("Single preprocessing applies MAC filtering to every aligned field", {
  geno <- cbind(
    common_a = c(0L, 1L, 2L, 0L, 1L, 2L),
    low_mac  = c(0L, 0L, 0L, 0L, 0L, 1L),
    common_b = c(2L, 1L, 0L, 2L, 1L, 0L)
  )
  metadata <- data.frame(
    chr = 1L,
    variant_id = c("a", "rare", "b"),
    pos = c(100, 200, 300),
    a1 = "A", a2 = "G", counted_allele = "A",
    stringsAsFactors = FALSE
  )

  # common_a/common_b are perfectly anticorrelated, so add a fourth independent
  # common variant to leave at least two columns after the LD representative
  # filter.
  geno <- cbind(geno, common_c = c(0L, 2L, 0L, 2L, 1L, 1L))
  metadata <- rbind(
    metadata,
    data.frame(chr = 1L, variant_id = "c", pos = 400,
               a1 = "A", a2 = "G", counted_allele = "A")
  )

  observed <- KnockoffPipeline:::Preprocess(
    geno, chr = 1L, window = 100L, variant_metadata = metadata,
    thres.ultrarare = 2
  )
  expect_false("rare" %in% observed$variant_metadata$variant_id)
  expect_equal(ncol(observed$G), nrow(observed$variant_metadata))
  expect_equal(as.numeric(colnames(observed$G)), observed$pos)
})
