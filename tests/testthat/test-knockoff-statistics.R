reference_mk_path <- function(kappa, tau, M) {
  eligible <- !is.na(kappa) & !is.na(tau) & tau > 0
  thresholds <- sort(unique(tau[eligible]), decreasing = TRUE)
  fdp <- vapply(thresholds, function(threshold) {
    included <- eligible & tau >= threshold
    n_original <- sum(kappa[included] == 0)
    n_knockoff <- sum(kappa[included] != 0)
    (1 + n_knockoff) / (M * max(1, n_original))
  }, numeric(1))
  list(threshold = thresholds, fdp = fdp)
}


reference_mk_q <- function(kappa, tau, M) {
  path <- reference_mk_path(kappa, tau, M)
  q <- rep(1, length(tau))
  original <- which(!is.na(kappa) & !is.na(tau) & kappa == 0 & tau > 0)
  for (i in original) {
    q[i] <- min(1, min(path$fdp[path$threshold <= tau[i]]))
  }
  q
}


reference_mk_threshold <- function(kappa, tau, M, fdr = 0.1) {
  path <- reference_mk_path(kappa, tau, M)
  accepted <- which(path$fdp <= fdr)
  if (length(accepted)) path$threshold[max(accepted)] else Inf
}


test_that("multiple-knockoff q values use complete positive-tau groups", {
  fixtures <- list(
    list(kappa = c(0, 0, 1, 0), tau = c(4, 3, 2, 1), M = 5L, fdr = 0.1),
    list(kappa = c(1, 0), tau = c(1, 1), M = 5L, fdr = 0.3),
    list(kappa = c(1, 1, 1, 1, 1, 0), tau = 6:1, M = 5L, fdr = 0.1),
    list(kappa = c(0, 1, 0), tau = c(3, 2, 1), M = 2L, fdr = 0.4),
    list(kappa = c(0, 1, 0), tau = c(0, 0, NA), M = 5L, fdr = 0.1)
  )

  for (x in fixtures) {
    expect_equal(
      KnockoffPipeline:::MK.q.byStat(x$kappa, x$tau, x$M),
      reference_mk_q(x$kappa, x$tau, x$M),
      tolerance = 0
    )
    expect_equal(
      KnockoffPipeline:::MK.threshold.byStat(
        x$kappa, x$tau, x$M, fdr = x$fdr
      ),
      reference_mk_threshold(x$kappa, x$tau, x$M, fdr = x$fdr),
      tolerance = 0
    )
  }

  # A tie is one threshold set, so its result cannot depend on row order.
  kappa <- c(1, 0, 1, 0)
  tau <- c(2, 2, 1, 1)
  permutation <- c(2, 1, 4, 3)
  q <- KnockoffPipeline:::MK.q.byStat(kappa, tau, M = 5L)
  q_permuted <- KnockoffPipeline:::MK.q.byStat(
    kappa[permutation], tau[permutation], M = 5L
  )
  expect_equal(q[permutation], q_permuted, tolerance = 0)

  # Raw FDP estimates can exceed one; a q value cannot.
  expect_equal(
    KnockoffPipeline:::MK.q.byStat(
      c(1, 1, 1, 1, 1, 0), 6:1, M = 5L
    )[6],
    1,
    tolerance = 0
  )
})


test_that("Rej.Bound does not split an equal-tau threshold group", {
  path <- KnockoffPipeline:::.MK.fdp.path(
    kappa = c(0, 1, 0, 1), tau = c(3, 2, 2, 1),
    M = 5L, Rej.Bound = 2L
  )
  expect_equal(path$index, 1:3)
  expect_equal(path$threshold, c(3, 2))
  expect_equal(path$fdp, c(0.2, 0.2), tolerance = 0)
})


test_that("multiple-knockoff summaries reject invalid inputs", {
  expect_error(
    KnockoffPipeline:::MK.q.byStat(c(0, 1), 1, M = 5L),
    "same length"
  )
  expect_error(
    KnockoffPipeline:::MK.q.byStat(0, 1, M = 0L),
    "positive integer"
  )
  expect_error(
    KnockoffPipeline:::MK.q.byStat(6, 1, M = 5L),
    "between 0 and M"
  )
  expect_error(
    KnockoffPipeline:::MK.q.byStat(0, -Inf, M = 5L),
    "negative infinity"
  )
  expect_error(
    KnockoffPipeline:::MK.threshold.byStat(0, 1, M = 5L, fdr = 1.1),
    "between 0 and 1"
  )
})


test_that("GeneScan3DKnock uses conservative ties and shared q formulas", {
  M <- 2L
  T <- rbind(
    c(4, 4, 0), # original/knockoff tie: original must not win
    c(3, 0, 0),
    c(0, 2, 0),
    c(1, 0, 0)
  )
  p <- 10^(-T)
  stat <- KnockoffPipeline:::MK.statistic(
    T[, 1L], T[, -1L, drop = FALSE], method = "median"
  )

  got <- KnockoffPipeline:::GeneScan3DKnock(
    M = M, p0 = p[, 1L], p_ko = p[, -1L, drop = FALSE],
    fdr = 0.4, gene_id = paste0("g", seq_len(nrow(p)))
  )
  expected_W <- (T[, 1L] - apply(T[, -1L, drop = FALSE], 1L, median)) *
    (stat[, "kappa"] == 0)

  expect_equal(unname(stat[1L, "kappa"]), 1)
  expect_equal(got$W[1L], 0)
  expect_equal(got$W, expected_W, tolerance = 0)
  expect_equal(
    got$Qvalue,
    reference_mk_q(stat[, "kappa"], stat[, "tau"], M),
    tolerance = 0
  )
  expect_equal(
    got$W.threshold,
    reference_mk_threshold(stat[, "kappa"], stat[, "tau"], M, fdr = 0.4),
    tolerance = 0
  )
})


test_that("missing tests are non-discoveries and zero p-values stay finite", {
  missing <- KnockoffPipeline:::GeneScan3DKnock(
    M = 2L, p0 = 1e-8, p_ko = matrix(c(0.5, NA), nrow = 1L),
    gene_id = "missing"
  )
  expect_equal(missing$W, 0)
  expect_equal(missing$Qvalue, 1)
  expect_equal(missing$W.threshold, Inf)

  underflow <- KnockoffPipeline:::GeneScan3DKnock(
    M = 2L, p0 = 0, p_ko = matrix(c(0.5, 0.5), nrow = 1L),
    fdr = 0.5, gene_id = "underflow"
  )
  expect_true(is.finite(underflow$W))
  expect_true(is.finite(underflow$W.threshold))
  expect_equal(underflow$Qvalue, 0.5)
  expect_true(is.finite(KnockoffPipeline:::MK.threshold.byStat(
    kappa = 0, tau = Inf, M = 2L, fdr = 0.5
  )))
})


test_that("KS_summary uses the supplied knockoff count", {
  columns <- c(
    "chr", "start", "end", "actual_start", "actual_end", "signal",
    "kappa", "tau", "W", "P_KS", "P_KS_k1", "P_KS_k2"
  )
  window <- matrix(
    c(1, 1, 1, 1, 1, 0, 0, 2, 2, 0.01, 0.5, 0.6),
    nrow = 1L, dimnames = list(NULL, columns)
  )
  single <- matrix(
    c(1, 2, 2, 2, 2, 0, 1, 1, 0, 0.5, 0.01, 0.6),
    nrow = 1L, dimnames = list(NULL, columns)
  )

  got <- KnockoffPipeline:::KS_summary(window, single, M = 2L, fdr = 0.4)
  expect_equal(
    got$Qvalue,
    reference_mk_q(c(0, 1), c(2, 1), M = 2L),
    tolerance = 0
  )
  expect_equal(
    unique(got$W_Threshold),
    reference_mk_threshold(c(0, 1), c(2, 1), M = 2L, fdr = 0.4),
    tolerance = 0
  )
})


test_that("GeneScan3D supports a one-variant binary enhancer", {
  set.seed(910)
  n <- 80L
  null <- KnockoffPipeline:::Fit_null_model(
    rep(c(0, 1), length.out = n), out_type = "D"
  )
  gene <- matrix(stats::rbinom(n * 4L, 2L, 0.3), nrow = n)
  enhancer <- matrix(stats::rbinom(n, 2L, 0.3), nrow = n, ncol = 1L)

  testthat::local_mocked_bindings(
    Get.p = function(X, result.null.model) matrix(0.5, ncol = ncol(X)),
    .package = "KnockoffPipeline"
  )
  got <- suppressWarnings(KnockoffPipeline:::GeneScan3D(
    G = gene, Z = NULL, G.promoter = NULL, Z.promoter = NULL,
    G.EnhancerAll = enhancer,
    Z.EnhancerAll = matrix(c(0.25, 0.75), nrow = 1L),
    R = 1L, p_Enhancer = 1L, window.size = 1000,
    pos = c(100, 300, 600, 900), Gsub.id = seq_len(n),
    result.null.model = null
  ))

  expect_length(got$GeneScan3D.Cauchy.pvalue, 3L)
  expect_true(is.finite(got$GeneScan3D.Cauchy.pvalue[1L]))
})
