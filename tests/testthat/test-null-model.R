test_that("GLM covariates distinguish continuous and categorical columns", {
  pheno <- data.frame(
    age = c(40, 50, 60, 70),
    centre = c("A", "B", "C", "A"),
    stringsAsFactors = FALSE
  )
  design <- KnockoffPipeline:::.build_glm_covariates(
    pheno, covar_cols = "age", cat_covar_cols = "centre"
  )

  expect_true("age" %in% colnames(design))
  expect_true(all(c("centreB", "centreC") %in% colnames(design)))
  expect_false("(Intercept)" %in% colnames(design))
  expect_equal(nrow(design), nrow(pheno))
})


test_that("SAIGE receives categorical columns in both required arguments", {
  args <- KnockoffPipeline:::.saige_covariate_args(
    covar_cols = "age", cat_covar_cols = c("centre", "array")
  )

  expect_identical(args$covarColList, c("age", "centre", "array"))
  expect_identical(args$qCovarCol, c("centre", "array"))
})


test_that("GLM covariate validation rejects ambiguous or invalid input", {
  pheno <- data.frame(age = c(40, 50), centre = c("A", "A"))
  expect_error(
    KnockoffPipeline:::.build_glm_covariates(
      pheno, covar_cols = "age", cat_covar_cols = "age"
    ),
    "both"
  )
  expect_error(
    KnockoffPipeline:::.build_glm_covariates(
      pheno, cat_covar_cols = "centre"
    ),
    "fewer than two"
  )
})


test_that("standard GLM null model drops linearly dependent covariate directions", {
  y <- c(1, 2, 3, 5, 8)
  x <- cbind(a = 1:5, duplicate = 2 * (1:5))
  fit <- KnockoffPipeline:::Fit_null_model(y, x, id = letters[1:5])

  expect_identical(as.character(fit$id), letters[1:5])
  expect_equal(ncol(fit$X0), 2L)
  expect_true(all(is.finite(fit$mu)))
})
