test_that("analysis checkpoint is created and an identical run can resume", {
  mid <- tempfile("checkpoint-")
  on.exit(unlink(mid, recursive = TRUE, force = TRUE), add = TRUE)
  context <- list(schema_version = 1L, phenotype = "Y", M = 5L)

  expect_identical(
    KnockoffPipeline:::.prepare_analysis_checkpoint(mid, context, TRUE),
    "new"
  )
  expect_true(file.exists(file.path(mid, "checkpoint_context.rds")))
  writeLines("done", file.path(mid, "done_0001.txt"))
  expect_identical(
    KnockoffPipeline:::.prepare_analysis_checkpoint(mid, context, TRUE),
    "resume"
  )
})


test_that("generated analysis contexts use the current algorithm schema", {
  context <- KnockoffPipeline:::.make_analysis_checkpoint_context(
    test_type = "Single_Window", phenotype = "Y", sample_ids = "sample-1",
    pheno_id = "IID", covar_cols = character(),
    cat_covar_cols = character(), out_type = "C", M = 5L, seed = 17L,
    genome_build = "hg19", sliding_window_length = c(1000, 5000),
    geno_missing_imputation = "fixed", sample_uncorrelated = TRUE,
    relatedness_cutoff = 0.125, n_markers_grm = 1000L, fdr = 0.1,
    sources = list()
  )

  expect_identical(context$schema_version, 2L)
})


test_that("incompatible and legacy checkpoints fail closed", {
  incompatible <- tempfile("checkpoint-incompatible-")
  legacy <- tempfile("checkpoint-legacy-")
  dir.create(incompatible)
  dir.create(legacy)
  on.exit(unlink(c(incompatible, legacy), recursive = TRUE, force = TRUE),
          add = TRUE)

  saveRDS(
    list(schema_version = 1L, phenotype = "Y", M = 5L),
    file.path(incompatible, "checkpoint_context.rds")
  )
  expect_error(
    KnockoffPipeline:::.prepare_analysis_checkpoint(
      incompatible,
      list(schema_version = 1L, phenotype = "Y", M = 3L),
      TRUE
    ),
    "mismatch: M"
  )

  writeLines("legacy", file.path(legacy, "done_0001.txt"))
  expect_error(
    KnockoffPipeline:::.prepare_analysis_checkpoint(
      legacy, list(schema_version = 1L), TRUE
    ),
    "lack compatibility metadata"
  )
})


test_that("read_mid_exist FALSE performs a genuinely fresh mid-directory run", {
  mid <- tempfile("checkpoint-fresh-")
  dir.create(file.path(mid, "blocks"), recursive = TRUE)
  writeLines("old", file.path(mid, "blocks", "old-result.txt"))
  writeLines("user-visible sentinel", file.path(mid, "old-chromosome.txt"))
  context <- list(schema_version = 1L, phenotype = "new")
  on.exit(unlink(mid, recursive = TRUE, force = TRUE), add = TRUE)

  expect_identical(
    KnockoffPipeline:::.prepare_analysis_checkpoint(mid, context, FALSE),
    "fresh"
  )
  expect_false(file.exists(file.path(mid, "blocks", "old-result.txt")))
  expect_false(file.exists(file.path(mid, "old-chromosome.txt")))
  expect_identical(
    readRDS(file.path(mid, "checkpoint_context.rds")), context
  )
})
