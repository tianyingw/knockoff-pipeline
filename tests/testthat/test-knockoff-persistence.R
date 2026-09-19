test_that("persistence plan honors full-run save requests", {
  transient <- KnockoffPipeline:::.resolve_knockoff_persistence(
    "full", save_knockoff = NULL, multi_pheno = FALSE
  )
  expect_false(transient$save_knockoff)
  expect_false(transient$write_first_pass)
  expect_false(transient$load_first_pass)
  expect_false(transient$cleanup_after_run)

  retained <- KnockoffPipeline:::.resolve_knockoff_persistence(
    "full", save_knockoff = TRUE, multi_pheno = FALSE
  )
  expect_true(retained$save_knockoff)
  expect_true(retained$write_first_pass)
  expect_false(retained$load_first_pass)
  expect_false(retained$cleanup_after_run)
})


test_that("multi-phenotype full runs save once and reuse later", {
  plan <- KnockoffPipeline:::.resolve_knockoff_persistence(
    "full", save_knockoff = FALSE, multi_pheno = TRUE
  )
  expect_true(plan$write_first_pass)
  expect_false(plan$load_first_pass)
  expect_true(plan$reuse_later_passes)
  expect_true(plan$cleanup_after_run)

  retained <- KnockoffPipeline:::.resolve_knockoff_persistence(
    "full", save_knockoff = TRUE, multi_pheno = TRUE
  )
  expect_true(retained$write_first_pass)
  expect_true(retained$reuse_later_passes)
  expect_false(retained$cleanup_after_run)

  first <- KnockoffPipeline:::.knockoff_pass_flags(plan, FALSE)
  later <- KnockoffPipeline:::.knockoff_pass_flags(plan, TRUE)
  expect_true(first$save_knockoff)
  expect_false(first$load_knockoff)
  expect_false(later$save_knockoff)
  expect_true(later$load_knockoff)
})


test_that("stage 1 persists and stage 2 always loads without cleanup", {
  stage1 <- KnockoffPipeline:::.resolve_knockoff_persistence(
    "stage1_knockoff", save_knockoff = NULL, multi_pheno = FALSE
  )
  expect_true(stage1$save_knockoff)
  expect_true(stage1$write_first_pass)
  expect_false(stage1$load_first_pass)
  expect_error(
    KnockoffPipeline:::.resolve_knockoff_persistence(
      "stage1_knockoff", save_knockoff = FALSE, multi_pheno = FALSE
    ),
    "requires save_knockoff"
  )

  for (multi in c(FALSE, TRUE)) {
    stage2 <- KnockoffPipeline:::.resolve_knockoff_persistence(
      "stage2_analysis", save_knockoff = NULL, multi_pheno = multi
    )
    expect_false(stage2$write_first_pass)
    expect_true(stage2$load_first_pass)
    expect_false(stage2$reuse_later_passes)
    expect_false(stage2$cleanup_after_run)
    flags <- KnockoffPipeline:::.knockoff_pass_flags(stage2, TRUE)
    expect_false(flags$save_knockoff)
    expect_true(flags$load_knockoff)
  }
})


test_that("stage 1 requires phenotype inputs to define complete cases", {
  expect_error(
    run_pipeline(
      outdir = tempfile("stage1-missing-phenotype-"),
      test_type = "Single_Window",
      geno_file = tempfile("missing-plink-"),
      pipeline_stage = "stage1_knockoff"
    ),
    "pheno_file.*required for every pipeline stage"
  )
})


test_that("directory preparation distinguishes new and existing paths", {
  new_dir <- tempfile("knockoff-new-")
  missing_dir <- tempfile("knockoff-missing-")
  on.exit(unlink(c(new_dir, missing_dir), recursive = TRUE, force = TRUE),
          add = TRUE)

  expect_true(KnockoffPipeline:::.prepare_knockoff_directory(
    new_dir, create = TRUE
  ))
  expect_true(dir.exists(new_dir))
  expect_false(KnockoffPipeline:::.prepare_knockoff_directory(
    new_dir, create = TRUE
  ))
  expect_false(KnockoffPipeline:::.prepare_knockoff_directory(
    new_dir, require_existing = TRUE
  ))
  expect_error(
    KnockoffPipeline:::.prepare_knockoff_directory(
      missing_dir, require_existing = TRUE
    ),
    "Knockoff directory not found"
  )
})


test_that("cleanup removes only directories owned by the current full run", {
  owned_dir <- tempfile("knockoff-owned-")
  user_dir <- tempfile("knockoff-user-")
  stage2_dir <- tempfile("knockoff-stage2-")
  dir.create(owned_dir)
  dir.create(user_dir)
  dir.create(stage2_dir)
  writeLines("keep", file.path(user_dir, "sentinel"))
  writeLines("keep", file.path(stage2_dir, "sentinel"))
  on.exit(unlink(c(owned_dir, user_dir, stage2_dir), recursive = TRUE,
                 force = TRUE), add = TRUE)

  expect_true(KnockoffPipeline:::.remove_owned_knockoff_directory(
    owned_dir, cleanup = TRUE, created_by_run = TRUE
  ))
  expect_false(dir.exists(owned_dir))

  expect_false(KnockoffPipeline:::.remove_owned_knockoff_directory(
    user_dir, cleanup = TRUE, created_by_run = FALSE
  ))
  expect_true(file.exists(file.path(user_dir, "sentinel")))

  stage2 <- KnockoffPipeline:::.resolve_knockoff_persistence(
    "stage2_analysis", save_knockoff = FALSE, multi_pheno = TRUE
  )
  expect_false(KnockoffPipeline:::.remove_owned_knockoff_directory(
    stage2_dir,
    cleanup = stage2$cleanup_after_run,
    created_by_run = TRUE
  ))
  expect_true(file.exists(file.path(stage2_dir, "sentinel")))
})


test_that("random genotype imputation requires an explicit base seed", {
  expect_error(
    run_pipeline(
      outdir = tempfile("random-imputation-"),
      test_type = "Single_Window",
      geno_file = tempfile("missing-plink-"),
      geno_missing_imputation = "random"
    ),
    "seed.*required"
  )
})


test_that("gene-centric mode rejects unsupported imputation choices", {
  expect_error(
    run_pipeline(
      outdir = tempfile("gene-imputation-"),
      test_type = "Gene_Centric",
      geno_file = tempfile("missing-plink-"),
      geno_missing_imputation = "bestguess"
    ),
    "supports geno_missing_imputation = 'fixed' only"
  )
})


test_that("logical control arguments reject missing values", {
  expect_error(
    run_pipeline(
      outdir = tempfile("logical-validation-"),
      test_type = "Single_Window",
      geno_file = tempfile("missing-plink-"),
      sample_uncorrelated = NA
    ),
    "sample_uncorrelated must be TRUE or FALSE"
  )
  expect_error(
    run_pipeline(
      outdir = tempfile("logical-validation-"),
      test_type = "Single_Window",
      geno_file = tempfile("missing-plink-"),
      read_mid_exist = NA
    ),
    "read_mid_exist must be TRUE or FALSE"
  )
})
