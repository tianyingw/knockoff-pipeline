test_that("SAIGE PLINK preparation combines sample filtering and marker thinning", {
  plan <- KnockoffPipeline:::.saige_plink_export_plan(
    plink_prefix = "/path with spaces/plink2",
    plink_file = "/data/common cohort",
    out_prefix = "/tmp/grm input",
    total_markers = 10000L,
    thin_target_markers = 5000L,
    plink_keep_file = "/tmp/analysis samples.fam",
    random_seed = 42L
  )

  expect_true(plan$thin_markers)
  expect_true(plan$subset_samples)
  expect_match(plan$command, "--keep", fixed = TRUE)
  expect_match(plan$command, "/tmp/analysis samples.fam", fixed = TRUE)
  expect_match(plan$command, "--thin 0.5", fixed = TRUE)
  expect_match(plan$command, "--seed 42", fixed = TRUE)
})


test_that("SAIGE still subsets samples when marker thinning is unnecessary", {
  plan <- KnockoffPipeline:::.saige_plink_export_plan(
    plink_prefix = "plink2",
    plink_file = "common",
    out_prefix = "prepared",
    total_markers = 100L,
    thin_target_markers = 5000L,
    plink_keep_file = "analysis.fam"
  )

  expect_false(plan$thin_markers)
  expect_true(plan$subset_samples)
  expect_match(plan$command, "--keep", fixed = TRUE)
  expect_match(plan$command, "analysis.fam", fixed = TRUE)
  expect_false(grepl("--thin", plan$command, fixed = TRUE))
  expect_false(grepl("--seed", plan$command, fixed = TRUE))
})


test_that("SAIGE reuses a small PLINK file when no sample subset is required", {
  plan <- KnockoffPipeline:::.saige_plink_export_plan(
    plink_prefix = "plink2",
    plink_file = "analysis",
    out_prefix = "prepared",
    total_markers = 100L,
    thin_target_markers = 5000L
  )

  expect_null(plan$command)
  expect_false(plan$thin_markers)
  expect_false(plan$subset_samples)
})
