test_that("Single_Window skips an LD block with no input variants", {
  blocks <- data.table::data.table(chr = 1L, start = 100L, stop = 200L)
  bim <- data.frame(
    chr = "1",
    variant_id = "rs-outside",
    cm = 0,
    pos = 500,
    a1 = "A",
    a2 = "G",
    stringsAsFactors = FALSE
  )

  expect_null(KnockoffPipeline:::run_single_block(
    blocks = blocks,
    kk = 1L,
    geno.file = tempfile("unused-genotype-prefix-"),
    obj_nullmodel = NULL,
    window_length = 100L,
    plink_prefix = tempfile("unused-plink-executable-"),
    impute.method = "fixed",
    M = 1L,
    Gsub.id = "sample-1",
    bim_metadata = bim,
    genome_build = "hg19",
    export_switch = "--export A",
    plink_threads = 1L
  ))
})


test_that("Gene_Centric skips a batch with no input variants", {
  genes <- data.table::data.table(
    chr = 1L, start = 100000L, end = 101000L, id = "GENE1"
  )
  bim <- data.frame(
    chr = "1",
    variant_id = "rs-outside",
    cm = 0,
    pos = 500000,
    a1 = "A",
    a2 = "G",
    stringsAsFactors = FALSE
  )

  expect_null(KnockoffPipeline:::run_batch_gene(
    genes = genes,
    b = 1L,
    batch_index = list(1L),
    geno.file = tempfile("unused-genotype-prefix-"),
    obj_nullmodel = NULL,
    window_length = 100L,
    plink_prefix = tempfile("unused-plink-executable-"),
    M = 1L,
    genome_build = "hg19",
    Gsub.id = "sample-1",
    bim_metadata = bim,
    abc_df = data.table::data.table(),
    gh_df = data.table::data.table(),
    export_switch = "--export A",
    plink_threads = 1L
  ))
})


test_that("an empty Single_Window result is marked complete", {
  mid_dir <- tempfile("single-empty-result-")
  on.exit(unlink(mid_dir, recursive = TRUE, force = TRUE), add = TRUE)

  status <- KnockoffPipeline:::.write_block_result(
    mid_dir, chr = 1L, kk = 2L, single_df = NULL, window_df = NULL
  )

  expect_s3_class(status, "knockoff_pipeline_block_status")
  expect_named(status, c("chr", "block", "empty"))
  expect_identical(status$chr, 1L)
  expect_identical(status$block, 2L)
  expect_true(status$empty)
  expect_identical(
    KnockoffPipeline:::.get_completed_blocks(mid_dir, chr = 1L), 2L
  )
  expect_false(any(grepl(
    "^(Single|Window)_block_", list.files(file.path(mid_dir, "blocks", "chr1"))
  )))
})


test_that("per-block merge returns the tables it checkpointed", {
  mid_dir <- tempfile("single-merge-result-")
  on.exit(unlink(mid_dir, recursive = TRUE, force = TRUE), add = TRUE)

  single <- data.table::data.table(chr = 1L, pos = 101L, score = 0.5)
  window <- data.table::data.table(chr = 1L, start = 100L, end = 200L)
  status <- KnockoffPipeline:::.write_block_result(
    mid_dir, chr = 1L, kk = 1L, single_df = single, window_df = window
  )
  merged <- KnockoffPipeline:::.merge_per_block_files(mid_dir, chr = 1L)

  expect_false(status$empty)
  expect_equal(merged$single, single)
  expect_equal(merged$window, window)
  expect_true(file.exists(file.path(
    mid_dir, "Single_mid_results_chr1.txt"
  )))
  expect_true(file.exists(file.path(
    mid_dir, "Window_mid_results_chr1.txt"
  )))
})
