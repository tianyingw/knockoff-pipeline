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
    genome_build = "hg19"
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
    gh_df = data.table::data.table()
  ))
})


test_that("an empty Single_Window result is marked complete", {
  mid_dir <- tempfile("single-empty-result-")
  on.exit(unlink(mid_dir, recursive = TRUE, force = TRUE), add = TRUE)

  KnockoffPipeline:::.write_block_result(
    mid_dir, chr = 1L, kk = 2L, single_df = NULL, window_df = NULL
  )

  expect_identical(
    KnockoffPipeline:::.get_completed_blocks(mid_dir, chr = 1L), 2L
  )
  expect_false(any(grepl(
    "^(Single|Window)_block_", list.files(file.path(mid_dir, "blocks", "chr1"))
  )))
})
