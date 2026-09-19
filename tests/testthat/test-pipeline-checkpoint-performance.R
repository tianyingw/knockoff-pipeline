test_that("per-batch merge returns the gene table it checkpointed", {
  mid_dir <- tempfile("gene-merge-result-")
  on.exit(unlink(mid_dir, recursive = TRUE, force = TRUE), add = TRUE)

  batch_dir <- KnockoffPipeline:::.ensure_batch_dir(mid_dir, chr = 1L)
  first <- data.table::data.table(
    chr = 1L, gene_id = "GENE1", gene_start = 100L, gene_end = 200L
  )
  second <- data.table::data.table(
    chr = 1L, gene_id = "GENE2", gene_start = 300L, gene_end = 400L
  )
  data.table::fwrite(
    first, file.path(batch_dir, "GeneCentric_batch_run_b0001.txt"),
    sep = "\t"
  )
  data.table::fwrite(
    second, file.path(batch_dir, "GeneCentric_batch_run_b0002.txt"),
    sep = "\t"
  )

  merged <- KnockoffPipeline:::.merge_batch_files(mid_dir, chr = 1L)

  expect_equal(merged, data.table::rbindlist(list(first, second)))
  expect_true(file.exists(file.path(
    mid_dir, "GeneCentric_mid_results_chr1.txt"
  )))
})
