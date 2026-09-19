test_that("PLINK rows are aligned by character IID without losing leading zeros", {
  raw <- data.frame(
    FID = c("fam-A", "fam-001", "fam-B"),
    IID = c("A-2", "001", "B03"),
    PAT = "0",
    MAT = "0",
    SEX = 1L,
    PHENOTYPE = -9L,
    "rs12345_A" = c(0, 1, 2),
    "rs9_C" = c(2, 0, 1),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  bim <- data.frame(
    chr = c("1", "1"),
    variant_id = c("rs9", "rs12345"),
    cm = c(0, 0),
    pos = c(90000000, 87654321),
    a1 = c("C", "A"),
    a2 = c("T", "G"),
    stringsAsFactors = FALSE
  )
  target_ids <- c("001", "B03", "A-2")

  prepared <- KnockoffPipeline:::.prepare_raw_genotypes(
    raw = raw, target_ids = target_ids, bim_metadata = bim
  )

  expect_identical(prepared$sample_ids, target_ids)
  expect_identical(rownames(prepared$geno), target_ids)
  expect_equal(unname(prepared$geno[, "rs12345_A"]), c(1, 2, 0))
  expect_equal(unname(prepared$geno[, "rs9_C"]), c(0, 1, 2))

  # Positions must come from BIM metadata: digits in an rsID are not base pairs.
  expect_identical(prepared$variant_metadata$variant_id,
                   c("rs12345", "rs9"))
  expect_equal(prepared$variant_metadata$pos, c(87654321, 90000000))
})


test_that("identity-ordered integer genotypes stay integer", {
  raw <- data.frame(
    FID = c("f1", "f2", "f3"), IID = c("001", "002", "003"),
    PAT = "0", MAT = "0", SEX = 1L, PHENOTYPE = -9L,
    "rs1_A" = c(0L, 1L, 2L), "rs2_C" = c(2L, 1L, 0L),
    "rs3_G" = c(NA, NA, NA),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  bim <- data.frame(
    chr = rep("1", 3L), variant_id = c("rs1", "rs2", "rs3"),
    pos = c(101, 202, 303), a1 = c("A", "C", "G"),
    a2 = c("G", "T", "A"),
    stringsAsFactors = FALSE
  )

  prepared <- KnockoffPipeline:::.prepare_raw_genotypes(
    raw, target_ids = raw$IID, bim_metadata = bim
  )

  expect_type(prepared$geno, "integer")
  expect_identical(
    unname(prepared$geno),
    matrix(c(0L, 1L, 2L, 2L, 1L, 0L, rep(NA_integer_, 3L)),
           nrow = 3L)
  )
})


test_that("reordered genotype rows equal an explicit identity-result subset", {
  raw <- data.frame(
    FID = c("f1", "f2", "f3"), IID = c("s1", "s2", "s3"),
    PAT = "0", MAT = "0", SEX = 1L, PHENOTYPE = -9L,
    "rs1_A" = c(0L, 1L, 2L), "rs2_C" = c(2L, 0L, 1L),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  bim <- data.frame(
    chr = c("1", "1"), variant_id = c("rs1", "rs2"),
    pos = c(101, 202), a1 = c("A", "C"), a2 = c("G", "T"),
    stringsAsFactors = FALSE
  )
  identity <- KnockoffPipeline:::.prepare_raw_genotypes(
    raw, target_ids = raw$IID, bim_metadata = bim
  )
  target_ids <- c("s3", "s1", "s2")
  reordered <- KnockoffPipeline:::.prepare_raw_genotypes(
    raw, target_ids = target_ids, bim_metadata = bim
  )

  expected_index <- match(target_ids, identity$sample_ids)
  expect_identical(
    reordered$geno,
    identity$geno[expected_index, , drop = FALSE]
  )
  expect_identical(reordered$variant_metadata, identity$variant_metadata)
})


test_that("non-numeric PLINK genotype columns fail explicitly", {
  raw <- data.frame(
    FID = c("f1", "f2"), IID = c("s1", "s2"),
    PAT = "0", MAT = "0", SEX = 1L, PHENOTYPE = -9L,
    "rs1_A" = c("0", "not-a-genotype"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  bim <- data.frame(
    chr = "1", variant_id = "rs1", pos = 101, a1 = "A", a2 = "G",
    stringsAsFactors = FALSE
  )

  expect_error(
    KnockoffPipeline:::.prepare_raw_genotypes(
      raw, target_ids = raw$IID, bim_metadata = bim
    ),
    "non-numeric genotype column.*rs1_A"
  )
})


test_that("all stages use phenotype and covariate complete cases in .fam order", {
  input_dir <- tempfile("analysis-samples-")
  dir.create(input_dir)
  on.exit(unlink(input_dir, recursive = TRUE, force = TRUE), add = TRUE)

  pheno_file <- file.path(input_dir, "phenotype.csv")
  fam_file <- file.path(input_dir, "genotype.fam")
  data.table::fwrite(
    data.table::data.table(
      IID = c("003", "001", "004", "002"),
      Y = c(1, 2, NA, 4),
      PC1 = c(0.3, NA, 0.4, 0.2),
      Batch = c("A", "A", "B", "B")
    ),
    pheno_file
  )
  data.table::fwrite(
    data.table::data.table(
      FID = c("f1", "f2", "f3", "f4", "f5"),
      IID = c("001", "002", "003", "004", "005"),
      PAT = "0", MAT = "0", SEX = "0", PHENO = "-9"
    ),
    fam_file, sep = "\t", col.names = FALSE
  )

  prepared <- KnockoffPipeline:::.prepare_analysis_samples(
    pheno_file = pheno_file,
    phenotype = "Y",
    pheno_id = "IID",
    covar_cols = "PC1",
    cat_covar_cols = "Batch",
    plink_fam = fam_file
  )

  expect_identical(prepared$sample_ids, c("002", "003"))
  expect_identical(as.character(prepared$pheno$IID), c("002", "003"))
  expect_identical(as.character(prepared$plink_keep_fam$IID),
                   c("002", "003"))
  expect_identical(prepared$all_covar_cols, c("PC1", "Batch"))
})


test_that("additive genotype export selects the installed PLINK dialect", {
  expect_identical(
    KnockoffPipeline:::.plink_additive_export_switch(
      "plink", version_text = "PLINK v1.90b6.21 64-bit"
    ),
    "--recode A"
  )
  expect_identical(
    KnockoffPipeline:::.plink_additive_export_switch(
      "plink", version_text = "PLINK v2.00a6.5LM 64-bit"
    ),
    "--export A"
  )
})


test_that("additive export accepts cached dialect and an explicit thread count", {
  skip_on_os("windows")
  tmpdir <- tempfile("fake-plink-")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)
  fake_plink <- file.path(tmpdir, "fake plink")
  args_file <- file.path(tmpdir, "args.txt")
  writeLines(
    c("#!/bin/sh", sprintf("printf '%%s\\n' \"$@\" > %s", shQuote(args_file))),
    fake_plink
  )
  Sys.chmod(fake_plink, mode = "0755")

  status <- KnockoffPipeline:::.run_plink_additive_export(
    plink_prefix = fake_plink,
    geno_file = file.path(tmpdir, "input prefix"),
    chr = 1L, start = 10L, stop = 20L, keep_arg = "",
    out_prefix = file.path(tmpdir, "output prefix"),
    export_switch = "--export A", plink_threads = 3L
  )
  args <- readLines(args_file)

  expect_identical(status, 0L)
  expect_false("--version" %in% args)
  expect_identical(args[match("--export", args) + 1L], "A")
  expect_identical(args[match("--threads", args) + 1L], "3")

  default_status <- KnockoffPipeline:::.run_plink_additive_export(
    fake_plink, "input", 1L, 10L, 20L, "", "output",
    export_switch = "--recode A"
  )
  default_args <- readLines(args_file)
  expect_identical(default_status, 0L)
  expect_false("--threads" %in% default_args)

  expect_error(
    KnockoffPipeline:::.run_plink_additive_export(
      fake_plink, "input", 1L, 10L, 20L, "", "output",
      export_switch = "--export A", plink_threads = 0L
    ),
    "positive integer"
  )
})


test_that("saved knockoff rows allow reordering but reject a different set", {
  saved_ids <- c("001", "A-2", "B03")
  target_ids <- c("B03", "001", "A-2")
  ko_obj <- list(
    sample_ids = saved_ids,
    G_k = list(matrix(seq_len(6), nrow = 3L,
                      dimnames = list(saved_ids, c("v1", "v2"))))
  )

  expect_false(KnockoffPipeline:::.need_regenerate_samples(
    target_ids, saved_ids
  ))
  aligned <- KnockoffPipeline:::.align_knockoff_samples(ko_obj, target_ids)
  expect_identical(aligned$sample_ids, target_ids)
  expect_equal(
    unname(aligned$G_k[[1L]]),
    unname(ko_obj$G_k[[1L]][c(3L, 1L, 2L), , drop = FALSE])
  )

  expect_true(KnockoffPipeline:::.need_regenerate_samples(
    target_ids[-1L], saved_ids
  ))
  expect_true(KnockoffPipeline:::.need_regenerate_samples(
    c(target_ids[-1L], "new-sample"), saved_ids
  ))
  expect_error(
    KnockoffPipeline:::.align_knockoff_samples(
      ko_obj, c(target_ids, "new-sample")
    ),
    "target sample"
  )
})


test_that("atomic RDS replacement preserves the previous file until commit", {
  checkpoint <- tempfile("atomic-checkpoint-", fileext = ".rds")
  on.exit(unlink(c(checkpoint, paste0(checkpoint, ".previous-", Sys.getpid())),
                 force = TRUE), add = TRUE)

  saveRDS(list(generation = 1L), checkpoint)
  KnockoffPipeline:::.atomic_save_rds(list(generation = 2L), checkpoint)
  expect_identical(readRDS(checkpoint)$generation, 2L)
  expect_false(file.exists(paste0(checkpoint, ".previous-", Sys.getpid())))
})


test_that("sample-list reconciliation is atomic and fails closed", {
  sample_file <- tempfile("knockoff-samples-", fileext = ".txt")
  on.exit(unlink(c(sample_file,
                   paste0(sample_file, ".previous-", Sys.getpid())),
                 force = TRUE), add = TRUE)

  created <- KnockoffPipeline:::.reconcile_knockoff_sample_file(
    c("001", "A-2", "B03"), sample_file, create = TRUE
  )
  expect_true(created$created)
  expect_identical(readLines(sample_file), c("001", "A-2", "B03"))

  reordered <- KnockoffPipeline:::.reconcile_knockoff_sample_file(
    c("B03", "001", "A-2"), sample_file, create = TRUE
  )
  expect_false(reordered$created)
  expect_identical(reordered$sample_ids, c("001", "A-2", "B03"))
  expect_identical(reordered$order, c(2L, 3L, 1L))

  expect_error(
    KnockoffPipeline:::.reconcile_knockoff_sample_file(
      c("001", "A-2", "new-sample"), sample_file, create = TRUE
    ),
    "must exactly match"
  )
  expect_identical(readLines(sample_file), c("001", "A-2", "B03"))
})


test_that("interrupt cleanup protects files in the installed manifest", {
  knockoff_dir <- tempfile("installed-manifest-")
  dir.create(knockoff_dir)
  on.exit(unlink(knockoff_dir, recursive = TRUE, force = TRUE), add = TRUE)

  manifest_file <- file.path(knockoff_dir, "block_0001_knockoff.rds")
  descriptor <- "new_matrix_1.desc"
  installed_files <- file.path(
    knockoff_dir, c(descriptor, sub("\\.desc$", ".bin", descriptor))
  )
  orphan <- file.path(knockoff_dir, "new_matrix_orphan.bin")
  saveRDS(
    list(storage = "bigmemory_filebacked", descriptor_files = descriptor),
    manifest_file
  )

  cleanup <- KnockoffPipeline:::.uncommitted_backing_cleanup(
    c(installed_files, orphan), manifest_file
  )
  expect_identical(cleanup, orphan)
})


test_that("random imputation is deterministic under a derived unit seed", {
  variants <- data.frame(
    chr = c("1", "1", "1"),
    variant_id = c("rs1", "rs2", "rs3"),
    pos = c(100, 200, 300),
    a1 = c("A", "C", "G"),
    a2 = c("G", "T", "A"),
    raw_name = c("rs1_A", "rs2_C", "rs3_G"),
    counted_allele = c("A", "C", "G"),
    stringsAsFactors = FALSE
  )
  geno <- matrix(
    c(0, 1, NA, 2, NA, 0, 1, 2, 1, 0, 1, 2),
    nrow = 4L, byrow = TRUE
  )
  seed <- KnockoffPipeline:::.derive_unit_seed(20260915L, "unit", 1L)
  run_once <- function() KnockoffPipeline:::.with_local_seed(
    KnockoffPipeline:::.derive_unit_seed(seed, "imputation"),
    function() KnockoffPipeline:::Preprocess(
      geno, chr = 1L, window = 100L, impute.method = "random",
      variant_metadata = variants, thres.ultrarare = 0
    )
  )

  first <- run_once()
  second <- run_once()
  expect_equal(as.matrix(first$G), as.matrix(second$G))
  expect_identical(first$variant_metadata, second$variant_metadata)
})


test_that("equal variant counts cannot bypass knockoff context validation", {
  variants <- data.frame(
    chr = c("1", "1"),
    variant_id = c("rs12345", "rs9"),
    pos = c(87654321, 90000000),
    a1 = c("A", "C"),
    a2 = c("G", "T"),
    coded_allele = c("G", "T"),
    stringsAsFactors = FALSE
  )
  saved <- KnockoffPipeline:::.make_knockoff_context(
    test_type = "Single_Window",
    M = 5L,
    genome_build = "hg19",
    variant_metadata = variants,
    reference_id = "blocks.bed:md5:reference-a"
  )
  expect_invisible(KnockoffPipeline:::.assert_knockoff_context(saved, saved))

  reordered <- KnockoffPipeline:::.make_knockoff_context(
    "Single_Window", 5L, "hg19", variants[c(2L, 1L), ],
    "blocks.bed:md5:reference-a"
  )
  expect_error(
    KnockoffPipeline:::.assert_knockoff_context(saved, reordered),
    "variant_fingerprint"
  )

  changed_alleles <- variants
  changed_alleles[1L, c("a1", "a2")] <- c("G", "A")
  allele_context <- KnockoffPipeline:::.make_knockoff_context(
    "Single_Window", 5L, "hg19", changed_alleles,
    "blocks.bed:md5:reference-a"
  )
  expect_error(
    KnockoffPipeline:::.assert_knockoff_context(saved, allele_context),
    "variant_fingerprint"
  )

  build_context <- KnockoffPipeline:::.make_knockoff_context(
    "Single_Window", 5L, "hg38", variants,
    "blocks.bed:md5:reference-a"
  )
  expect_error(
    KnockoffPipeline:::.assert_knockoff_context(saved, build_context),
    "genome_build"
  )

  m_context <- KnockoffPipeline:::.make_knockoff_context(
    "Single_Window", 3L, "hg19", variants,
    "blocks.bed:md5:reference-a"
  )
  expect_error(
    KnockoffPipeline:::.assert_knockoff_context(saved, m_context),
    "M"
  )

  reference_context <- KnockoffPipeline:::.make_knockoff_context(
    "Single_Window", 5L, "hg19", variants,
    "blocks.bed:md5:reference-b"
  )
  expect_error(
    KnockoffPipeline:::.assert_knockoff_context(saved, reference_context),
    "reference_id"
  )

  constructor_context <- KnockoffPipeline:::.make_knockoff_context(
    "Single_Window", 5L, "hg19", variants,
    "blocks.bed:md5:reference-a", construction_id = "changed-settings"
  )
  expect_error(
    KnockoffPipeline:::.assert_knockoff_context(saved, constructor_context),
    "construction_id"
  )

  seeded_context <- KnockoffPipeline:::.make_knockoff_context(
    "Single_Window", 5L, "hg19", variants,
    "blocks.bed:md5:reference-a", random_seed = 123L
  )
  different_seed_context <- KnockoffPipeline:::.make_knockoff_context(
    "Single_Window", 5L, "hg19", variants,
    "blocks.bed:md5:reference-a", random_seed = 456L
  )
  expect_error(
    KnockoffPipeline:::.assert_knockoff_context(
      seeded_context, different_seed_context
    ),
    "random_seed"
  )
})


test_that("reference identity changes when file content changes", {
  reference_file <- tempfile("ld-blocks-", fileext = ".bed")
  on.exit(unlink(reference_file, force = TRUE), add = TRUE)

  writeLines(c("chr\tstart\tstop", "1\t1\t100"), reference_file)
  first <- KnockoffPipeline:::.reference_file_id(reference_file)
  writeLines(c("chr\tstart\tstop", "1\t1\t101"), reference_file)
  second <- KnockoffPipeline:::.reference_file_id(reference_file)

  expect_false(identical(first, second))
})


test_that("unit-specific seeds are deterministic", {
  block_seed <- KnockoffPipeline:::.derive_unit_seed(
    123L, "Single_Window", 1L, 2L
  )

  expect_identical(
    block_seed,
    KnockoffPipeline:::.derive_unit_seed(123L, "Single_Window", 1L, 2L)
  )
  expect_identical(block_seed, 643501276L)
  expect_false(identical(
    block_seed,
    KnockoffPipeline:::.derive_unit_seed(123L, "Single_Window", 1L, 3L)
  ))
  expect_null(KnockoffPipeline:::.derive_unit_seed(
    NULL, "Single_Window", 1L, 2L
  ))
})


test_that("sparse GRM is reordered and subset to the SAIGE model IDs", {
  grm <- Matrix::Matrix(
    matrix(c(
      1, 0.1, 0.2,
      0.1, 1, 0.3,
      0.2, 0.3, 1
    ), nrow = 3L, byrow = TRUE),
    sparse = TRUE
  )
  aligned <- KnockoffPipeline:::.align_sparse_grm(
    grm, grm_ids = c("001", "A-2", "B03"),
    model_ids = c("B03", "001")
  )

  expect_equal(as.matrix(aligned), as.matrix(grm[c(3L, 1L), c(3L, 1L)]))
  expect_error(
    KnockoffPipeline:::.align_sparse_grm(
      grm, c("001", "A-2", "B03"), c("missing", "001")
    ),
    "absent from the sparse GRM"
  )
  expect_error(
    KnockoffPipeline:::.align_sparse_grm(
      grm, c("001", "A-2"), c("001", "A-2")
    ),
    "dimensions"
  )
})


test_that("a file-backed knockoff manifest survives an R process boundary", {
  skip_if_not_installed("bigmemory")
  rscript <- file.path(R.home("bin"), "Rscript")
  skip_if(!file.exists(rscript), "Rscript is required for this process-boundary test")

  knockoff_dir <- tempfile("knockoff-manifest-")
  dir.create(knockoff_dir)
  on.exit(unlink(knockoff_dir, recursive = TRUE, force = TRUE), add = TRUE)

  manifest_file <- file.path(knockoff_dir, "block_0001_knockoff.rds")
  backing_file <- "block_0001_knockoff_matrix_1.bin"
  descriptor_file <- "block_0001_knockoff_matrix_1.desc"
  matrix_values <- matrix(c(0, 1, 2, 2, 1, 0), nrow = 3L)
  variant_metadata <- data.frame(
    chr = c("1", "1"),
    variant_id = c("rs1", "rs2"),
    pos = c(100, 200),
    a1 = c("A", "C"),
    a2 = c("G", "T"),
    coded_allele = c("G", "T"),
    stringsAsFactors = FALSE
  )

  file_backed <- bigmemory::filebacked.big.matrix(
    nrow = 3L,
    ncol = 2L,
    init = 0,
    backingfile = backing_file,
    descriptorfile = descriptor_file,
    backingpath = knockoff_dir
  )
  file_backed[,] <- matrix_values
  saveRDS(
    list(
      storage = "bigmemory_filebacked",
      G_k = NULL,
      descriptor_files = descriptor_file,
      sample_ids = c("001", "A-2", "B03"),
      snp_pos = variant_metadata$pos,
      context = KnockoffPipeline:::.make_knockoff_context(
        "Single_Window", 1L, "hg19", variant_metadata,
        "blocks.bed:md5:reference-a"
      )
    ),
    manifest_file
  )
  rm(file_backed)
  invisible(gc())

  child_script <- file.path(knockoff_dir, "reattach.R")
  child_result <- file.path(knockoff_dir, "reattached.rds")
  writeLines(
    c(
      "args <- commandArgs(trailingOnly = TRUE)",
      "manifest <- readRDS(args[[1L]])",
      "stopifnot(identical(manifest$storage, 'bigmemory_filebacked'))",
      "stopifnot(identical(manifest$context$schema_version, 4L))",
      "x <- bigmemory::attach.big.matrix(",
      "  manifest$descriptor_files[[1L]], path = dirname(args[[1L]])",
      ")",
      "saveRDS(x[, ], args[[2L]])"
    ),
    child_script
  )
  status <- system2(
    rscript,
    c("--vanilla", shQuote(child_script), shQuote(manifest_file),
      shQuote(child_result))
  )

  expect_identical(status, 0L)
  expect_equal(readRDS(child_result), matrix_values)
  expect_true(all(file.exists(file.path(
    knockoff_dir, c(manifest_file = basename(manifest_file),
                    descriptor_file, backing_file)
  ))))
})
