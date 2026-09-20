make_gene_region_fixture <- function() {
  samples <- c("s1", "s2", "s3")
  enhancer_start <- 4000000L
  enhancer_end <- 4000500L
  positions <- c(
    895000L, 945000L, 995100L, 999900L, 1000200L, 1005900L,
    1051000L, 1106000L, 3950000L, 3995000L, 4000050L, 4000100L,
    4000150L, 4000200L, 4000250L, 4000300L, 4005000L, 4050500L
  )
  variant_id <- paste0("rs", seq_along(positions))
  bim <- data.frame(
    chr = "1", variant_id = variant_id, cm = 0, pos = positions,
    a1 = "A", a2 = "G", stringsAsFactors = FALSE
  )
  attr(bim, "kp_pos_sorted") <- TRUE

  genotype <- vapply(seq_along(positions), function(j) {
    as.integer((seq_along(samples) + j) %% 3L)
  }, integer(length(samples)))
  dimnames(genotype) <- list(samples, variant_id)

  list(
    samples = samples,
    positions = positions,
    bim = bim,
    genotype = genotype,
    genes = data.table::data.table(
      chr = 1L, start = 1000000L, end = 1001000L, id = "GENE1"
    ),
    abc = data.table::data.table(
      TargetGene = "GENE1", start = enhancer_start, end = enhancer_end
    ),
    gh = data.table::data.table(
      gene = character(), GH_start = integer(), GH_end = integer()
    ),
    enhancer_start = enhancer_start,
    enhancer_end = enhancer_end,
    enhancer_target_positions = positions[
      positions >= enhancer_start & positions <= enhancer_end
    ]
  )
}


write_mock_plink_raw <- function(fixture, variant_index, out_prefix,
                                 sample_order = fixture$samples) {
  sample_index <- match(sample_order, fixture$samples)
  raw <- data.frame(
    FID = sample_order, IID = sample_order,
    PAT = 0L, MAT = 0L, SEX = 1L, PHENOTYPE = -9L,
    check.names = FALSE
  )
  for (j in variant_index) {
    raw[[paste0(fixture$bim$variant_id[j], "_", fixture$bim$a1[j])]] <-
      fixture$genotype[sample_index, j]
  }
  data.table::fwrite(raw, paste0(out_prefix, ".raw"))
  invisible(NULL)
}


make_two_gene_region_fixture <- function() {
  samples <- c("s1", "s2", "s3")
  positions <- sort(c(
    945000L, 995100L, 999900L, 1000200L, 1005900L, 1051000L,
    1945000L, 1995100L, 1999900L, 2000200L, 2005900L, 2051000L,
    3995000L, seq(4000050L, 4000300L, by = 50L), 4005000L,
    4995000L, seq(5000050L, 5000300L, by = 50L), 5005000L
  ))
  variant_id <- paste0("rs", seq_along(positions))
  bim <- data.frame(
    chr = "1", variant_id = variant_id, cm = 0, pos = positions,
    a1 = "A", a2 = "G", stringsAsFactors = FALSE
  )
  attr(bim, "kp_pos_sorted") <- TRUE

  genotype <- vapply(seq_along(positions), function(j) {
    as.integer((seq_along(samples) + j) %% 3L)
  }, integer(length(samples)))
  dimnames(genotype) <- list(samples, variant_id)

  list(
    samples = samples,
    positions = positions,
    bim = bim,
    genotype = genotype,
    genes = data.table::data.table(
      chr = 1L,
      start = c(1000000L, 2000000L),
      end = c(1001000L, 2001000L),
      id = c("GENE1", "GENE2")
    ),
    abc = data.table::data.table(
      TargetGene = "GENE1", start = 4000000L, end = 4000500L
    ),
    gh = data.table::data.table(
      gene = "GENE2", GH_start = 5000000L, GH_end = 5000500L
    )
  )
}


test_that("gene-region specifications retain method-specific source flanks", {
  glm <- KnockoffPipeline:::.gene_region_spec(use_glmm = FALSE)
  glmm <- KnockoffPipeline:::.gene_region_spec(use_glmm = TRUE)

  expect_identical(
    glm,
    list(
      gene_buffer_bp = 5000L,
      gene_neighbor_bp = 10000L,
      gene_source_flank_bp = 55000L,
      enhancer_source_flank_bp = 10000L
    )
  )
  expect_identical(
    glmm,
    list(
      gene_buffer_bp = 5000L,
      gene_neighbor_bp = 100000L,
      gene_source_flank_bp = 105000L,
      enhancer_source_flank_bp = 50000L
    )
  )
})


test_that("ABC and GH enhancers are combined in stable order and de-duplicated", {
  abc <- data.table::data.table(
    TargetGene = c("GENE1", "GENE1", "GENE1", "OTHER"),
    start = c(100L, 100L, 300L, 700L),
    end = c(200L, 200L, 400L, 800L)
  )
  gh <- data.table::data.table(
    gene = c("GENE1", "GENE1", "GENE1", "OTHER"),
    GH_start = c(100L, 500L, 500L, 900L),
    GH_end = c(200L, 600L, 600L, 1000L)
  )

  got <- KnockoffPipeline:::.gene_enhancers("GENE1", abc, gh)
  expect_named(got, c("start", "end"))
  expect_equal(
    unname(as.matrix(got)),
    matrix(c(100L, 300L, 500L, 200L, 400L, 600L), ncol = 2L)
  )

  empty <- KnockoffPipeline:::.gene_enhancers("MISSING", abc, gh)
  expect_named(empty, c("start", "end"))
  expect_equal(nrow(empty), 0L)
})


test_that("malformed enhancer references fail instead of dropping regions", {
  expect_error(
    KnockoffPipeline:::.gene_enhancers(
      "GENE1", data.frame(start = 1, end = 2), NULL
    ),
    "ABC enhancer table is missing: TargetGene"
  )
  expect_error(
    KnockoffPipeline:::.gene_enhancers(
      "GENE1", NULL, data.frame(gene = "GENE1", GH_start = 1)
    ),
    "GeneHancer table is missing: GH_end"
  )
  expect_error(
    KnockoffPipeline:::.gene_enhancers(
      "GENE1",
      data.frame(TargetGene = "GENE1", start = "bad", end = 2),
      NULL
    ),
    "invalid coordinates"
  )
})


test_that("interval BIM subsetting preserves source order without repeated rows", {
  bim <- data.frame(
    chr = "1",
    variant_id = c("rs1", "rs2", "rs3a", "rs3b", "rs4", "rs5", "rs6"),
    cm = 0,
    pos = c(50L, 100L, 125L, 125L, 160L, 210L, 300L),
    a1 = "A", a2 = "G",
    stringsAsFactors = FALSE
  )
  attr(bim, "kp_pos_sorted") <- TRUE
  intervals <- data.frame(
    start = c(140L, 90L, 90L),
    end = c(220L, 150L, 150L)
  )

  got <- KnockoffPipeline:::.subset_bim_intervals(bim, intervals)
  expect_equal(got, bim[c(2L, 3L, 4L, 5L, 6L), , drop = FALSE],
               ignore_attr = TRUE)
  expect_false(anyDuplicated(got$variant_id) > 0L)

  empty <- KnockoffPipeline:::.subset_bim_intervals(
    bim, data.frame(start = integer(), end = integer())
  )
  expect_equal(empty, bim[0, , drop = FALSE], ignore_attr = TRUE)
})


test_that("a distal enhancer is loaded independently of the gene batch range", {
  fixture <- make_gene_region_fixture()
  observed <- new.env(parent = emptyenv())
  observed$range_calls <- list()
  observed$extract_calls <- list()
  observed$analysis <- NULL

  testthat::local_mocked_bindings(
    .run_plink_additive_export = function(..., start, stop, out_prefix) {
      observed$range_calls[[length(observed$range_calls) + 1L]] <-
        c(start = start, stop = stop)
      index <- which(fixture$positions >= start & fixture$positions <= stop)
      write_mock_plink_raw(fixture, index, out_prefix)
      0L
    },
    .run_plink_additive_extract = function(..., extract_file, out_prefix) {
      requested <- scan(extract_file, what = character(), quiet = TRUE)
      requested <- unique(requested[requested %in% fixture$bim$variant_id])
      observed$extract_calls[[length(observed$extract_calls) + 1L]] <- requested
      index <- rev(match(requested, fixture$bim$variant_id))
      # Deliberately reverse both rows and variant columns.  The independent
      # import must restore sample order and chromosome-BIM variant order.
      write_mock_plink_raw(
        fixture, index, out_prefix,
        sample_order = rev(fixture$samples)
      )
      0L
    },
    GeneScan3D.KnockoffGeneration = function(...) {
      dots <- list(...)
      observed$analysis <- dots
      list(
        GeneScan3D.Cauchy = c(0.5, 0.5, 0.5),
        GeneScan3D.Cauchy_knockoff = matrix(0.5, nrow = dots$M, ncol = 3L)
      )
    },
    .package = "KnockoffPipeline"
  )

  got <- KnockoffPipeline:::run_batch_gene(
    genes = fixture$genes, b = 1L, batch_index = list(1L),
    geno.file = "unused", obj_nullmodel = list(), window_length = 1000L,
    plink_prefix = "unused", M = 1L, genome_build = "hg19",
    Gsub.id = fixture$samples, bim_metadata = fixture$bim,
    abc_df = fixture$abc, gh_df = fixture$gh,
    use_glmm = FALSE, user_cores = 1L,
    export_switch = "--export A", plink_threads = 1L
  )

  expect_equal(nrow(got), 1L)
  expect_length(observed$range_calls, 1L)
  expect_true(length(observed$extract_calls) >= 1L)
  expect_true(all(vapply(observed$range_calls, function(x) {
    x[["stop"]] < fixture$enhancer_start || x[["start"]] > fixture$enhancer_end
  }, logical(1))))

  enhancer_ids <- fixture$bim$variant_id[
    fixture$bim$pos %in% fixture$enhancer_target_positions
  ]
  expect_true(all(enhancer_ids %in% unlist(observed$extract_calls)))

  analysis <- observed$analysis
  expect_identical(analysis$R, 1L)
  expect_equal(
    as.numeric(analysis$Enhancer.pos),
    c(fixture$enhancer_start, fixture$enhancer_end)
  )
  expect_equal(analysis$p.EnhancerAll,
               length(fixture$enhancer_target_positions))
  expect_true(all(
    fixture$enhancer_target_positions %in%
      analysis$variants_EnhancerAll_surround
  ))
  enhancer_surround_index <- which(
    fixture$positions >= fixture$enhancer_start - 10000L &
      fixture$positions <= fixture$enhancer_end + 10000L
  )
  expect_equal(
    analysis$variants_EnhancerAll_surround,
    fixture$positions[enhancer_surround_index],
    tolerance = 0
  )
  expect_false(any(
    fixture$enhancer_target_positions %in%
      analysis$variants_gene_buffer_surround
  ))
  expect_identical(
    rownames(analysis$G_EnhancerAll_surround), fixture$samples
  )
  expect_equal(
    ncol(analysis$G_EnhancerAll_surround),
    length(analysis$variants_EnhancerAll_surround)
  )
  expect_equal(
    unname(analysis$G_EnhancerAll_surround),
    unname(fixture$genotype[, enhancer_surround_index, drop = FALSE])
  )
})


test_that("BIGKnock dispatch uses its wider independent source regions", {
  fixture <- make_gene_region_fixture()
  observed <- new.env(parent = emptyenv())
  observed$gene_range <- NULL
  observed$requested <- character()
  observed$analysis <- NULL

  testthat::local_mocked_bindings(
    .run_plink_additive_export = function(..., start, stop, out_prefix) {
      observed$gene_range <- c(start = start, stop = stop)
      index <- which(fixture$positions >= start & fixture$positions <= stop)
      write_mock_plink_raw(fixture, index, out_prefix)
      0L
    },
    .run_plink_additive_extract = function(..., extract_file, out_prefix) {
      requested <- scan(extract_file, what = character(), quiet = TRUE)
      observed$requested <- requested
      index <- match(requested, fixture$bim$variant_id)
      write_mock_plink_raw(fixture, index, out_prefix)
      0L
    },
    GeneScan3D.UKB.GLMM.KnockoffGeneration = function(...) {
      dots <- list(...)
      observed$analysis <- dots
      list(
        GeneScan3D.Cauchy = c(0.5, 0.5, 0.5),
        GeneScan3D.Cauchy_knockoff = matrix(0.5, nrow = dots$M, ncol = 3L)
      )
    },
    .package = "KnockoffPipeline"
  )

  got <- KnockoffPipeline:::run_batch_gene(
    genes = fixture$genes, b = 1L, batch_index = list(1L),
    geno.file = "unused", obj_nullmodel = list(), window_length = 1000L,
    plink_prefix = "unused", M = 1L, genome_build = "hg19",
    Gsub.id = fixture$samples, bim_metadata = fixture$bim,
    abc_df = fixture$abc, gh_df = fixture$gh,
    use_glmm = TRUE, user_cores = 1L,
    export_switch = "--export A", plink_threads = 1L
  )

  expect_equal(nrow(got), 1L)
  expect_equal(
    observed$gene_range,
    c(start = fixture$genes$start - 105000L,
      stop = fixture$genes$end + 105000L),
    tolerance = 0
  )
  expected_enhancer <- which(
    fixture$positions >= fixture$enhancer_start - 50000L &
      fixture$positions <= fixture$enhancer_end + 50000L
  )
  expect_identical(
    observed$requested, fixture$bim$variant_id[expected_enhancer]
  )
  expect_equal(
    observed$analysis$variants_gene_buffer_surround,
    fixture$positions[fixture$positions >= fixture$genes$start - 105000L &
      fixture$positions <= fixture$genes$end + 105000L],
    tolerance = 0
  )
  expect_equal(
    observed$analysis$variants_EnhancerAll_surround,
    fixture$positions[expected_enhancer], tolerance = 0
  )
})


test_that("gene inputs are invariant to joined or split batches", {
  fixture <- make_two_gene_region_fixture()

  run_layout <- function(batch_index) {
    observed <- new.env(parent = emptyenv())
    observed$by_gene <- list()

    testthat::local_mocked_bindings(
      .run_plink_additive_export = function(..., start, stop, out_prefix) {
        index <- which(fixture$positions >= start & fixture$positions <= stop)
        write_mock_plink_raw(fixture, index, out_prefix)
        0L
      },
      .run_plink_additive_extract = function(..., extract_file, out_prefix) {
        requested <- scan(extract_file, what = character(), quiet = TRUE)
        index <- match(requested, fixture$bim$variant_id)
        # PLINK is not required to emit columns or samples in extract-file
        # order.  Reverse both here so invariance relies on explicit identity
        # alignment rather than incidental ordering.
        write_mock_plink_raw(
          fixture, rev(index), out_prefix,
          sample_order = rev(fixture$samples)
        )
        0L
      },
      GeneScan3D.KnockoffGeneration = function(...) {
        dots <- list(...)
        gene <- if (mean(dots$gene_buffer.pos) < 1500000) {
          "GENE1"
        } else {
          "GENE2"
        }
        observed$by_gene[[gene]] <- list(
          gene_positions = dots$variants_gene_buffer_surround,
          enhancer_positions = dots$variants_EnhancerAll_surround,
          enhancer_coordinates = unname(as.matrix(dots$Enhancer.pos)),
          enhancer_surround_sizes = dots$p_EnhancerAll_surround,
          enhancer_target_sizes = dots$p.EnhancerAll,
          gene_matrix = unname(as.matrix(dots$G_gene_buffer_surround)),
          enhancer_matrix = unname(as.matrix(dots$G_EnhancerAll_surround))
        )
        list(
          GeneScan3D.Cauchy = c(0.5, 0.5, 0.5),
          GeneScan3D.Cauchy_knockoff =
            matrix(0.5, nrow = dots$M, ncol = 3L)
        )
      },
      .package = "KnockoffPipeline"
    )

    results <- lapply(seq_along(batch_index), function(b) {
      KnockoffPipeline:::run_batch_gene(
        genes = fixture$genes, b = b, batch_index = batch_index,
        geno.file = "unused", obj_nullmodel = list(),
        window_length = 1000L, plink_prefix = "unused", M = 1L,
        genome_build = "hg19", Gsub.id = fixture$samples,
        bim_metadata = fixture$bim, abc_df = fixture$abc,
        gh_df = fixture$gh, use_glmm = FALSE, user_cores = 1L,
        export_switch = "--export A", plink_threads = 1L
      )
    })

    list(
      inputs = observed$by_gene,
      results = data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
    )
  }

  joined <- run_layout(list(1:2))
  split <- run_layout(list(1L, 2L))

  expect_setequal(names(joined$inputs), c("GENE1", "GENE2"))
  expect_setequal(names(split$inputs), c("GENE1", "GENE2"))
  for (gene in c("GENE1", "GENE2")) {
    expect_equal(joined$inputs[[gene]], split$inputs[[gene]], tolerance = 0)
  }
  expect_equal(
    joined$results[order(gene_id)],
    split$results[order(gene_id)],
    tolerance = 0
  )
})


test_that("stage 1 does not inspect or extract enhancer regions", {
  fixture <- make_gene_region_fixture()
  observed <- new.env(parent = emptyenv())
  observed$generation_called <- FALSE

  testthat::local_mocked_bindings(
    .run_plink_additive_export = function(..., start, stop, out_prefix) {
      index <- which(fixture$positions >= start & fixture$positions <= stop)
      write_mock_plink_raw(fixture, index, out_prefix)
      0L
    },
    .run_plink_additive_extract = function(...) {
      stop("stage 1 must not extract enhancer genotypes")
    },
    .gene_enhancers = function(...) {
      stop("stage 1 must not inspect enhancer maps")
    },
    GeneScan3D.KnockoffGeneration = function(..., stage1_only) {
      expect_true(stage1_only)
      observed$generation_called <- TRUE
      NULL
    },
    .package = "KnockoffPipeline"
  )

  expect_null(KnockoffPipeline:::run_batch_gene(
    genes = fixture$genes, b = 1L, batch_index = list(1L),
    geno.file = "unused", obj_nullmodel = NULL, window_length = 1000L,
    plink_prefix = "unused", M = 1L, genome_build = "hg19",
    Gsub.id = fixture$samples, bim_metadata = fixture$bim,
    abc_df = NULL, gh_df = NULL,
    use_glmm = FALSE, user_cores = 1L,
    save_knockoff = TRUE, stage1_only = TRUE,
    export_switch = "--export A", plink_threads = 1L
  ))
  expect_true(observed$generation_called)
})
