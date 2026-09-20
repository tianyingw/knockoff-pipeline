.make_gene_bundle_fixture <- function(n = 60L) {
  samples <- sprintf("sample_%03d", seq_len(n))
  genes <- data.table::data.table(
    chr = 1L,
    start = c(1000000L, 2000000L),
    end = c(1001000L, 2001000L),
    id = c("GENE1", "GENE2")
  )
  enhancer_start <- c(4000000L, 5000000L)
  enhancer_end <- enhancer_start + 500L
  positions <- sort(c(
    945000L, 996000L, 998000L, 1000000L, 1000500L, 1001000L,
    1003000L, 1005000L, 1056000L,
    1945000L, 1996000L, 1998000L, 2000000L, 2000500L, 2001000L,
    2003000L, 2005000L, 2056000L,
    3990000L, seq(4000050L, 4000300L, by = 50L), 4010000L,
    4990000L, seq(5000050L, 5000300L, by = 50L), 5010000L
  ))
  variant_id <- sprintf("bundle_rs%03d", seq_along(positions))
  bim <- data.frame(
    chr = "1", variant_id = variant_id, cm = 0, pos = positions,
    a1 = "A", a2 = "G", stringsAsFactors = FALSE
  )
  attr(bim, "kp_pos_sorted") <- TRUE

  genotype <- vapply(seq_along(positions), function(j) {
    as.integer((seq_len(n) + 2L * j + seq_len(n) %/% 7L) %% 3L)
  }, integer(n))
  dimnames(genotype) <- list(samples, variant_id)

  list(
    samples = samples,
    positions = positions,
    bim = bim,
    genotype = genotype,
    genes = genes,
    abc = data.table::data.table(
      TargetGene = genes$id,
      start = enhancer_start,
      end = enhancer_end
    ),
    gh = data.table::data.table(
      gene = character(), GH_start = integer(), GH_end = integer()
    )
  )
}


.write_gene_bundle_raw <- function(fixture, variant_index, out_prefix,
                                   sample_order = rev(fixture$samples)) {
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


.gene_bundle_array <- function(x, positions, lower, upper, M) {
  index <- which(positions >= lower & positions <= upper)
  original <- as.matrix(x[, index, drop = FALSE])
  out <- array(NA_real_, dim = c(M, nrow(original), ncol(original)))
  for (k in seq_len(M)) out[k, , ] <- original + k / 100
  out
}


.gene_bundle_score <- function(G, G_enhancer, ids) {
  combined <- as.matrix(G)
  if (!is.null(G_enhancer)) combined <- cbind(combined, G_enhancer)
  row_weight <- match(as.character(ids), sort(as.character(ids)))
  value <- sum(combined * row_weight) / 1e7
  c(value, value + 0.01, value + 0.02)
}


.record_gene_bundle_call <- function(observed, G, G_enhancer, ids) {
  if (is.null(observed)) return(invisible(NULL))
  if (is.null(observed$calls)) observed$calls <- list()
  observed$calls[[length(observed$calls) + 1L]] <- list(
    gene = unname(as.matrix(G)),
    enhancer = if (is.null(G_enhancer)) NULL else
      unname(as.matrix(G_enhancer)),
    ids = as.character(ids)
  )
  invisible(NULL)
}


.gene_bundle_null <- function(ids, use_glmm) {
  y <- seq_along(ids) / length(ids)
  if (isTRUE(use_glmm)) {
    return(list(
      fitted.values = rep(mean(y), length(y)),
      residuals = y - mean(y),
      sampleID = ids,
      X = matrix(1, nrow = length(y), ncol = 1L),
      traitType = "C",
      theta = 1
    ))
  }
  list(
    nullglm = list(fitted.values = rep(mean(y), length(y))),
    Y = y,
    re.Y.res = y - mean(y),
    X0 = matrix(1, nrow = length(y), ncol = 1L),
    out_type = "C",
    id = ids
  )
}


.gene_bundle_has_knockoff <- function(entry) {
  if (!is.list(entry) || is.null(names(entry))) return(FALSE)
  candidates <- grep("knockoff", names(entry), ignore.case = TRUE)
  any(vapply(candidates, function(i) {
    value <- entry[[i]]
    is.array(value) && length(dim(value)) == 3L && length(value) > 0L
  }, logical(1)))
}


.run_gene_bundle_layout <- function(fixture, batch_index, ids,
                                    use_glmm = FALSE,
                                    stage1_only = FALSE,
                                    save_knockoff = FALSE,
                                    load_knockoff = FALSE,
                                    knockoff_dir = NULL,
                                    read_mid_exist = FALSE,
                                    forbid_generation = FALSE,
                                    observed = NULL) {
  range_export <- function(..., start, stop, out_prefix) {
    index <- which(fixture$positions >= start & fixture$positions <= stop)
    .write_gene_bundle_raw(fixture, index, out_prefix)
    0L
  }
  interval_export <- function(..., extract_file, out_prefix) {
    requested <- scan(extract_file, what = character(), quiet = TRUE)
    index <- match(requested, fixture$bim$variant_id)
    .write_gene_bundle_raw(fixture, rev(index), out_prefix)
    0L
  }

  gene_generator <- function(X, pos, gene_buffer_start, gene_buffer_end,
                             M, ...) {
    if (isTRUE(forbid_generation)) stop("gene generation must not run")
    .gene_bundle_array(
      X, pos, gene_buffer_start, gene_buffer_end, M
    )
  }
  enhancer_generator <- function(X, pos, Enhancer_start, Enhancer_end,
                                 M, ...) {
    if (isTRUE(forbid_generation)) stop("enhancer generation must not run")
    .gene_bundle_array(X, pos, Enhancer_start, Enhancer_end, M)
  }
  big_generator <- function(X, positions, lower, upper, M) {
    if (isTRUE(forbid_generation)) stop("BIGKnock generation must not run")
    index <- which(positions >= lower & positions <= upper)
    list(
      knockoff = .gene_bundle_array(X, positions, lower, upper, M),
      selected_input_index = index,
      positions = positions[index]
    )
  }

  null_model <- if (isTRUE(stage1_only)) NULL else
    .gene_bundle_null(ids, use_glmm)
  common <- list(
    genes = fixture$genes,
    geno.file = "unused",
    obj_nullmodel = null_model,
    window_length = 1000L,
    plink_prefix = "unused",
    M = 2L,
    genome_build = "hg19",
    Gsub.id = ids,
    bim_metadata = fixture$bim,
    reference_id = "gene-reference-v1",
    enhancer_reference_id = "enhancer-reference-v1",
    seed = 41L,
    abc_df = fixture$abc,
    gh_df = fixture$gh,
    use_glmm = use_glmm,
    user_cores = 1L,
    save_knockoff = save_knockoff,
    load_knockoff = load_knockoff,
    knockoff_dir = knockoff_dir,
    stage1_only = stage1_only,
    read_mid_exist = read_mid_exist,
    export_switch = "--export A",
    plink_threads = 1L
  )
  if (isTRUE(use_glmm)) {
    common$sparseSigma <- Matrix::Diagonal(length(ids))
    common$ratio <- 1
    common$glmm_precomputed <- list(
      C = matrix(1, 1L, 1L), inv_vX = matrix(1, 1L, 1L), outcome = "C"
    )
  }

  run_batches <- function() {
    out <- lapply(seq_along(batch_index), function(b) {
      do.call(
        KnockoffPipeline:::run_batch_gene,
        c(common, list(b = b, batch_index = batch_index))
      )
    })
    out <- Filter(Negate(is.null), out)
    if (length(out) == 0L) return(NULL)
    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  }

  if (!isTRUE(use_glmm)) {
    testthat::local_mocked_bindings(
      .run_plink_additive_export = range_export,
      .run_plink_additive_extract = interval_export,
      create.MK.AL_gene_buffer = gene_generator,
      create.MK.AL_Enhancer = enhancer_generator,
      GeneScan3D = function(G, G.EnhancerAll = NULL, Gsub.id, ...) {
        .record_gene_bundle_call(observed, G, G.EnhancerAll, Gsub.id)
        list(GeneScan3D.Cauchy.pvalue =
          .gene_bundle_score(G, G.EnhancerAll, Gsub.id))
      },
      .package = "KnockoffPipeline"
    )
    return(run_batches())
  }

  testthat::local_mocked_bindings(
    .run_plink_additive_export = range_export,
    .run_plink_additive_extract = interval_export,
    Knockoffgeneration.gene.buffer = function(
      G_gene_buffer_surround, positions, gene_buffer_start,
      gene_buffer_end, M, ...
    ) {
      big_generator(
        G_gene_buffer_surround, positions, gene_buffer_start,
        gene_buffer_end, M
      )
    },
    Knockoffgeneration.enhancer = function(
      G_enhancer_surround, positions, enhancer_start, enhancer_end, M, ...
    ) {
      big_generator(
        G_enhancer_surround, positions, enhancer_start, enhancer_end, M
      )
    },
    GeneScan3D.UKB.GLMM = function(
      G, G.EnhancerAll = NULL, Gsub.id, ...
    ) {
      .record_gene_bundle_call(observed, G, G.EnhancerAll, Gsub.id)
      list(GeneScan3D.Cauchy.pvalue =
        .gene_bundle_score(G, G.EnhancerAll, Gsub.id))
    },
    .package = "KnockoffPipeline"
  )
  run_batches()
}


test_that("GeneScan bundles persist enhancers across batch and sample order", {
  fixture <- .make_gene_bundle_fixture()
  knockoff_dir <- tempfile("genescan-bundles-")
  dir.create(knockoff_dir)
  on.exit(unlink(knockoff_dir, recursive = TRUE, force = TRUE), add = TRUE)

  expect_null(.run_gene_bundle_layout(
    fixture, batch_index = list(1:2), ids = fixture$samples,
    stage1_only = TRUE, save_knockoff = TRUE,
    knockoff_dir = knockoff_dir
  ))

  for (gene in fixture$genes$id) {
    path <- file.path(knockoff_dir, paste0("gene_", gene, "_ko.rds"))
    expect_true(file.exists(path))
    bundle <- readRDS(path)
    expect_identical(bundle$bundle_schema_version, 2L)
    expect_true(is.array(bundle$G_gene_buffer_knockoff))
    expect_identical(as.character(bundle$sample_ids), fixture$samples)
    expect_true(length(bundle$enhancers) > 0L)
    expect_true(any(vapply(
      bundle$enhancers, .gene_bundle_has_knockoff, logical(1)
    )))
  }

  # A resumed stage-1 run must reuse each complete per-gene bundle even when
  # genes are regrouped into different batches.
  expect_null(.run_gene_bundle_layout(
    fixture, batch_index = list(1L, 2L), ids = fixture$samples,
    stage1_only = TRUE, save_knockoff = TRUE,
    knockoff_dir = knockoff_dir, read_mid_exist = TRUE,
    forbid_generation = TRUE
  ))

  reversed_ids <- rev(fixture$samples)
  loaded_calls <- new.env(parent = emptyenv())
  fresh_calls <- new.env(parent = emptyenv())
  loaded <- .run_gene_bundle_layout(
    fixture, batch_index = list(1L, 2L), ids = reversed_ids,
    load_knockoff = TRUE, knockoff_dir = knockoff_dir,
    forbid_generation = TRUE, observed = loaded_calls
  )
  fresh <- .run_gene_bundle_layout(
    fixture, batch_index = list(1L, 2L), ids = reversed_ids,
    observed = fresh_calls
  )

  expect_equal(loaded_calls$calls, fresh_calls$calls, tolerance = 0)
  expect_equal(
    loaded[order(gene_id)], fresh[order(gene_id)], tolerance = 0
  )
})


test_that("GeneScan stage 2 rejects missing and obsolete gene bundles", {
  fixture <- .make_gene_bundle_fixture()
  fixture$genes <- fixture$genes[1L]
  fixture$abc <- fixture$abc[TargetGene == "GENE1"]
  knockoff_dir <- tempfile("genescan-invalid-bundle-")
  dir.create(knockoff_dir)
  on.exit(unlink(knockoff_dir, recursive = TRUE, force = TRUE), add = TRUE)

  .run_gene_bundle_layout(
    fixture, batch_index = list(1L), ids = fixture$samples,
    stage1_only = TRUE, save_knockoff = TRUE,
    knockoff_dir = knockoff_dir
  )
  path <- file.path(knockoff_dir, "gene_GENE1_ko.rds")
  valid <- readRDS(path)

  obsolete <- valid
  obsolete$bundle_schema_version <- NULL
  saveRDS(obsolete, path)
  expect_error(
    .run_gene_bundle_layout(
      fixture, list(1L), fixture$samples,
      load_knockoff = TRUE, knockoff_dir = knockoff_dir,
      forbid_generation = TRUE
    ),
    "obsolete schema"
  )

  incomplete <- valid
  incomplete$enhancers <- NULL
  saveRDS(incomplete, path)
  expect_error(
    .run_gene_bundle_layout(
      fixture, list(1L), fixture$samples,
      load_knockoff = TRUE, knockoff_dir = knockoff_dir,
      forbid_generation = TRUE
    ),
    "enhancer"
  )

  unlink(path)
  expect_error(
    .run_gene_bundle_layout(
      fixture, list(1L), fixture$samples,
      load_knockoff = TRUE, knockoff_dir = knockoff_dir,
      forbid_generation = TRUE
    ),
    "Required saved knockoff file not found"
  )
})


test_that("BIGKnock bundles load both gene and enhancer payloads", {
  fixture <- .make_gene_bundle_fixture()
  fixture$genes <- fixture$genes[1L]
  fixture$abc <- fixture$abc[TargetGene == "GENE1"]
  knockoff_dir <- tempfile("bigknock-bundles-")
  dir.create(knockoff_dir)
  on.exit(unlink(knockoff_dir, recursive = TRUE, force = TRUE), add = TRUE)

  expect_null(.run_gene_bundle_layout(
    fixture, list(1L), fixture$samples, use_glmm = TRUE,
    stage1_only = TRUE, save_knockoff = TRUE,
    knockoff_dir = knockoff_dir
  ))
  path <- file.path(knockoff_dir, "gene_GENE1_ko.rds")
  bundle <- readRDS(path)
  expect_identical(bundle$bundle_schema_version, 2L)
  expect_true(is.array(bundle$G_gene_buffer_knockoff))
  expect_true(length(bundle$enhancers) > 0L)
  expect_true(any(vapply(
    bundle$enhancers, .gene_bundle_has_knockoff, logical(1)
  )))

  reversed_ids <- rev(fixture$samples)
  loaded_calls <- new.env(parent = emptyenv())
  fresh_calls <- new.env(parent = emptyenv())
  loaded <- .run_gene_bundle_layout(
    fixture, list(1L), reversed_ids, use_glmm = TRUE,
    load_knockoff = TRUE, knockoff_dir = knockoff_dir,
    forbid_generation = TRUE, observed = loaded_calls
  )
  fresh <- .run_gene_bundle_layout(
    fixture, list(1L), reversed_ids, use_glmm = TRUE,
    observed = fresh_calls
  )
  expect_equal(loaded_calls$calls, fresh_calls$calls, tolerance = 0)
  expect_equal(loaded, fresh, tolerance = 0)
})
