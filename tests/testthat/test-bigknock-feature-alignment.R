make_bigknock_alignment_toy <- function() {
  n <- 300L
  ld_a <- rep(c(0, 1, 2, 1), length.out = n)
  x <- cbind(
    ld_a = ld_a,
    ld_b = ld_a,
    independent = rep(c(2, 0, 1, 0, 2, 1), length.out = n),
    outside = rep(c(0, 2, 0, 1, 2, 1, 1, 0, 2, 1), length.out = n)
  )
  positions <- c(100, 110, 180, 400)
  metadata <- data.frame(
    chr = 1,
    variant_id = paste0("rs", seq_along(positions)),
    pos = positions,
    a1 = "A",
    a2 = "G",
    coded_allele = "A",
    stringsAsFactors = FALSE
  )
  list(x = x, positions = positions, metadata = metadata)
}


test_that("BIGKnock leverage sampling handles one retained variant", {
  x <- Matrix::Matrix(matrix(rep(c(0, 1, 2, 1), 75L), ncol = 1L),
                      sparse = TRUE)
  prob <- KnockoffPipeline:::.bigknock_shrinkage_prob(x)

  expect_length(prob, nrow(x))
  expect_equal(prob, rep(1 / nrow(x), nrow(x)))
  expect_equal(sum(prob), 1)
})


test_that("BIGKnock high-LD gene-buffer columns remain aligned", {
  toy <- make_bigknock_alignment_toy()
  set.seed(11)
  selected <- KnockoffPipeline:::.bigknock_prepare_region(
    toy$x, toy$positions, region_start = 90, region_end = 250,
    LD_filter = 0.75, label = "gene buffer"
  )

  # The perfectly correlated pair is represented once, while the independent
  # target feature remains.  Returned indices refer to the original columns.
  expect_equal(length(intersect(selected$target_input_index, 1:2)), 1L)
  expect_true(3L %in% selected$target_input_index)
  expect_equal(
    selected$target_positions,
    toy$positions[selected$target_input_index]
  )

  M <- 2L
  set.seed(11)
  generated <- KnockoffPipeline:::Knockoffgeneration.gene.buffer(
    toy$x, positions = toy$positions,
    gene_buffer_start = 90, gene_buffer_end = 250,
    M = M, LD.filter = 0.75, return_details = TRUE
  )
  expect_equal(generated$selected_input_index, selected$target_input_index)
  aligned <- KnockoffPipeline:::.bigknock_align_generated_features(
    generated, toy$x, toy$positions, M,
    input_metadata = toy$metadata, label = "gene buffer"
  )

  expect_equal(aligned$feature_index, selected$target_input_index)
  expect_equal(aligned$positions, toy$positions[aligned$feature_index])
  expect_equal(aligned$original, toy$x[, aligned$feature_index, drop = FALSE])
  expect_equal(
    aligned$metadata$variant_id,
    toy$metadata$variant_id[aligned$feature_index]
  )
  expect_identical(
    dim(aligned$knockoff),
    c(M, nrow(toy$x), length(aligned$feature_index))
  )

  wrong_dim <- generated
  wrong_dim$knockoff <- array(0, dim = c(M, nrow(toy$x), 1L))
  expect_error(
    KnockoffPipeline:::.bigknock_align_generated_features(
      wrong_dim, toy$x, toy$positions, M, input_metadata = toy$metadata
    ),
    "dimensions do not match"
  )
  wrong_identity <- generated
  wrong_identity$positions[1L] <- wrong_identity$positions[1L] + 1
  expect_error(
    KnockoffPipeline:::.bigknock_align_generated_features(
      wrong_identity, toy$x, toy$positions, M,
      input_metadata = toy$metadata
    ),
    "positions do not match"
  )
})


test_that("BIGKnock preserves two target columns when LD reduction would leave one", {
  toy <- make_bigknock_alignment_toy()
  set.seed(17)
  selected <- KnockoffPipeline:::.bigknock_prepare_region(
    toy$x, toy$positions, region_start = 90, region_end = 150,
    LD_filter = 0.75, label = "enhancer"
  )

  expect_equal(selected$target_input_index, 1:2)
  set.seed(17)
  generated <- KnockoffPipeline:::Knockoffgeneration.enhancer(
    toy$x, positions = toy$positions,
    enhancer_start = 90, enhancer_end = 150,
    M = 3L, LD.filter = 0.75, return_details = TRUE
  )
  expect_equal(generated$selected_input_index, selected$target_input_index)
  aligned <- KnockoffPipeline:::.bigknock_align_generated_features(
    generated, toy$x, toy$positions, 3L, label = "enhancer"
  )

  expect_equal(aligned$positions, toy$positions[aligned$feature_index])
  expect_equal(aligned$original, toy$x[, aligned$feature_index, drop = FALSE])
  expect_identical(dim(aligned$knockoff), c(3L, nrow(toy$x), 2L))
})


test_that("BIGKnock validates and persists selected gene-buffer identity", {
  toy <- make_bigknock_alignment_toy()
  set.seed(23)
  selected <- KnockoffPipeline:::.bigknock_prepare_region(
    toy$x, toy$positions, region_start = 90, region_end = 250,
    LD_filter = 0.75, label = "gene buffer"
  )
  M <- 2L
  generated <- list(
    knockoff = array(
      seq_len(M * nrow(toy$x) * length(selected$target_input_index)),
      dim = c(M, nrow(toy$x), length(selected$target_input_index))
    ),
    selected_input_index = selected$target_input_index,
    positions = selected$target_positions
  )
  context <- KnockoffPipeline:::.make_knockoff_context(
    test_type = "Gene_Centric_GLMM", M = M, genome_build = "hg19",
    variant_metadata = toy$metadata, reference_id = "toy",
    construction_id = "BIGKnock-alignment-test"
  )
  ids <- sprintf("id%03d", seq_len(nrow(toy$x)))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  first <- KnockoffPipeline:::.bigknock_gene_ko_load_or_gen(
    load_knockoff = FALSE, save_knockoff = TRUE,
    knockoff_file = path, matched_ids = ids,
    original_matrix = toy$x, input_positions = toy$positions,
    input_metadata = toy$metadata, gen_fun = function() generated,
    context = context
  )
  saved <- readRDS(path)
  expect_identical(saved$feature_schema_version, 1L)
  expect_equal(saved$selected_input_index, first$feature_index)
  expect_identical(saved$feature_fingerprint, first$feature_fingerprint)

  loaded <- KnockoffPipeline:::.bigknock_gene_ko_load_or_gen(
    load_knockoff = TRUE, save_knockoff = FALSE,
    knockoff_file = path, matched_ids = rev(ids),
    original_matrix = toy$x[nrow(toy$x):1L, , drop = FALSE],
    input_positions = toy$positions, input_metadata = toy$metadata,
    gen_fun = function() stop("must not regenerate"), context = context
  )
  expect_equal(loaded$feature_index, first$feature_index)
  expect_equal(loaded$positions, first$positions)
  expect_equal(
    loaded$knockoff,
    first$knockoff[, nrow(toy$x):1L, , drop = FALSE]
  )

  saved$feature_fingerprint[1L] <- "corrupt"
  saveRDS(saved, path)
  expect_error(
    KnockoffPipeline:::.bigknock_gene_ko_load_or_gen(
      load_knockoff = TRUE, save_knockoff = FALSE,
      knockoff_file = path, matched_ids = ids,
      original_matrix = toy$x, input_positions = toy$positions,
      input_metadata = toy$metadata, gen_fun = function() NULL,
      context = context
    ),
    "feature identity is incompatible"
  )
})
