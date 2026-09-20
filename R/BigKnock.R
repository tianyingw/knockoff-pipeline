utils::globalVariables(c('G_gene_buffer_surround','LD.filter',
                         'surround.region','G_gene_buffer','G_EnhancerAll','p_EnhancerAll',
                         'pos_gene_buffer','G_Enhancer','n','G_enhancer_surround','pos_enhancer'))


# BIGKnock performs an additional MAC/LD reduction before constructing a
# knockoff.  Keep that reduction in one helper and carry the original column
# indices through it so the observed matrix and its metadata can be reduced in
# exactly the same way as the generated knockoff.
.bigknock_pick_representatives <- function(indices, clusters) {
  if (length(indices) == 0L) return(integer(0))
  unname(vapply(unique(clusters[indices]), function(cluster_id) {
    candidates <- indices[clusters[indices] == cluster_id]
    candidates[sample.int(length(candidates), 1L)]
  }, integer(1)))
}


.bigknock_shrinkage_prob <- function(X) {
  nr <- nrow(X)
  nc <- ncol(X)
  if (nr < 2L || nc < 1L)
    stop("BIGKnock leverage sampling requires at least two samples and one variant.",
         call. = FALSE)

  # irlba requires a positive truncated rank strictly below the smaller
  # matrix dimension.  LD reduction can legitimately leave one surrounding
  # variant; in that case its leverage is constant and the shrinkage mixture
  # reduces to uniform sampling.
  if (nc == 1L) return(rep(1 / nr, nr))
  nv <- min(
    floor(sqrt(nc * log(nc))),
    nr - 1L,
    nc - 1L
  )
  if (nv < 1L) return(rep(1 / nr, nr))

  smaller_dim <- min(nr, nc)
  u <- if (nv >= smaller_dim / 2) {
    # irlba warns and offers no computational advantage when most singular
    # vectors are requested; this branch is reached only for a small smaller
    # dimension in normal gene-region inputs.
    base::svd(as.matrix(X), nu = nv, nv = 0L)$u
  } else {
    irlba(X, nv = nv)$u
  }
  leverage <- rowSums(u^2)
  leverage_total <- sum(leverage)
  if (!is.finite(leverage_total) || leverage_total <= 0)
    stop("BIGKnock leverage scores are not finite and positive.", call. = FALSE)
  0.5 * (leverage / leverage_total) + 0.5 * rep(1 / nr, nr)
}


.bigknock_prepare_region <- function(
  X, positions, region_start, region_end, LD_filter = 0.75,
  min_mac = 25, label = "region", min_target_variants = 2L
) {
  # Keep the full surround sparse; converting n-by-p biobank matrices to base
  # matrices before MAC filtering can require many unnecessary gigabytes.
  X <- Matrix::Matrix(X, sparse = TRUE)
  positions <- as.numeric(positions)
  if (ncol(X) != length(positions))
    stop("BIGKnock ", label, " positions do not match genotype columns.",
         call. = FALSE)
  if (ncol(X) == 0L || nrow(X) == 0L)
    stop("BIGKnock ", label, " genotype matrix is empty.", call. = FALSE)
  if (anyNA(positions))
    stop("BIGKnock ", label,
         " positions must be non-missing for feature alignment.",
         call. = FALSE)

  input_index <- seq_len(ncol(X))
  X[X < 0 | X > 2] <- NA_real_
  if (anyNA(X)) {
    means <- colMeans(X, na.rm = TRUE)
    missing <- which(is.na(X), arr.ind = TRUE)
    X[missing] <- means[missing[, 2L]]
  }

  raw_maf <- colMeans(X) / 2
  minor_maf <- pmin(raw_maf, 1 - raw_maf)
  minor_mac <- 2 * nrow(X) * minor_maf
  variances <- colMeans(X^2) - colMeans(X)^2
  keep <- which(
    is.finite(minor_maf) & minor_maf > 0 &
      is.finite(minor_mac) & minor_mac >= min_mac &
      is.finite(variances) & variances != 0
  )
  if (length(keep) <= 1L) {
    stop("BIGKnock ", label,
         " has <=1 variant after MAC and variance filtering.", call. = FALSE)
  }

  X <- X[, keep, drop = FALSE]
  positions <- positions[keep]
  input_index <- input_index[keep]
  ord <- order(positions)
  X <- X[, ord, drop = FALSE]
  positions <- positions[ord]
  input_index <- input_index[ord]

  maf <- colMeans(X) / 2
  flip <- is.finite(maf) & maf > 0.5
  if (any(flip)) X[, flip] <- 2 - X[, flip, drop = FALSE]

  initial_target <- positions >= region_start & positions <= region_end
  if (sum(initial_target) < min_target_variants) {
    stop("BIGKnock ", label,
         " has fewer than ", min_target_variants,
         " target variant(s) after MAC and variance filtering.",
         call. = FALSE)
  }

  repeat {
    if (ncol(X) <= 1L) break
    cor_X <- .kp_sparse_cov_cor(
      Matrix::Matrix(X, sparse = TRUE), need_cov = FALSE, need_cor = TRUE
    )$cor
    diag(cor_X) <- 0
    max_corr <- suppressWarnings(max(abs(cor_X), na.rm = TRUE))
    if (!is.finite(max_corr) || max_corr < LD_filter) break

    clusters <- stats::cutree(
      stats::hclust(stats::as.dist(1 - abs(cor_X)), method = "complete"),
      h = 1 - LD_filter
    )
    in_region <- positions >= region_start & positions <= region_end
    region_reps <- .bigknock_pick_representatives(which(in_region), clusters)
    if (length(region_reps) == 0L)
      stop("BIGKnock ", label,
           " has no target variant after LD filtering.", call. = FALSE)
    # The downstream gene statistic requires at least two target variants.
    # Match the original BIGKnock rule: if representative filtering would
    # collapse the target to one column, retain the current pre-filter matrix.
    if (length(region_reps) <= 1L) break

    outside <- which(!in_region)
    outside <- outside[!clusters[outside] %in% clusters[region_reps]]
    outside_reps <- .bigknock_pick_representatives(outside, clusters)
    selected <- sort(c(region_reps, outside_reps))
    selected <- selected[order(positions[selected])]
    if (length(selected) >= ncol(X)) break

    X <- X[, selected, drop = FALSE]
    positions <- positions[selected]
    input_index <- input_index[selected]
  }

  target <- positions >= region_start & positions <= region_end
  if (!any(target))
    stop("BIGKnock ", label,
         " has no target variant after feature selection.", call. = FALSE)

  X <- Matrix::Matrix(X, sparse = TRUE)
  colnames(X) <- as.character(positions)
  list(
    surround_matrix = X,
    surround_positions = positions,
    surround_input_index = input_index,
    target_matrix = X[, target, drop = FALSE],
    target_positions = positions[target],
    target_input_index = input_index[target]
  )
}


.bigknock_align_generated_features <- function(
  generated, original_matrix, input_positions, M,
  input_metadata = NULL, label = "region"
) {
  required <- c("knockoff", "selected_input_index", "positions")
  if (!is.list(generated) || !all(required %in% names(generated)))
    stop("Generated BIGKnock ", label,
         " object is missing feature-selection metadata.", call. = FALSE)

  idx <- generated$selected_input_index
  if (!is.numeric(idx) || length(idx) == 0L || anyNA(idx) ||
      any(idx != as.integer(idx)) || anyDuplicated(idx) ||
      any(idx < 1L | idx > ncol(original_matrix))) {
    stop("Generated BIGKnock ", label,
         " object has invalid selected column indices.", call. = FALSE)
  }
  idx <- as.integer(idx)
  expected_positions <- as.numeric(input_positions[idx])
  reported_positions <- as.numeric(generated$positions)
  if (!identical(expected_positions, reported_positions))
    stop("Generated BIGKnock ", label,
         " feature positions do not match the selected original columns.",
         call. = FALSE)

  arr <- generated$knockoff
  expected_dim <- c(as.integer(M), nrow(original_matrix), length(idx))
  if (length(dim(arr)) != 3L || !identical(as.integer(dim(arr)), expected_dim))
    stop("Generated BIGKnock ", label,
         " knockoff dimensions do not match M, samples, and selected features.",
         call. = FALSE)

  metadata <- NULL
  fingerprint <- NULL
  if (!is.null(input_metadata)) {
    if (nrow(input_metadata) != ncol(original_matrix))
      stop("BIGKnock ", label,
           " metadata do not match genotype columns.", call. = FALSE)
    metadata <- input_metadata[idx, , drop = FALSE]
    fingerprint <- .variant_fingerprint(metadata)
  }

  list(
    knockoff = arr,
    original = original_matrix[, idx, drop = FALSE],
    positions = expected_positions,
    metadata = metadata,
    feature_index = idx,
    feature_fingerprint = fingerprint
  )
}


.bigknock_gene_ko_load_or_gen <- function(
  load_knockoff, save_knockoff, knockoff_file, matched_ids,
  original_matrix, input_positions, input_metadata, gen_fun, context,
  knockoff_object = NULL
) {
  if (isTRUE(load_knockoff)) {
    if (is.null(knockoff_file) || !file.exists(knockoff_file))
      stop("Required saved knockoff file not found: ", knockoff_file)
    ko_obj <- if (is.null(knockoff_object))
      readRDS(knockoff_file) else knockoff_object
    .assert_knockoff_context(ko_obj$context, context, knockoff_file)
    if (!identical(ko_obj$feature_schema_version, 1L) ||
        is.null(ko_obj$selected_input_index) ||
        is.null(ko_obj$feature_fingerprint)) {
      stop("Saved BIGKnock feature-selection metadata are missing or obsolete: ",
           knockoff_file, ". Regenerate the knockoff.", call. = FALSE)
    }

    current_ids <- as.character(matched_ids)
    saved_ids <- as.character(ko_obj$sample_ids)
    if (.need_regenerate_samples(current_ids, saved_ids)) {
      stop(
        "Saved knockoff cannot be reused because its sample-ID set differs ",
        "from the current run: ", knockoff_file, ". Regenerate knockoffs for ",
        "the exact analysis sample set."
      )
    }
    row_map <- match(current_ids, saved_ids)
    saved_arr <- ko_obj$G_gene_buffer_knockoff
    if (length(dim(saved_arr)) != 3L ||
        dim(saved_arr)[1L] != context$M ||
        dim(saved_arr)[2L] != length(saved_ids)) {
      stop("Saved BIGKnock gene-buffer knockoff dimensions are incompatible: ",
           knockoff_file, ". Regenerate the knockoff.", call. = FALSE)
    }
    if (!identical(row_map, seq_along(current_ids))) {
      saved_arr <- saved_arr[, row_map, , drop = FALSE]
    }
    generated <- list(
      knockoff = saved_arr,
      selected_input_index = ko_obj$selected_input_index,
      positions = ko_obj$snp_pos
    )
    aligned <- .bigknock_align_generated_features(
      generated, original_matrix, input_positions, context$M,
      input_metadata = input_metadata, label = "gene buffer"
    )
    if (!identical(ko_obj$feature_fingerprint,
                   aligned$feature_fingerprint)) {
      stop("Saved BIGKnock gene-buffer feature identity is incompatible: ",
           knockoff_file, ". Regenerate the knockoff.", call. = FALSE)
    }
    return(aligned)
  }

  aligned <- .bigknock_align_generated_features(
    gen_fun(), original_matrix, input_positions, context$M,
    input_metadata = input_metadata, label = "gene buffer"
  )
  if (isTRUE(save_knockoff) && !is.null(knockoff_file)) {
    dir.create(dirname(knockoff_file), recursive = TRUE, showWarnings = FALSE)
    .atomic_save_rds(
      list(
        feature_schema_version = 1L,
        G_gene_buffer_knockoff = aligned$knockoff,
        sample_ids = matched_ids,
        snp_pos = aligned$positions,
        selected_input_index = aligned$feature_index,
        feature_fingerprint = aligned$feature_fingerprint,
        context = context
      ),
      path = knockoff_file
    )
  }
  aligned
}

# Shared BIGKnock GLMM matrices ---------------------------------------------
# X, fitted values, sparseSigma and theta are phenotype-level quantities.
# Computing these two small generalized inverses for every gene needlessly
# repeats the expensive sparseSigma solve.  Keep the calculation here so the
# pipeline can do it once while direct callers still have the same fallback.
.bigknock_svd_inverse <- function(x, label, ridge = 1e-2) {
  tryCatch({
    s <- svd(x)
    keep <- s$d > 1e-10 * s$d[1L]
    s$v[, keep, drop = FALSE] %*%
      (s$d[keep]^(-1) * t(s$u[, keep, drop = FALSE]))
  }, error = function(e) {
    message("  [BigKnock] ", label, " SVD failed, using ridge fallback")
    solve(x + diag(ridge, nrow(x)))
  })
}


.bigknock_glmm_precompute <- function(result.null.model, sparseSigma) {
  if (is.null(result.null.model) || is.null(result.null.model$X))
    stop("BIGKnock GLMM precomputation requires a fitted null model with X.",
         call. = FALSE)
  if (is.null(sparseSigma))
    stop("BIGKnock GLMM precomputation requires sparseSigma.", call. = FALSE)

  X <- result.null.model$X
  if (nrow(sparseSigma) != nrow(X) || ncol(sparseSigma) != nrow(X))
    stop("sparseSigma dimensions do not match the null-model samples.",
         call. = FALSE)

  invSigma_X <- .safe_solve_sparse(sparseSigma, X)
  C <- .bigknock_svd_inverse(
    t(X) %*% invSigma_X, label = "C matrix"
  )

  outcome <- as.character(result.null.model$traitType)[1L]
  if (!outcome %in% c("C", "D"))
    stop("BIGKnock null-model traitType must be 'C' or 'D'.", call. = FALSE)
  if (outcome == "D") {
    mu <- as.vector(result.null.model$fitted.values)
    if (length(mu) != nrow(X))
      stop("Null-model fitted values do not match X.", call. = FALSE)
    v <- mu * (1 - mu)
  } else {
    theta <- as.numeric(result.null.model$theta)[1L]
    if (!is.finite(theta) || theta <= 0)
      stop("Continuous-trait BIGKnock requires a positive theta.",
           call. = FALSE)
    v <- 1 / theta
  }
  inv_vX <- .bigknock_svd_inverse(
    t(X) %*% (v * X), label = "weighted-X matrix"
  )

  list(C = C, inv_vX = inv_vX, outcome = outcome)
}


.validate_bigknock_glmm_precompute <- function(precomputed, result.null.model) {
  required <- c("C", "inv_vX", "outcome")
  if (!is.list(precomputed) ||
      !all(required %in% names(precomputed))) {
    stop("Invalid BIGKnock GLMM precomputation object.", call. = FALSE)
  }
  p <- ncol(result.null.model$X)
  expected_outcome <- as.character(result.null.model$traitType)[1L]
  if (!identical(dim(precomputed$C), c(p, p)) ||
      !identical(dim(precomputed$inv_vX), c(p, p)) ||
      !identical(as.character(precomputed$outcome)[1L], expected_outcome)) {
    stop("BIGKnock GLMM precomputation is incompatible with the null model.",
         call. = FALSE)
  }
  precomputed
}


# GeneScan3D.UKB.GLMM.KnockoffGeneration ------------------------------------
# Changes vs original:
#  * save/load one atomic gene-buffer + enhancer knockoff bundle.
#  * Knockoff validation: check both row count AND column count.
#  * stage1_only=TRUE: finish and save the complete bundle, then return NULL.
#  * NULL null model supported when stage1_only=TRUE.
# ----------------------------------------------------------------------------
GeneScan3D.UKB.GLMM.KnockoffGeneration <- function(
  G_gene_buffer_surround,
  variants_gene_buffer_surround,
  gene_buffer.pos,
  R                             = 0,
  G_EnhancerAll_surround        = NULL,
  variants_EnhancerAll_surround = NULL,
  p_EnhancerAll_surround        = NULL,
  Enhancer.pos                  = NULL,
  p.EnhancerAll                 = NULL,
  window.size                   = NULL,
  result.null.model             = NULL,
  M                             = 5,
  sparseSigma                   = NULL,
  ratio                         = NULL,
  MAC.threshold                 = 10,
  MAF.threshold                 = 0.01,
  Gsub.id                       = NULL,
  variant_metadata_gene_buffer_surround = NULL,
  variant_metadata_EnhancerAll_surround = NULL,
  genome_build                  = NULL,
  reference_id                  = NULL,
  enhancer_reference_id         = NULL,
  knockoff_seed                 = NULL,
  save_knockoff                 = FALSE,
  load_knockoff                 = FALSE,
  knockoff_file                 = NULL,
  knockoff_sample_ids           = NULL,
  stage1_only                   = FALSE,
  glmm_precomputed              = NULL
) {
  impute.method <- "fixed"

  if (is.null(variant_metadata_gene_buffer_surround) ||
      nrow(variant_metadata_gene_buffer_surround) != ncol(G_gene_buffer_surround)) {
    stop("Variant metadata from the matching PLINK .bim rows must be supplied for every gene-surround column.")
  }

  persistent_knockoff <- isTRUE(save_knockoff) || isTRUE(load_knockoff)
  if (persistent_knockoff && R != 0L &&
      (is.null(variant_metadata_EnhancerAll_surround) ||
       is.null(G_EnhancerAll_surround) ||
       nrow(variant_metadata_EnhancerAll_surround) !=
         ncol(G_EnhancerAll_surround))) {
    stop(paste0(
      "Variant metadata from the matching PLINK .bim rows must be supplied ",
      "for every enhancer-surround column when saving or loading knockoffs."
    ))
  }

  # ---- Sample matching (supports NULL null model for stage1) ---------------
  if (is.null(result.null.model)) {
    if (!isTRUE(stage1_only))
      stop("result.null.model is NULL but stage1_only is not TRUE.")
    # Stage 1: use all rows in Gsub.id / sequential order
    n           <- nrow(G_gene_buffer_surround)
    match.index <- seq_len(n)
  } else {
    mu    <- as.vector(result.null.model$fitted.values)
    Y.res <- as.vector(result.null.model$residuals)
    n     <- length(mu)
    if (length(Gsub.id) == 0) {
      match.index <- match(result.null.model$sampleID,
                           seq_len(nrow(G_gene_buffer_surround)))
    } else {
      match.index <- match(result.null.model$sampleID, Gsub.id)
    }
    if (mean(is.na(match.index)) > 0)
      warning(sprintf("Some individuals not matched with genotype. Rate = %f",
                      mean(is.na(match.index))), call. = FALSE)
  }

  # IDs in matched row order (used as saved sample_ids in the knockoff file)
  matched_ids <- if (!is.null(Gsub.id)) Gsub.id[match.index] else match.index

  # ---- QC: gene buffer surround -------------------------------------------
  if (!identical(match.index, seq_len(nrow(G_gene_buffer_surround)))) {
    G_gene_buffer_surround <-
      G_gene_buffer_surround[match.index, , drop = FALSE]
  }
  G_gene_buffer_surround <- Matrix::Matrix(G_gene_buffer_surround)
  G_gene_buffer_surround[G_gene_buffer_surround == -9 |
                         G_gene_buffer_surround ==  9] <- NA
  missing_mask <- is.na(G_gene_buffer_surround)
  N_MISS    <- sum(missing_mask)
  MISS.freq <- .kp_col_means(missing_mask)
  rm(missing_mask)
  if (N_MISS > 0) {
    warning(sprintf("Missing genotype rate = %f. Imputation applied.",
                    N_MISS / nrow(G_gene_buffer_surround) / ncol(G_gene_buffer_surround)),
            call. = FALSE)
    G_gene_buffer_surround <- Impute(G_gene_buffer_surround, impute.method)
  }
  MAF       <- .kp_col_means(G_gene_buffer_surround) / 2
  flip_to_minor <- MAF > 0.5 & !is.na(MAF)
  G_gene_buffer_surround[, flip_to_minor] <-
    2 - G_gene_buffer_surround[, flip_to_minor, drop = FALSE]
  MAF       <- .kp_col_means(G_gene_buffer_surround) / 2
  variance  <- .kp_col_means(G_gene_buffer_surround^2) -
    .kp_col_means(G_gene_buffer_surround)^2
  minor_mac <- 2 * nrow(G_gene_buffer_surround) * MAF
  SNP.index <- which(MAF > 0 & variance > 0 & !is.na(MAF) & MISS.freq < 0.1)
  if (length(SNP.index) <= 1) {
    warning("Number of variants passing QC in gene buffer surround is <=1", call. = FALSE)
    return(NULL)
  }
  G_gene_buffer_surround               <- Matrix::Matrix(G_gene_buffer_surround[, SNP.index])
  variants_gene_buffer_surround_filter <- variants_gene_buffer_surround[SNP.index]
  variant_metadata_filter <- as.data.frame(
    variant_metadata_gene_buffer_surround[SNP.index, , drop = FALSE],
    stringsAsFactors = FALSE
  )
  counted <- as.character(variant_metadata_filter$counted_allele)
  opposite <- ifelse(
    counted == as.character(variant_metadata_filter$a1),
    as.character(variant_metadata_filter$a2),
    as.character(variant_metadata_filter$a1)
  )
  variant_metadata_filter$coded_allele <- ifelse(
    flip_to_minor[SNP.index], opposite, counted
  )
  colnames(G_gene_buffer_surround)     <-
    extract_position_universal(colnames(G_gene_buffer_surround))

  # Positions within the gene buffer (subset of surround)
  positions_gene_buffer <- variants_gene_buffer_surround_filter[
    variants_gene_buffer_surround_filter <= gene_buffer.pos[2] &
    variants_gene_buffer_surround_filter >= gene_buffer.pos[1]
  ]
  if (length(positions_gene_buffer) == 0) return(NULL)
  in_gene_buffer <- variants_gene_buffer_surround_filter <= gene_buffer.pos[2] &
    variants_gene_buffer_surround_filter >= gene_buffer.pos[1]
  target_mac <- minor_mac[SNP.index][in_gene_buffer]
  if (sum(is.finite(target_mac) & target_mac >= 25) <= 1L) {
    warning("Gene buffer has <=1 target variant after BIGKnock MAC and variance filtering; skipping gene.",
            call. = FALSE)
    return(NULL)
  }
  current_context <- if (persistent_knockoff) {
    .make_knockoff_context(
      test_type = "Gene_Centric_GLMM",
      M = M,
      genome_build = genome_build,
      # Construction uses the complete post-QC surround matrix, so safe reuse
      # must fingerprint every predictor, not only the returned buffer columns.
      variant_metadata = variant_metadata_filter,
      reference_id = reference_id,
      construction_id = "BIGKnock-gene-buffer-v5;corrected_skip_index;impute=fixed;gene_target_flank=5000;neighbor_bp=100000;source_flank=105000;MAC_min=25;LD_filter=0.75;corr_base=0.05;thres_ultrarare=25;retain_if_target_reps_le_1",
      random_seed = knockoff_seed
    )
  } else {
    list(M = as.integer(M))
  }

  # A persistent BIGKnock file is a single per-gene bundle.  Read it once so
  # the gene and every enhancer are validated against the same manifest.
  knockoff_bundle <- NULL
  enhancer_row_map <- NULL
  current_enhancer_reference_id <-
    .normalize_knockoff_reference_id(enhancer_reference_id)
  if (isTRUE(load_knockoff)) {
    loaded_bundle <- .load_gene_knockoff_bundle(
      knockoff_file = knockoff_file,
      matched_ids = matched_ids,
      enhancer_reference_id = enhancer_reference_id,
      n_enhancers = as.integer(R),
      method_label = "BIGKnock"
    )
    knockoff_bundle <- loaded_bundle$object
    enhancer_row_map <- loaded_bundle$row_map
  }

  # ---- Gene buffer knockoff: save / load / generate -----------------------
  # BIGKnock performs an additional MAC/LD selection.  Its selected column
  # indices are saved and used for both the original matrix and the knockoff.
  gene_ko <- .bigknock_gene_ko_load_or_gen(
    load_knockoff     = load_knockoff,
    # The complete bundle is written atomically after all enhancer entries
    # have succeeded; never leave a gene-only partial checkpoint behind.
    save_knockoff     = FALSE,
    knockoff_file     = knockoff_file,
    matched_ids       = matched_ids,
    original_matrix   = G_gene_buffer_surround,
    input_positions   = variants_gene_buffer_surround_filter,
    input_metadata    = if (persistent_knockoff) variant_metadata_filter else NULL,
    gen_fun           = function() .with_local_seed(
      .derive_unit_seed(knockoff_seed, "gene_buffer"),
      function() {
        ko <- NULL
        invisible(capture.output(
          ko <- Knockoffgeneration.gene.buffer(
            G_gene_buffer_surround = G_gene_buffer_surround,  # surround matrix
            positions              = variants_gene_buffer_surround_filter,
            gene_buffer_start      = gene_buffer.pos[1],
            gene_buffer_end        = gene_buffer.pos[2],
            M                      = M,
            return_details         = TRUE
          )
        ))
        ko
      }
    ),
    context           = current_context,
    knockoff_object   = knockoff_bundle
  )
  G_gene_buffer_knockoff <- gene_ko$knockoff
  G_gene_buffer <- gene_ko$original
  positions_gene_buffer <- gene_ko$positions

  # ---- R enhancers (generated in stage 1, loaded in stage 2) --------------
  G_EnhancerAll          <- NULL
  p_EnhancerAll_out      <- integer(0)
  G_EnhancerAll_knockoff <- NULL
  R_input                <- R
  R                      <- 0L
  enhancer_entries <- vector("list", as.integer(R_input))
  names(enhancer_entries) <- if (R_input > 0L) {
    sprintf("enhancer_%04d", seq_len(as.integer(R_input)))
  } else {
    character(0)
  }

  if (R_input != 0) {
    if (is.null(G_EnhancerAll_surround) ||
        length(p_EnhancerAll_surround) != R_input ||
        length(variants_EnhancerAll_surround) !=
          ncol(G_EnhancerAll_surround) ||
        sum(as.integer(p_EnhancerAll_surround)) !=
          ncol(G_EnhancerAll_surround) ||
        is.null(Enhancer.pos) || nrow(Enhancer.pos) != R_input) {
      stop("Enhancer-surround inputs are inconsistent with R.", call. = FALSE)
    }
    enhancer_column_ends <- cumsum(as.integer(p_EnhancerAll_surround))

    for (r in seq_len(R_input)) {
      # Slice enhancer surround columns from the batch matrix
      enhancer_column_start <- if (r == 1L) 1L else
        enhancer_column_ends[r - 1L] + 1L
      enhancer_column_index <- seq.int(
        enhancer_column_start, enhancer_column_ends[r]
      )
      G_Enh_surround <- G_EnhancerAll_surround[,
        enhancer_column_index, drop = FALSE]
      pos_Enh_surround <- variants_EnhancerAll_surround[
        enhancer_column_index]
      metadata_Enh_surround <- if (persistent_knockoff) {
        as.data.frame(
          variant_metadata_EnhancerAll_surround[
            enhancer_column_index, , drop = FALSE
          ],
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
      enhancer_coords <- c(
        start = as.numeric(Enhancer.pos[r, 1L]),
        end = as.numeric(Enhancer.pos[r, 2L])
      )

      # QC: enhancer surround
      if (!identical(match.index, seq_len(nrow(G_Enh_surround)))) {
        G_Enh_surround <- G_Enh_surround[match.index, , drop = FALSE]
      }
      G_Enh_surround <- Matrix::Matrix(G_Enh_surround)
      G_Enh_surround[G_Enh_surround == -9 | G_Enh_surround == 9] <- NA
      missing_mask <- is.na(G_Enh_surround)
      N_MISS    <- sum(missing_mask)
      MISS.freq <- .kp_col_means(missing_mask)
      rm(missing_mask)
      if (N_MISS > 0) {
        warning(sprintf("Enhancer %d: missing rate = %f. Imputation applied.", r,
                        N_MISS / nrow(G_Enh_surround) / ncol(G_Enh_surround)),
                call. = FALSE)
        G_Enh_surround <- Impute(G_Enh_surround, impute.method)
      }
      MAF <- .kp_col_means(G_Enh_surround) / 2
      flip_to_minor_enhancer <- MAF > 0.5 & !is.na(MAF)
      if (any(flip_to_minor_enhancer)) {
        G_Enh_surround[, flip_to_minor_enhancer] <-
          2 - G_Enh_surround[, flip_to_minor_enhancer, drop = FALSE]
      }
      MAF       <- .kp_col_means(G_Enh_surround) / 2
      variance  <- .kp_col_means(G_Enh_surround^2) -
        .kp_col_means(G_Enh_surround)^2
      minor_mac <- 2 * nrow(G_Enh_surround) * MAF
      SNP.index <- which(MAF > 0 & variance > 0 & !is.na(MAF) & MISS.freq < 0.1)
      G_Enh_surround <- Matrix::Matrix(
        G_Enh_surround[, SNP.index, drop = FALSE]
      )
      pos_Enh_filter <- pos_Enh_surround[SNP.index]
      if (ncol(G_Enh_surround) > 0L) {
        colnames(G_Enh_surround) <-
          extract_position_universal(colnames(G_Enh_surround))
      }
      variant_metadata_Enh_filter <- if (persistent_knockoff) {
        metadata <- metadata_Enh_surround[SNP.index, , drop = FALSE]
        counted <- as.character(metadata$counted_allele)
        opposite <- ifelse(
          counted == as.character(metadata$a1),
          as.character(metadata$a2),
          as.character(metadata$a1)
        )
        metadata$coded_allele <- ifelse(
          flip_to_minor_enhancer[SNP.index], opposite, counted
        )
        metadata
      } else {
        NULL
      }

      in_enhancer <- pos_Enh_filter >= enhancer_coords[["start"]] &
        pos_Enh_filter <= enhancer_coords[["end"]]
      mac_pass <- is.finite(minor_mac[SNP.index]) & minor_mac[SNP.index] >= 25
      enhancer_skip_reason <- if (length(SNP.index) <= 1L) {
        "variants_passing_qc_le_1"
      } else if (sum(mac_pass) <= 1L ||
                 sum(mac_pass & in_enhancer) == 0L) {
        "no_testable_target_or_total_mac_variants_le_1"
      } else {
        NA_character_
      }
      enhancer_status <- if (is.na(enhancer_skip_reason)) "used" else "skipped"
      enhancer_seed <- .derive_unit_seed(knockoff_seed, "enhancer", r)
      enhancer_context <- if (persistent_knockoff) {
        .make_knockoff_context(
          test_type = "Gene_Centric_GLMM_Enhancer",
          M = M,
          genome_build = genome_build,
          variant_metadata = variant_metadata_Enh_filter,
          reference_id = enhancer_reference_id,
          construction_id = "BIGKnock-enhancer-v2;impute=fixed;neighbor_bp=50000;source_flank=50000;MAC_min=25;LD_filter=0.75;corr_base=0.05;thres_ultrarare=25;retain_if_target_reps_le_1",
          random_seed = enhancer_seed
        )
      } else {
        list(M = as.integer(M))
      }
      enhancer_entry <- list(
        coords = enhancer_coords,
        used = identical(enhancer_status, "used"),
        skipped = identical(enhancer_status, "skipped"),
        status = enhancer_status,
        skip_reason = enhancer_skip_reason,
        context = enhancer_context,
        G_Enhancer_knockoff = NULL,
        selected_input_index = integer(0),
        target_positions = numeric(0),
        feature_fingerprint = character(0)
      )

      saved_enhancer_entry <- NULL
      if (isTRUE(load_knockoff)) {
        entry_name <- sprintf("enhancer_%04d", r)
        saved_enhancer_entry <- knockoff_bundle$enhancers[[entry_name]]
        required_entry_fields <- c(
          "coords", "used", "skipped", "status", "skip_reason", "context",
          "G_Enhancer_knockoff", "selected_input_index", "target_positions",
          "feature_fingerprint"
        )
        if (!is.list(saved_enhancer_entry) ||
            !all(required_entry_fields %in% names(saved_enhancer_entry))) {
          stop(
            "Saved BIGKnock enhancer entry is incomplete: ", entry_name,
            " in ", knockoff_file, ". Regenerate the knockoff bundle.",
            call. = FALSE
          )
        }
        if (!identical(saved_enhancer_entry$coords, enhancer_coords) ||
            !identical(saved_enhancer_entry$used, enhancer_entry$used) ||
            !identical(saved_enhancer_entry$skipped,
                       enhancer_entry$skipped) ||
            !identical(saved_enhancer_entry$status, enhancer_status) ||
            !identical(saved_enhancer_entry$skip_reason,
                       enhancer_skip_reason)) {
          stop(
            "Saved BIGKnock enhancer coordinates or status are incompatible: ",
            entry_name, " in ", knockoff_file,
            ". Regenerate the knockoff bundle.", call. = FALSE
          )
        }
        .assert_knockoff_context(
          saved_enhancer_entry$context, enhancer_context,
          paste0(knockoff_file, " [", entry_name, "]")
        )
      }

      if (!is.na(enhancer_skip_reason)) {
        if (isTRUE(load_knockoff) &&
            (!is.null(saved_enhancer_entry$G_Enhancer_knockoff) ||
             !identical(saved_enhancer_entry$selected_input_index,
                        integer(0)) ||
             !identical(saved_enhancer_entry$target_positions,
                        numeric(0)) ||
             !identical(saved_enhancer_entry$feature_fingerprint,
                        character(0)))) {
          stop(
            "Saved skipped BIGKnock enhancer contains generated features: ",
            sprintf("enhancer_%04d", r), " in ", knockoff_file,
            ". Regenerate the knockoff bundle.", call. = FALSE
          )
        }
        if (identical(enhancer_skip_reason, "variants_passing_qc_le_1")) {
          warning(
            sprintf("Enhancer %d: variants passing QC <=1; skipping.", r),
            call. = FALSE
          )
        } else {
          warning(sprintf(
            paste0(
              "Enhancer %d: no testable target or <=1 total variant after ",
              "BIGKnock MAC and variance filtering; skipping."
            ),
            r
          ), call. = FALSE)
        }
        if (persistent_knockoff) {
          enhancer_entries[[r]] <- if (isTRUE(load_knockoff)) {
            saved_enhancer_entry
          } else {
            enhancer_entry
          }
        }
        next
      }

      if (isTRUE(load_knockoff)) {
        saved_arr <- saved_enhancer_entry$G_Enhancer_knockoff
        saved_ids <- as.character(knockoff_bundle$sample_ids)
        if (length(dim(saved_arr)) != 3L ||
            dim(saved_arr)[1L] != as.integer(M) ||
            dim(saved_arr)[2L] != length(saved_ids)) {
          stop(
            "Saved BIGKnock enhancer knockoff dimensions are incompatible: ",
            sprintf("enhancer_%04d", r), " in ", knockoff_file,
            ". Regenerate the knockoff bundle.", call. = FALSE
          )
        }
        if (!identical(enhancer_row_map, seq_along(enhancer_row_map))) {
          saved_arr <- saved_arr[, enhancer_row_map, , drop = FALSE]
        }
        enhancer_aligned <- .bigknock_align_generated_features(
          list(
            knockoff = saved_arr,
            selected_input_index =
              saved_enhancer_entry$selected_input_index,
            positions = saved_enhancer_entry$target_positions
          ),
          G_Enh_surround, pos_Enh_filter, M,
          input_metadata = variant_metadata_Enh_filter,
          label = paste0("enhancer ", r)
        )
        if (!identical(saved_enhancer_entry$feature_fingerprint,
                       enhancer_aligned$feature_fingerprint)) {
          stop(
            "Saved BIGKnock enhancer feature identity is incompatible: ",
            sprintf("enhancer_%04d", r), " in ", knockoff_file,
            ". Regenerate the knockoff bundle.", call. = FALSE
          )
        }
        enhancer_entries[[r]] <- saved_enhancer_entry
      } else {
        # Keep the original deterministic seed and feature-selection path.
        enhancer_generated <- .with_local_seed(
          enhancer_seed,
          function() {
            ko <- NULL
            invisible(capture.output(
              ko <- Knockoffgeneration.enhancer(
                G_enhancer_surround = G_Enh_surround,
                positions = pos_Enh_filter,
                enhancer_start = enhancer_coords[["start"]],
                enhancer_end = enhancer_coords[["end"]],
                M = M,
                return_details = TRUE
              )
            ))
            ko
          }
        )
        enhancer_aligned <- .bigknock_align_generated_features(
          enhancer_generated, G_Enh_surround, pos_Enh_filter, M,
          input_metadata = if (persistent_knockoff) {
            variant_metadata_Enh_filter
          } else {
            NULL
          },
          label = paste0("enhancer ", r)
        )
        if (persistent_knockoff) {
          enhancer_entry$G_Enhancer_knockoff <- enhancer_aligned$knockoff
          enhancer_entry$selected_input_index <- enhancer_aligned$feature_index
          enhancer_entry$target_positions <- enhancer_aligned$positions
          enhancer_entry$feature_fingerprint <-
            enhancer_aligned$feature_fingerprint
          enhancer_entries[[r]] <- enhancer_entry
        }
      }

      # Stage 1 only needs the aligned arrays in the bundle; avoid building a
      # second concatenated copy solely for an association test that will not
      # run in this stage.
      if (isTRUE(stage1_only)) next

      G_Enh_knockoff <- enhancer_aligned$knockoff
      G_enhancer <- enhancer_aligned$original
      G_EnhancerAll <- if (is.null(G_EnhancerAll)) G_enhancer else
        cbind(G_EnhancerAll, G_enhancer)
      p_EnhancerAll_out <- c(p_EnhancerAll_out, ncol(G_enhancer))
      G_EnhancerAll_knockoff <- if (is.null(G_EnhancerAll_knockoff)) {
        G_Enh_knockoff
      } else {
        abind::abind(G_EnhancerAll_knockoff, G_Enh_knockoff, along = 3L)
      }
      R <- R + 1L
    }
  }

  # Write only after the gene and every enhancer entry are complete.  This is
  # the sole write for a gene, so readers never observe a gene-only bundle.
  if (isTRUE(save_knockoff) && !isTRUE(load_knockoff) &&
      !is.null(knockoff_file)) {
    .atomic_save_rds(
      list(
        bundle_schema_version = 2L,
        feature_schema_version = 1L,
        G_gene_buffer_knockoff = gene_ko$knockoff,
        sample_ids = matched_ids,
        snp_pos = gene_ko$positions,
        selected_input_index = gene_ko$feature_index,
        feature_fingerprint = gene_ko$feature_fingerprint,
        context = current_context,
        enhancer_reference_id = current_enhancer_reference_id,
        enhancers = enhancer_entries
      ),
      path = knockoff_file
    )
  }

  # Stage 1 includes enhancer construction and the atomic bundle save above.
  if (isTRUE(stage1_only)) return(invisible(NULL))

  # Direct callers retain a local fallback; run_pipeline supplies the same
  # phenotype-level object once to every gene and batch.
  if (is.null(glmm_precomputed)) {
    glmm_precomputed <- .bigknock_glmm_precompute(
      result.null.model, sparseSigma
    )
  }
  glmm_precomputed <- .validate_bigknock_glmm_precompute(
    glmm_precomputed, result.null.model
  )
  C_pre       <- glmm_precomputed$C
  inv_vX_pre  <- glmm_precomputed$inv_vX
  outcome_pre <- glmm_precomputed$outcome

  # ---- Association tests ---------------------------------------------------
  tmp <- GeneScan3D.UKB.GLMM(
    G                    = G_gene_buffer,
    G.EnhancerAll        = G_EnhancerAll,
    R                    = R,
    p_Enhancer           = p_EnhancerAll_out,
    window.size          = window.size,
    pos                  = positions_gene_buffer,
    MAC.threshold        = MAC.threshold,
    MAF.threshold        = MAF.threshold,
    Gsub.id              = Gsub.id[match.index],
    result.null.model.GLMM = result.null.model,
    outcome              = outcome_pre,
    sparseSigma          = sparseSigma,
    ratio                = ratio,
    C_precomputed        = C_pre,
    inv_vX_precomputed   = inv_vX_pre
  )$GeneScan3D.Cauchy.pvalue
  GeneScan3D.Cauchy <- tmp

  GeneScan3D.Cauchy_knockoff <- matrix(NA, nrow = M, ncol = 3)
  for (k in seq_len(M)) {
    G_gbk <- matrix(
      G_gene_buffer_knockoff[k, , ],
      nrow = nrow(G_gene_buffer), ncol = ncol(G_gene_buffer)
    )
    invisible(capture.output(
      tmp <- GeneScan3D.UKB.GLMM(
        G                    = G_gbk,
        G.EnhancerAll        = if (R > 0 && length(G_EnhancerAll_knockoff) > 0)
          matrix(G_EnhancerAll_knockoff[k, , ],
                 nrow = nrow(G_gene_buffer),
                 ncol = sum(p_EnhancerAll_out)) else NULL,
        R                    = R,
        p_Enhancer           = p_EnhancerAll_out,
        window.size          = window.size,
        pos                  = positions_gene_buffer,
        MAC.threshold        = MAC.threshold,
        MAF.threshold        = MAF.threshold,
        Gsub.id              = Gsub.id[match.index],
        result.null.model.GLMM = result.null.model,
        outcome              = outcome_pre,
        sparseSigma          = sparseSigma,
        ratio                = ratio,
        C_precomputed        = C_pre,
        inv_vX_precomputed   = inv_vX_pre
      )$GeneScan3D.Cauchy.pvalue
    ))
    GeneScan3D.Cauchy_knockoff[k, ] <- tmp
  }

  return(list(
    GeneScan3D.Cauchy          = GeneScan3D.Cauchy,
    GeneScan3D.Cauchy_knockoff = GeneScan3D.Cauchy_knockoff
  ))
}


Knockoffgeneration.gene.buffer <- function(
  G_gene_buffer_surround = G_gene_buffer_surround,
  positions = NULL,
  gene_buffer_start = gene_buffer_start,
  gene_buffer_end = gene_buffer_end,
  M = 5, surround.region = 100000, LD.filter = 0.75,
  return_details = FALSE
) {
  if (is.null(positions))
    positions <- extract_position_universal(colnames(G_gene_buffer_surround))
  prepared <- .bigknock_prepare_region(
    G_gene_buffer_surround, positions,
    gene_buffer_start, gene_buffer_end, LD.filter,
    min_mac = 25, label = "gene buffer"
  )
  n <- nrow(prepared$surround_matrix)
  knockoff <- create.MK.AL_gene_buffer_bigknock(
    X = prepared$surround_matrix,
    pos = prepared$surround_positions,
    gene_buffer_start = gene_buffer_start,
    gene_buffer_end = gene_buffer_end,
    M = M, corr_max = LD.filter, maxN.neighbor = Inf,
    maxBP.neighbor = surround.region, corr_base = 0.05,
    n.AL = floor(10 * n^(1/3) * log(n)),
    thres.ultrarare = 25, R2.thres = LD.filter
  )
  details <- list(
    knockoff = knockoff,
    selected_input_index = prepared$target_input_index,
    positions = prepared$target_positions
  )
  if (isTRUE(return_details)) details else knockoff
}


Knockoffgeneration.enhancer <- function(
  G_enhancer_surround = G_enhancer_surround,
  positions = NULL,
  enhancer_start = enhancer_start,
  enhancer_end = enhancer_start,
  M = 5, surround.region = 50000, LD.filter = 0.75,
  return_details = FALSE
) {
  if (is.null(positions))
    positions <- extract_position_universal(colnames(G_enhancer_surround))
  prepared <- .bigknock_prepare_region(
    G_enhancer_surround, positions,
    enhancer_start, enhancer_end, LD.filter,
    min_mac = 25, label = "enhancer", min_target_variants = 1L
  )
  n <- nrow(prepared$surround_matrix)
  knockoff <- create.MK.AL_enhancer(
    X = prepared$surround_matrix,
    pos = prepared$surround_positions,
    enhancer_start = enhancer_start,
    enhancer_end = enhancer_end,
    M = M, corr_max = LD.filter, maxN.neighbor = Inf,
    maxBP.neighbor = surround.region, corr_base = 0.05,
    n.AL = floor(10 * n^(1/3) * log(n)),
    thres.ultrarare = 25, R2.thres = LD.filter
  )
  details <- list(
    knockoff = knockoff,
    selected_input_index = prepared$target_input_index,
    positions = prepared$target_positions
  )
  if (isTRUE(return_details)) details else knockoff
}


######### Other functions #########
#Optimize create.MK.AL function provided by Zihuai
#Knockoff generation for gene buffer regions
create.MK.AL_gene_buffer_bigknock <- function(X=G_gene_buffer_surround,pos,gene_buffer_start,gene_buffer_end,M,corr_max=LD.filter,maxN.neighbor=Inf,
                                              maxBP.neighbor=surround.region,corr_base=0.05,n.AL=floor(10*n^(1/3)*log(n)),
                                              thres.ultrarare=25,R2.thres=LD.filter) {

  method='shrinkage'
  cor.X <- .kp_sparse_cov_cor(
    X, need_cov = FALSE, need_cor = TRUE
  )$cor

  #svd to get leverage score, can be optimized;update: tried fast leveraging, but the R matrix is singular possibly because X is sparse.
  #Fast Truncated Singular Value Decomposition
  if(method=='shrinkage'){
    prob <- .bigknock_shrinkage_prob(X)
  }

  index.AL<-sample(1:nrow(X),min(n.AL,nrow(X)),replace = FALSE,prob=prob) #sampling r samples from n samples, using shrinkage leveraging estimator
  w<-1/sqrt(n.AL*prob[index.AL])
  X.AL<-w*X[index.AL, , drop = FALSE] #n.AL samples

  cov.X.AL <- .kp_sparse_cov_cor(
    X.AL, need_cov = TRUE, need_cor = FALSE
  )$cov
  skip.index <- which(colSums(X.AL != 0) <= thres.ultrarare)

  Sigma.distance = as.dist(1 - abs(cor.X))
  if(ncol(X)>1){
    fit = hclust(Sigma.distance, method="single") #hierarchical clustering
    corr_max = corr_max
    clusters = cutree(fit, h=1-corr_max)  #variants from two different clusters do not have a correlation greater than 0.75.
  }else{clusters<-1}

  X_k<-list()
  for(k in 1:M){
    X_k[[k]]<-matrix(0,nrow=nrow(X),ncol=ncol(X))
    #X_k[[k]]<-big.matrix(nrow=nrow(X),ncol=ncol(X),init=0,shared=FALSE)
  }

  ##only run snps within gene buffer
  snps_ind=which(pos<=gene_buffer_end&pos>=gene_buffer_start)

  index.exist<-c()
  for (k in unique(clusters[snps_ind])){
    #print(paste0('cluster',k))
    cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
    for(i in which(clusters==k)[which(clusters==k)%in%snps_ind]){
      #print(i)
      rate<-1;R2<-1;temp.maxN.neighbor<-maxN.neighbor
      while(R2>=R2.thres){ #avoid over-fitting
        temp.maxN.neighbor<-floor(temp.maxN.neighbor/rate)
        snp.pos=as.numeric(gsub("^.*\\:","",names(clusters[i])))
        #+-100kb surrounding region
        index.pos<-which(pos>=max(snp.pos-maxBP.neighbor,pos[1]) & pos<=min(snp.pos+maxBP.neighbor,pos[length(pos)]))
        #correlation between this snp with other snps in +-100kb surrounding region
        temp<-abs(cor.X[i,])
        temp[which(clusters==k)]<-0 #exclude variants if they are in the same cluster as the target variant
        temp[-index.pos]<-0 #only focus on +-100kb surrounding region
        temp[which(temp<=corr_base)]<-0
        index<-order(temp,decreasing=T)
        if(sum(temp!=0,na.rm=T)==0 | temp.maxN.neighbor==0){index<-NULL}else{
          index<-setdiff(index[1:min(length(index),floor((nrow(X))^(1/3)),temp.maxN.neighbor,sum(temp!=0,na.rm=T))],i)
        } #top K snps up to K=n^1/3=75

        y<-X[,i] #n samples
        if(length(index)==0){fitted.values<-0}
        if(i %in% skip.index){fitted.values<-0}
        if(!(i %in% skip.index |length(index)==0)){
          x.AL<-X.AL[,index,drop=F]; #n.AL by K
          n.exist<-length(intersect(index,index.exist))
          x.exist.AL<-matrix(0,nrow=nrow(X.AL),ncol=n.exist*M)
          if(length(intersect(index,index.exist))!=0){
            for(j in 1:M){ # this is the most time-consuming part
              x.exist.AL[,((j-1)*n.exist+1):(j*n.exist)]<-w*X_k[[j]][index.AL,intersect(index,index.exist),drop=F]
            }
          }
          y.AL<-w*X[index.AL,i]; #n.AL

          temp.xy<-rbind(mean(y.AL),crossprod(x.AL,y.AL)/length(y.AL)-colMeans(x.AL)*mean(y.AL))
          temp.xy<-rbind(temp.xy,crossprod(x.exist.AL,y.AL)/length(y.AL)-colMeans(x.exist.AL)*mean(y.AL))
          temp.cov.cross <- .kp_sparse_cross_cov(x.AL, x.exist.AL)
          temp.cov <- .kp_sparse_cov_cor(
            x.exist.AL, need_cov = TRUE, need_cor = FALSE
          )$cov
          temp.xx<-cov.X.AL[index,index]
          temp.xx<-rbind(cbind(temp.xx,temp.cov.cross),cbind(t(temp.cov.cross),temp.cov))
          temp.xx<-cbind(0,temp.xx)
          temp.xx<-rbind(c(1,rep(0,ncol(temp.xx)-1)),temp.xx)

          svd.fit<-svd(temp.xx)
          v<-svd.fit$v
          cump<-cumsum(svd.fit$d)/sum(svd.fit$d)
          n.svd<-which(cump>=0.999)[1]
          svd.index<-intersect(1:n.svd,which(svd.fit$d!=0))
          temp.inv<-v[,svd.index,drop=F]%*%(svd.fit$d[svd.index]^(-1)*t(v[,svd.index,drop=F]))
          temp.beta<-temp.inv%*%temp.xy #least square estimate for regression coefficient, alpha and beta_k

          x<-X[,index,drop=F]
          temp.j<-1
          fitted.values<-temp.beta[1]+x%*%temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F]-sum(colMeans(x)*temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])

          if(length(intersect(index,index.exist))!=0){
            temp.j<-temp.j+ncol(x)
            for(j in 1:M){
              temp.x<-X_k[[j]][,intersect(index,index.exist),drop=F]
              if(ncol(temp.x)>=1){
                fitted.values<-fitted.values+temp.x%*%temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F]-sum(colMeans(temp.x)*temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])
              }
              temp.j<-temp.j+ncol(temp.x)
            }
          }
        }
        residuals<-as.numeric(y-fitted.values)
        #overfitted model
        R2<-1-var(residuals,na.rm=T)/var(y,na.rm=T)
        rate<-rate*2;temp.maxN.neighbor<-length(index)
      }
      cluster.fitted[,match(i,which(clusters==k))]<-as.vector(fitted.values)
      cluster.residuals[,match(i,which(clusters==k))]<-as.vector(residuals)
      index.exist<-c(index.exist,i)
    }
    #sample mutiple knockoffs
    cluster.sample.index<-sapply(1:M,function(x)sample(1:nrow(X)))
    for(j in 1:M){
      X_k[[j]][,which(clusters==k)]<-round(cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F],digits=1)
    }
  }

  #save knockoffs of gene buffer region
#   print('saving knockoffs of gene buffer region')
  #G_gene_buffer=X[,snps_ind]
  G_gene_buffer_knockoff <- array(0, dim = c(M, nrow(X), length(snps_ind)))
  for (j in 1:M) {
    G_gene_buffer_knockoff[j, ,] <-X_k[[j]][,snps_ind]
  }
  rm(X_k)

  #G_gene_buffer_knockoff=list(G_gene_buffer=G_gene_buffer,G_gene_buffer_knockoff=G_gene_buffer_knockoff)
  return(G_gene_buffer_knockoff)
}

create.MK.AL_enhancer <- function(X=G_enhancer_surround,pos,enhancer_start,enhancer_end,M,corr_max=0.75,maxN.neighbor=Inf,
                                  maxBP.neighbor=50000,corr_base=0.05,n.AL=floor(10*n^(1/3)*log(n)),
                                  thres.ultrarare=25,R2.thres=0.75) {

  method='shrinkage'
  cor.X <- .kp_sparse_cov_cor(
    X, need_cov = FALSE, need_cor = TRUE
  )$cor

  #svd to get leverage score, can be optimized;update: tried fast leveraging, but the R matrix is singular possibly because X is sparse.
  if(method=='shrinkage'){
    prob <- .bigknock_shrinkage_prob(X)
  }

  index.AL<-sample(1:nrow(X),min(n.AL,nrow(X)),replace = FALSE,prob=prob)
  w<-1/sqrt(n.AL*prob[index.AL])
  X.AL<-w*X[index.AL, , drop = FALSE]
  cov.X.AL <- .kp_sparse_cov_cor(
    X.AL, need_cov = TRUE, need_cor = FALSE
  )$cov
  skip.index <- which(colSums(X.AL != 0) <= thres.ultrarare)

  Sigma.distance = as.dist(1 - abs(cor.X))
  if(ncol(X)>1){
    fit = hclust(Sigma.distance, method="single")
    corr_max = corr_max
    clusters = cutree(fit, h=1-corr_max)
  }else{clusters<-1}

  X_k<-list()
  ##only focus on snps within gene buffer
  for(k in 1:M){
    #X_k[[k]]<-big.matrix(nrow=nrow(X),ncol=ncol(X),init=0,shared=FALSE)
    X_k[[k]]<-matrix(0,nrow=nrow(X),ncol=ncol(X))
  }

  snps_ind=which(pos<=enhancer_end&pos>=enhancer_start)

  index.exist<-c()
  for (k in unique(clusters[snps_ind])){
    #print(paste0('cluster',k))
    cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
    for(i in which(clusters==k)[which(clusters==k)%in%snps_ind]){
      #print(i)
      rate<-1;R2<-1;temp.maxN.neighbor<-maxN.neighbor

      while(R2>=R2.thres){

        temp.maxN.neighbor<-floor(temp.maxN.neighbor/rate)
        snp.pos=as.numeric(gsub("^.*\\:","",names(clusters[i])))
        index.pos<-which(pos>=max(snp.pos-maxBP.neighbor,pos[1]) & pos<=min(snp.pos+maxBP.neighbor,pos[length(pos)]))

        temp<-abs(cor.X[i,]);temp[which(clusters==k)]<-0;temp[-index.pos]<-0
        temp[which(temp<=corr_base)]<-0

        index<-order(temp,decreasing=T)
        if(sum(temp!=0,na.rm=T)==0 | temp.maxN.neighbor==0){index<-NULL}else{
          index<-setdiff(index[1:min(length(index),floor((nrow(X))^(1/3)),temp.maxN.neighbor,sum(temp!=0,na.rm=T))],i)
        }

        y<-X[,i]
        if(length(index)==0){fitted.values<-0}
        if(i %in% skip.index){fitted.values<-0}
        if(!(i %in% skip.index |length(index)==0)){

          x.AL<-X.AL[,index,drop=F];
          n.exist<-length(intersect(index,index.exist))
          x.exist.AL<-matrix(0,nrow=nrow(X.AL),ncol=n.exist*M)
          if(length(intersect(index,index.exist))!=0){
            for(j in 1:M){ # this is the most time-consuming part
              x.exist.AL[,((j-1)*n.exist+1):(j*n.exist)]<-w*X_k[[j]][index.AL,intersect(index,index.exist),drop=F]
            }
          }
          y.AL<-w*X[index.AL,i];

          temp.xy<-rbind(mean(y.AL),crossprod(x.AL,y.AL)/length(y.AL)-colMeans(x.AL)*mean(y.AL))
          temp.xy<-rbind(temp.xy,crossprod(x.exist.AL,y.AL)/length(y.AL)-colMeans(x.exist.AL)*mean(y.AL))
          temp.cov.cross <- .kp_sparse_cross_cov(x.AL, x.exist.AL)
          temp.cov <- .kp_sparse_cov_cor(
            x.exist.AL, need_cov = TRUE, need_cor = FALSE
          )$cov
          temp.xx<-cov.X.AL[index,index]
          temp.xx<-rbind(cbind(temp.xx,temp.cov.cross),cbind(t(temp.cov.cross),temp.cov))
          temp.xx<-cbind(0,temp.xx)
          temp.xx<-rbind(c(1,rep(0,ncol(temp.xx)-1)),temp.xx)

          svd.fit<-svd(temp.xx)
          v<-svd.fit$v
          cump<-cumsum(svd.fit$d)/sum(svd.fit$d)
          n.svd<-which(cump>=0.999)[1]
          svd.index<-intersect(1:n.svd,which(svd.fit$d!=0))
          temp.inv<-v[,svd.index,drop=F]%*%(svd.fit$d[svd.index]^(-1)*t(v[,svd.index,drop=F]))
          temp.beta<-temp.inv%*%temp.xy

          x<-X[,index,drop=F]
          temp.j<-1
          fitted.values<-temp.beta[1]+x%*%temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F]-sum(colMeans(x)*temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])

          if(length(intersect(index,index.exist))!=0){
            temp.j<-temp.j+ncol(x)
            for(j in 1:M){
              temp.x<-X_k[[j]][,intersect(index,index.exist),drop=F]
              if(ncol(temp.x)>=1){
                fitted.values<-fitted.values+temp.x%*%temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F]-sum(colMeans(temp.x)*temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])
              }
              temp.j<-temp.j+ncol(temp.x)
            }
          }
        }
        residuals<-as.numeric(y-fitted.values)
        #overfitted model
        R2<-1-var(residuals,na.rm=T)/var(y,na.rm=T)
        rate<-rate*2;temp.maxN.neighbor<-length(index)
      }
      cluster.fitted[,match(i,which(clusters==k))]<-as.vector(fitted.values)
      cluster.residuals[,match(i,which(clusters==k))]<-as.vector(residuals)
      index.exist<-c(index.exist,i)
    }
    #sample mutiple knockoffs
    cluster.sample.index<-sapply(1:M,function(x)sample(1:nrow(X)))
    for(j in 1:M){
      X_k[[j]][,which(clusters==k)]<-round(cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F],digits=1)
    }
  }

  #save knockoffs of enhancer
#   print('saving knockoffs of enhancer')
  G_enhancer_knockoff <- array(0, dim = c(M, nrow(X), length(snps_ind)))
  for (j in 1:M) {
    G_enhancer_knockoff[j, ,] <-X_k[[j]][,snps_ind]
  }
  rm(X_k)

  return(G_enhancer_knockoff)
}



# Robust sparse solve: try Cholesky, fall back to ridge-regularized LU
.safe_solve_sparse <- function(A, B, ridge = 1e-4) {
  res <- tryCatch(
    Matrix::solve(A, B, sparse = TRUE),
    error = function(e) NULL
  )
  if (!is.null(res)) return(res)
  # Ridge fallback: A + ridge * I
  n <- nrow(A)
  A_ridge <- A + Matrix::Diagonal(n, ridge)
  tryCatch(
    Matrix::solve(A_ridge, B, sparse = TRUE),
    error = function(e) stop("Sparse solve failed even with ridge: ", conditionMessage(e))
  )
}

GeneScan3D.UKB.GLMM<-function(G=G_gene_buffer,G.EnhancerAll=G_EnhancerAll,R=length(p_EnhancerAll),
                              p_Enhancer=p_EnhancerAll,window.size=c(1000,5000,10000),pos=pos_gene_buffer,
                              MAC.threshold=10,MAF.threshold=0.01,Gsub.id=Gsub.id,
                              result.null.model.GLMM=result.null.model.GLMM,outcome='C',
                              sparseSigma=sparseSigma,ratio=ratio,
                              C_precomputed=NULL, inv_vX_precomputed=NULL){
  #load preliminary features
  mu<-as.vector(result.null.model.GLMM$fitted.values)
  Y.res<-as.vector(result.null.model.GLMM$residuals)
  X<-result.null.model.GLMM$X #covariates include intercept

  # Direct calls can omit either cache; build the common pair once rather
  # than maintaining a second implementation of the same sparse solve/SVD.
  shared_precomputed <- NULL
  if (is.null(C_precomputed) || is.null(inv_vX_precomputed)) {
    shared_precomputed <- .bigknock_glmm_precompute(
      result.null.model.GLMM, sparseSigma
    )
  }
  C <- if (is.null(C_precomputed)) shared_precomputed$C else C_precomputed
  #genotype filtering/checking/missing values imputation
  G_filter <- tryCatch(
    Genotype_filter(G, pos, impute.method = 'fixed'),
    error = function(e) NULL
  )
  if (is.null(G_filter) || ncol(G_filter$G) <= 1L) return(NULL)
  G   <- G_filter$G
  pos <- G_filter$pos

  #match phenotype id (phecode) and genotype id
  if(length(Gsub.id)==0){match.index<-match(as.numeric(result.null.model.GLMM$sampleID),1:nrow(G))}else{
    match.index<-match(result.null.model.GLMM$sampleID,Gsub.id)
  }
  if(mean(is.na(match.index))>0){
    msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
    warning(msg,call.=F)
  }
  
  #individuals ids are matched with genotype
  G=Matrix(G[match.index,])
#   print("match")
  #generate window matrix to specify the variants in each window
  window.matrix0_gene_buffer<-c()
  for(size in window.size){
    if (size==1){next}
    pos.tag<-seq(min(pos),max(pos),by=size*1/2)
    pos.tag<-sapply(pos.tag,function(x)pos[which.min(abs(x-pos))])
    window.matrix0_gene_buffer<-cbind(window.matrix0_gene_buffer,sapply(pos.tag,function(x)as.numeric(pos>=x & pos<x+size)))
  }

  window.string_gene_buffer<-apply(window.matrix0_gene_buffer,2,function(x)paste(as.character(x),collapse = ""))
  window.matrix_gene_buffer<-Matrix(window.matrix0_gene_buffer[,match(unique(window.string_gene_buffer),window.string_gene_buffer)])
  #Number of 1-D windows to scan the gene buffer region
  M_gene_buffer=dim(window.matrix_gene_buffer)[2]
#   print("window matrix")
  ##single variant score tests, related samples using SAIGE null GLMM
  if(outcome=='D'){v=as.numeric((mu*(1-mu)))}
  if(outcome=='C'){v=1/result.null.model.GLMM$theta[1]} #phi is residual variance

  # pseudoinverse of t(X) %*% (v*X) (use precomputed if provided)
  inv_vX <- if (is.null(inv_vX_precomputed)) {
    shared_precomputed$inv_vX
  } else {
    inv_vX_precomputed
  }

  #covariate adjusted genotypes
  G_tilde=G-X%*%inv_vX%*%(t(X)%*%(v*G))

  #variance-adjusted score statistics
  #as.vector(t(G_tilde)%*%Y.res)==as.vector(t(G)%*%Y.res)
  S=as.vector(t(G_tilde)%*%Y.res)/result.null.model.GLMM$theta[1]

  ##GLMM
  #adjusted score statistics, without SPA
  if(outcome=='C'){
    invSigma_G_tilde<-.safe_solve_sparse(sparseSigma, G_tilde)
    V=t(G_tilde)%*%invSigma_G_tilde
    p.single=pchisq(S^2/(ratio*diag(V)),df=1,lower.tail=F)
  }

  #with SPA
  if(outcome=='D'){
    qtilde =S/sqrt(ratio) +as.vector(t(G_tilde)%*%mu)
    #The term as.vector(t(G_tilde)%*%mu) would be removed in SPAtest::Saddle_Prob
    #keep the ratio to estimate variance of scores in Saddle_Prob
    p.single=rep(NA,ncol(G))
    for (p in 1:ncol(G)){
      p.single[p]=SPAtest::Saddle_Prob(q=as.vector(qtilde)[p], mu = mu, g = G_tilde[,p])$p.value
    }
  }
#   print("GLMM")
  GeneScan1D.Cauchy.window=matrix(NA,nrow=M_gene_buffer,ncol=3)
  #Burden test: for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically; for binary traits, use SPA gene- or region-based score test
  #SKAT test: for continuous traits, compute p-value use Davies; for binary traits, use SPA gene- or region-based score test

  for (m in 1:M_gene_buffer){
    # print(paste0('1D-window',m))

    #Create index for each window
    index.window<-(window.matrix_gene_buffer[,m]==1)
    G.window=G[,index.window]
    G.window=Matrix(G.window)
    #if there is no variant in this window, then do not conduct combined test in this window, move to the next one
    if(dim(G.window)[2]<=1){
      next
    }

    MAF.window<-apply(G.window,2,mean)/2
    MAC.window<-apply(G.window,2,sum)
    weight.beta_125<-dbeta(MAF.window,1,25)
    weight.beta_1<-dbeta(MAF.window,1,1)

    weight.matrix<-cbind(MAC.window<MAC.threshold,(MAF.window<MAF.threshold&MAC.window>=MAC.threshold)*weight.beta_125,(MAF.window>=MAF.threshold)*weight.beta_1)
    #ultra-rare variants, rare and common variants
    colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
    weight.matrix<-Matrix(weight.matrix)

    #Single variant score test for all variants in the window, SPA p-values for binary traits
    p.single.window<-p.single[index.window]

    #approximation the covariance matrix for GLMM: t(G) P_S G
    #G.window=G.window-X%*%solve(t(X)%*%(v*X))%*%(t(X)%*%(v*G.window))
    invSigma_G.window<-.safe_solve_sparse(sparseSigma, G.window)

    A<-t(G.window)%*%invSigma_G.window
    B<-t(X)%*%invSigma_G.window
    K_S=A-t(B)%*%C%*%B
    #adjusted covariance matrix
    K=K_S*ratio

    #SPA gene-based tests
    if(outcome=='D'){
      V=diag(K)
      #adjusted variance
      v_tilde=as.vector(S^2)[index.window]/qchisq(p.single.window,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
      #adjusted covariance matrix
      K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
    }

    #Burden test: for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically
    #for binary traits, calculate the SPA gene-based p-value of Burden
    p.burden<-matrix(NA,1,ncol(weight.matrix))
    for (j in 1:ncol(weight.matrix)){
      if (sum(weight.matrix[,j]!=0)>1){
        #only conduct Burden test for at least 1 variants
        temp.window.matrix<-weight.matrix[,j]
        G.window2<-as.matrix(G.window%*%temp.window.matrix)
        weights=as.vector(weight.matrix[,j])
        if(outcome=='D'){ #SPA-adjusted
          p.burden[,j]<-pchisq(as.numeric((t(G.window2)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F) ;
        }else{
          #continuous
          p.burden[,j]<-pchisq(as.numeric((t(G.window2)%*%Y.res/result.null.model.GLMM$theta[1])^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F) ;
        }
      }
    }

    score<-as.vector(S)[index.window]
    p.dispersion<-matrix(NA,1,ncol(weight.matrix))
    #For extremely rare variants, do not conduct SKAT, change MAC.threshold to 10, do not apply resampling based moment matching
    weight.matrix0=(MAC.window>=MAC.threshold)*weight.matrix
    for (j in 2:ncol(weight.matrix)){
      if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
        if(outcome=='D'){
          #binary
          p.dispersion[,j]<-Get.p.SKAT_noMA(score,K=K_tilde,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j])
        }else{
          #continuous
          p.dispersion[,j]<-Get.p.SKAT_noMA(score,K=K,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j])
        }
      }
    }

    p.individual1<-Get.cauchy.scan(p.single.window,as.matrix((MAC.window>=MAC.threshold & MAF.window<MAF.threshold))) #rare variants
    p.individual2<-Get.cauchy.scan(p.single.window,as.matrix((MAF.window>=MAF.threshold))) #common and low frequency variants
    p.individual<-cbind(p.burden,p.dispersion,p.individual1,p.individual2);
    colnames(p.individual)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')

    p.Cauchy<-as.matrix(apply(p.individual,1,Get.cauchy))
    #aggregated Cauchy association test
    test.common<-grep('MAF>=MAF.threshold',colnames(p.individual))
    p.Cauchy.common<-as.matrix(apply(p.individual[,test.common,drop=FALSE],1,Get.cauchy))
    p.Cauchy.rare<-as.matrix(apply(p.individual[,-test.common,drop=FALSE],1,Get.cauchy))
    GeneScan1D.Cauchy.window[m,]=c(p.Cauchy,p.Cauchy.common,p.Cauchy.rare)
  }

  GeneScan1D.Cauchy=c(Get.cauchy(GeneScan1D.Cauchy.window[,1]),Get.cauchy(GeneScan1D.Cauchy.window[,2]),Get.cauchy(GeneScan1D.Cauchy.window[,3]))
#   print("1d scan")
  ###Obtain p-values for R enhancers
  GeneScan3D.Cauchy.EnhancerAll=c()
  if(R!=0){
    for (r in 1:R){ #Loop for each enhancer
    #   print(paste0('Enhancer',r))
      if (r==1){
        G.Enhancer=as.matrix(G.EnhancerAll[,1:cumsum(p_Enhancer)[r]])
      }else{
        G.Enhancer=as.matrix(G.EnhancerAll[,(cumsum(p_Enhancer)[r-1]+1):cumsum(p_Enhancer)[r]])
      }

      G.Enhancer=Genotype_filter_Enhancer(G.Enhancer=G.Enhancer,impute.method='fixed')

      #individuals ids are matched with genotype
      G.window.Enhancer=Matrix(G.Enhancer[match.index,])
      MAF.window.Enhancer<-apply(G.window.Enhancer,2,mean)/2
      MAC.window.Enhancer<-apply(G.window.Enhancer,2,sum)

      weight.beta_125<-dbeta(MAF.window.Enhancer,1,25)
      weight.beta_1<-dbeta(MAF.window.Enhancer,1,1)
      weight.matrix<-cbind(MAC.window.Enhancer<MAC.threshold,(MAF.window.Enhancer<MAF.threshold&MAC.window.Enhancer>=MAC.threshold)*weight.beta_125,(MAF.window.Enhancer>=MAF.threshold)*weight.beta_1)
      colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
      weight.matrix<-Matrix(weight.matrix)

      #Single variant score test for all variants in the enhancer
      G_tilde.Enhancer=G.window.Enhancer-X%*%inv_vX%*%(t(X)%*%(v*G.window.Enhancer))

      S.Enhancer=as.vector(t(G_tilde.Enhancer)%*%Y.res)/result.null.model.GLMM$theta[1]

      ##GLMM
      #adjusted score statistics, without SPA
      if(outcome=='C'){
        invSigma_G_tilde.Enhancer<-.safe_solve_sparse(sparseSigma, G_tilde.Enhancer)
        V.Enhancer=t(G_tilde.Enhancer)%*%invSigma_G_tilde.Enhancer
        p.single.Enhancer=pchisq(S.Enhancer^2/(ratio*diag(V.Enhancer)),df=1,lower.tail=F)
      }
      #with SPA
      if(outcome=='D'){
        #Observed test statistic
        qtilde.Enhancer =as.vector(S.Enhancer)/sqrt(ratio) +as.vector(t(G_tilde.Enhancer)%*%mu)
        #The term as.vector(t(G_tilde.Enhancer)%*%mu) would be removed in SPAtest::Saddle_Prob
        #keep the ratio to estimate variance of scores in Saddle_Prob
        p.single.Enhancer=rep(NA,ncol(G_tilde.Enhancer))
        for (p in 1:ncol(G_tilde.Enhancer)){
          p.single.Enhancer[p]=SPAtest::Saddle_Prob(q=as.vector(qtilde.Enhancer)[p], mu = mu, g = G_tilde.Enhancer[,p])$p.value
        }
      }

      p.burden.Enhancer<-matrix(NA,1,ncol(weight.matrix))
      p.dispersion.Enhancer<-matrix(NA,1,ncol(weight.matrix))

      if(length(p.single.Enhancer)>1){
        #enhancer have more than 1 variant, then conduct SKAT and burden; otherwise only conduct single variant score test
        #approximation the covariance matrix for GLMM
        #t(G) P_S G
        invSigma_G.window.Enhancer<-.safe_solve_sparse(sparseSigma, G.window.Enhancer)
        A<-t(G.window.Enhancer)%*%invSigma_G.window.Enhancer
        B<-t(X)%*%invSigma_G.window.Enhancer
        K_S=A-t(B)%*%C%*%B
        #adjusted covariance matrix
        K=K_S*ratio

        #SPA gene-based tests
        if(outcome=='D'){
          V=diag(K)
          #adjusted variance
          v_tilde=as.vector(S.Enhancer^2)/qchisq(p.single.Enhancer,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
          #adjusted covariance matrix
          K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
        }
      }

      #Burden
      for (j in 1:ncol(weight.matrix)){
        if (sum(weight.matrix[,j]!=0)>1){
          #only conduct Burden test for at least 1 variants
          temp.window.matrix<-weight.matrix[,j]
          G.window.Enhancer2<-as.matrix(G.window.Enhancer%*%temp.window.matrix)
          weights=as.vector(weight.matrix[,j])
          if(outcome=='D'){ #SPA-adjusted
            p.burden.Enhancer[,j]<-pchisq(as.numeric((t(G.window.Enhancer2)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F)
          }else{
            #continuous
            p.burden.Enhancer[,j]<-pchisq(as.numeric((t(G.window.Enhancer2)%*%Y.res/result.null.model.GLMM$theta[1])^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F)
          }
        }
      }

      #SKAT
      #For extremely rare variants, do not conduct SKAT
      for (j in 2:ncol(weight.matrix)){
        if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
          if(outcome=='D'){
            #binary
            p.dispersion.Enhancer[,j]<-Get.p.SKAT_noMA(S.Enhancer,K=K_tilde,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j])
          }else{
            #continuous
            p.dispersion.Enhancer[,j]<-Get.p.SKAT_noMA(S.Enhancer,K=K,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j])
          }
        }
      }

      p.individual1.Enhancer<-Get.cauchy.scan(p.single.Enhancer,as.matrix((MAC.window.Enhancer>=MAC.threshold & MAF.window.Enhancer<MAF.threshold))) #rare variants
      p.individual2.Enhancer<-Get.cauchy.scan(p.single.Enhancer,as.matrix((MAF.window.Enhancer>=MAF.threshold))) #common and low frequency variants
      p.individual.Enhancer<-cbind(p.burden.Enhancer ,p.dispersion.Enhancer,p.individual1.Enhancer,p.individual2.Enhancer);
      colnames(p.individual.Enhancer)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')

      #aggregated Cauchy association test
      p.Cauchy.Enhancer<-as.matrix(apply(p.individual.Enhancer,1,Get.cauchy))
      test.common<-grep('MAF>=MAF.threshold',colnames(p.individual.Enhancer))
      p.Cauchy.common.Enhancer<-as.matrix(apply(p.individual.Enhancer[,test.common,drop=FALSE],1,Get.cauchy))
      p.Cauchy.rare.Enhancer<-as.matrix(apply(p.individual.Enhancer[,-test.common,drop=FALSE],1,Get.cauchy))
      GeneScan3D.Cauchy.Enhancer=c(p.Cauchy.Enhancer,p.Cauchy.common.Enhancer,p.Cauchy.rare.Enhancer)

      GeneScan3D.Cauchy.EnhancerAll=rbind(GeneScan3D.Cauchy.EnhancerAll,GeneScan3D.Cauchy.Enhancer)
    }  #end of the loop of R enhancers
  }
#   print("enhancer scan")
  ##Obtain 3D windows and p-values
  #do not add promoter
  #M 1D windows + Enhancer r, r=1, ..., R
  GeneScan3D.window.EnhancerAll=c()
  if(R!=0){
    for (r in 1:dim(GeneScan3D.Cauchy.EnhancerAll)[1]){
      GeneScan3D.window.enhancer=data.frame(apply(cbind(GeneScan1D.Cauchy.window[,1],GeneScan3D.Cauchy.EnhancerAll[r,1]),1,Get.cauchy),
                                            apply(cbind(GeneScan1D.Cauchy.window[,2],GeneScan3D.Cauchy.EnhancerAll[r,2]),1,Get.cauchy),
                                            apply(cbind(GeneScan1D.Cauchy.window[,3],GeneScan3D.Cauchy.EnhancerAll[r,3]),1,Get.cauchy))
      colnames(GeneScan3D.window.enhancer)=c('all','common','rare')
      GeneScan3D.window.EnhancerAll=rbind(GeneScan3D.window.EnhancerAll,GeneScan3D.window.enhancer)
    }
  }else{
    GeneScan3D.window.enhancer=data.frame(Get.cauchy(GeneScan1D.Cauchy.window[,1]),
                                          Get.cauchy(GeneScan1D.Cauchy.window[,2]),
                                          Get.cauchy(GeneScan1D.Cauchy.window[,3]))
    colnames(GeneScan3D.window.enhancer)=c('all','common','rare')
    GeneScan3D.window.EnhancerAll=rbind(GeneScan3D.window.EnhancerAll,GeneScan3D.window.enhancer)
  }

  GeneScan3D.Cauchy.RE=GeneScan3D.window.EnhancerAll
  GeneScan3D.Cauchy=c(Get.cauchy(GeneScan3D.Cauchy.RE[,1]), Get.cauchy(GeneScan3D.Cauchy.RE[,2]), Get.cauchy(GeneScan3D.Cauchy.RE[,3]))

  ###min-p and RE with min-p
  RE_minp.all=NA;RE_minp.common=NA;RE_minp.rare=NA
  if(R!=0){
    RE.indicator=c(rep(1:R,each=M_gene_buffer))

    if(!is.infinite(min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE))){
      RE_minp.all=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,1]==min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE))])
    }

    if(!is.infinite(min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE))){
      RE_minp.common=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,2]==min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE))])
    }

    if(!is.infinite(min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE))){
      RE_minp.rare=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,3]==min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE))])
    }
  }

  #best enhancer
  return(list(GeneScan3D.Cauchy.pvalue=GeneScan3D.Cauchy,M=M_gene_buffer,R=R,
              minp=c(min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE)),
              RE_minp=cbind(RE_minp.all,RE_minp.common,RE_minp.rare)))  #GeneScan1D.Cauchy.pvalue=GeneScan1D.Cauchy,
}


######### Other functions #########
Get.p.SKAT_noMA<-function(score,K,window.matrix,weight){

  Q<-as.vector(t(score^2)%*%(weight*window.matrix)^2) #SKAT statistics
  K.temp<-weight*t(weight*K)

  temp<-K.temp[window.matrix[,1]!=0,window.matrix[,1]!=0]
  if(sum(temp^2)==0){p<-NA}else{
    lambda=eigen(temp,symmetric=T,only.values=T)$values #eigenvalues, mixture of chi-square
    temp.p<-SKAT_davies(Q,lambda,acc=10^(-6))$Qq

    if(length(temp.p)==0 || temp.p > 1 || temp.p <= 0){
      temp.p<-Get_Liu_PVal.MOD.Lambda(Q,lambda)
    }
    p<-temp.p
  }
  return(p)
}
SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
  r <- length(lambda)
  if (length(h) != r) warning("lambda and h should have the same length!")
  if (length(delta) != r) warning("lambda and delta should have the same length!")
  #out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0),PACKAGE="SKAT")
  out=davies(q, lambda, h = rep(1, length(lambda)), delta = rep(0,length(lambda)), sigma = 0, lim = 10000, acc = 0.0001)
  out$res <- 1 - out$res
  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
}
Get_Liu_PVal.MOD.Lambda<-function(Q.all, lambda, log.p=FALSE){
  param<-Get_Liu_Params_Mod_Lambda(lambda)
  Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
  Q.Norm1<-Q.Norm * param$sigmaX + param$muX
  p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)
  return(p.value)
}
Get_Liu_Params_Mod_Lambda<-function(lambda){
  ## Helper function for getting the parameters for the null approximation

  c1<-rep(0,4)
  for(i in 1:4){
    c1[i]<-sum(lambda^i)
  }

  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0

  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a

  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}

Get.cauchy.scan<-function(p,window.matrix){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16)
  temp<-rep(0,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[!is.small]<-as.numeric(tan((0.5-p[!is.small])*pi))

  cct.stat<-as.numeric(t(temp)%*%window.matrix/apply(window.matrix,2,sum))
  is.large<-cct.stat>1e+15 & !is.na(cct.stat)
  is.regular<-cct.stat<=1e+15 & !is.na(cct.stat)
  pval<-rep(NA,length(cct.stat))
  pval[is.large]<-(1/cct.stat[is.large])/pi
  pval[is.regular]<-1-pcauchy(cct.stat[is.regular])
  return(pval)
}
Get.cauchy<-function(p){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16) & !is.na(p)
  is.regular<-(p>=1e-16) & !is.na(p)
  temp<-rep(NA,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))

  cct.stat<-mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}
Genotype_filter=function(G,pos,impute.method='fixed'){

  if(ncol(G)==0|ncol(G)==1){
    stop('Number of variants in the gene buffer region is 0 or 1')
  }

  #missing genotype imputation
  G <- Matrix::Matrix(G, sparse = TRUE)
  G[G==-9 | G==9]=NA
  missing_mask <- is.na(G)
  N_MISS=sum(missing_mask)
  MISS.freq=.kp_col_means(missing_mask)
  rm(missing_mask)

  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
    warning(msg,call.=F)
    G=Impute(G,impute.method)
  }

  #MAF filtering
  MAF<-.kp_col_means(G)/2 #MAF of nonfiltered variants
  G[,MAF>0.5 & !is.na(MAF)]<-2-G[,MAF>0.5 & !is.na(MAF)]
  MAF<-.kp_col_means(G)/2
  variance <- .kp_col_means(G^2) - .kp_col_means(G)^2
  SNP.index<-which(MAF>0 & variance>0 & !is.na(MAF))

  check.index<-which(MAF>0 & variance>0 & !is.na(MAF)  & MISS.freq<0.1)
  if(length(check.index)<=1 ){
    stop('Number of variants with missing rate <=10% in the gene plus buffer region is <=1')
  }

  G<-Matrix::Matrix(G[,SNP.index,drop=FALSE])
  pos=pos[SNP.index]
  genotype_filter=list(G=G,pos=pos)
  return(genotype_filter)
}
Genotype_filter_Enhancer=function(G.Enhancer,impute.method='fixed'){

  #missing genotype imputation
  G.Enhancer <- Matrix::Matrix(G.Enhancer, sparse = TRUE)
  G.Enhancer[G.Enhancer==-9 | G.Enhancer==9]=NA
  missing_mask <- is.na(G.Enhancer)
  N_MISS.Enhancer=sum(missing_mask)
  MISS.freq.Enhancer=.kp_col_means(missing_mask)
  rm(missing_mask)
  if(N_MISS.Enhancer>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS.Enhancer/nrow(G.Enhancer)/ncol(G.Enhancer))
    warning(msg,call.=F)
    G.Enhancer=Impute(G.Enhancer,impute.method)
  }

  #MAF filtering
  MAF.Enhancer<-.kp_col_means(G.Enhancer)/2 #MAF of nonfiltered variants
  G.Enhancer[,MAF.Enhancer>0.5 & !is.na(MAF.Enhancer)]<-2-G.Enhancer[,MAF.Enhancer>0.5 & !is.na(MAF.Enhancer)]
  MAF.Enhancer<-.kp_col_means(G.Enhancer)/2
  variance.Enhancer <- .kp_col_means(G.Enhancer^2) -
    .kp_col_means(G.Enhancer)^2
  SNP.index.Enhancer<-which(MAF.Enhancer>0 & variance.Enhancer>0 & !is.na(MAF.Enhancer))

  G.Enhancer<-Matrix::Matrix(G.Enhancer[,SNP.index.Enhancer,drop=FALSE])
  return(G.Enhancer)
}

#####knockoff AL functions
max_nth<-function(x,n){return(sort(x,partial=length(x)-(n-1))[length(x)-(n-1)])}
