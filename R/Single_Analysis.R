# run_single_block -----------------------------------------------------------
# New params vs original: save_knockoff, load_knockoff, knockoff_file,
#   knockoff_sample_ids, stage1_only.
# When stage1_only=TRUE: saves knockoff then returns NULL immediately.
# ----------------------------------------------------------------------------
run_single_block <- function(
  blocks,
  kk,
  geno.file,
  obj_nullmodel,
  window_length,
  plink_prefix,
  impute.method,
  M,
  Gsub.id,
  save_knockoff       = FALSE,
  load_knockoff       = FALSE,
  knockoff_file       = NULL,
  knockoff_sample_ids = NULL,
  stage1_only         = FALSE
) {
  chr   <- blocks[kk, chr]
  start <- blocks[kk, start]
  stop  <- blocks[kk, stop]

  tmpdir       <- tempdir()
  block_prefix <- file.path(tmpdir, sprintf("temp_chr%d_block%d", chr, kk))
  cmd <- sprintf(
    "%s --bfile %s --chr %s --from-bp %d --to-bp %d --recode A --out %s --silent",
    plink_prefix, geno.file, chr, start, stop, block_prefix
  )
  system(paste(cmd, "2>/dev/null"), ignore.stdout = TRUE)

  raw_file <- paste0(block_prefix, ".raw")
  if (!file.exists(raw_file)) return(NULL)

  raw <- data.table::fread(raw_file, data.table = FALSE)
  unlink(paste0(block_prefix, c(".raw", ".log", ".nosex")), force = TRUE)
  if (ncol(raw) <= 6) return(NULL)

  df <- as.matrix(raw[, -(1:6), drop = FALSE])
  rm(raw)
  cat(sprintf("chr: %s, start: %d, end: %d, snp count: %d\n", chr, start, stop, ncol(df)))

  results <- Single_Window_Analysis(
      nullobj             = obj_nullmodel,
      geno                = df,
      chr                 = chr,
      window_length       = window_length,
      M                   = M,
      impute.method       = impute.method,
      Gsub.id             = Gsub.id,
      save_knockoff       = save_knockoff,
      load_knockoff       = load_knockoff,
      knockoff_file       = knockoff_file,
      knockoff_sample_ids = knockoff_sample_ids,
      stage1_only         = stage1_only
    )
  rm(df); gc()
  return(results)
}


# Single_Window_Analysis -----------------------------------------------------
# New params: save_knockoff, load_knockoff, knockoff_file,
#   knockoff_sample_ids, stage1_only.
#
# Knockoff save format (RDS list):
#   $G_k        — list of M sparse matrices [n × p], one per knockoff copy
#   $sample_ids — numeric IIDs in row order (length n)
#   $snp_pos    — numeric SNP positions (length p); used for column validation
#
# Load validation:
#   1. Column count (p): if saved p ≠ current p after QC → warn + regenerate
#   2. Row reindexing : align saved rows to current Gsub.id via .align_knockoff_samples()
# ----------------------------------------------------------------------------
Single_Window_Analysis <- function(
  nullobj,
  geno,
  chr,
  window_length       = NULL,
  M                   = 5,
  thres.single        = 0.01,
  thres.ultrarare     = 25,
  thres.missing       = 0.10,
  midout.dir          = NULL,
  jobtitle            = NULL,
  impute.method       = "fixed",
  Gsub.id             = NULL,
  bigmemory           = TRUE,
  leveraging          = TRUE,
  LD.filter           = NULL,
  save_knockoff       = FALSE,
  load_knockoff       = FALSE,
  knockoff_file       = NULL,
  knockoff_sample_ids = NULL,
  stage1_only         = FALSE
) {
  preprocess <- Preprocess(geno = geno, chr = chr, window = window_length,
                           impute.method = impute.method)
  if (is.null(preprocess)) return(NULL)

  G          <- preprocess$G
  pos        <- preprocess$pos
  window.bed <- preprocess$window.bed

  # ---- Knockoff: load or generate -----------------------------------------
  G_k          <- NULL
  need_generate <- TRUE

  if (isTRUE(load_knockoff) && !is.null(knockoff_file) && file.exists(knockoff_file)) {
    tryCatch({
      ko_obj    <- readRDS(knockoff_file)
      saved_p   <- length(ko_obj$snp_pos)
      current_p <- ncol(G)

      if (saved_p != current_p) {
        # Column count mismatch: QC filtering likely changed with the new sample set
        warning(sprintf(
          "Block chr%s kk=%s: saved knockoff has %d SNP columns, current QC gives %d; regenerating.",
          chr, "?", saved_p, current_p
        ))
        # need_generate remains TRUE → fall through to generation below
      } else {
        # Column counts match — reindex rows to current Gsub.id order
        target_ids <- if (!is.null(Gsub.id)) Gsub.id else knockoff_sample_ids
        if (!is.null(target_ids) && !is.null(ko_obj$sample_ids)) {
          ko_obj <- .align_knockoff_samples(ko_obj, target_ids)
        }
        G_k           <- ko_obj$G_k
        need_generate <- FALSE
      }
    }, error = function(e) {
      warning("Failed to load knockoff file ", knockoff_file, ": ", conditionMessage(e),
              "; regenerating.")
    })
  }

  if (need_generate) {
    G_k <- create.KS(X = G, pos = pos, M = M)

    if (isTRUE(save_knockoff) && !is.null(knockoff_file)) {
      saveRDS(
        list(
          G_k        = G_k,
          sample_ids = if (!is.null(Gsub.id)) Gsub.id else knockoff_sample_ids,
          snp_pos    = pos
        ),
        file = knockoff_file
      )
    }
  }

  # ---- Stage 1: return after saving knockoff (no association test) --------
  if (isTRUE(stage1_only)) return(invisible(NULL))

  # ---- Association test ---------------------------------------------------
  fit <- KS.chr(
    result.prelim = nullobj,
    input.X       = G,
    window.bed    = window.bed,
    input.G_k     = G_k,
    M             = M,
    thres.single  = thres.single,
    Gsub.id       = Gsub.id
  )
  rm(G_k, G, pos, window.bed); gc()
  return(fit)
}


# Preprocess -----------------------------------------------------------------
Preprocess <- function(geno, chr, window = NULL, thres.maf = 0,
                       thres.missing = 0.1, impute.method = "fixed") {

  G <- as.matrix(geno)[, -c(1:6)]
  G <- 2 - G

  if (length(G) == 0 || ncol(G) == 0) {
    warning("Number of variants in the specified range is 0", call. = FALSE)
    return(NULL)   # FIX: was `next`
  }
  if (ncol(G) == 1) {
    warning("Number of variants in the specified range is 1", call. = FALSE)
    return(NULL)   # FIX: was `next`
  }

  G <- G[, match(unique(colnames(G)), colnames(G)), drop = FALSE]

  # Missing imputation
  G[G < 0 | G > 2] <- NA
  G <- Impute(G, "fixed")

  # Filter constant variants
  s <- apply(G, 2, sd)
  G <- G[, s != 0, drop = FALSE]
  if (ncol(G) < 2) return(NULL)

  # Reorder by position
  pos <- extract_position_universal(colnames(G))
  G   <- G[, order(pos), drop = FALSE]
  MAF <- colMeans(G) / 2
  G   <- as.matrix(G)
  G[, MAF > 0.5 & !is.na(MAF)] <- 2 - G[, MAF > 0.5 & !is.na(MAF)]
  MAF <- colMeans(G) / 2
  MAC <- colSums(G)
  G   <- Matrix::Matrix(G, sparse = TRUE)

  pos <- extract_position_universal(colnames(G))
  colnames(G) <- pos
  start <- min(pos); end <- max(pos)

  # Clustering to remove highly correlated SNPs
  cor.X      <- sparse.cor(Matrix::Matrix(G))$cor
  Sigma.dist <- as.dist(1 - abs(cor.X))
  fit_clust  <- hclust(Sigma.dist, method = "single")
  clusters   <- cutree(fit_clust, h = 1 - 0.75)

  cluster.idx <- match(unique(clusters), clusters)
  G   <- G[, cluster.idx, drop = FALSE]
  MAF <- MAF[cluster.idx]; MAC <- MAC[cluster.idx]; pos <- pos[cluster.idx]

  unique.idx <- match(unique(pos), pos)
  G   <- G[, unique.idx, drop = FALSE]
  MAF <- MAF[unique.idx]; MAC <- MAC[unique.idx]; pos <- pos[unique.idx]

  if (ncol(G) < 2) return(NULL)

  # Window bed
  if (length(window) != 0) {
    window.bed <- c()
    for (size in window) {
      pos.tag    <- seq(start, end, by = size * 0.5)
      window.bed <- rbind(window.bed, cbind(chr, pos.tag, pos.tag + size))
    }
    window.bed <- window.bed[order(as.numeric(window.bed[, 2])), ]
    return(list(G = G, chr = chr, pos = pos, window.bed = window.bed))
  } else {
    return(list(G = G, chr = chr, pos = pos, window.bed = NULL))
  }
}


extract_position_universal <- function(col_names) {
  positions <- sapply(col_names, function(col) {
    numbers <- regmatches(col, gregexpr("\\d+", col))[[1]]
    if (length(numbers) == 0) {
      return(NA)
    } else if (length(numbers) == 1) {
      return(as.numeric(numbers))
    } else {
      candidate <- numbers[nchar(numbers) >= 3 & nchar(numbers) <= 9]
      if (length(candidate) > 0)
        return(as.numeric(utils::tail(candidate, 1)))
      return(as.numeric(max(numbers)))
    }
  })
  return(as.numeric(positions))
}