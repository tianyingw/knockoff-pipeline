# Internal helpers for reading PLINK exports without losing sample or variant
# identity.  PLINK .raw files contain six leading sample columns; callers must
# remove those columns exactly once and retain IID for explicit row alignment.

utils::globalVariables(c(
  ".", "chr", "FID", "IID", "start", "end", "id", "TargetGene",
  "gene", "GH_start", "GH_end"
))

.as_sample_id <- function(x, label = "sample IDs") {
  ids <- trimws(as.character(x))
  if (anyNA(ids) || any(!nzchar(ids)))
    stop(label, " contain missing or empty values.")
  if (anyDuplicated(ids))
    stop(label, " contain duplicate values; sample alignment is ambiguous.")
  ids
}

.plink_additive_export_switch <- function(plink_prefix, version_text = NULL) {
  if (is.null(version_text)) {
    version_text <- tryCatch(
      suppressWarnings(system2(
        plink_prefix, "--version", stdout = TRUE, stderr = TRUE
      )),
      error = function(e) character(0)
    )
  }
  is_v2 <- any(grepl("PLINK[[:space:]]+v?2([.]|$)", version_text,
                     ignore.case = TRUE)) ||
    grepl("plink2", basename(plink_prefix), ignore.case = TRUE)
  if (is_v2) "--export A" else "--recode A"
}

.run_plink_additive_export <- function(
  plink_prefix, geno_file, chr, start, stop, keep_arg, out_prefix
) {
  export_switch <- .plink_additive_export_switch(plink_prefix)
  cmd <- sprintf(
    "%s --bfile %s --chr %s --from-bp %d --to-bp %d %s %s --out %s --silent",
    shQuote(plink_prefix), shQuote(geno_file), as.integer(chr),
    as.integer(start), as.integer(stop), keep_arg, export_switch,
    shQuote(out_prefix)
  )
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
}

.read_plink_bim_chr <- function(geno_file, chr, plink_prefix = "plink2") {
  source_bim <- paste0(geno_file, ".bim")
  if (!file.exists(source_bim)) stop("PLINK .bim not found: ", source_bim)

  tmp_prefix <- tempfile(sprintf("KnockoffPipeline_chr%s_", chr))
  out_bim <- paste0(tmp_prefix, ".bim")
  on.exit(unlink(paste0(tmp_prefix, c(".bim", ".log", ".nosex")), force = TRUE),
          add = TRUE)

  cmd <- sprintf(
    "%s --bfile %s --chr %s --make-just-bim --out %s --silent",
    shQuote(plink_prefix), shQuote(geno_file), as.integer(chr),
    shQuote(tmp_prefix)
  )
  status <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (!identical(status, 0L) || !file.exists(out_bim)) {
    stop(
      "PLINK failed to export chromosome ", chr,
      " variant metadata with --make-just-bim. Check 'plink_path' and the input .bim file."
    )
  }

  bim <- data.table::fread(
    out_bim,
    header = FALSE,
    col.names = c("chr", "variant_id", "cm", "pos", "a1", "a2"),
    colClasses = c("character", "character", "numeric", "numeric",
                   "character", "character"),
    showProgress = FALSE
  )
  bim <- as.data.frame(bim, stringsAsFactors = FALSE)
  if (nrow(bim) == 0L) return(bim)
  if (anyNA(bim$pos)) stop("Non-numeric positions found in exported .bim metadata.")
  bim
}

.match_raw_variants <- function(raw_names, bim_metadata) {
  raw_names <- as.character(raw_names)
  if (length(raw_names) == 0L) return(bim_metadata[0, , drop = FALSE])
  if (anyDuplicated(raw_names))
    stop("PLINK .raw contains duplicate genotype column names.")
  if (is.null(bim_metadata) || nrow(bim_metadata) == 0L)
    stop("No .bim metadata are available for the exported variants.")

  id <- as.character(bim_metadata$variant_id)
  a1 <- as.character(bim_metadata$a1)
  a2 <- as.character(bim_metadata$a2)
  key_a1 <- paste0(id, "_", a1)
  key_a2 <- paste0(id, "_", a2)

  idx <- match(raw_names, key_a1)
  counted <- rep(NA_character_, length(raw_names))
  counted[!is.na(idx)] <- a1[idx[!is.na(idx)]]

  missing <- is.na(idx)
  idx_a2 <- match(raw_names[missing], key_a2)
  idx[missing] <- idx_a2
  counted[missing & !is.na(idx)] <- a2[idx[missing & !is.na(idx)]]

  missing <- is.na(idx)
  idx_id <- match(raw_names[missing], id)
  idx[missing] <- idx_id
  counted[missing & !is.na(idx)] <- a1[idx[missing & !is.na(idx)]]

  if (anyNA(idx)) {
    bad <- raw_names[is.na(idx)]
    stop(
      "Could not match ", length(bad), " PLINK .raw column(s) to .bim variant IDs/alleles. Examples: ",
      paste(utils::head(bad, 5L), collapse = ", ")
    )
  }
  if (anyDuplicated(idx))
    stop("Multiple PLINK .raw columns mapped to the same .bim row.")

  out <- bim_metadata[idx, c("chr", "variant_id", "pos", "a1", "a2"),
                      drop = FALSE]
  out$raw_name <- raw_names
  out$counted_allele <- counted
  rownames(out) <- NULL
  out
}

.prepare_raw_genotypes <- function(raw, target_ids, bim_metadata) {
  raw <- as.data.frame(raw, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(raw) <= 6L)
    stop("PLINK .raw contains no genotype columns.")

  normalized_names <- toupper(sub("^#", "", names(raw)))
  iid_col <- which(normalized_names == "IID")
  if (length(iid_col) != 1L)
    stop("PLINK .raw must contain exactly one IID column among its first six columns.")
  if (iid_col > 6L)
    stop("PLINK .raw IID column was not found among the six sample columns.")

  raw_ids <- .as_sample_id(raw[[iid_col]], "PLINK .raw IIDs")
  if (is.null(target_ids)) {
    target_ids <- raw_ids
  } else {
    target_ids <- .as_sample_id(target_ids, "target sample IDs")
  }
  row_index <- match(target_ids, raw_ids)
  if (anyNA(row_index)) {
    missing_ids <- target_ids[is.na(row_index)]
    stop(
      length(missing_ids), " target sample(s) are absent from the PLINK .raw export. Examples: ",
      paste(utils::head(missing_ids, 5L), collapse = ", ")
    )
  }

  genotype_names <- names(raw)[-(1:6)]
  variant_metadata <- .match_raw_variants(genotype_names, bim_metadata)
  geno <- as.matrix(raw[row_index, -(1:6), drop = FALSE])
  storage.mode(geno) <- "numeric"
  rownames(geno) <- target_ids
  colnames(geno) <- genotype_names

  list(geno = geno, sample_ids = target_ids,
       variant_metadata = variant_metadata)
}

.variant_fingerprint <- function(metadata) {
  required <- c("chr", "variant_id", "pos", "a1", "a2", "coded_allele")
  missing <- setdiff(required, names(metadata))
  if (length(missing) > 0L)
    stop("Variant metadata are missing: ", paste(missing, collapse = ", "))
  apply(metadata[, required, drop = FALSE], 1L, paste, collapse = "\t")
}

.make_knockoff_context <- function(test_type, M, genome_build,
                                   variant_metadata, reference_id = NULL,
                                   construction_id = NULL,
                                   random_seed = NULL) {
  list(
    schema_version = 4L,
    test_type = as.character(test_type),
    M = as.integer(M),
    genome_build = as.character(genome_build),
    reference_id = if (is.null(reference_id)) NA_character_ else as.character(reference_id),
    construction_id = if (is.null(construction_id)) NA_character_ else as.character(construction_id),
    random_seed = if (is.null(random_seed)) NA_integer_ else as.integer(random_seed),
    variant_fingerprint = .variant_fingerprint(variant_metadata)
  )
}

.assert_knockoff_context <- function(saved_context, current_context,
                                     knockoff_file = NULL) {
  where <- if (is.null(knockoff_file)) "saved knockoff" else knockoff_file
  if (is.null(saved_context) || !identical(saved_context$schema_version, 4L)) {
    stop(
      "Saved knockoff metadata are missing or use an obsolete schema: ", where,
      ". Regenerate knockoffs with the current package version."
    )
  }
  fields <- c("test_type", "M", "genome_build", "reference_id",
              "construction_id", "random_seed",
              "variant_fingerprint")
  mismatched <- fields[!vapply(fields, function(field) {
    identical(saved_context[[field]], current_context[[field]])
  }, logical(1))]
  if (length(mismatched) > 0L) {
    stop(
      "Saved knockoff is incompatible with the current run (mismatch: ",
      paste(mismatched, collapse = ", "), "): ", where,
      ". Regenerate rather than reusing it."
    )
  }
  invisible(TRUE)
}

.derive_unit_seed <- function(seed, ...) {
  if (is.null(seed)) return(NULL)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      seed < 0 || seed > .Machine$integer.max)
    stop("'seed' must be NULL or one integer between 0 and .Machine$integer.max.")
  key <- paste(c(as.integer(seed), ...), collapse = "|")
  value <- as.double(as.integer(seed)) %% 2147483646
  for (code in utf8ToInt(enc2utf8(key)))
    value <- (value * 131 + code) %% 2147483646
  as.integer(value + 1)
}

.with_local_seed <- function(seed, fun) {
  if (is.null(seed)) return(fun())
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  fun()
}

.atomic_save_rds <- function(object, path) {
  parent <- dirname(path)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(
    pattern = paste0(".", basename(path), ".tmp-"),
    tmpdir = parent
  )
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  saveRDS(object, tmp)

  # On POSIX, rename(2) atomically replaces an existing file. Some platforms
  # refuse to rename over an existing destination, so fall back to a guarded
  # two-rename sequence and restore the previous manifest if the swap fails.
  if (!file.rename(tmp, path)) {
    backup <- paste0(path, ".previous-", Sys.getpid())
    if (file.exists(backup))
      stop("Refusing to overwrite an existing checkpoint backup: ", backup)
    had_previous <- file.exists(path)
    if (had_previous && !file.rename(path, backup))
      stop("Could not preserve the previous checkpoint manifest: ", path)
    replaced <- file.rename(tmp, path)
    if (!replaced) {
      if (had_previous) file.rename(backup, path)
      stop("Could not install the new checkpoint manifest: ", path)
    }
    if (had_previous) unlink(backup, force = TRUE)
  }
  invisible(path)
}

.atomic_write_lines <- function(text, path) {
  parent <- dirname(path)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(
    pattern = paste0(".", basename(path), ".tmp-"),
    tmpdir = parent
  )
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  writeLines(as.character(text), tmp, useBytes = TRUE)

  if (!file.rename(tmp, path)) {
    backup <- paste0(path, ".previous-", Sys.getpid())
    if (file.exists(backup))
      stop("Refusing to overwrite an existing sample-list backup: ", backup)
    had_previous <- file.exists(path)
    if (had_previous && !file.rename(path, backup))
      stop("Could not preserve the previous sample list: ", path)
    replaced <- file.rename(tmp, path)
    if (!replaced) {
      if (had_previous) file.rename(backup, path)
      stop("Could not install the new sample list: ", path)
    }
    if (had_previous) unlink(backup, force = TRUE)
  }
  invisible(path)
}

.reconcile_knockoff_sample_file <- function(current_ids, sample_file,
                                             create = FALSE) {
  current_ids <- .as_sample_id(current_ids, "current knockoff sample IDs")
  if (!file.exists(sample_file)) {
    if (!isTRUE(create))
      stop("Knockoff sample list not found: ", sample_file,
           "\nRun pipeline_stage = 'stage1_knockoff' first.")
    .atomic_write_lines(current_ids, sample_file)
    return(list(
      sample_ids = current_ids,
      order = seq_along(current_ids),
      created = TRUE
    ))
  }

  saved_ids <- .as_sample_id(
    readLines(sample_file, warn = FALSE), "saved knockoff sample IDs"
  )
  n_only_saved <- length(setdiff(saved_ids, current_ids))
  n_only_current <- length(setdiff(current_ids, saved_ids))
  if (n_only_saved > 0L || n_only_current > 0L) {
    stop(
      "Current sample IDs must exactly match the saved knockoff set (apart from row order). ",
      n_only_saved, " saved-only and ", n_only_current,
      " current-only sample(s) were found. Use a new or empty knockoff_dir ",
      "to generate knockoffs for a different analysis sample set."
    )
  }

  list(
    sample_ids = saved_ids,
    order = match(saved_ids, current_ids),
    created = FALSE
  )
}

.manifest_backing_files <- function(manifest, manifest_path) {
  if (is.null(manifest) ||
      !identical(manifest$storage, "bigmemory_filebacked") ||
      is.null(manifest$descriptor_files)) return(character(0))
  desc <- file.path(dirname(manifest_path), manifest$descriptor_files)
  unique(c(desc, sub("\\.desc$", ".bin", desc)))
}

.uncommitted_backing_cleanup <- function(new_backing_files, manifest_path) {
  if (length(new_backing_files) == 0L) return(character(0))
  installed <- if (file.exists(manifest_path)) {
    tryCatch(readRDS(manifest_path), error = function(e) NULL)
  } else {
    NULL
  }
  protected <- .manifest_backing_files(installed, manifest_path)
  setdiff(new_backing_files, protected)
}
