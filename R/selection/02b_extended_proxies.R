#!/usr/bin/env Rscript
# Phase 2b: Extended Proxy Acquisition
# Adds three new proxy dimensions beyond the original research-attention proxies:
#   1. Historical priority (year of valid description → species age)
#   2. Accessibility (culture collection deposits)
#   3. Genomic characterization (assembly quality, CheckM completeness)
#
# Sources: NCBI Taxonomy XML, RefSeq assembly_summary.txt, CheckM report
# All downloads are cached to data/*.rds for re-runs.

library(dplyr)
library(rentrez)
library(xml2)
library(data.table)

# ---- Paths & config ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg)) else getwd()
proj_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))
source(file.path(proj_dir, "R", "selection", "00_config.R"))

paths <- resolve_paths()
data_dir <- paths$data_dir

df <- readRDS(file.path(data_dir, "bugphyzz_with_proxies.rds"))
cat("Loaded:", nrow(df), "rows\n")

ncbi_ids <- unique(df$NCBI_ID[!is.na(df$NCBI_ID)])
cat("Unique NCBI_IDs:", length(ncbi_ids), "\n")

# ========================================================================
# 1. NCBI Taxonomy: authority year + culture collections
# ========================================================================
cache_taxonomy <- file.path(data_dir, "taxonomy_meta.rds")

if (file.exists(cache_taxonomy)) {
  cat("\n[Taxonomy] Loading cached taxonomy_meta.rds\n")
  tax_meta <- readRDS(cache_taxonomy)
  cat("  Cached entries:", nrow(tax_meta), "\n")
} else {
  cat("\n[Taxonomy] Fetching from NCBI (batch size 200)...\n")

  batch_size <- 200
  batches <- split(ncbi_ids, ceiling(seq_along(ncbi_ids) / batch_size))
  cat("  Batches:", length(batches), "\n")

  results_list <- list()

  for (i in seq_along(batches)) {
    batch_ids <- batches[[i]]

    xml_text <- tryCatch({
      entrez_fetch(db = "taxonomy", id = batch_ids, rettype = "xml")
    }, error = function(e) {
      cat(sprintf("  Batch %d/%d FAILED: %s\n", i, length(batches), conditionMessage(e)))
      NULL
    })

    if (is.null(xml_text)) {
      Sys.sleep(1)
      next
    }

    doc <- read_xml(xml_text)
    # Only top-level Taxon elements (not lineage ancestors nested inside)
    taxons <- xml_find_all(doc, "/TaxaSet/Taxon")

    for (taxon in taxons) {
      tid <- xml_text(xml_find_first(taxon, "TaxId"))

      # Authority and type material are inside OtherNames/Name elements,
      # distinguished by ClassCDE ("authority" or "type material").
      name_nodes <- xml_find_all(taxon, ".//OtherNames/Name")
      authority_strs <- character(0)
      type_strs <- character(0)
      for (nn in name_nodes) {
        class_cde <- xml_text(xml_find_first(nn, "ClassCDE"))
        disp_name <- xml_text(xml_find_first(nn, "DispName"))
        if (class_cde == "authority") {
          authority_strs <- c(authority_strs, disp_name)
        } else if (class_cde == "type material") {
          type_strs <- c(type_strs, disp_name)
        }
      }

      # Parse year from authority strings — take current accepted name (last entry)
      year_valid <- NA_integer_
      if (length(authority_strs) > 0) {
        # The last authority string is typically the current accepted name
        auth_text <- paste(authority_strs, collapse = " ")
        year_matches <- regmatches(auth_text, gregexpr("\\b(1[7-9]\\d{2}|20[0-2]\\d)\\b", auth_text))[[1]]
        if (length(year_matches) > 0) {
          year_valid <- max(as.integer(year_matches))  # latest year = current accepted
        }
      }

      # Type material — count distinct culture collection prefixes
      collection_prefixes <- c("ATCC", "DSM", "DSMZ", "JCM", "NCTC", "CCUG",
                               "NBRC", "CIP", "LMG", "KCTC", "NCIMB", "CECT",
                               "BCRC", "NRRL", "VKM", "CGMCC", "IAM", "IFO")
      n_collections <- 0L
      has_type <- FALSE
      if (length(type_strs) > 0) {
        has_type <- TRUE
        all_type_text <- paste(type_strs, collapse = " ")
        found_prefixes <- sapply(collection_prefixes, function(p) {
          grepl(paste0("\\b", p, "\\b"), all_type_text, ignore.case = TRUE)
        })
        n_collections <- sum(found_prefixes)
      }

      auth_combined <- if (length(authority_strs) > 0) paste(authority_strs, collapse = " | ") else NA_character_
      results_list[[tid]] <- data.frame(
        NCBI_ID = as.integer(tid),
        authority = auth_combined,
        year_valid_pub = year_valid,
        n_culture_collections = n_collections,
        has_type_material = has_type,
        stringsAsFactors = FALSE
      )
    }

    if (i %% 10 == 0 || i == length(batches)) {
      cat(sprintf("  Batch %d/%d done (%d records)\n", i, length(batches), length(results_list)))
    }
    Sys.sleep(0.4)  # respect NCBI rate limits
  }

  tax_meta <- do.call(rbind, results_list)
  rownames(tax_meta) <- NULL
  saveRDS(tax_meta, cache_taxonomy)
  cat("  Saved:", cache_taxonomy, "(", nrow(tax_meta), "rows)\n")
}

# Summary
cat("\n[Taxonomy] Summary:\n")
cat("  year_valid_pub: ", sum(!is.na(tax_meta$year_valid_pub)), "/", nrow(tax_meta), " non-NA\n", sep = "")
cat("  year range: ", min(tax_meta$year_valid_pub, na.rm = TRUE), "-",
    max(tax_meta$year_valid_pub, na.rm = TRUE), "\n")
cat("  has_type_material: ", sum(tax_meta$has_type_material), "/", nrow(tax_meta), "\n", sep = "")
cat("  n_culture_collections: median =", median(tax_meta$n_culture_collections[tax_meta$has_type_material]),
    ", max =", max(tax_meta$n_culture_collections), "\n")

# ========================================================================
# 2. RefSeq assembly_summary.txt — genomic quality metrics
# ========================================================================
cache_assembly <- file.path(data_dir, "assembly_quality.rds")
asm_file <- file.path(data_dir, "assembly_summary.txt")

if (file.exists(cache_assembly)) {
  cat("\n[Assembly] Loading cached assembly_quality.rds\n")
  asm_quality <- readRDS(cache_assembly)
  cat("  Cached entries:", nrow(asm_quality), "\n")
} else {
  # Download bacteria + archaea assembly summaries
  asm_bact_file <- file.path(data_dir, "assembly_summary_bacteria.txt")
  asm_arch_file <- file.path(data_dir, "assembly_summary_archaea.txt")
  if (!file.exists(asm_bact_file)) {
    cat("\n[Assembly] Downloading bacteria assembly_summary.txt...\n")
    download.file(
      "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",
      asm_bact_file, method = "auto", quiet = FALSE
    )
  }
  if (!file.exists(asm_arch_file)) {
    cat("[Assembly] Downloading archaea assembly_summary.txt...\n")
    download.file(
      "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt",
      asm_arch_file, method = "auto", quiet = FALSE
    )
  }
  cat("[Assembly] Reading assembly summaries (bacteria + archaea)...\n")
  read_asm <- function(f) {
    dt <- fread(f, skip = 1, sep = "\t", header = TRUE, quote = "")
    old_names <- names(dt)
    old_names[1] <- sub("^#[ ]?", "", old_names[1])
    setnames(dt, old_names)
    dt
  }
  asm <- rbind(read_asm(asm_bact_file), read_asm(asm_arch_file))

  cat("  Raw rows:", nrow(asm), "\n")
  cat("  Columns:", paste(names(asm), collapse = ", "), "\n")

  # Filter to latest versions only
  if ("version_status" %in% names(asm)) {
    asm <- asm[version_status == "latest"]
    cat("  After version_status=latest:", nrow(asm), "\n")
  }

  # Aggregate per species_taxid
  asm[, seq_year := as.integer(substr(seq_rel_date, 1, 4))]

  asm_agg <- asm[, .(
    best_contig_count = min(contig_count, na.rm = TRUE),
    has_complete_genome = any(assembly_level == "Complete Genome"),
    year_first_seq = suppressWarnings(min(seq_year, na.rm = TRUE)),
    n_assemblies = .N,
    has_type_assembly = any(relation_to_type_material != "na")
  ), by = species_taxid]

  # Clean up infinite values from min/max on empty groups
  asm_agg[is.infinite(best_contig_count), best_contig_count := NA_integer_]
  asm_agg[is.infinite(year_first_seq), year_first_seq := NA_integer_]

  setnames(asm_agg, "species_taxid", "NCBI_ID")
  asm_quality <- as.data.frame(asm_agg)

  saveRDS(asm_quality, cache_assembly)
  cat("  Saved:", cache_assembly, "(", nrow(asm_quality), "rows)\n")
}

cat("\n[Assembly] Summary:\n")
cat("  Species with assemblies:", nrow(asm_quality), "\n")
cat("  has_complete_genome:", sum(asm_quality$has_complete_genome), "\n")
cat("  best_contig_count: median =", median(asm_quality$best_contig_count, na.rm = TRUE), "\n")
cat("  year_first_seq range:", min(asm_quality$year_first_seq, na.rm = TRUE), "-",
    max(asm_quality$year_first_seq, na.rm = TRUE), "\n")

# ========================================================================
# 3. CheckM report — genome completeness/contamination
# ========================================================================
cache_checkm <- file.path(data_dir, "checkm_quality.rds")
checkm_file <- file.path(data_dir, "checkm_report.txt")

if (file.exists(cache_checkm)) {
  cat("\n[CheckM] Loading cached checkm_quality.rds\n")
  checkm_quality <- readRDS(cache_checkm)
  cat("  Cached entries:", nrow(checkm_quality), "\n")
} else {
  if (!file.exists(checkm_file)) {
    cat("\n[CheckM] Downloading CheckM_report_prokaryotes.txt...\n")
    download.file(
      "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/CheckM_report_prokaryotes.txt",
      checkm_file, method = "auto", quiet = FALSE
    )
  }
  cat("[CheckM] Reading CheckM report...\n")
  ckm <- fread(checkm_file, sep = "\t", header = TRUE, quote = "")
  # Clean header (may have # prefix)
  old_names <- names(ckm)
  old_names[1] <- sub("^# ", "", old_names[1])
  setnames(ckm, old_names)

  cat("  Raw rows:", nrow(ckm), "\n")
  cat("  Columns:", paste(names(ckm), collapse = ", "), "\n")

  # Identify the species-taxid and completeness columns
  # Column names may vary; look for likely matches
  taxid_col <- grep("species.taxid|species_taxid|taxid", names(ckm), value = TRUE, ignore.case = TRUE)
  compl_col <- grep("completeness", names(ckm), value = TRUE, ignore.case = TRUE)
  contam_col <- grep("contamination", names(ckm), value = TRUE, ignore.case = TRUE)

  cat("  taxid col:", taxid_col, "\n")
  cat("  completeness col:", compl_col, "\n")
  cat("  contamination col:", contam_col, "\n")

  if (length(taxid_col) > 0 && length(compl_col) > 0) {
    # Use first matching column if multiple
    tid_c <- taxid_col[1]
    cmp_c <- compl_col[1]
    cnt_c <- if (length(contam_col) > 0) contam_col[1] else NULL

    ckm_agg <- ckm[, .(
      best_checkm_completeness = max(get(cmp_c), na.rm = TRUE)
    ), by = c(tid_c)]

    if (!is.null(cnt_c)) {
      ckm_contam <- ckm[, .(
        best_checkm_contamination = min(get(cnt_c), na.rm = TRUE)
      ), by = c(tid_c)]
      ckm_agg <- merge(ckm_agg, ckm_contam, by = tid_c)
    }

    # Clean infinities
    ckm_agg[is.infinite(best_checkm_completeness), best_checkm_completeness := NA_real_]
    if ("best_checkm_contamination" %in% names(ckm_agg)) {
      ckm_agg[is.infinite(best_checkm_contamination), best_checkm_contamination := NA_real_]
    }

    setnames(ckm_agg, tid_c, "NCBI_ID")
    checkm_quality <- as.data.frame(ckm_agg)
  } else {
    cat("  WARNING: Could not identify columns in CheckM report. Skipping.\n")
    checkm_quality <- data.frame(NCBI_ID = integer(0),
                                  best_checkm_completeness = numeric(0))
  }

  saveRDS(checkm_quality, cache_checkm)
  cat("  Saved:", cache_checkm, "(", nrow(checkm_quality), "rows)\n")
}

cat("\n[CheckM] Summary:\n")
cat("  Species with CheckM data:", nrow(checkm_quality), "\n")
if (nrow(checkm_quality) > 0) {
  cat("  best_checkm_completeness: median =",
      median(checkm_quality$best_checkm_completeness, na.rm = TRUE), "\n")
}

# ========================================================================
# 4. Join all new proxies to main data
# ========================================================================
cat("\n=== Joining extended proxies ===\n")

# Remove any pre-existing extended proxy columns (from previous runs)
extended_cols <- c("authority", "year_valid_pub", "n_culture_collections", "has_type_material",
                   "best_contig_count", "has_complete_genome", "year_first_seq",
                   "has_type_assembly", "n_assemblies",
                   "best_checkm_completeness", "best_checkm_contamination",
                   "species_age", "log_species_age", "log_best_contig_count",
                   "has_assembly", "has_checkm_eval",
                   "contig_imp_median", "contig_imp_p95", "checkm_imp")
existing_ext <- intersect(names(df), extended_cols)
if (length(existing_ext) > 0) {
  cat("  Removing pre-existing extended columns:", paste(existing_ext, collapse = ", "), "\n")
  df <- df[, setdiff(names(df), existing_ext), drop = FALSE]
}

n_before <- ncol(df)

# Join taxonomy metadata
df <- merge(df, tax_meta[, c("NCBI_ID", "year_valid_pub", "n_culture_collections",
                              "has_type_material")],
            by = "NCBI_ID", all.x = TRUE)
cat("  After taxonomy join:", ncol(df), "cols\n")

# Join assembly quality
df <- merge(df, asm_quality[, c("NCBI_ID", "best_contig_count", "has_complete_genome",
                                 "year_first_seq", "has_type_assembly")],
            by = "NCBI_ID", all.x = TRUE)
cat("  After assembly quality join:", ncol(df), "cols\n")

# Join CheckM quality
if (nrow(checkm_quality) > 0 && "best_checkm_completeness" %in% names(checkm_quality)) {
  merge_cols <- c("NCBI_ID", "best_checkm_completeness")
  if ("best_checkm_contamination" %in% names(checkm_quality)) {
    merge_cols <- c(merge_cols, "best_checkm_contamination")
  }
  df <- merge(df, checkm_quality[, merge_cols], by = "NCBI_ID", all.x = TRUE)
  cat("  After CheckM join:", ncol(df), "cols\n")
}

# ---- Binary indicators: has_assembly, has_checkm_eval ----
# These capture "has this species been sequenced / evaluated at all?"
df$has_assembly   <- as.integer(!is.na(df$best_contig_count))
df$has_checkm_eval <- as.integer(df$NCBI_ID %in% checkm_quality$NCBI_ID)

cat("\n=== Assembly / CheckM cross-tabulation (before imputation) ===\n")
cat("                       Has CheckM eval   No CheckM eval\n")
cat(sprintf("  Has assembly:         %5d              %5d\n",
    sum(df$has_assembly == 1 & df$has_checkm_eval == 1),
    sum(df$has_assembly == 1 & df$has_checkm_eval == 0)))
cat(sprintf("  No assembly:          %5d              %5d\n",
    sum(df$has_assembly == 0 & df$has_checkm_eval == 1),
    sum(df$has_assembly == 0 & df$has_checkm_eval == 0)))

# ---- 0-imputation for columns where NA → 0 is semantically clear ----
impute_zero <- c("has_complete_genome", "has_type_assembly",
                 "has_type_material", "n_culture_collections")
cat("\n=== 0-imputation (NA → 0 where semantically valid) ===\n")
for (col in impute_zero) {
  if (col %in% names(df)) {
    n_imp <- sum(is.na(df[[col]]))
    if (n_imp > 0) {
      df[[col]][is.na(df[[col]])] <- 0
      cat(sprintf("  %-30s %5d NAs → 0\n", col, n_imp))
    }
  }
}

# ---- Contig count imputation (median + P95 variants) ----
# Species with no assembly get imputed contig count; has_assembly captures the binary signal.
contig_obs <- df$best_contig_count[!is.na(df$best_contig_count)]
contig_median <- median(contig_obs)
contig_p95    <- quantile(contig_obs, 0.95)
cat(sprintf("\n=== Contig count imputation reference ===\n"))
cat(sprintf("  Observed: N=%d, median=%g, P95=%g, max=%g\n",
    length(contig_obs), contig_median, contig_p95, max(contig_obs)))

# contig_imp_median: observed value if available, else median
df$contig_imp_median <- df$best_contig_count
df$contig_imp_median[is.na(df$contig_imp_median)] <- contig_median
# contig_imp_p95: observed value if available, else P95 (high = poor quality → worse than any real assembly)
df$contig_imp_p95 <- df$best_contig_count
df$contig_imp_p95[is.na(df$contig_imp_p95)] <- contig_p95
cat(sprintf("  Imputed %d NAs with median=%g and P95=%g\n",
    sum(is.na(df$best_contig_count)), contig_median, contig_p95))

# ---- CheckM imputation (corrected) ----
# Three situations:
#   1. Has CheckM evaluation → use real value (already present)
#   2. Has assembly but no CheckM → median-impute (genome exists, just not evaluated)
#   3. No assembly, no CheckM → 0 (no genome data at all)
checkm_obs <- df$best_checkm_completeness[df$has_checkm_eval == 1 &
                                            !is.na(df$best_checkm_completeness)]
checkm_median <- median(checkm_obs, na.rm = TRUE)

n_asm_no_ckm <- sum(df$has_assembly == 1 & df$has_checkm_eval == 0)
n_no_asm     <- sum(df$has_assembly == 0 & is.na(df$best_checkm_completeness))

cat(sprintf("\n=== CheckM imputation (corrected) ===\n"))
cat(sprintf("  Observed CheckM: N=%d, median=%.1f\n", length(checkm_obs), checkm_median))
cat(sprintf("  Has assembly, no CheckM: %d → median-impute (%.1f)\n", n_asm_no_ckm, checkm_median))
cat(sprintf("  No assembly, no CheckM:  %d → 0\n", n_no_asm))

# checkm_imp: corrected imputation
df$checkm_imp <- df$best_checkm_completeness
# Case 2: has assembly but no CheckM → median
idx_asm_no_ckm <- df$has_assembly == 1 & df$has_checkm_eval == 0 & is.na(df$checkm_imp)
df$checkm_imp[idx_asm_no_ckm] <- checkm_median
# Case 3: no assembly → 0
idx_no_asm <- df$has_assembly == 0 & is.na(df$checkm_imp)
df$checkm_imp[idx_no_asm] <- 0
# Any remaining NAs (e.g., no assembly but has CheckM with NA) → 0
df$checkm_imp[is.na(df$checkm_imp)] <- 0

cat(sprintf("  checkm_imp: %d non-NA (should be %d)\n",
    sum(!is.na(df$checkm_imp)), nrow(df)))

# ---- Derived variables ----
df$species_age <- 2026L - df$year_valid_pub
df$log_species_age <- log1p(df$species_age)
df$log_best_contig_count <- log1p(df$best_contig_count)
df$has_type_material <- as.integer(df$has_type_material)
df$has_complete_genome <- as.integer(df$has_complete_genome)
df$has_type_assembly <- as.integer(df$has_type_assembly)

cat("\n=== Coverage rates (after imputation) ===\n")
report_cols <- c("year_valid_pub", "species_age",
                 "n_culture_collections", "has_type_material",
                 "has_assembly", "has_checkm_eval",
                 "best_contig_count", "contig_imp_median", "contig_imp_p95",
                 "has_complete_genome",
                 "best_checkm_completeness", "checkm_imp")
for (col in report_cols) {
  if (col %in% names(df)) {
    n_ok <- sum(!is.na(df[[col]]))
    cat(sprintf("  %-30s %5d / %d (%.1f%%)\n", col, n_ok, nrow(df),
                100 * n_ok / nrow(df)))
  }
}

cat(sprintf("\nColumns: %d → %d (+%d new)\n", n_before, ncol(df), ncol(df) - n_before))

# ---- Save ----
out_path <- file.path(data_dir, "bugphyzz_with_proxies.rds")
saveRDS(df, out_path)
cat("Saved:", out_path, "\n")

cat("\nPhase 2b complete.\n")
