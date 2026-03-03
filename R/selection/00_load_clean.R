#!/usr/bin/env Rscript
# Phase 0: Data Loading & Cleaning
# Reads bugphyzz.xlsx, converts types, creates derived columns, saves cleaned RDS.

library(readxl)
library(dplyr)
library(stringr)

# ---- Paths ----
# Resolve project root: parent of R/ directory where this script lives
this_script <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error = function(e) {
    # Fallback: use commandArgs to find script path
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg))
    else getwd()
  }
)
proj_dir <- normalizePath(file.path(dirname(this_script), "..", ".."))
data_dir  <- file.path(proj_dir, "data")
input_file <- file.path(proj_dir, "bugphyzz.xlsx")

cat("Project dir:", proj_dir, "\n")
cat("Reading:", input_file, "\n")

# ---- 1. Read Excel ----
df <- read_xlsx(input_file)
cat("Raw dimensions:", nrow(df), "x", ncol(df), "\n")

# ---- 2. Define column groups ----
binary_phenotypes <- c(
  "Motility", "Extreme environment tolerance", "Biofilm formation",
  "Animal pathogenicity", "Health association", "Host association",
  "Plant pathogenicity", "Spore formation"
)

categorical_phenotypes <- c(
  "Gram staining", "Aerophilicity", "Biosafety level",
  "Hemolysis", "Cell shape"
)

all_phenotypes <- c(binary_phenotypes, categorical_phenotypes)

taxonomy_cols <- c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus")

# ---- 3. Convert string "NA" → proper NA across phenotype columns ----
for (col in all_phenotypes) {
  df[[col]] <- ifelse(df[[col]] == "NA", NA_character_, df[[col]])
}

# Verify: no residual string "NA"
for (col in all_phenotypes) {
  n_string_na <- sum(df[[col]] == "NA", na.rm = TRUE)
  if (n_string_na > 0) stop("Residual string 'NA' found in column: ", col)
}
cat("String 'NA' conversion: OK\n")

# ---- 4. Convert binary phenotypes: "TRUE"/"FALSE" → logical ----
for (col in binary_phenotypes) {
  df[[col]] <- as.logical(df[[col]])
}

# ---- 5. Convert categorical phenotypes → factor ----
for (col in categorical_phenotypes) {
  df[[col]] <- factor(df[[col]])
}

# ---- 6. Convert NCBI_ID to integer ----
df$NCBI_ID <- as.integer(df$NCBI_ID)

# ---- 6b. Resolve merged NCBI taxids ----
# Some taxids in bugphyzz have been merged into other species by NCBI since
# the dataset was compiled. We fetch the current taxid for each, cache the
# remap table, and either update the ID (simple rename) or merge phenotype
# annotations into the surviving row (duplicate).
remap_file <- file.path(data_dir, "taxid_remap.rds")

if (!file.exists(remap_file)) {
  cat("\n=== Building taxid remap table (one-time NCBI fetch) ===\n")
  library(rentrez)
  library(xml2)

  all_ids <- unique(df$NCBI_ID[!is.na(df$NCBI_ID)])
  remap_rows <- list()
  batch_size <- 200

  for (i in seq(1, length(all_ids), by = batch_size)) {
    batch <- all_ids[i:min(i + batch_size - 1, length(all_ids))]
    xml_text <- tryCatch(
      entrez_fetch(db = "taxonomy", id = batch, rettype = "xml"),
      error = function(e) NULL
    )
    if (is.null(xml_text)) { Sys.sleep(1); next }
    doc <- read_xml(xml_text)
    taxons <- xml_find_all(doc, "/TaxaSet/Taxon")
    for (taxon in taxons) {
      current_id <- as.integer(xml_text(xml_find_first(taxon, "TaxId")))
      sci_name <- xml_text(xml_find_first(taxon, "ScientificName"))
      aka_nodes <- xml_find_all(taxon, ".//AkaTaxIds/TaxId")
      aka_ids <- as.integer(xml_text(aka_nodes))
      # Only record entries where an old ID differs from the current one
      old_in_data <- intersect(aka_ids, all_ids)
      for (old_id in old_in_data) {
        remap_rows[[as.character(old_id)]] <- data.frame(
          old_id = old_id, new_id = current_id, new_name = sci_name,
          stringsAsFactors = FALSE
        )
      }
    }
    if ((i %/% batch_size) %% 20 == 0)
      cat(sprintf("  Scanned %d / %d IDs\n", min(i + batch_size - 1, length(all_ids)), length(all_ids)))
    Sys.sleep(0.4)
  }

  if (length(remap_rows) > 0) {
    remap <- do.call(rbind, remap_rows)
    rownames(remap) <- NULL
  } else {
    remap <- data.frame(old_id = integer(0), new_id = integer(0),
                        new_name = character(0), stringsAsFactors = FALSE)
  }
  saveRDS(remap, remap_file)
  cat("  Saved:", remap_file, "(", nrow(remap), "merged taxids)\n")
} else {
  remap <- readRDS(remap_file)
  cat("\nLoaded taxid remap:", nrow(remap), "merged taxids\n")
}

if (nrow(remap) > 0) {
  # Classify: does the new_id already exist in df?
  remap$is_dup <- remap$new_id %in% df$NCBI_ID

  # ---- Simple renames (new_id not already in df) ----
  renames <- remap[!remap$is_dup, ]
  if (nrow(renames) > 0) {
    for (i in seq_len(nrow(renames))) {
      idx <- which(df$NCBI_ID == renames$old_id[i])
      if (length(idx) == 1) {
        df$NCBI_ID[idx] <- renames$new_id[i]
        df[["Binomial name"]][idx] <- renames$new_name[i]
      }
    }
    cat(sprintf("  Renamed %d taxids (no duplicate)\n", nrow(renames)))
  }

  # ---- Duplicates: merge phenotype annotations into surviving row, drop old ----
  dupes <- remap[remap$is_dup, ]
  if (nrow(dupes) > 0) {
    rows_to_drop <- integer(0)
    fills <- 0L
    for (i in seq_len(nrow(dupes))) {
      old_idx <- which(df$NCBI_ID == dupes$old_id[i])
      new_idx <- which(df$NCBI_ID == dupes$new_id[i])
      if (length(old_idx) != 1 || length(new_idx) != 1) next
      # Fill NAs in the surviving row from the old row (no overwrite of existing values)
      for (p in all_phenotypes) {
        if (is.na(df[[p]][new_idx]) && !is.na(df[[p]][old_idx])) {
          df[[p]][new_idx] <- df[[p]][old_idx]
          fills <- fills + 1L
        }
      }
      rows_to_drop <- c(rows_to_drop, old_idx)
    }
    if (length(rows_to_drop) > 0) {
      df <- df[-rows_to_drop, , drop = FALSE]
      rownames(df) <- NULL
    }
    cat(sprintf("  Merged %d duplicate rows (%d phenotype values filled)\n",
                length(rows_to_drop), fills))
  }

  cat(sprintf("  Rows after taxid cleanup: %d\n", nrow(df)))
}

# ---- 7. Convert taxonomy columns → factor (ordered by frequency) ----
for (col in taxonomy_cols) {
  freq_order <- names(sort(table(df[[col]]), decreasing = TRUE))
  df[[col]] <- factor(df[[col]], levels = freq_order)
}

# ---- 8. Convert "Member of WA subset" to logical ----
df[["Member of WA subset"]] <- as.logical(df[["Member of WA subset"]])

# ---- 9. Derived columns ----
# Phylum:Order interaction
df$phylum_order <- paste(df$Phylum, df$Order, sep = ":")

# Number of annotated phenotypes per species
df$n_phenotypes_annotated <- rowSums(!is.na(df[, all_phenotypes]))


# ---- 10. Validation checks ----
cat("\n=== Validation ===\n")

# Row count
cat("Total rows:", nrow(df), "\n")
stopifnot(nrow(df) > 15000)

# Spore formation has TRUE, FALSE, and NA
spore_vals <- unique(df[["Spore formation"]])
stopifnot(TRUE %in% spore_vals)
stopifnot(FALSE %in% spore_vals)
stopifnot(any(is.na(spore_vals)))
cat("Spore formation values: TRUE/FALSE/NA present\n")

# Taxonomy NA check
for (col in taxonomy_cols) {
  n_na <- sum(is.na(df[[col]]))
  cat(sprintf("  %s: %d NAs (%.1f%%)\n", col, n_na, 100 * n_na / nrow(df)))
}

# Cross-tabulate: Phylum × Spore formation annotation rate
cat("\nPhylum × Spore formation annotation rate (top 15 phyla):\n")
phylum_spore <- df %>%
  group_by(Phylum) %>%
  summarise(
    n = n(),
    n_annotated = sum(!is.na(`Spore formation`)),
    pct_annotated = round(100 * n_annotated / n, 1),
    n_true = sum(`Spore formation` == TRUE, na.rm = TRUE),
    pct_true_among_annotated = round(100 * n_true / pmax(n_annotated, 1), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n)) %>%
  head(15)
print(as.data.frame(phylum_spore), row.names = FALSE)

# Annotation depth summary
cat("\nAnnotation depth (n_phenotypes_annotated) summary:\n")
print(summary(df$n_phenotypes_annotated))

# ---- 11. Save ----
out_path <- file.path(data_dir, "bugphyzz_clean.rds")
saveRDS(df, out_path)
cat("\nSaved:", out_path, "\n")
cat("Done.\n")
