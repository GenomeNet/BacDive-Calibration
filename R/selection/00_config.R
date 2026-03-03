# 00_config.R — Shared configuration for the MNAR pipeline.
# Source this file from any phase script:
#   source(file.path(proj_dir, "R", "selection", "00_config.R"))

# ---- Column groups ----
BINARY_PHENOTYPES <- c(
  "Motility", "Extreme environment tolerance", "Biofilm formation",
  "Animal pathogenicity", "Health association", "Host association",
  "Plant pathogenicity", "Spore formation"
)

CATEGORICAL_PHENOTYPES <- c(
  "Gram staining", "Aerophilicity", "Biosafety level",
  "Hemolysis", "Cell shape"
)

ALL_PHENOTYPES <- c(BINARY_PHENOTYPES, CATEGORICAL_PHENOTYPES)

TAXONOMY_COLS <- c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Short names for plotting
SHORT_NAMES <- c(
  "Motility" = "Motil", "Extreme environment tolerance" = "ExtEnv",
  "Biofilm formation" = "Biofilm", "Animal pathogenicity" = "AnimPath",
  "Health association" = "Health", "Host association" = "Host",
  "Plant pathogenicity" = "PlantPath", "Spore formation" = "Spore",
  "Gram staining" = "Gram", "Aerophilicity" = "Aero",
  "Biosafety level" = "BSL", "Hemolysis" = "Hemol", "Cell shape" = "Shape"
)

# Minimum annotation rate for a phenotype to be analysed in Phases 3/5/6
MIN_ANNOTATION_RATE <- 0.02  # skip if <2% annotated

# IPW (inverse probability weighting) settings
IPW_TRIM_BOUNDS <- c(0.01, 0.99)  # propensity score trimming
IPW_BOOT_REPS   <- 2000           # bootstrap replications for IPW CIs

# ---- Helpers ----
resolve_paths <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg)) else getwd()
  proj <- normalizePath(file.path(dirname(script_path), "..", ".."))
  list(
    proj_dir = proj,
    data_dir = file.path(proj, "data"),
    fig_dir  = file.path(proj, "figures"),
    out_dir  = file.path(proj, "output")
  )
}

collapse_rare <- function(x, min_count = 50, other_label = "Other") {
  counts <- table(x)
  rare <- names(counts[counts < min_count])
  x_out <- as.character(x)
  x_out[x_out %in% rare] <- other_label
  factor(x_out)
}

expit <- function(x) 1 / (1 + exp(-x))

# Parse --phenotype= CLI argument.
# pool: "binary" (default for phases 5/6), "all" (phases 1/3), or a character vector
# Returns phenotype names above the annotation threshold.
parse_phenotype_args <- function(df, pool = "binary") {
  candidate <- if (is.character(pool) && length(pool) == 1) {
    switch(pool, binary = BINARY_PHENOTYPES, all = ALL_PHENOTYPES, BINARY_PHENOTYPES)
  } else pool

  args <- commandArgs(trailingOnly = TRUE)
  pheno_arg <- grep("^--phenotype=", args, value = TRUE)
  if (length(pheno_arg)) {
    phenos <- unlist(strsplit(sub("^--phenotype=", "", pheno_arg), ","))
    phenos <- trimws(phenos)
    bad <- setdiff(phenos, ALL_PHENOTYPES)
    if (length(bad)) warning("Unknown phenotype(s): ", paste(bad, collapse = ", "))
    candidate <- intersect(phenos, candidate)
  }
  rates <- sapply(candidate, function(p) mean(!is.na(df[[p]])))
  candidate[rates >= MIN_ANNOTATION_RATE]
}

# ---- Proxy columns ----
# Original Phase 2 proxies (research attention axis)
PROXY_COLS_ORIGINAL <- c("log_pub_count", "log_assembly_count", "log_wiki_size")

# Extended Phase 2b proxies (historical, accessibility, genomic quality axes)
# log_species_age: NA for 209 species with unparseable authority year (dropped)
# log_best_contig_count: NA for 2066 species with no assembly (use imputed variants)
PROXY_COLS_EXTENDED <- c("log_species_age", "n_culture_collections",
                         "has_complete_genome", "log_best_contig_count",
                         "best_checkm_completeness")

# All proxy columns for M2 selection model (complete-case version)
PROXY_COLS <- c(PROXY_COLS_ORIGINAL, PROXY_COLS_EXTENDED)

# Imputed proxy set: replaces contig + checkm with imputed versions + indicators.
# No NAs except log_species_age (209 species), so modeling uses ~18,827 rows.
PROXY_COLS_IMPUTED <- c("log_pub_count", "log_assembly_count", "log_wiki_size",
                        "log_species_age", "n_culture_collections",
                        "has_complete_genome", "has_assembly", "has_checkm_eval",
                        "checkm_imp")

# Convenience: safe filename / column name from phenotype name
safe_pheno_name <- function(pheno) gsub(" ", "_", tolower(pheno))

# Named mapping: original name → safe name (for mice column renaming)
SAFE_PHENO <- setNames(safe_pheno_name(ALL_PHENOTYPES), ALL_PHENOTYPES)
# Reverse mapping: safe name → original name
ORIG_PHENO <- setNames(ALL_PHENOTYPES, safe_pheno_name(ALL_PHENOTYPES))
