#!/usr/bin/env Rscript
# Phase 2: Selection Proxy Acquisition
# PubMed publication counts (by species name), assembly counts (by NCBI_ID),
# Wikipedia article size (by species name).
# All queries are cached incrementally to survive interruptions.

library(dplyr)
library(rentrez)
library(jsonlite)
library(httr)

# ---- Paths ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg)) else getwd()
proj_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))
data_dir <- file.path(proj_dir, "data")

df <- readRDS(file.path(data_dir, "bugphyzz_clean.rds"))
cat("Loaded:", nrow(df), "rows\n")

batch_size <- 500  # save cache every 500 single-key queries

# ---- Helper: single-key query with caching ----
query_with_cache <- function(keys, cache_file, query_fn, desc, sleep = 0.35) {
  if (file.exists(cache_file)) {
    cache <- readRDS(cache_file)
    cat(sprintf("[%s] Loaded cache: %d entries\n", desc, length(cache)))
  } else {
    cache <- integer(0)
    names(cache) <- character(0)
  }

  # Retry any previously-failed queries (NA values in cache)
  na_keys <- names(cache[is.na(cache)])
  cache <- cache[!is.na(cache)]
  remaining <- setdiff(keys, names(cache))
  n_retry <- length(intersect(na_keys, remaining))
  cat(sprintf("[%s] Total: %d, Cached: %d, Remaining: %d (incl. %d retries)\n",
              desc, length(keys), length(keys) - length(remaining), length(remaining), n_retry))

  if (length(remaining) == 0) return(cache)

  n_batches <- ceiling(length(remaining) / batch_size)
  for (b in seq_len(n_batches)) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, length(remaining))
    batch <- remaining[start_idx:end_idx]

    for (key in batch) {
      cache[key] <- tryCatch(query_fn(key), error = function(e) NA_integer_)
      Sys.sleep(sleep)
    }

    saveRDS(cache, cache_file)
    cat(sprintf("  [%s] Batch %d/%d done (%d cached)\n", desc, b, n_batches, length(cache)))
  }
  cache
}

# ---- Helper: batched Wikipedia query with caching ----
query_wikipedia_batched <- function(species_names, cache_file, batch_size_wp = 50) {
  if (file.exists(cache_file)) {
    cache <- readRDS(cache_file)
    cat(sprintf("[Wikipedia] Loaded cache: %d entries\n", length(cache)))
  } else {
    cache <- integer(0)
    names(cache) <- character(0)
  }

  remaining <- setdiff(species_names, names(cache))
  cat(sprintf("[Wikipedia] Total: %d, Cached: %d, Remaining: %d\n",
              length(species_names), length(species_names) - length(remaining),
              length(remaining)))

  if (length(remaining) == 0) return(cache)

  n_batches <- ceiling(length(remaining) / batch_size_wp)
  save_every <- 10  # save cache every 10 API batches (~500 species)

  for (b in seq_len(n_batches)) {
    start_idx <- (b - 1) * batch_size_wp + 1
    end_idx <- min(b * batch_size_wp, length(remaining))
    batch <- remaining[start_idx:end_idx]

    titles_str <- paste(batch, collapse = "|")

    result <- tryCatch({
      resp <- GET("https://en.wikipedia.org/w/api.php",
                  query = list(
                    action = "query",
                    titles = titles_str,
                    prop = "revisions",
                    rvprop = "size",
                    format = "json",
                    formatversion = "2"
                  ))
      parsed <- fromJSON(content(resp, "text", encoding = "UTF-8"))
      pages <- parsed$query$pages

      # Extract title → size mapping
      for (i in seq_len(nrow(pages))) {
        title <- pages$title[i]
        missing <- if ("missing" %in% names(pages)) pages$missing[i] else FALSE
        if (isTRUE(missing) || is.null(pages$revisions[[i]])) {
          cache[title] <- 0L  # no article = 0 bytes
        } else {
          cache[title] <- as.integer(pages$revisions[[i]]$size[1])
        }
      }
      TRUE
    }, error = function(e) {
      # On error, mark entire batch as 0 (conservative: no data = no article)
      for (sp in batch) cache[sp] <<- 0L
      cat(sprintf("  [Wikipedia] Batch %d error: %s\n", b, conditionMessage(e)))
      FALSE
    })

    Sys.sleep(1)  # respect rate limits

    if (b %% save_every == 0 || b == n_batches) {
      saveRDS(cache, cache_file)
      cat(sprintf("  [Wikipedia] Batch %d/%d done (%d cached)\n",
                  b, n_batches, length(cache)))
    }
  }
  cache
}

# ---- Test connectivity ----
cat("\nTesting NCBI API... ")
ncbi_ok <- tryCatch({
  res <- entrez_search(db = "pubmed",
                       term = "\"Escherichia coli\"[Organism]",
                       retmax = 0)
  cat("OK (E. coli:", res$count, "papers)\n")
  TRUE
}, error = function(e) {
  cat("FAILED:", conditionMessage(e), "\n")
  FALSE
})

cat("Testing Wikipedia API... ")
wiki_ok <- tryCatch({
  resp <- GET("https://en.wikipedia.org/w/api.php",
              query = list(action = "query", titles = "Escherichia coli",
                           prop = "revisions", rvprop = "size",
                           format = "json", formatversion = "2"))
  parsed <- fromJSON(content(resp, "text", encoding = "UTF-8"))
  sz <- parsed$query$pages$revisions[[1]]$size[1]
  cat("OK (E. coli article:", sz, "bytes)\n")
  TRUE
}, error = function(e) {
  cat("FAILED:", conditionMessage(e), "\n")
  FALSE
})

# ---- Unique species names ----
species_names <- unique(df[["Binomial name"]])
species_names <- species_names[!is.na(species_names)]

# ================================================================
# Proxy 1: PubMed publication count (keyed by Binomial name)
# ================================================================
if (ncbi_ok) {
  pub_cache <- query_with_cache(
    keys = species_names,
    cache_file = file.path(data_dir, "pubmed_counts.rds"),
    query_fn = function(name) {
      res <- entrez_search(db = "pubmed",
                           term = paste0("\"", name, "\"[Organism]"),
                           retmax = 0)
      as.integer(res$count)
    },
    desc = "PubMed"
  )
  df$pub_count <- pub_cache[df[["Binomial name"]]]
  df$log_pub_count <- log1p(df$pub_count)
} else {
  cat("\nNCBI unavailable. Skipping PubMed + assembly proxies.\n")
  df$pub_count <- NA_integer_
  df$log_pub_count <- NA_real_
}

# ================================================================
# Proxy 2: Genome assembly count (keyed by NCBI_ID, txid format)
# ================================================================
if (ncbi_ok) {
  all_ids <- unique(as.character(df$NCBI_ID))
  all_ids <- all_ids[!is.na(all_ids) & all_ids != "NA"]

  asm_cache <- query_with_cache(
    keys = all_ids,
    cache_file = file.path(data_dir, "assembly_counts.rds"),
    query_fn = function(tid) {
      res <- entrez_search(db = "assembly",
                           term = paste0("txid", tid, "[Organism]"),
                           retmax = 0)
      as.integer(res$count)
    },
    desc = "Assembly"
  )
  df$assembly_count <- asm_cache[as.character(df$NCBI_ID)]
  df$log_assembly_count <- log1p(df$assembly_count)
} else {
  df$assembly_count <- NA_integer_
  df$log_assembly_count <- NA_real_
}

# ================================================================
# Proxy 3: Wikipedia article size (keyed by Binomial name)
# ================================================================
if (wiki_ok) {
  wiki_cache <- query_wikipedia_batched(
    species_names = species_names,
    cache_file = file.path(data_dir, "wikipedia_sizes.rds"),
    batch_size_wp = 50
  )
  df$wiki_size <- wiki_cache[df[["Binomial name"]]]
  df$log_wiki_size <- log1p(as.numeric(df$wiki_size))
} else {
  cat("\nWikipedia unavailable. Skipping Wikipedia proxy.\n")
  df$wiki_size <- NA_integer_
  df$log_wiki_size <- NA_real_
}

# ================================================================
# Summary & save
# ================================================================
cat("\n=== Proxy summary ===\n")
proxy_cols <- c("pub_count", "assembly_count", "wiki_size")
for (col in proxy_cols) {
  vals <- df[[col]]
  cat(sprintf("%-18s  %d non-NA, median = %s, mean = %s\n",
              paste0(col, ":"),
              sum(!is.na(vals)),
              ifelse(all(is.na(vals)), "NA", as.character(median(vals, na.rm = TRUE))),
              ifelse(all(is.na(vals)), "NA", round(mean(vals, na.rm = TRUE), 1))))
}
cat("n_phenotypes:      median =", median(df$n_phenotypes_annotated),
    ", mean =", round(mean(df$n_phenotypes_annotated), 1), "\n")

out_path <- file.path(data_dir, "bugphyzz_with_proxies.rds")
saveRDS(df, out_path)
cat("\nSaved:", out_path, "\n")

# ========================================================================
# Figure 02b: Documentation Bias Landscape
# ========================================================================
cat("\n=== Documentation bias landscape ===\n")

fig_dir <- file.path(proj_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Load config for ALL_PHENOTYPES
source(file.path(proj_dir, "R", "selection", "00_config.R"))

# Annotation completeness: fraction of phenotypes annotated per species
df$annotation_completeness <- df$n_phenotypes_annotated / length(ALL_PHENOTYPES)

# Only plot if at least two proxies are available
has_pub <- !all(is.na(df$log_pub_count))
has_asm <- !all(is.na(df$log_assembly_count))

if (has_pub && has_asm) {
  library(ggplot2)
  library(cowplot)

  plot_df <- df[!is.na(df$log_pub_count) & !is.na(df$log_assembly_count), ]

  panel_a <- ggplot(plot_df, aes(x = log_pub_count, y = log_assembly_count,
                                  color = annotation_completeness)) +
    geom_point(alpha = 0.3, size = 0.6) +
    scale_color_viridis_c(option = "plasma", name = "Annotation\ncompleteness") +
    labs(x = "log(PubMed count + 1)", y = "log(Assembly count + 1)",
         title = "A) Proxy space") +
    theme_minimal(base_size = 11)

  panel_b <- ggplot(plot_df, aes(x = log_pub_count, y = annotation_completeness)) +
    geom_point(alpha = 0.1, size = 0.4, color = "grey50") +
    geom_smooth(method = "loess", color = "steelblue", fill = "steelblue",
                alpha = 0.2, se = TRUE) +
    labs(x = "log(PubMed count + 1)", y = "Annotation completeness",
         title = "B) PubMed vs completeness") +
    theme_minimal(base_size = 11)

  panel_c <- ggplot(plot_df, aes(x = log_assembly_count, y = annotation_completeness)) +
    geom_point(alpha = 0.1, size = 0.4, color = "grey50") +
    geom_smooth(method = "loess", color = "tomato", fill = "tomato",
                alpha = 0.2, se = TRUE) +
    labs(x = "log(Assembly count + 1)", y = "Annotation completeness",
         title = "C) Assembly vs completeness") +
    theme_minimal(base_size = 11)

  fig02b <- plot_grid(panel_a, panel_b, panel_c, ncol = 3, rel_widths = c(1.3, 1, 1))
  ggsave(file.path(fig_dir, "fig02b_documentation_bias_landscape.pdf"),
         fig02b, width = 16, height = 5)
  cat("Saved: fig02b_documentation_bias_landscape.pdf\n")
} else {
  cat("Skipping fig02b (need both PubMed and assembly proxies)\n")
}

cat("Done.\n")
