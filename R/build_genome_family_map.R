# ── Build genome → genus → family mapping ────────────────────────────────────
# Extracts genus from FASTA headers (bulk bash), joins to microbe.cards S1
# for Family. Saves genome_family_map.csv for reuse.
#
# Usage: conda run -n genome Rscript R/build_genome_family_map.R

library(dplyr)
library(readxl)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg)) else normalizePath(getwd())
script_dir <- if (dir.exists(script_path)) script_path else dirname(script_path)
proj_dir <- if (basename(script_dir) == "R") dirname(script_dir) else script_dir
source(file.path(proj_dir, "R", "config_paths.R"))
cfg <- load_bacdive_config(
  proj_dir,
  required = c("fasta_dir", "bd_splits_csv", "microbe_cards_xlsx")
)
fasta_dir <- cfg$fasta_dir

# All genomes in the dataset
splits <- data.table::fread(
  cfg$bd_splits_csv,
  data.table = FALSE)
all_files <- splits$file
cat("Total genomes:", length(all_files), "\n")

# Bulk extract FASTA headers via bash
tmp_list <- tempfile(fileext = ".txt")
tmp_hdr  <- tempfile(fileext = ".tsv")
writeLines(all_files, tmp_list)

cat("Extracting FASTA headers...\n")
t0 <- proc.time()
cmd <- sprintf(
  "xargs -a %s -P 8 -I{} sh -c 'f=\"%s/{}\"; [ -f \"$f\" ] && printf \"%%s\\t%%s\\n\" \"{}\" \"$(head -1 \"$f\")\"' > %s 2>/dev/null",
  tmp_list, fasta_dir, tmp_hdr)
system(cmd)
cat("Headers extracted in", round((proc.time() - t0)[3], 1), "s\n")

hdr <- data.table::fread(tmp_hdr, header = FALSE, col.names = c("file", "header"),
                          sep = "\t", quote = "", data.table = FALSE)
cat("Headers loaded:", nrow(hdr), "\n")

# Parse genus
parse_genus_vec <- function(files, headers) {
  n <- length(files)
  genus <- character(n)
  for (i in seq_len(n)) {
    f <- files[i]; h <- headers[i]
    if (grepl("^GC[AF]_", f)) {
      sp <- sub("^>[^ ]+ ", "", h)
    } else if (grepl("\\[", h)) {
      sp <- sub(".*\\[", "", h)
      sp <- sub("\\]", "", sp)
      sp <- sub(" \\|.*", "", sp)
    } else {
      sp <- sub("^>[^ ]+ ", "", h)
      sp <- sub(" :.*", "", sp)
    }
    sp <- sub(" strain .*", "", sp)
    sp <- sub(" chromosome.*", "", sp)
    sp <- sub(" plasmid.*", "", sp)
    sp <- sub(" contig.*", "", sp)
    sp <- sub(",.*", "", sp)
    sp <- trimws(sp)
    g <- strsplit(sp, "\\s+")[[1]][1]
    genus[i] <- gsub("[\\[\\]]", "", g)
  }
  genus
}

hdr$genus <- parse_genus_vec(hdr$file, hdr$header)

# Join to xlsx for Family
s1 <- read_xlsx(cfg$microbe_cards_xlsx)
genus_tax <- s1 %>% select(Genus, Family, Order, Phylum) %>% distinct(Genus, .keep_all = TRUE)

hdr <- hdr %>% left_join(genus_tax, by = c("genus" = "Genus"))

cat("Family matched:", sum(!is.na(hdr$Family)), "/", nrow(hdr),
    "(", round(100 * sum(!is.na(hdr$Family)) / nrow(hdr), 1), "%)\n")
cat("Order matched:", sum(!is.na(hdr$Order)), "/", nrow(hdr),
    "(", round(100 * sum(!is.na(hdr$Order)) / nrow(hdr), 1), "%)\n")
cat("Unique genera:", length(unique(hdr$genus)), "\n")
cat("Unique families:", length(unique(hdr$Family[!is.na(hdr$Family)])), "\n")
cat("Unique orders:", length(unique(hdr$Order[!is.na(hdr$Order)])), "\n")
cat("Unique phyla:", length(unique(hdr$Phylum[!is.na(hdr$Phylum)])), "\n")

# Save
out <- hdr[, c("file", "genus", "Family", "Order", "Phylum")]
write.csv(out, cfg$genome_family_map_csv, row.names = FALSE)
cat("Saved:", cfg$genome_family_map_csv, "\n")
