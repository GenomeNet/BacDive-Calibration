read_env_file <- function(path) {
  if (!file.exists(path)) {
    return(list())
  }

  lines <- readLines(path, warn = FALSE)
  out <- list()

  for (line in lines) {
    line <- trimws(line)
    if (line == "" || startsWith(line, "#")) {
      next
    }
    line <- sub("^export\\s+", "", line)
    if (!grepl("=", line, fixed = TRUE)) {
      next
    }

    kv <- strsplit(line, "=", fixed = TRUE)[[1]]
    key <- trimws(kv[1])
    val <- trimws(paste(kv[-1], collapse = "="))
    val <- sub("^\"(.*)\"$", "\\1", val)
    val <- sub("^'(.*)'$", "\\1", val)
    out[[key]] <- val
  }

  out
}

resolve_path <- function(x, proj_dir) {
  if (is.null(x) || x == "") {
    return(x)
  }
  if (grepl("^/", x)) {
    return(x)
  }
  normalizePath(file.path(proj_dir, x), mustWork = FALSE)
}

load_bacdive_config <- function(proj_dir,
                                required = character(),
                                require_selection = FALSE) {
  proj_dir <- normalizePath(proj_dir, mustWork = FALSE)
  env_file <- Sys.getenv("BACDIVE_ENV_FILE", unset = file.path(proj_dir, ".env"))
  env_map <- read_env_file(env_file)

  get_value <- function(key, default = NULL, required = FALSE) {
    val <- Sys.getenv(key, unset = "")
    if (val == "" && !is.null(env_map[[key]])) {
      val <- env_map[[key]]
    }
    if (val == "" && !is.null(default)) {
      val <- default
    }
    if (required && (is.null(val) || val == "")) {
      stop(sprintf("Missing required config key: %s", key), call. = FALSE)
    }
    val
  }

  reports_dir <- resolve_path(get_value("REPORTS_DIR", "."), proj_dir)

  cfg <- list(
    project_dir = proj_dir,
    reports_dir = reports_dir,
    output_dir = resolve_path(get_value("OUTPUT_DIR", file.path(reports_dir, "output")), proj_dir),
    figures_dir = resolve_path(get_value("FIGURES_DIR", file.path(reports_dir, "figures")), proj_dir),
    slurm_output_dir = resolve_path(get_value("SLURM_OUTPUT_DIR", file.path(reports_dir, "slurm_output")), proj_dir),
    bd_pred_csv = resolve_path(get_value("BD_PRED_CSV"), proj_dir),
    bd_labels_csv = resolve_path(get_value("BD_LABELS_CSV"), proj_dir),
    bd_splits_csv = resolve_path(get_value("BD_SPLITS_CSV"), proj_dir),
    genome_family_map_csv = resolve_path(get_value("GENOME_FAMILY_MAP_CSV",
                                                   file.path(reports_dir, "genome_family_map.csv")), proj_dir),
    selection_probs_rds = resolve_path(get_value("SELECTION_PROBS_RDS"), proj_dir),
    bugphyzz_proxies_rds = resolve_path(get_value("BUGPHYZZ_PROXIES_RDS"), proj_dir),
    count_list_rds = resolve_path(get_value("COUNT_LIST_RDS"), proj_dir),
    bacdive_complete_rds = resolve_path(get_value("BACDIVE_COMPLETE_RDS"), proj_dir),
    bacdive_test_rds = resolve_path(get_value("BACDIVE_TEST_RDS",
                                              file.path(reports_dir, "bacdive_test.rds")), proj_dir),
    fasta_dir = resolve_path(get_value("FASTA_DIR"), proj_dir),
    microbe_cards_xlsx = resolve_path(get_value("MICROBE_CARDS_XLSX"), proj_dir),
    checkpoints_dir = resolve_path(get_value("CHECKPOINTS_DIR"), proj_dir)
  )

  if (isTRUE(require_selection)) {
    required <- unique(c(required, "selection_probs_rds", "bugphyzz_proxies_rds"))
  }
  required <- unique(required)
  known_keys <- names(cfg)
  unknown <- setdiff(required, known_keys)
  if (length(unknown) > 0) {
    stop(sprintf("Unknown config keys in required=: %s", paste(unknown, collapse = ", ")),
         call. = FALSE)
  }

  key_kinds <- c(
    reports_dir = "dir",
    output_dir = "dir",
    figures_dir = "dir",
    slurm_output_dir = "dir",
    bd_pred_csv = "file",
    bd_labels_csv = "file",
    bd_splits_csv = "file",
    genome_family_map_csv = "file",
    selection_probs_rds = "file",
    bugphyzz_proxies_rds = "file",
    count_list_rds = "file",
    bacdive_complete_rds = "file",
    bacdive_test_rds = "path",
    fasta_dir = "dir",
    microbe_cards_xlsx = "file",
    checkpoints_dir = "dir"
  )

  for (k in c("reports_dir", "output_dir", "figures_dir", "slurm_output_dir")) {
    if (!is.null(cfg[[k]]) && cfg[[k]] != "") {
      dir.create(cfg[[k]], recursive = TRUE, showWarnings = FALSE)
    }
  }

  for (k in required) {
    val <- cfg[[k]]
    if (is.null(val) || val == "") {
      stop(sprintf("Missing required config key: %s", k), call. = FALSE)
    }
    kind <- key_kinds[[k]]
    if (identical(kind, "file") && !file.exists(val)) {
      stop(sprintf("Configured file does not exist: %s -> %s", k, val), call. = FALSE)
    }
    if (identical(kind, "dir") && !dir.exists(val)) {
      stop(sprintf("Configured directory does not exist: %s -> %s", k, val), call. = FALSE)
    }
  }

  cfg
}
