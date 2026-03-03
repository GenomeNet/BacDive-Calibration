# Selection Artifacts (Vendored from `phenomnar`)

This folder contains the minimal scripts needed to reproduce the two artifacts
required by IPW calibration in this repo:

- `data/selection/bugphyzz_with_proxies.rds`
- `data/selection/selection_probs.rds`

Copied scripts:

- `00_config.R`
- `00_load_clean.R`
- `02_selection_proxies.R`
- `02b_extended_proxies.R`
- `03_selection_model.R`

These were adapted to run from `R/selection/` and write into this repo's
`data/`, `output/`, and `figures/` directories.
