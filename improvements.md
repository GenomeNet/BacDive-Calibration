• Current state and a reliable implementation plan:

  What I Found

  1. /home/yhan/bacdive is not under git yet (.git missing), while legacy BacDiveComplete is its own git repo (/home/yhan/bacdive/
     BacDiveComplete/.git).
  2. bacdive scripts currently contain 49 absolute-path literals (highest density in calibration_ipw_platt.R:27).
  3. phenomnar is already mostly portable (relative path resolution via config), but includes broad MNAR phases not needed for calibration-
     only workflows (README:26).
  4. Calibration currently depends on two external phenomnar artifacts (selection_probs.rds read:59, bugphyzz_with_proxies.rds read:60).
  5. Deployment logic in legacy code is centralized in eval_non_train_self_pad.R:149, which is the right insertion point for calibrated
     probabilities.

  Recommended Execution Plan

  1. Initialize version control first.
     Deliverables: init git in /home/yhan/bacdive, add .gitignore for generated PDFs/RDS/CSV, first baseline commit.
     Exit criterion: clean reproducible baseline before large refactors.
  2. Centralize all paths into one env-driven config layer.
     Deliverables: conf/.env.example + R/config_paths.R with strict validation (must_exist, must_writable), and scripts reading only config
     keys.
     Exit criterion: zero absolute paths left in executable scripts.
     Suggested env keys: BD_PRED_CSV, BD_LABELS_CSV, BD_SPLITS_CSV, BD_FASTA_DIR, BD_CHECKPOINT_DIR, SELECTION_PROBS_RDS,
     BUGPHYZZ_PROXIES_RDS, MICROBE_CARDS_XLSX, OUTPUT_DIR, FIGURES_DIR.
  3. Split “core calibration” from “peripheral analysis”.
     Deliverables: folders scripts/ (entrypoints), R/ (shared functions), analysis/ (standalone exploratory work), docs/.
     Exit criterion: calibration pipeline runnable without touching peripheral scripts.
  4. Make selection-model dependency end-to-end but minimal.
     Deliverables: add a selection/ module with only required steps to produce the two artifacts needed by calibration; keep MNAR phases
     optional.
     Exit criterion: one command can produce/validate selection_probs.rds and bugphyzz_with_proxies.rds without manual copying.
     Reliability guard: keep “fast mode” (consume frozen artifacts) and “full rebuild mode” (recompute proxies/models), both checksum-
     validated.
  5. Implement multiclass calibration properly.
     Deliverables: add multiclass calibrators for gram/cellshape/flagellum in LOFO pipeline with model selection by held-out NLL/ECE.
     Recommended order: temperature baseline -> vector/matrix scaling -> Dirichlet calibration.
     Exit criterion: per-head method chosen by cross-validated metric, stored as deployable parameter files.
     Research basis: Guo et al. (temperature scaling) and Kull et al. (Dirichlet calibration).
     Sources: https://proceedings.mlr.press/v70/guo17a.html and

  https://papers.nips.cc/paper/9397-beyond-temperature-scaling-obtaining-well-calibrated-multi-class-probabilities-with-dirichlet-calibration
  6. Deploy calibration into BacDiveComplete with backward compatibility.
     Deliverables: patch eval_non_train_self_pad.R:149 to load calibration params and emit calibrated outputs (new columns), keep legacy
     _prob temporarily.
     Exit criterion: side-by-side output parity checks + feature flag for rollback.
  7. Add test and release gates.
     Deliverables: small fixture dataset, unit tests for transforms, integration smoke test for full pipeline, schema checks for required
     columns.
     Exit criterion: CI/local gate passes before any deploy artifact generation.

  Pragmatic order to implement

  1. Git bootstrap.
  2. Env/config refactor (vertical slice on one script, then roll out).
  3. Repo structure cleanup.
  4. Minimal selection module integration.
  5. Multiclass calibration implementation.
  6. BacDiveComplete deployment patch.
  7. Freeze release with runbook.