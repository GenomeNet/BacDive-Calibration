# Implementation Plan

## Recommended Repo Name

`bacdive-calibration`

Why: it is specific, short, and consistent with the production objective (calibration + deployment readiness).

## Recommended Repository Split

Use **two repos**, not one:

1. `bacdive-calibration` (production)
- Selection artifact builder (minimal, deterministic subset only)
- Calibration methods (LOFO/IPW + future multiclass calibration)
- Deployment adapters for BacDive integration output

2. `phenomnar` (research)
- MNAR sensitivity and exploratory analyses
- Broader proxy/model experiments not required for production calibration

Rationale:
- Selection probabilities are required for end-to-end production calibration.
- MNAR exploratory phases are valuable but should not block production runs.
- This keeps production reproducible while preserving research flexibility.

## Current Status (Completed)

- Git initialized at `/home/yhan/bacdive`
- `.gitignore` added, including `BacDiveComplete/`, `bacdivecomplete/`, and generated artifacts
- Env/config scaffold finalized:
  - `.env.example` expanded for all workflows
  - `R/config_paths.R` supports strict per-script required keys
  - local `.env` refreshed from template (ignored)
- Path migration completed for all active scripts:
  - `calibration_logo.R`
  - `calibration_ipw_platt.R`
  - `calibration_val_test.R`
  - `deployed_calibration_audit.R`
  - `temperature_scaling.R`
  - `build_genome_family_map.R`
  - `eval_test.R`
  - `smoke_test.R`
  - `run_eval_test.sh`
- Shared helper extraction started:
  - `R/calibration_common.R` created and used by calibration scripts
- Validation completed:
  - Parse checks passed for all edited scripts
  - Runtime checks succeeded: `calibration_logo.R`, `calibration_ipw_platt.R`,
    `calibration_val_test.R`, `deployed_calibration_audit.R`, `temperature_scaling.R`
  - Config-only checks passed for build/eval/smoke required keys

## Next Phases

### Phase 2: Minimal Selection Module (Production)

- Keep `phenomnar` as research repo, but vendor only required production steps:
  - `00_load_clean.R`
  - `02_selection_proxies.R`
  - `02b_extended_proxies.R`
  - `03_selection_model.R`
- Production deliverables remain:
  - `data/bugphyzz_with_proxies.rds`
  - `data/selection_probs.rds`
- Add deterministic runner in this repo (e.g., `selection/`) to either:
  - rebuild these artifacts, or
  - validate frozen artifact versions/checksums.

Exit criteria:
- One command validates or rebuilds both required selection artifacts.
- Calibration can run without manual cross-repo path hunting.

### Phase 3: Multiclass Calibration

- Add multiclass calibrators:
  - temperature baseline (already present),
  - vector/matrix scaling,
  - Dirichlet calibration.
- Use LOFO CV for model selection by NLL/ECE.
- Persist chosen calibrator parameters per head.

Exit criteria:
- `gram`, `cellshape`, `flagellum` have explicit calibrated parameter outputs.

### Phase 4: Deployment Integration

- Patch legacy deployment script in `BacDiveComplete` to optionally emit calibrated probabilities.
- Keep legacy behavior behind a feature flag for safe rollback.

Exit criteria:
- Backward-compatible outputs + calibrated outputs available in same run.
