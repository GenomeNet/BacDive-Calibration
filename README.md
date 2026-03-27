# BacDive Calibration

Post-hoc calibration and deployment packaging for BacDive phenotype predictions from the deepG multi-head model.

Detailed internal documentation was moved to [docs/documentation.md](docs/documentation.md).

## Repo Structure

```text
.
├── R/
│   ├── calibration_logo.R               # LOFO CV: fit Platt/temp per phylum fold, pick best method by event ECE
│   ├── calibration_ipw_platt.R          # IPW-weighted Platt scaling to correct for selection bias
│   ├── calibration_val_test.R           # Val→test and within-test 70/30 calibration comparison
│   ├── deploy_calibrated_csv.R          # Build date-stamped calibrated CSVs + ECE summary for deployment
│   ├── deployed_ece_curves_all_methods.R  # Audit PDF: reliability curves + ECE table (all methods)
│   ├── deployed_calibration_audit.R     # Audit of legacy _prob confidence columns from BacDive website
│   ├── build_genome_family_map.R        # Extract genus from FASTA headers, join to Family via microbe.cards
│   ├── temperature_scaling.R            # Standalone temperature scaling exploration (early prototype)
│   ├── config_paths.R                   # .env loader and path resolution helpers
│   ├── calibration_common.R             # Shared helpers: sigmoid, logit, ECE, Platt/Dirichlet fitting
│   └── selection/                       # Vendored phenomnar selection-model code for IPW weights
├── data/
│   ├── selection/             # Selection probability artifacts from phenomnar
│   └── deployment/            # Date-stamped calibrated CSVs + ECE summaries
├── output/                    # LOFO results (.rds), audit PDFs, comparison tables
├── docs/
│   ├── documentation.md       # Full technical documentation
└── .env.example               # Public-safe environment template
```

## Usage

1. Configure paths:

```bash
cp .env.example .env
# edit .env with local/project paths
```

2. Run calibration workflows:

```bash
Rscript R/calibration_logo.R
Rscript R/calibration_ipw_platt.R
Rscript R/calibration_val_test.R
```

3. Build deployment CSVs (calibrated + ECE summary):

```bash
Rscript R/deploy_calibrated_csv.R
```

4. Build submission audit PDF/CSV (all methods, confidence + event):

```bash
Rscript R/deployed_ece_curves_all_methods.R
```

Deployment outputs are written to `data/deployment/`.

## Methodology

- Base model: deepG multi-head phenotype predictor (binary, 2-class, multiclass, and linear heads).
- Calibration selection: per-head method chosen using event-probability ECE from LOFO cross-validation, with an explicit `none` fallback.
- Reporting: every deployment summary now reports both confidence-ECE and event-probability ECE.
- Deployment parameter source: LOFO median parameters (`median_lofo_*`) by default.
- Methods used:
  - Platt scaling (binary and 2-class heads)
  - Linear-to-Platt for pathogenicity heads
  - Temperature scaling for multiclass heads
- Selection-bias analysis: optional IPW-Platt sensitivity analysis using external selection probabilities.
- Deployment packaging: emits calibrated full/integration CSVs plus per-head ECE summary CSV.

### `_val` vs `_prob` (deployment semantics)

- `_val` stores the raw model output per head (sigmoid/softmax vector or linear score).
- `_val_cal` stores calibrated event outputs in calibrated full CSVs.
- `_prob` is deployment confidence for the predicted label.
- Calibration keeps all predicted labels unchanged.

## Main Results (Confidence + Event)

Source: `data/deployment/bd_calibration_ece_YYYYMMDD.csv` (date-stamped at generation).

| Head | Quality | Method | ECE conf (none→sel) | ECE event (none→sel) |
|---|---|---|---:|---:|
| pathogenicity_human | PASS | platt_linear | 0.5171 → 0.0050 | 0.4539 → 0.0027 |
| pathogenicity_animal | PASS | platt_linear | 0.2180 → 0.0024 | 0.4011 → 0.0001 |
| pathogenicity_plant | PASS | platt_linear | 0.1058 → 0.0010 | 0.4255 → 0.0000 |
| oxygen_facultative | PASS | platt | 0.0516 → 0.0222 | 0.3208 → 0.0222 |
| oxygen_microaerophile | PASS | platt | 0.0424 → 0.0151 | 0.2028 → 0.0150 |
| oxygen_growth | PASS | platt | 0.0198 → 0.0100 | 0.0258 → 0.0126 |
| spore | PASS | platt | 0.0177 → 0.0132 | 0.0521 → 0.0179 |
| motility | PASS | platt | 0.0119 → 0.0179 | 0.0740 → 0.0276 |
| oxygen_obligate | PASS | none | 0.0505 → 0.0505 | 0.0767 → 0.0767 |
| gram | PASS | temp | 0.0727 → 0.0704 | 0.0435 → 0.0462 |
| cellshape | FAIL | temp | 0.1594 → 0.1126 | 0.0277 → 0.0227 |
| flagellum | FAIL | temp | 0.1244 → 0.1166 | 0.0324 → 0.0294 |

`biosafety` and `cellsize` are listed as `EXCLUDED` in the CSV summary.

