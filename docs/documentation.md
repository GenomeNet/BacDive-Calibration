# BacDive Phenotype Model Calibration

Post-hoc probability calibration for the BacDive multi-head bacterial phenotype
prediction model, with selection-bias correction via inverse probability weighting
from the phenomnar MNAR analysis pipeline.

## Model

- **Checkpoint**: `bacdive_attn_lstm_sum_28/Ep.023.hdf5`
- **Architecture**: CNN-LSTM with time-distributed set learning (deepG), 14 output heads
- **Input**: 200 fragments of 10kb each per genome (2Mb total), one-hot encoded DNA
- **Predictions**: 43,086 genomes; val + test pool = 4,686 genomes (after Phylum filter)

## Phenotype Heads

| Head | Activation | Calibration | IPW-Platt |
|------|-----------|-------------|-----------|
| motility | sigmoid | Platt | sel_prob_motility |
| spore | sigmoid | Platt | sel_prob_spore_formation |
| oxygen_microaerophile | sigmoid | Platt | sel_prob_aerophilicity |
| oxygen_growth | softmax(2) | Platt | sel_prob_aerophilicity |
| oxygen_facultative | softmax(2) | Platt | sel_prob_aerophilicity |
| oxygen_obligate | softmax(2) | Platt | sel_prob_aerophilicity |
| pathogenicity_human | linear | Platt | sel_prob_animal_pathogenicity |
| pathogenicity_animal | linear | Platt | sel_prob_animal_pathogenicity |
| pathogenicity_plant | linear | Platt | sel_prob_plant_pathogenicity |
| gram | softmax(3) | Temperature | -- (multiclass, skip) |
| cellshape | softmax(11) | Temperature | -- (multiclass, skip) |
| flagellum | softmax(6) | Temperature | -- (multiclass, skip) |
| biosafety | linear | Excluded | -- |
| cellsize | linear(4) | Excluded | -- |

### Per-head accuracy (val vs test)

Val shares 93.7% of genera with train (in-distribution); test is a genus-level
holdout (12.7% overlap). Balanced accuracy = mean per-class recall.

| Head | Type | N val | N test | Acc val | Acc test | Bal acc val | Bal acc test | Delta |
|------|------|------:|-------:|--------:|---------:|------------:|-------------:|------:|
| spore | binary | 573 | 690 | 0.900 | 0.878 | 0.905 | 0.890 | -0.015 |
| oxygen_obligate | softmax2 | 42 | 80 | 0.881 | 0.875 | 0.831 | 0.840 | +0.009 |
| oxygen_growth | softmax2 | 1055 | 1326 | 0.906 | 0.807 | 0.883 | 0.779 | -0.104 |
| oxygen_microaerophile | binary | 1313 | 1635 | 0.828 | 0.805 | 0.808 | 0.767 | -0.041 |
| motility | binary | 677 | 953 | 0.749 | 0.682 | 0.758 | 0.685 | -0.073 |
| gram | multiclass(3) | 834 | 1112 | 0.866 | 0.772 | 0.854 | 0.599 | -0.254 |
| oxygen_facultative | softmax2 | 111 | 233 | 0.676 | 0.631 | 0.672 | 0.633 | -0.038 |
| flagellum | multiclass(6) | 78 | 145 | 0.423 | 0.531 | 0.331 | 0.458 | +0.127 |
| cellshape | multiclass(11) | 709 | 996 | 0.165 | 0.167 | 0.262 | 0.151 | -0.111 |
| pathogenicity_human | linear | 2215 | 2575 | 0.276 | 0.312 | 0.376 | 0.413 | +0.037 |
| pathogenicity_animal | linear | 2215 | 2575 | 0.573 | 0.689 | 0.460 | 0.606 | +0.146 |
| pathogenicity_plant | linear | 2215 | 2575 | 0.932 | 0.935 | 0.509 | 0.474 | -0.035 |

**Notes:**

- Gram drops 25 points in balanced accuracy on test because the minority class
  (gram-variable, ~3% of labels) collapses to near-zero recall on out-of-genus
  genomes — the model learned genus-specific shortcuts for a rare class.
- Cellshape (11 classes) is near-random on both splits (bal acc 0.15-0.26);
  the head has no useful discriminative power beyond the majority class.
- Pathogenicity heads have low accuracy + low balanced accuracy because the
  linear head's raw scores cluster near the decision threshold (see calibration
  discussion below). Accuracy appears inflated for plant (93%) only because
  prevalence is 1.3% — predicting "no" for everyone scores 98.7%.
- Flagellum improves on test (+13 points) possibly because test genera have
  simpler flagellum distributions; however N is very small (78/145).

---

## Calibration Methods

### 1. Temperature Scaling

Single scalar T per head (Guo et al., 2017).

- Binary/2-class: `p_cal = sigmoid(logit(p) / T)`
- Multiclass: `p_cal = softmax(log(p) / T)`

Largely ineffective (T ~ 1.0) because miscalibration was predominantly a
**bias problem**, not a sharpness problem.

### 2. Platt Scaling

Two parameters (a, b) per head (Platt, 1999).

- Sigmoid/softmax heads: `p_cal = sigmoid(a * logit(p) + b)`
- Linear heads (pathogenicity): `p_cal = sigmoid(a * raw_output + b)`

Resolves bias. Fitted via NLL minimization (BFGS).

### 3. Leave-One-Phylum-Out (LOFO) Cross-Validation

For each Phylum: hold out all genomes in that group, fit calibration on the rest,
predict on the held-out group. Aggregate held-out predictions for ECE. This avoids
phylogenetic leakage where closely related organisms inflate calibration quality.

Taxonomy mapping: genus extracted from FASTA headers (three formats: PATRIC brackets,
GCA/GCF, IMG-style), joined to microbe.cards S1.xlsx for Family/Order/Phylum.
Saved as `genome_family_map.csv` (12,823 genomes, 97.5% coverage).

### 4. IPW-Platt Scaling (Selection-Bias Correction)

Standard Platt scaling optimizes on labeled genomes, which are biased toward
well-studied species. IPW-Platt reweights the NLL to correct for this:

```
nll_ipw = sum(w_i * nll_i) / sum(w_i)
```

**Weights**: Stabilized Hajek weights from phenomnar selection probabilities:

```
w_i = P(R=1) / pi_hat(genus_i)
```

where `pi_hat` is the genus-level median selection probability from the phenomnar
Phase 3 selection model (see below). Weights are percentile-trimmed and normalized.

**Linkage**: BacDive genomes are joined to phenomnar via genus name. 97.7% of
genomes match; unmapped genomes receive w = 1 (no correction).

---

## Phenomnar Selection Model (Phases 2-3)

The IPW weights come from the phenomnar MNAR analysis pipeline, which models the
probability that a species has a given phenotype annotated as a function of
research attention, historical priority, accessibility, and genomic characterization.

### Phase 2b: Proxy Variable Collection

| Proxy | Source | Imputation | Axis |
|-------|--------|------------|------|
| `log_pub_count` | PubMed `"Binomial name"[Organism]` | -- (from Phase 2a) | Research attention |
| `log_assembly_count` | NCBI Assembly search | -- (from Phase 2a) | Research attention |
| `log_wiki_size` | Wikipedia MediaWiki API | -- (from Phase 2a) | Research attention |
| `log_species_age` | NCBI Taxonomy authority field (year parsed) | 209 species without year dropped (not imputed) | Historical priority |
| `n_culture_collections` | NCBI Taxonomy type material (18 collection prefixes) | 0 if no type material | Accessibility |
| `has_complete_genome` | RefSeq assembly_summary.txt | 0 if no assembly | Genomic quality |
| `has_assembly` | RefSeq assembly_summary.txt | binary indicator | Genomic quality |
| `has_checkm_eval` | NCBI CheckM report | binary indicator | Genomic quality |
| `checkm_imp` | NCBI CheckM report | observed / median (98.7) if assembly but no CheckM / 0 if no assembly | Genomic quality |
| `log_contig_imp_median` | RefSeq assembly_summary.txt (min contig count) | median (14) if no assembly | Genomic quality |

**CheckM imputation rationale**: 1,092 species have assemblies but no CheckM
evaluation (329 with complete genomes). These cannot be 0-imputed since they
have real genomes that simply weren't evaluated. They receive the median observed
completeness (98.7). Species with no assembly get 0.

**Coverage**: 19,036 total species in bugphyzz. 209 lack a parseable valid
description year and are dropped (not imputable from NCBI taxonomy), leaving
18,827 species (98.9%). Within this filtered set, all other proxy columns
achieve 100% coverage via the imputation strategies above.

### Phase 3: Selection Model Fitting

Per-phenotype logistic regression predicting annotation status:

```
R_j ~ Phylum_collapsed + proxy_1 + proxy_2 + ... + proxy_k
```

where `R_j = 1` if species has phenotype j annotated (non-NA). Phyla with <50
species are collapsed to "Other" to prevent overfitting.

**Model variants** (for comparison/sensitivity):

| Model | Proxies | Coverage |
|-------|---------|----------|
| M1 | Taxonomy only | 18,827 |
| M2_old | log_pub + log_asm + log_wiki (complete-case) | ~15,700 |
| M2_cc | All extended proxies (complete-case) | ~15,700 |
| M2_med | All imputed (median contig) | 18,827 |
| M2_p95 | All imputed (P95 contig) | 18,827 |

**M2_med** is used for final selection probabilities (best coverage, identical
AUC to M2_p95 since `has_assembly` binary dominates contig count value).

**Diagnostics**: Correlation matrix, PCA, VIF, 5-fold CV-AUC, McFadden R^2,
propensity density plots, KS test + density overlap (annotated vs unannotated).

**Selection probability output**: `selection_probs.rds` (18,827 species x 12
phenotypes). Joined to BacDive genomes at the genus level (median aggregation).

---

## IPW Diagnostics

### Weight sanity

| Head | n_eff/n (P99) | w range | Issue |
|------|---------------|---------|-------|
| motility | 0.501 | [0.42, 9.05] | -- |
| oxygen (all 4) | 0.692 | [0.25, 3.38] | -- |
| pathogenicity_human/animal | 0.705 | [0.10, 3.02] | -- |
| pathogenicity_plant | 0.897 | [0.76, 4.03] | -- |
| **spore** | **0.010** | [0.00, 99.7] | Extreme weights; use P90 trim |

### Trimming sensitivity

Aggressive trimming (P90) rescues extreme cases while preserving the correction
for heads that benefit:

| Head | ECE plain | P90 | P95 | P99 |
|------|-----------|-----|-----|-----|
| motility | 0.0376 | 0.0354 | 0.0345 | 0.0308 |
| pathogenicity_human | 0.0334 | 0.0320 | 0.0311 | 0.0308 |
| pathogenicity_plant | 0.0011 | 0.0002 | 0.0002 | 0.0011 |
| spore | 0.0231 | 0.0271 | 0.0295 | 0.0266 |

### Parameter shift (selection bias effect on calibration)

The b parameter (intercept) shift quantifies selection bias:

| Head | b_plain | b_ipw | delta_b | Interpretation |
|------|---------|-------|---------|----------------|
| pathogenicity_animal | -3.341 | -4.028 | -0.687 | Strong bias correction |
| pathogenicity_human | 0.632 | 0.319 | -0.313 | Moderate bias correction |
| oxygen_growth | -0.199 | -0.421 | -0.223 | Moderate bias correction |
| motility | -0.394 | -0.354 | +0.040 | Minimal shift |

### Stratified ECE (low-propensity tertile = understudied species)

IPW helps most where selection bias is strongest:

| Head | ECE plain | ECE IPW | Improvement |
|------|-----------|---------|-------------|
| pathogenicity_human | 0.0378 | 0.0138 | 63% |
| pathogenicity_animal | 0.0365 | 0.0158 | 57% |
| oxygen_growth | 0.0415 | 0.0345 | 17% |

---

## Files

```
bacdive/
  # Executable code
  R/
    calibration_logo.R
    calibration_ipw_platt.R
    calibration_val_test.R
    deploy_calibrated_csv.R      # builds calibrated deployment CSVs + ECE summary
    deployed_ece_curves_all_methods.R  # audit: all-methods reliability curves (incl. pathogenicity)
    deployed_calibration_audit.R # legacy single-method audit
    temperature_scaling.R
    build_genome_family_map.R
    config_paths.R
    calibration_common.R
    selection/                   # vendored minimal selection-model scripts

  # Outputs — IPW calibration
  output/
    calibration_logo.pdf/rds
    calibration_val_test.pdf/rds
    calibration_results.pdf/rds
    deployed_calibration_audit.pdf
    deployed_ece_curves_all_methods.pdf  # all-methods reliability curves
    ipw_platt_comparison.csv     # Per-head: (a,b) vs (a_ipw,b_ipw), ECE, Brier
    ipw_weight_diagnostics.csv   # Weight summary stats, n_eff per head
    ipw_smd_balance.csv          # SMD before/after weighting per covariate
    ipw_tipping_point.csv        # Gamma sweep: a(g), b(g), ECE(g) per head
    ipw_trimming_sensitivity.csv # P90/P95/P99 trim: ECE, n_eff per head
    ipw_stratified_ece.csv       # ECE by propensity tertile
    ipw_platt_results.rds        # Full results object
    fig_ipw_weights.pdf          # Weight distributions per head
    fig_ipw_calibration.pdf      # Reliability diagrams: plain vs IPW
    fig_ipw_trimming.pdf         # ECE + n_eff by trim level
    fig_ipw_tipping.pdf          # Parameter trajectories vs gamma

  # Local selection artifacts (no cross-repo dependency)
  data/selection/
    selection_probs.rds          # 18,827 species x 12 phenotypes
    bugphyzz_with_proxies.rds

  # Legacy
  BacDiveComplete/               # Legacy repo (training, eval, data generation)
```

## Usage

```bash
# 1) Configure paths
cp .env.example .env
# edit .env for your environment

# 2) Run scripts from repo root
RSCRIPT=Rscript

# One-time: build taxonomy mapping (requires conda genome env for FASTA access)
conda run -n genome Rscript R/build_genome_family_map.R

# Plain LOFO calibration
$RSCRIPT R/calibration_logo.R

# IPW-Platt (uses local data/selection artifacts)
$RSCRIPT R/calibration_ipw_platt.R

# Build deployment CSVs (date-stamped: bd_calibrated_*_YYYYMMDD.csv)
$RSCRIPT R/deploy_calibrated_csv.R

# All-methods audit (reliability curves for all heads including pathogenicity)
$RSCRIPT R/deployed_ece_curves_all_methods.R

# Additional analyses
$RSCRIPT R/calibration_val_test.R
$RSCRIPT R/deployed_calibration_audit.R
$RSCRIPT R/temperature_scaling.R
```

### Environment keys

`.env.example` is the source of truth. Main groups:

- Core data: `BD_PRED_CSV`, `BD_LABELS_CSV`, `BD_SPLITS_CSV`, `GENOME_FAMILY_MAP_CSV`
- Selection/IPW dependency: `SELECTION_PROBS_RDS`, `BUGPHYZZ_PROXIES_RDS`
- Inference inputs: `COUNT_LIST_RDS`, `BACDIVE_COMPLETE_RDS`, `FASTA_DIR`, `CHECKPOINTS_DIR`, `MICROBE_CARDS_XLSX`
- Outputs: `REPORTS_DIR`, `OUTPUT_DIR`, `FIGURES_DIR`, `SLURM_OUTPUT_DIR`, `BACDIVE_TEST_RDS`
- SLURM helper runtime: `BACDIVE_CONDA_SH`, `BACDIVE_TF_ENV`, `EVAL_SCRIPT_PATH`
- Multiclass option: `ENABLE_DIRICHLET_MULTICLASS` (default `false`)

## External Data

### Data flow

```
Raw FASTA (51,340 genomes, e.g. 1000562.3.fasta)
  /vol/projects/pmuench/bacdive_new/bacdive_references/
  Header: >contig_id  Species description  [Species | PATRIC_ID]
       │
       ├──[training only]──> one-hot RDS (pre-encoded for training loop)
       │    /vol/projects/rmreches/bacdive_time_dist_2m_rds/{train,validation}/
       │
       └──[val+test+unlabeled]──> legacy FASTA inference (reference: eval_non_train_self_pad.R)
            ├── per-window: eval_bd_new.rds (43,086 genomes × 14 heads)
            ├── aggregated raw: bd_pred_new_full.csv
            └── deployment calibration: R/deploy_calibrated_csv.R -> data/deployment/bd_calibrated_*_YYYYMMDD.csv
```

FASTA filenames match the `file` column in labels/splits CSVs.
`build_genome_family_map.R` also reads FASTA headers to extract genus names.

There is **no test/ RDS directory**. Train genomes were pre-encoded to RDS for the
training loop; val/test/unlabeled genomes are encoded on-the-fly during inference
(`deepG::seq_encoding_label()` on raw FASTA). The train/val/test split is defined
solely by the `type` column in the splits CSV.

### Labels & splits

- **Ground truth labels** (12,823 genomes × 14 heads, phenotype values):
  `/vol/projects/rmreches/bacdive_labels/training_morphcomplete_norm_na_rm.csv`
- **Train/val/test splits** (12,823 genomes, `file` + `type` column only):
  `/vol/projects/rmreches/bacdive_labels/training_morphcomplete_labels_2024-01-31.csv`
  - train: 8,030 | validation: 2,216 | test: 2,577
  - Calibration pool = val + test = 4,793 (4,686 after Phylum filter)
- **Other label files** (all in `/vol/projects/rmreches/bacdive_labels/`):
  - `training_morphcomplete_norm.csv` — normalized labels (before NA removal)
  - `training_morphcomplete_2024-01-31.csv` — raw labels (before normalization)
  - Per-task label files: `training_cell_morph_*.csv`, `training_cell_oxy_*.csv`,
    `training_cell_size_*.csv`, `training_cell_spore_*.csv`, `training_patho_*.csv`
  - `morphcomplete_mean_sd.rds` — normalization parameters (for cellsize head)

### Data split provenance

**No generating script was preserved.** The split was created by pmuench (Jan 31
2024), likely interactively, following the same pattern as `process_abiogram_data.rmd`
in the GenomeNet/BacDive GitHub repo. The original file is at
`/vol/projects/pmuench/bacdive_new/labels/` (columns `filename`, `label`); the copy
at `/vol/projects/rmreches/bacdive_labels/` renames columns to `file`, `type`.

**Split procedure** (inferred from abiogram template + observed overlap statistics):

1. **Test set (~20%)**: genus-level holdout — randomly shuffle genera, greedily
   accumulate entire genera until ~20% of genomes reached. This makes test a
   hard out-of-genus generalization set.
2. **Train/val split of remainder**: within-genus stratified split via
   `caret::createDataPartition(genus, p=0.7)` — 70% of each genus's genomes
   go to train, 30% to val. Val and train share almost all genera.

**Taxonomic overlap between splits:**

|  | Val | Test |
|--|-----|------|
| Genera | 541 | 679 |
| Genus overlap with train | 93.7% (507/541) | 12.7% (86/679) |
| Genome-level genus overlap with train | 98.4% | 28.9% |
| Exclusive genera (not in either other set) | 30 | 589 |
| Family overlap with train | 99.5% | 73.5% |

**Implication**: val is an in-distribution split (same genera as train), while test
is a genuine out-of-genus holdout. Val performance overestimates generalization to
unseen genera. This asymmetry motivated LOFO calibration (hold out entire Phyla)
rather than fitting calibration on val alone.

### Training RDS (pre-encoded, training loop only)

- `/vol/projects/rmreches/bacdive_time_dist_2m_rds/train/`
- `/vol/projects/rmreches/bacdive_time_dist_2m_rds/validation/`
- Created by `rds_time_dist.R` from FASTA; not used for calibration

### Model

- **Checkpoint**: `/vol/projects/BIFO/genomenet/states/bacdive_attn_lstm_sum_28/Ep.023.hdf5`
- **Checkpoint (alt path)**: `/vol/projects/BIFO/genomenet/checkpoints/bacdive_attn_lstm_sum_28/`

### Predictions

Raw prediction sources:

- **Per-window predictions** (200 fragments per genome, raw per-head matrices):
  `/vol/projects/BIFO/genomenet/states/eval_bd_new.rds` (38 MB, 43,086 genomes)
- **Full aggregated raw CSV**:
  `/vol/projects/BIFO/genomenet/states/bd_pred_new_full.csv` (50 MB)
- **Integration raw CSV**:
  `/vol/projects/BIFO/genomenet/states/bd_pred_integrate_202502.csv` (12 MB)

Current deployment-calibrated outputs are generated locally via
`R/deploy_calibrated_csv.R` into `data/deployment/`.

### Taxonomy

- **Microbe.cards**: `/home/yhan/GenomeInterpretation/sporulation/microbe.cards table S1.xlsx`

### Selection model artifacts (local)

- **Selection probabilities**: `data/selection/selection_probs.rds` (18,827 species × 12 phenotypes)
- **Bugphyzz with proxies**: `data/selection/bugphyzz_with_proxies.rds`

### `_prob` column semantics

- In raw source CSVs, `_prob` is legacy confidence:
  - Binary: `max(p, 1-p)`
  - 2-class softmax: `max(p1, p2)`
  - 3+ class softmax: `0.5 + 0.5 * log(K * max(p)) / log(K)`
  - Pathogenicity: `(tanh(5 * |raw - threshold|) + 1) / 2`
- In calibrated deployment CSVs (`data/deployment/bd_calibrated_*_YYYYMMDD.csv`), `_prob` is
  confidence-calibrated while predicted label columns are kept unchanged.
- `_val` remains raw model output; `_val_cal` stores calibrated event outputs.

---

## Note on ECE Formulations

Two ECE formulations are used and should not be mixed:

- **Probability ECE** (`calibration_logo.R`, `calibration_val_test.R`, `calibration_ipw_platt.R`):
  calibrates event probabilities (`P(y=1)` or class probabilities).
- **Confidence ECE** (`deployed_calibration_audit.R`, `deployed_ece_curves_all_methods.R`,
  `deploy_calibrated_csv.R`): calibrates deployed confidence for fixed predicted labels.

Deployment packaging selects calibration method per head using **event-probability ECE
from LOFO cross-validation** (with a `none` fallback). The quality gate (PASS/FAIL)
no longer blocks method selection — all heads are eligible for calibration if ECE
improves. Both confidence-ECE and event-ECE are reported for transparency.

Submitted deployment artifacts now include both views simultaneously:

- `data/deployment/bd_calibration_ece_YYYYMMDD.csv` with `ECE_conf_*` and `ECE_event_*` columns.
- `output/deployed_ece_all_methods.csv` with per-method `ECE_conf` and `ECE_event`.
- `output/deployed_ece_curves_all_methods.pdf` with confidence and event reliability curves
  side-by-side (all 12 heads including pathogenicity).
