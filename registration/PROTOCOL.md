# Closed-form floors versus learned models in spatial omics: confirmatory protocol

Amirhossein Shahsavand, Shayan Forghani. Hematologic Malignancies Research Center, Research Institute for Oncology,
Hematology and Cell Therapy, Tehran University of Medical Sciences, Tehran, Iran.

Version 1, 2026-11-04. Registration DOI: 10.5281/zenodo.23931264. Code: https://github.com/AH-Shah/spatialfloors, tag
`registration-v1`; release commit and code-archive hash in `SHA256SUMS`.

## 1. Aim and hypotheses

Aim: on units not used for any design choice, test whether pre-specified closed-form floors are non-inferior to published
learned models in five spatial-omics task families, and estimate the difference where they are not.

For a test pair p (family, learned model, primary metric) and unit u, D_pu = M(learned) − M(floor) on the same cells.
μ_p is the population mean of D_pu and δ_f the family margin (§4).

- H1, per pair (primary): the floor is non-inferior, μ_p < δ_f.
- H1, joint: the floors are non-inferior on every pair (intersection-union; each pair at one-sided α = .05, no correction).
- H2, per pair (secondary): the learned model is superior, μ_p > 0 (Holm over the 12 pairs).
- Reported without decision weight: equivalence (90% CI inside (−δ_f, δ_f)) and floor superiority (upper bound < 0).

## 2. Task families, floors and comparators

| Family | Primary metric | Floor (locked) | Learned comparators (locked configuration) | δ_f |
|---|---|---|---|---|
| F1 cell identity | count-split kNN Poisson score | Pearson-residual PCA, 50 PCs | scVI (n_latent 50 human, 100 mouse; NB); Nicheformer (zero-shot, final layer) | .0195 |
| F2 spatial domains | balanced ARI at the true number of regions | GMM on [expression PCA, neighbour mean] | CellCharter (scVI n_latent 10 human / 30 mouse, 3 layers); Novae (zero-shot; mouse units) | .05 |
| F3 spatial context → receiver state | held-out Spearman, quadrant CV | ridge on [type mix, bleed-safe context, same-lineage rings] | LightGBM (15 leaves, lr .05); NCEM (authors' linear interactions model; 19 µm human, 46 µm mouse) | .0043 |
| F4 held-out-type counterfactual | Pearson of predicted vs observed log fold change, official and spike-free | L6 additive type + domain WLS on pseudobulks | Cellina, CPA (authors' configurations) | .03 |
| F5 label transfer | held-out-marker AUROC | InSituType ML + spatial prior | scANVI (joint with the reference; n_latent 30, NB); Tangram (cluster mode, uniform density) | .0255 |

Learned models run with their authors' architectures; only the development grids (Appendix A) and the input adaptations
needed to run inside the shared design were applied. Floors are non-learned, closed-form or fixed-step procedures; they
were revised on development units (Appendix A) and the v2 floors are reported next to v3. Configuration files:
`configs/locked_v3.json`, `configs/locked_v3_ncem.json`; metric definitions with their tests: `SPEC.md` in the code archive.

**F1.** X = X_a + X_b with X_a ~ Binomial(X, .5) (seed 0). Gene selection (seurat_v3, 2,000 genes on whole-transcriptome
panels; all filtered genes otherwise) and every fit use X_a. For each embedding, k = 15 cosine neighbours (self excluded);
μ̂_b,ig = L_b,i · ((1 − λ) · mean_{j∈N(i)} X_a,jg / L_a,j + λ · p̄_a,g), λ = .5; score = 1 − ΣDev(X_b, μ̂) / ΣDev(X_b, L_b p̄_a)
(Poisson deviance over selected genes and cells; L = totals over filtered genes). Floor: analytic Pearson residuals
(θ 100, clipped at √n) of X_a, centred, PCA 50. Nicheformer: mean over non-padding token outputs, per-platform inputs.

**F2.** Scored cells: regions present in the crop at ≥ 2% of the slide's labelled cells; k = number of present regions;
every method clusters every cell at k. Balanced ARI = mean ARI over 10 equal-size resamples (200 cells of each region
with ≥ 200 cells), averaged over 10 clustering seeds (CellCharter: 3 scVI seeds × 10 GMM seeds). Floor: full-covariance
GMM on [log1p-CP10k PCA 50, its Gaussian-weighted mean over the 200 nearest cells of the same FOV (nominal ℓ 400 µm)],
each block centred and scaled to unit total SD, weights √.2 and √.8. Secondary floor: multi-scale context across FOVs
(ℓ 50 / 200 / 800 µm plus composition at 200 µm).

**F3.** Receivers: label-free Leiden identities with ≥ 1,000 cells, fixed once per crop and shared by all models. The crop
is split into four quadrants at the median x and y; each is held out in turn; graphs, kernels and rings stay within a
quadrant; cells within 100 µm of a split line are not scored. Target of receiver r: its bleed-safe genes renormalised on
themselves and centred within each FOV over all of r's cells there (one target for every model; a MERFISH section is one
FOV; on development units 21–40% of scored CosMx cells lie in FOVs that cross a split line); then depth removed,
z-scored and projected on its top-10 state PCs, these three fitted on training quadrants. Score = mean Spearman over the
10 PCs, averaged over folds and receivers. Floor: ridge (α 1000) on the type mix and bleed-safe context of other lineages (human: Gaussian ℓ 25 µm,
10 context components; mouse: context only, ℓ 100 µm, 50 components), plus two same-lineage rings of the receiver's own
genes (10 components fitted outside the test quadrant, FOV-centred): cells 50–150 µm away among the 600 nearest and
150–400 µm away among the 1,500 nearest cells of any lineage, same FOV and quadrant. The far ring reaches a median
220–330 µm on CosMx and the full 400 µm on MERFISH (the design note in `locked_v3.json` gives the CosMx reach,
150–250 µm). LightGBM receives every floor feature group plus sender-lineage-resolved context. The rings carry
tissue-scale gradients of the receiver's own state, so F3 tests prediction of receiver state from spatial context, not
heterotypic signalling alone; the decomposition below separates the parts.

**F4.** Cellina and CPA at the authors' configurations (control REF, held-out CRC; GSE1 222: held-out TVA; MERFISH:
control Thalamus, held-out Isocortex and Fiber tracts), five held-out cell types per section. Official metric: Pearson over the top-50
|LFC| ground-truth genes with the authors' pseudocount; spike-free: pseudocount .1 and its own top 50 (co-primary).
Prospective units: 3 seeds per learned model, averaged per task, then over tasks.

**F5.** Classes: epithelial, immune non-myeloid, myeloid. Held-out markers (15 for CRC; for each new cancer the top 5 per
class by one-vs-rest AUROC in its reference) are removed from every model's input. Per class, AUROC of the log1p-CP10k
sum of its held-out markers between cells predicted that class and all other cells ("other" included); classes with
< 50 cells on either side score .5; mean over classes. Reference: Pelka et al. CRC atlas; for new cancers one
CELLxGENE dataset per tissue chosen by a fixed rule (`configs/refs_*.json`); learned models use a stratified subsample
(≤ 80,000 cells, ≤ 4,000 per type). Floor: InSituType ML log-likelihoods against the reference's mean profiles, scaled by
their median top-1 − top-2 margin, then a mean-field prior over the 15 nearest cells (β 1, 5 sweeps). F5 therefore tests
coarse-lineage transfer; finer types have no held-out validation markers and are not scored.

**Secondary outputs** (no decision weight): v2 floors; F2 multi-scale floor and CellCharter on PCA input; F3 ED-NCEM and a
decomposition of the NCEM pair: NCEM vs a ridge on exactly NCEM's inputs (binarised neighbour types, same graph and
folds), NCEM vs the v2 floor (other lineages only, no rings), v3 vs v2 floor; F5 InSituType without the prior, logistic
regression, IF-gated macro-F1 (GSE1), marker-positive recall and class shares; F2 geometry reference (k-means on xy).

## 3. Data and roles

Unit = patient, donor or animal. Each slide contributes one contiguous ~50,000-cell crop chosen from coordinates and
region labels before any model ran (every present region ≥ 20% of the crop, else the crop maximising the smallest share;
Bruker slides: densest square). Gene panels: protein-coding and annotated genes after removing probes, mitochondrial and
cytoplasmic ribosomal genes, pseudogenes and unannotated lncRNAs. `units/units.tsv` lists every crop with its role,
unit, platform stratum, cell count and the SHA-256 of its cell-ID manifest (`units/crops/`); panels in `units/panels/`.
The estimand is region-scale: in F1–F3 and F5 every model is fitted and scored per crop, or per whole sample where the
sample is smaller (at least 14 of the 29 confirmation crops are whole samples or sections; the GSE1 and SPATCH crops
cover 7–36% of their slides); F4 follows the authors' whole-slide protocol except GSE1 110 (table below). Whole slides are not
part of confirmation in F1–F3 and F5; the validation check is in Appendix A.

- **Development** (every choice): GSE1 232 (CosMx whole-transcriptome CRC, the Cellina CRC data; whole slide for F4);
  GSE0 S0 (CosMx whole-transcriptome colon cancer, public release); GSE303070 slide 148191A samples 118N, 278N, 278T,
  46T (CosMx 1k CRC; donors 118, 278, 46); Allen Brain Cell Atlas MERFISH (Zhuang) animal 4 sections 4.001, 4.002, 4.003.
- **Validation** (locked configurations run once; pre-specified gate): GSE1 231 (same patient as 232) and its whole
  slide; Bruker CosMx whole-transcriptome breast IDC.
- **Confirmation**:

| Family | Units (N) | Strata |
|---|---|---|
| F1 | GSE1 110, 120, 210, 221, 222, 242; GSE303070 donors 7, 64, 73, 74, 92, 130, 137, 142, 150, 170, 280; SPATCH Xenium 5K; Allen animals 1, 2, 3 and Yao 638850; Bruker WTx breast multiomic, skin SCC, kidney RCC (25) | CosMx WTx, CosMx 1k, Xenium 5K, MERFISH |
| F2 | GSE1 (6); Allen (4) (10) | CosMx WTx, MERFISH |
| F3 | GSE1 (6); GSE303070 (11); Allen (4); Bruker (3) (24) | CosMx WTx, CosMx 1k, MERFISH |
| F4 | prospective: GSE1 110 (crop), 222, Allen animals 1, 3, Yao (3 sections each); retrospective: GSE1 120, 210, 221, 242, Allen animal 2 (sections 2.036, 2.039, 2.041) (10) | CosMx WTx, MERFISH |
| F5 | GSE1 (6); GSE303070 (11); Bruker (3) (20) | CosMx WTx, CosMx 1k |

F4 retrospective units carry the authors' published Cellina and CPA results; their official-metric comparison with L6
was computed before registration. No other confirmation metric has been computed. Confirmation data were otherwise
touched only by metric-blind steps: crop selection, panels, format conversion, reference preparation, and a 5,000-cell
smoke window outside each Bruker crop (loading, class coverage, marker detection, IF QC). The kidney reference classes
were taken from the authors' cluster names after the ontology rule placed macrophages in the epithelial class (decided
before any kidney score).

## 4. Analysis (tier 1)

- Crops are averaged into units (GSE303070 donor 280 = samples 280N and 280T; F4 MERFISH animal = mean over its
  sections). Model D_u ~ 1 + mean-centred platform-stratum indicators; the intercept is the mean of D; exact t with
  df = N − number of strata; one-sided 95% bounds.
- Pairs (12): F1 scVI, Nicheformer; F2 CellCharter, Novae (mouse units, N = 4); F3 LightGBM, NCEM; F4 Cellina, CPA ×
  {official (N = 10), spike-free (prospective, N = 5)}; F5 scANVI, Tangram.
- Decisions: non-inferior if upper bound < δ_f; learned superior if lower bound > 0 after Holm over 12 pairs.
- Margins: δ_f = min(cap, 10% of the median above-chance development score of the floor), caps F1 .03, F2 .05, F3 .03,
  F5 .05 (chance: F1 permutation, F2 0, F3 0, F5 .5); F4 keeps .03. Computed once at the first lock; not changed since.
- Sensitivity: one-sided Wilcoxon signed-rank on D (learned superior) and on D − δ_f (floor non-inferior);
  leave-one-cancer-out (F1, F3, F5); F1 Nicheformer without the Yao animal (in its pretraining material); per-stratum
  means. Novae collapse rule: effective rank < 5, or > 50% of cells in one domain; collapsed runs are reported as such.
- Scope: verdicts refer to region scale (§3). On the validation slide, D stayed below 0 in every family on the whole
  slide, with the floor's lead 21–91% smaller than on its crop (Appendix A); floor superiority is not claimed beyond
  region scale.
- Power (`configs/power_v3.json`): sd of D from development and validation units, confirmation N, true mean D equal to
  the development mean, or 0.

| Pair | Development mean D | Power at that D | Power at D = 0 |
|---|---|---|---|
| F1 scVI / Nicheformer | −.013 / −.049 | 1.00 / 1.00 | 1.00 / .99 |
| F2 CellCharter / Novae | −.184 (3 units) / one unit | .92 / – | .16 / – |
| F3 LightGBM / NCEM | −.015 / −.033 | 1.00 / 1.00 | .93 / .18 |
| F4 Cellina / CPA, official | +.011 / −.164 | .99 / 1.00 | 1.00 / .22 |
| F5 scANVI / Tangram | −.027 / −.143 | 1.00 / 1.00 | .66 / .35 |

  Joint power for the intersection-union claim is .90 at the development effects and ≈ 0 if every learned model equals
  its floor; per-pair verdicts are therefore the main result. Spike-free F4 has no development estimate.
- Code: `scripts/p5_analyze.py` at the release commit, run once on the sealed outputs.

## 5. Tier 2: tuning-transfer and capacity diagnostic

Purpose: separate how much a model loses because its configuration was chosen on development units from the difference
between models when each is at its best on the unit, for floors and learned models alike. Tier 2 is computed only after
the tier-1 report is written, never changes a tier-1 verdict, and carries no multiplicity control.

Grids (`configs/tier2.json`; the locked configuration is a member of each):

| Model | Grid |
|---|---|
| F1 floor | Pearson residuals × scaling {centred, z-scored} × PCs {50, 100} |
| scVI | n_latent {10, 30, 50, 100}, NB |
| Nicheformer | layer {final, penultimate} |
| F2 floor | {ℓ 400 µm, λ .8}, {ℓ 400, λ 1}, {ℓ 200, λ .8}, multi-scale secondary floor |
| CellCharter | scVI n_latent {10, 30} × layers {1, 3} (3 scVI seeds) |
| Novae (mouse) | {zero-shot, fine-tuned} |
| F3 floor | ridge α {300, 1000, 3000, 10000} |
| LightGBM | leaves {15, 63} × learning rate {.05, .1} |
| NCEM | {ED-NCEM, linear interactions} × radius {19, 46 µm} |
| F5 floor | InSituType ML without the prior; spatial prior β {.5, 1, 2} |
| scANVI | n_latent {10, 30} × likelihood {NB, ZINB} |
| Tangram | density prior {uniform, RNA count} |

F4 is excluded: the learned models have no grid (authors' configurations) and the floor is locked on both sides.

Units: floors, LightGBM, Nicheformer and Novae (the 4 mouse units; a subset of 2 would leave one degree of freedom) on
every confirmation unit of the family. The other models trained on the unit (scVI, CellCharter, NCEM, scANVI, Tangram)
run on a seeded subset: within each dataset, ⌈n/2⌉
units drawn without replacement by `numpy.random.default_rng(seed)` from the sorted unit list (datasets in sorted
order); seed = first 32 bits of SHA-256(`locked_v3.json`) = 545993150. Drawn: GSE1 210, 221, 222; GSE303070 donors 7, 64,
92, 130, 170, 280; Bruker kidney, breast multiomic; SPATCH Xenium 5K; Allen animal 1, Yao 638850.

Quantities, per model m and unit u: S*_mu = highest primary over the grid (ties: grid order); transfer gap
T_mu = S*_mu − S_mu(locked) ≥ 0. Per pair: capacity gap C_pu = S*_learned − S*_floor = D_pu + T_learned,u − T_floor,u,
on the units where both are available. Unit means with one-sided 95% bounds from the tier-1 model and margin.

Reading (fixed in `scripts/p5_tier2.py`):
- T: upper bound < δ_f → the locked configuration transfers; lower bound > δ_f → development tuning under-represents
  the model; otherwise inconclusive.
- C: upper bound < δ_f → the floor holds at best tuning; lower bound > 0 → the learned model is ahead at best tuning;
  otherwise inconclusive. The decomposition attributes any change between D and C to the learned or the floor side.
- Selection on the scored unit makes S* optimistic; the bias grows with grid size and run-to-run noise. Floors are
  deterministic on CPU and no floor grid is larger than the largest comparator grid in its family.

Order: tier-2 jobs sit behind gate T2, which requires the tier-1 report (`P5.md`, `P5.json`) to exist; `p5_tier2.py`
refuses to run without it and records its SHA-256 in the tier-2 output.

## 6. Execution, blinding and deviations

- The confirmation job list is sealed (`p4_jobs.txt`, SHA-256 in `SHA256SUMS`). Learned-model jobs open at gate G7,
  floors at P4LOCK, tier 2 at T2. Outputs and logs are sealed; a failing job shows only its exit code and traceback.
- After registration and before unblinding, code changes are limited to metric-blind crash fixes (a job runs without a
  change in what it computes), each logged with commit, reason and affected jobs, and checked on a development unit.
- Any other change is a deviation: logged with date, reason and whether any confirmation output had been seen, and
  reported with the results; where possible the registered analysis is reported alongside.
- A model that fails on a unit after the logged fixes (crash, non-finite output) is reported as failed and its pair is
  analysed on the remaining units. No unit is excluded after scoring.
- Exploratory analyses (development stress panels, protein shift-null on IF stains, F1 depth-matched reference, moderators
  such as unit size) run after the tier-1 report and are labelled as such.

## 7. Compute

One Apple M5 laptop (10-core GPU via MPS, 24 GB unified memory), one GPU job at a time, CPU floors bit-identical across
runs. Tier 1: ≈ 80 GPU h (dominated by the F4 Cellina tasks) plus CPU jobs in parallel. Tier 2: ≈ 24 GPU h (scANVI 8.5,
scVI 8, CellCharter 6, Tangram 1.5) and ≈ 25 CPU h (NCEM 12, LightGBM 4, F1 re-scoring 2, the rest < 2 each), about one
day of wall time after the tier-1 report.

## Appendix A. Development history

All choices below used development units unless stated; dates are commit dates of the development repository of
`spatialfloors` (tags `p0-lock`, `p3-lock` there); the public repository starts at the registered state.

| Date | Version | Change and the rule that made it |
|---|---|---|
| 2026-10-28 | metrics (tag `p0-lock`) | Metrics, split, statistics and crop rule fixed; F1 shrinkage λ .5 set before the lock (1/(k+1) gave negative scores on sparse data). |
| 2026-10-28 | review | Margins set to the 10%-of-floor-signal rule (the F3 cap exceeded the floor's whole signal); F2 floor clusterer k-means → full-covariance GMM (ablation located the development gap in the clusterer); F5 primary → held-out-marker AUROC (GSE303070 IF stains failed a model-free QC on every unit); Allen animal 4 added as the mouse development unit. |
| 2026-10-29 | landscape | Each floor characterised over parameters, input perturbations (thinning, subsampling), seeds and development units; choice rule = centre of the broadest plateau that holds on every stress-panel angle. F1 floor → Pearson residuals; F2 primary → balanced ARI. |
| 2026-10-30 | v1 (tag `p3-lock`) | Floors at plateau centres: F1 Pearson residuals, z-scored, 100 PCs; F2 GMM ℓ 200 µm, λ .8; F3 ridge (type mix + context, ℓ 25 µm, 50 components, α by CV); F4 L6; F5 InSituType ML. Learned configurations = highest development mean in each ≤ 4-configuration grid. Margins computed. |
| 2026-10-30 – 11-01 | validation gate | Locked configurations run once on the validation units under pre-specified gates (finite runs, floor within .02 of the unit's landscape best, stress panels, CPU determinism, label QC, power, crop vs whole slide). F2 failed the plateau gate on GSE1 231. Crop check (GSE1 231 whole slide, 298k cells, vs its 50k crop; pass = abs(D_whole − D_crop) ≤ δ_f): F1 scVI failed (.0198), F2 CellCharter failed (.0803), F3 LightGBM (.0042) and F5 scANVI (.0004) passed; F4 compared descriptively (the whole-slide learned values are the authors' runs). |
| 2026-11-02 | v2 | F2 re-centred over development + validation units by the plateau rule (ℓ 400 µm); F1 → 50 PCs, centred (the plateau rule extended to slide scale; 100 PCs lost on the validation whole slide); F3 per species with plateau tolerance = margin (human ℓ 25 µm, 10 components; mouse context only, ℓ 100 µm, 50 components; α 1000). scVI grid extended to n_latent 50 / 100; CellCharter averaged over 3 scVI seeds; GAT removed as a test pair (failed its planted nonlinear control); NCEM added as the published F3 model (authors' models; radius chosen on development; positive control on the authors' data). Crop check after v2: F1 passed (.0067); F2 still failed (.0616), but CellCharter trained and clustered on the whole slide and scored on the crop's cells gave D −.103 against −.097 on the crop, so the difference comes from the scored cells, not from training data; the registered consequence ("crops grow") was amended to keep the crops and report the whole slide descriptively. F3 LightGBM under v2 gave .0081 (not recorded at the time). |
| 2026-11-03 | v3 | Candidates and acceptance rules fixed before any v3 run; development units only; validation run once afterwards. Accepted: F3 same-lineage rings (human .052 → .061, 6/6 units; mouse .125 → .153); F5 spatial prior (AUROC up on 5/5 units, recall within .03). F1, F2, F4 kept v2 (no candidate met every rule). LightGBM re-run on the v3 features. Confirmation jobs built. |
| 2026-11-03 | registration | NCEM test configuration named by the lock rule (linear interactions NCEM, highest development mean); tier-2 design and subset; power at v3. Crop check under v3: F3 LightGBM .0039 < .0043, F5 scANVI .0031 < .0255; on the whole slide D < 0 in every family, with the floor's lead 21–91% smaller than on the crop. A design review before sealing (no confirmation output existed) added: F3 floor jobs ordered after the frozen receiver ids; Wilcoxon for non-inferiority; unit-weighted descriptive means; Novae fine-tuning on every mouse unit and a three-way reading of T; the F3 decomposition; the region-scale estimand; exact wording of the F3 target and rings. Secondary outputs without a confirmation job (GAT, residual boosting, F1 depth-matched reference) were removed from the list or moved to exploratory. This protocol. |

Margins have not changed since v1; learned-model configurations changed only as listed above.
