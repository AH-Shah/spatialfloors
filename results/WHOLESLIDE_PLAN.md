# Whole-slide arm: analysis plan (exploratory)

Written 2026-11-09, after `results-v1` and before any whole-slide output of a confirmation slide existed. It is
committed before the runs start; the commit time is the record. This is an exploratory analysis in the sense of
PROTOCOL §6. It carries no decision weight and cannot change a registered verdict.

## 1. Question

The registered estimand is region scale: ~50,000-cell crops, or whole samples where a sample is smaller
(PROTOCOL §3–4). On the one whole slide checked before registration (validation GSE1 231, 298k cells, same patient as
development unit 232), D stayed below 0 in every family, but the floor's lead was 21–91% smaller than on the crop
(PROTOCOL Appendix A). This arm asks, on the confirmation slides themselves: do the four pairs that train on the unit
stay non-inferior on whole slides, and how much does D move between crop and whole slide?

## 2. Units

The six GSE1 confirmation slides (110, 120, 210, 221, 222, 242), every cell of each (141k–741k cells; 3.04M in all).
All six are run; none is selected or dropped after scoring. Their registered crops cover 7–36% of each slide. GSE1 is
the only confirmation dataset where every family applies (F2 needs region labels) and the crop is a small part of the
slide; the GSE303070 samples and the Allen sections are already whole. Unit = slide (= patient). One stratum
(CosMx WTx), df = 5.

## 3. Pairs and configurations

Locked configurations (`manifests/locked_v3.json`), registered metrics and scripts:

| Pair | Floor | Learned model |
|---|---|---|
| F1 cell identity | Pearson-residual PCA, centred, 50 PCs | scVI, n_latent 50, NB |
| F2 spatial domains (balanced ARI) | GMM on [PCA 50, 200-NN mean], ℓ 400 µm, λ .8 | CellCharter, n_latent 10, 3 layers, scVI seeds 0–2 × GMM seeds 0–9 |
| F3 context → receiver state | v3 floor (ridge, α 1000, with rings); v2 floor (no rings) as secondary | LightGBM, 15 leaves, learning rate .05, v3 features |
| F5 label transfer (held-out-marker AUROC) | InSituType ML + spatial prior | scANVI, n_latent 30, NB |

Not run:
- Nicheformer: about 50 h of MPS inference for 3M cells. Its per-cell embedding does not depend on the other cells;
  only the k-NN population would change.
- NCEM, Tangram and Novae: human Novae is not a test pair.
- F4: already follows the authors' whole-slide protocol, except GSE1 110.

Receivers (F3), the panel and the count split are recomputed on each whole slide by the registered rules.

## 4. Pipeline and memory

The registered scripts run with `SF_LEAN=1`. Without it nothing changes.
- **Cache** (`scripts/xw_prepare.py`). The registered loader needs ~23 GB on slide 120 (741k cells, 968M non-zeros),
  so the raw file is streamed instead. The script reproduces the registered whole-slide definitions exactly:
  - every cell, in on-disk order;
  - X_a = Binomial(X, .5) from `default_rng(0)` over the non-zeros in CSR order;
  - seurat_v3 HVG 2,000 on X_a;
  - obs as `load_crop`.

  The cache keeps the panel, the F5 held-out markers, the negative probes, and one column with the counts of every
  other filtered gene, so library sizes are unchanged. The script was checked on GSE1 231 against the registered
  `gse1_231w` cache and panel: panel, X, X_a, the rest column and obs are all identical.
- **PCA** (F1 and F2 floors). The registered code builds a dense cells × 2,000 float64 matrix (11.9 GB on 741k cells)
  and runs randomized PCA. The lean path computes the same rows per chunk and runs exact PCA from the accumulated
  covariance.
- **F3 rings.** The registered code queries a dense table of up to 1,500 neighbours per cell (17.8 GB transient on
  741k cells). The lean path queries in chunks and stores a sparse matrix with the same selection.
- **InSituType.** The same matrix goes to R as binary CSC arrays instead of Matrix Market text. InSituType itself is
  not chunked, because its background fit and cluster frequencies use all cells.
- **Equivalence check** (`scripts/xw_check.py`), with pass rules fixed before the run. The lean paths run on GSE1 231
  as a whole slide and are compared with that slide's registered-path outputs (`runs/p3v2`, `runs/g7v3`):
  - F1: |Δ score| < .001
  - F2: |Δ balanced ARI| < .02
  - F3 receiver ids: identical
  - F3 ridge: |Δ| < .001
  - F3 LightGBM: |Δ| < .002
  - F5: InSituType log-likelihoods identical, and |Δ AUROC| < 1e-9

  Results are in §7.
- **Execution** (amended 2026-11-09, before any learned-model output; see the feasibility review below).
  - `runs/sched/xw_jobs.txt` (from `scripts/xw_jobs.py`) runs through `scripts/sched.sh`, with one GPU job at a time
    and at most two CPU jobs.
  - **System limits** (the machine's own): RAM 23 GB, which is every job's `mem_guard` footprint cap; and swap at
    most 2 GB for at most 15 minutes.
  - The largest slide (120, 741k cells) runs first, all of its jobs. Every job of the other five slides waits until
    every slide-120 job has finished (gate `xw_gate120`).
  - **Stop rule** (`scripts/xw_watch.sh` writes `runs/xw/ARM_FAILED`). The arm is declared failed and stopped, and
    reported as infeasible at whole-slide scale on this machine, with nothing rerun, if either:
    - a slide-120 job fails (cap kill or crash); or
    - swap stays above 2 GB for more than 15 minutes at any point.
  - On the smaller slides, a model that fails is reported as failed on that slide.
- **Compute estimate.** From the feasibility review below: 30–60 min per learned model on slide 120, about 10–12 h of
  wall time for the arm on the M5.

### Feasibility review and amendment (2026-11-09)

The first start of the pool was stopped by the investigators after 45 min, before any learned model had produced
an output, to check whether each learned model can run on a whole slide as-is.

- **What had finished by then:** the caches of slides 120, 221, 210 and 242; the F3 receiver ids of 120, 221 and 242;
  and the F2 floor on 120 and 221. These outputs were kept unread.
- **Review method:** `scripts/xw_scaling.py` read the peak memory footprint and active run time of 1,369 earlier
  jobs (crops of 13.6k–54k cells and the 298k-cell GSE1 231 slide) and related them to the number of cells
  (`runs/xw/XW_SCALING.md`).
- **The four learned models of this arm:** scVI, CellCharter, LightGBM and scANVI run unchanged (no change to their
  code since the results). Their memory grows ~0.5–1.2 GB per 100k cells, giving a forecast of ≤ ~10 GB and 30–60 min
  each at 741k cells.
- **Three learned models that cannot run on a whole slide on this machine:**
  - Nicheformer: per-cell inference at ~16 cells/s, 7–23 h per slide;
  - NCEM: 2–4 h per 50k crop, with images padded to the largest;
  - Tangram: a dense cells × genes matrix on the GPU, ~12 GB per 100k cells.

  They stay outside the arm (§3), and this limit is reported.
- **Slide 120 is the worst case.** No other slide is larger, hence the stop rule.

## 5. Analysis (`scripts/xw_analyze.py`, run once after every job is done)

Per pair, over the six slides:
- **D_whole** = learned − floor on the whole slide, with one-sided 95% bounds from the exact t and the registered
  family margin δ_f.
- **Non-inferior at slide scale:** upper bound < δ_f.
- **Learned ahead at slide scale:** lower bound > 0 after Holm over the four pairs.
- **Scale effect:** Δ = D_whole − D_crop per slide, where D_crop is the registered tier-1 value of the same unit
  (`results/reports/P5_corrected.json`). Reported as the mean with a two-sided 90% CI, plus the change in the
  floor's lead, 1 − mean D_whole / mean D_crop.
- **Secondary (descriptive):**
  - F3 LightGBM − v2 floor (no rings);
  - F5 scANVI − floor at scANVI's per-class call rate (as `scripts/x_f5_matched.py`).

## 6. Reporting commitments

- All four pairs are reported with the numbers above, whatever they show, in `RESULTS.md` and in the paper.
- If a pair is not non-inferior at slide scale, or a learned model is ahead after Holm, the abstract says so next to
  the region-scale verdicts.
- The registered verdicts and their wording are unchanged. Superiority of a floor is not claimed at slide scale
  unless its upper bound is < 0 there.

## 7. Equivalence check results

Run on GSE1 231 as a whole slide (298k cells) before any confirmation-slide job, and compared with that slide's
registered-path outputs. Source: `runs/xw/check/XW_CHECK.json` in the development repository.

| Check | Registered path | Lean path | abs(Δ) | Rule | Result |
|---|---|---|---|---|---|
| Cache: panel, X, X_a, rest column, obs vs the `gse1_231w` cache | – | – | identical | identical | pass |
| F1 floor score | .081621 | .081621 | 4.5e-8 | < .001 | pass |
| F2 floor balanced ARI | .190874 | .191205 | 3.3e-4 | < .02 | pass |
| F3 receiver ids | – | – | identical | identical | pass |
| F3 v3 floor (ridge) | .041036 | .041036 | 0 | < .001 | pass |
| F3 LightGBM | .026449 | .026449 | 0 | < .002 | pass |
| F5 InSituType log-likelihoods | – | – | 0 | < 1e-9 | pass |
| F5 floor AUROC | .633973 | .633973 | 0 | < 1e-9 | pass |

- **First attempt.** The lean PCA was first an exact PCA from the covariance. It failed the F1 rule (abs(Δ) .0017):
  the registered floor runs sklearn's randomized PCA, whose trailing components differ from the exact ones. The lean
  PCA now runs sklearn's randomized algorithm itself, over row chunks, and the check above is the second run. The
  first run's outputs are kept in `runs/xw/check/exact_pca`.
- **Remaining differences.** F2 differs at 3e-4 because the registered F2 PCA runs in float32, so summation order
  matters there.
- **Unit tests** (`tests/test_lean.py`):
  - lean randomized PCA = registered to 1e-6 (float64);
  - signal components to 1e-3 (float32);
  - chunked k-NN kernels bitwise identical;
  - sparse rings equal to the dense rings.
- **Full suite:** 36 passed, 1 skipped (the registered 33 passed + 1 skipped, plus the 3 new tests). Registered paths are unchanged
  without `SF_LEAN`.
- **Peak memory at 298k cells, lean path** (registered path in brackets):

  | Job | Lean | Registered |
  |---|---|---|
  | Cache build | ~2.5 GB | ~23 GB at 741k |
  | F1 floor | 2.6 GB | 9.8 GB |
  | F2 floor | 1.9 GB | 5.0 GB |
  | F3 ridge | 4.6 GB | – |
  | F3 LightGBM | 4.9 GB | – |
  | F5 floor | 2.6 GB | 4.6 GB |

- **Memory plan for slide 120** (741k cells): job estimates are 4.5–9.5 GB, with kill caps 10–16 GB.
