# Results: closed-form floors versus learned models in spatial omics

Confirmatory results of the protocol registered as DOI [10.5281/zenodo.23931264](https://doi.org/10.5281/zenodo.23931264)
(`registration/PROTOCOL.md`, tag `registration-v1`). Release `results-v1.2`, 2026-11-11, deposited on Zenodo as version 4
of the registration record, DOI [10.5281/zenodo.24135678](https://doi.org/10.5281/zenodo.24135678). Earlier releases:
`results-v1.1` (2026-11-10, version 3, DOI [10.5281/zenodo.24109807](https://doi.org/10.5281/zenodo.24109807)) and
`results-v1` (2026-11-08, version 2, DOI [10.5281/zenodo.24044400](https://doi.org/10.5281/zenodo.24044400)); the
registration (version 1) is unchanged. The registered verdicts and every tier-1 and tier-2 number are the same in
all releases (`CHANGELOG.md`). Every number below is read from
the reports in [`results/reports/`](results/reports/); events and deviations after registration are in
[`results/DEVIATIONS.md`](results/DEVIATIONS.md).

D = M(learned) − M(floor) on the same cells, averaged per unit (patient, donor or animal); estimate = intercept of
D ~ 1 + centred platform strata (= mean D); one-sided 95% bounds from the exact t. Non-inferior: upper bound < δ_f.
Learned superior: lower bound > 0 after Holm over the 12 pairs. Verdicts refer to region scale (~50,000-cell crops or
whole small samples, PROTOCOL §3).

## Tier 1: confirmatory verdicts

| Pair | Units | Floor | Learned | D | One-sided 95% bounds | δ_f | Non-inferior | Learned superior (Holm) | Equivalent | Floor superior |
|---|---|---|---|---|---|---|---|---|---|---|
| F1 cell identity, scVI | 25 | .163 | .155 | −.0087 | −.0112, −.0062 | .0195 | yes | no | yes | yes |
| F1 cell identity, Nicheformer | 25 | .163 | .126 | −.0379 | −.0432, −.0327 | .0195 | yes | no | no | yes |
| F2 spatial domains, CellCharter | 10 | .461 | .360 | −.1003 | −.1814, −.0193 | .05 | yes | no | no | yes |
| F2 spatial domains, Novae (mouse) | 4 | .439 | .324 | −.1147 | −.1437, −.0857 | .05 | yes | no | no | yes |
| F3 context → receiver state, LightGBM | 22 | .073 | .057 | −.0155 | −.0182, −.0127 | .0043 | yes | no | no | yes |
| F3 context → receiver state, NCEM | 22 | .073 | .037 | −.0355 | −.0413, −.0298 | .0043 | yes | no | no | yes |
| F4 counterfactual, Cellina, official | 10 | .873 | .882 | +.0095 | +.0011, +.0180 | .03 | yes | no | yes | no |
| F4 counterfactual, CPA, official | 10 | .873 | .768 | −.1045 | −.1501, −.0588 | .03 | yes | no | no | yes |
| F4 counterfactual, Cellina, spike-free | 5 | .887 | .896 | +.0086 | −.0096, +.0268 | .03 | yes | no | yes | no |
| F4 counterfactual, CPA, spike-free | 5 | .887 | .820 | −.0677 | −.1102, −.0253 | .03 | yes | no | no | yes |
| F5 label transfer, scANVI | 20 | .719 | .704 | −.0143 | −.0225, −.0062 | .0255 | yes | no | yes | yes |
| F5 label transfer, Tangram | 20 | .719 | .605 | −.1132 | −.1431, −.0833 | .0255 | yes | no | no | yes |

- **H1, per pair:** the floor is non-inferior on all 12 pairs.
- **H1, joint (intersection-union):** met.
- **H2:** no learned model is superior after Holm. Before correction, Cellina's lower bound on the official metric is
  above 0 (one-sided p = .035).
- **Floor superiority** (upper bound < 0; no decision weight) holds on 10 of 12 pairs.
- **Units:** F3 lost 3 of 25 crops before scoring. Their label-free receiver identities form a single lineage, so no
  model has an other-lineage context (`results/DEVIATIONS.md`). Donors 7 and 170 are therefore missing, and donor 280
  is scored on sample 280N only.
- **F4 retrospective units:** 5 of the 10 official-metric units carry the authors' published Cellina and CPA results,
  and their comparison with L6 was computed before registration (PROTOCOL §3). The spike-free pairs use the 5
  prospective units only.

### Registered sensitivity analyses

| Pair | Wilcoxon, D − δ_f < 0 (p) | Wilcoxon, D > 0 (p) | Leave-one-cancer-out: non-inferior in every split | Mean D per platform stratum |
|---|---|---|---|---|
| F1 scVI | <.001 | 1.00 | yes | CosMx 1k −.007; CosMx WTx −.014; MERFISH −.002; Xenium 5K −.007 |
| F1 Nicheformer | <.001 | 1.00 | yes | CosMx 1k −.022; CosMx WTx −.073; MERFISH −.012; Xenium 5K −.002 |
| F2 CellCharter | .001 | .995 | yes | CosMx WTx −.140; MERFISH −.041 |
| F2 Novae | .063 | 1.00 | (one cancer) | MERFISH −.115 |
| F3 LightGBM | <.001 | 1.00 | yes | CosMx 1k −.013; CosMx WTx −.013; MERFISH −.027 |
| F3 NCEM | <.001 | 1.00 | yes | CosMx 1k −.028; CosMx WTx −.020; MERFISH −.087 |
| F4 Cellina, official | .005 | .097 | no: CRC only, upper +.0314 | CosMx WTx +.017; MERFISH −.001 |
| F4 CPA, official | .001 | 1.00 | yes | CosMx WTx −.151; MERFISH −.035 |
| F4 Cellina, spike-free | .063 | .50 | no: CRC only (2 units), upper +.155 | CosMx WTx +.023; MERFISH −.001 |
| F4 CPA, spike-free | .031 | 1.00 | no: CRC only (2 units), upper +.178 | CosMx WTx −.078; MERFISH −.061 |
| F5 scANVI | <.001 | .993 | no: without CRC (3 units), upper +.0525 | CosMx 1k −.026; CosMx WTx −.000 |
| F5 Tangram | <.001 | 1.00 | no: without CRC (3 units), upper +.0298 | CosMx 1k −.125; CosMx WTx −.098 |

- **Wilcoxon:** with 4 units the smallest attainable one-sided signed-rank p is .0625, so the Novae entry cannot reach
  .05. With 5 units the smallest is .031, and Cellina spike-free gives .063.
- **Leave-one-cancer-out:** SPATCH Xenium 5K (F1 only) is colon adenocarcinoma and counts as CRC. The registered
  analysis code labelled it "tissue unverified" and so left it out as a group of its own; the label is corrected
  (`results/DEVIATIONS.md`, `results/reports/P5_loco_crc.md`). This changes only F1 without CRC (7 units; upper
  bound scVI −.0061, Nicheformer −.0269; with the registered label 8 units, −.0060, −.0232) and removes the SPATCH
  split. Every F1 split is non-inferior with either label.
- **Nicheformer without the Yao animal** (in its pretraining material): D −.0388, upper −.0333.
- **Novae collapse rule** (effective rank < 5, or > 50% of cells in one domain): no mouse run collapsed (effective rank
  7.6–8.7, largest domain ≤ 39%). 5 of the 6 human zero-shot runs met the rule; human Novae is not a test pair.

## Tier 2: tuning transfer and capacity (diagnostic, no decision weight)

The transfer gap T = best configuration on the unit − locked configuration (≥ 0). The capacity gap C = best learned − best
floor. Models trained on the unit ran on the seeded subset (PROTOCOL §5).

| Model | Units | T | Upper | Locked = best | Reading |
|---|---|---|---|---|---|
| F1 floor | 25 | +.0043 | +.0057 | 19% | transfers |
| F1 scVI | 14 | +.0021 | +.0031 | 20% | transfers |
| F1 Nicheformer | 25 | +.0001 | +.0003 | 92% | transfers |
| F2 floor | 10 | +.0835 | +.1356 | 0% | inconclusive |
| F2 CellCharter | 5 | +.0241 | +.0766 | 60% | inconclusive |
| F2 Novae (mouse) | 4 | .0000 | .0000 | 100% | transfers |
| F3 floor | 22 | +.0014 | +.0020 | 14% | transfers |
| F3 LightGBM | 22 | +.0004 | +.0008 | 91% | transfers |
| F3 NCEM | 11 | +.0022 | +.0049 | 73% | inconclusive |
| F5 floor | 20 | +.0025 | +.0031 | 5% | transfers |
| F5 scANVI | 11 | +.0024 | +.0043 | 33% | transfers |
| F5 Tangram | 11 | +.0034 | +.0070 | 75% | transfers |

| Pair | Units | D (same units) | C | One-sided 95% bounds | T learned | T floor | Reading |
|---|---|---|---|---|---|---|---|
| F1 scVI | 14 | −.0093 | −.0122 | −.0162, −.0082 | +.0019 | +.0047 | floor holds at best tuning |
| F1 Nicheformer | 25 | −.0379 | −.0421 | −.0471, −.0371 | +.0001 | +.0042 | floor holds at best tuning |
| F2 CellCharter | 5 | −.1242 | −.2464 | −.4262, −.0666 | +.0241 | +.1463 | floor holds at best tuning |
| F2 Novae (mouse) | 4 | −.1147 | −.1402 | −.1724, −.1080 | .0000 | +.0255 | floor holds at best tuning |
| F3 LightGBM | 22 | −.0155 | −.0166 | −.0194, −.0137 | +.0004 | +.0014 | floor holds at best tuning |
| F3 NCEM | 11 | −.0352 | −.0347 | −.0429, −.0264 | +.0022 | +.0017 | floor holds at best tuning |
| F5 scANVI | 11 | −.0170 | −.0166 | −.0297, −.0036 | +.0022 | +.0021 | floor holds at best tuning |
| F5 Tangram | 11 | −.1196 | −.1183 | −.1575, −.0791 | +.0045 | +.0021 | floor holds at best tuning |

- F2 floor: the multi-scale secondary floor is the best configuration on 8 of 10 units.
- Fine-tuned Novae, balanced ARI .08–.20 on the four animals; zero-shot .29–.41.

## Secondary outputs (no decision weight)

From `results/reports/P5_secondary.md` and the F3 decomposition in `P5_corrected.md`. D is the first named minus the
second, with one-sided 95% bounds.

| Comparison | Units | D | Bounds |
|---|---|---|---|
| F2: multi-scale floor − locked floor | 10 | +.052 | −.038, +.141 |
| F2: CellCharter − multi-scale floor | 10 | −.152 | −.269, −.035 |
| F2: CellCharter on PCA input − locked floor | 10 | −.115 | −.197, −.033 |
| F2: CellCharter (scVI input) − CellCharter (PCA input) | 10 | +.015 | +.000, +.029 |
| F3: v3 floor − v2 floor (same-lineage rings) | 22 | +.0150 | +.0120, +.0179 |
| F3: LightGBM − v2 floor | 22 | −.0005 | −.0034, +.0024 |
| F3: NCEM − v2 floor | 22 | −.0206 | −.0246, −.0166 |
| F3: NCEM − ridge on NCEM's own inputs | 22 | −.0051 | −.0066, −.0037 |
| F3: ED-NCEM − floor | 22 | −.0412 | −.0466, −.0359 |
| F5: floor − InSituType without the spatial prior (v2 floor) | 20 | +.0043 | +.0032, +.0054 |
| F5: scANVI − v2 floor | 20 | −.0100 | −.0177, −.0023 |
| F5: logistic regression (C .1, reference scaling) − floor | 20 | −.0888 | −.1024, −.0752 |
| F5: IF-gated macro-F1 on GSE1, scANVI − floor | 6 | −.053 | −.088, −.018 |
| F5: IF-gated macro-F1 on GSE1, Tangram − floor | 6 | −.397 | −.443, −.352 |

| F5 model | Held-out-marker AUROC | IF-gated macro-F1 (GSE1) | Marker-positive recall, epi / imm / mye | Predicted share, epi / imm / mye / other |
|---|---|---|---|---|
| Floor (InSituType ML + spatial prior) | .723 | .657 | .84 / .26 / .47 | .51 / .10 / .13 / .26 |
| InSituType ML without the prior | .719 | .654 | .83 / .27 / .48 | .49 / .11 / .13 / .27 |
| Logistic regression | .635 | .528 | .75 / .36 / .41 | .48 / .25 / .14 / .13 |
| scANVI | .709 | .605 | .83 / .45 / .47 | .47 / .24 / .11 / .18 |
| Tangram | .605 | .260 | .20 / .62 / .32 | .11 / .51 / .16 / .22 |

- F3: the v2 floor uses other lineages only; v3 adds same-lineage rings (PROTOCOL §2, F3).
- F5 marker-positive recall is undefined on crops where another class's 90th percentile is 0: it is defined on 19, 21
  and 19 of the 21 F5 crops.
- F2 geometry reference: k-means on x, y alone, mean balanced ARI .317 (floor .461); at or above the floor on GSE1 210
  and 221.

## Exploratory analyses (run after the tier-1 report; no decision weight)

D = learned − floor with one-sided 95% bounds unless stated; the report of each analysis is in `results/reports/`.

- **Interior-only re-score** (`X_interior`): cells within R µm of a crop's cut edge left unscored.
  - CellCharter: −.100 (R = 0), −.093 (100 µm), −.095 (250 µm), −.073 (500 µm).
  - Novae: −.115 to −.130. scANVI: −.014 to −.015. Tangram: −.111 to −.113.
- **Unit size** (`X_moderators`): Spearman ρ of D against log cells scored; F5 scANVI ρ = +.53 (p = .017, uncorrected);
  other pairs |ρ| ≤ .36, p ≥ .12.
- **F5 at a matched call rate** (`X_f5_matched`; added after `results-v1`): the floor's calls per class re-made at the
  comparator's call rate (its top cells by class posterior); nothing refitted.
  - scANVI − floor at scANVI's call rate: −.0040 (−.0075, −.0005); registered −.0143.
  - By class: immune non-myeloid +.025 (+.019, +.031), epithelial −.012, myeloid −.025.
  - By stratum: CosMx 1k +.002, CosMx WTx −.011.
  - Tangram at its own call rate: −.035 (−.048, −.022).
  - Floor ranked by class posterior at its own call rate: AUROC .741 (registered .723).
- **F5 with depth held fixed** (`X_f5_depth`; F5 stress panel, PROTOCOL §6): the held-out-marker AUROC computed within
  each library-size quintile, mean over quintiles; registered scores recomputed exactly on every crop.
  - scANVI − floor: +.0015 (−.0044, +.0074); registered −.0143.
  - Tangram − floor: −.1526 (−.1798, −.1253).
  - Share of epithelial calls in the deepest quintile: Tangram 70%, scANVI 34%, floor 31%.
- **F1 with Euclidean k-NN** (`X_f1_euclid`; added 2026-11-09): every saved locked embedding re-scored with Euclidean
  instead of cosine neighbours.
  - Crop means, cosine → Euclidean: scVI .154 → .111; floor .162 → .149.
  - Euclidean for both: scVI −.0376, Nicheformer −.0256 (registered, cosine: −.0087, −.0379).
- **Whole-slide arm** (`X_wholeslide`; plan `results/WHOLESLIDE_PLAN.md`, committed before any whole-slide output):
  the six GSE1 confirmation slides whole (141k–741k cells), for the four pairs whose learned model trains on the unit.
  68 jobs, none failed.
  - All four pairs non-inferior; no learned model superior after Holm.
  - D whole slide vs crop of the same slides: scVI −.0076 vs −.0118; CellCharter −.155 vs −.140; LightGBM −.0071 vs
    −.0126; scANVI −.0639 vs −.0015.
  - LightGBM − v2 floor: +.0167 (+.0036, +.0297).
- **F5 at slide scale** (`X_f5_scale`; diagnostic after the arm). Query-to-reference cell ratio on the whole slides:
  6–11 to 1.
  - GSE1 221, the crop's number of training epochs on every cell: AUROC .588 (registered recipe .588).
  - scANVI trained on the reference plus 50,000 random cells of the slide, then labelling every cell: D −.0106
    (−.0253, +.0041) over the six slides (registered recipe −.0639; crops −.0015).
- **F4 by cancer and provenance** (`X_f4_split`; subsets chosen after tier 1).
  - Cellina, official metric, CRC (6 slides): +.0165 (+.0017, +.0314); D > 0 on 5 of 6 slides, > .03 on 2 (110
    +.044, run here; 120 +.031, the authors' run).
  - Cellina, mouse (4 animals): −.0010 (−.0051, +.0031).
  - CPA: D < 0 in every subset.
- **Replication levels** (`X_levels`): tier-1 D re-aggregated per slide (GSE303070's 11 donors sit on 5 slides) and
  per study (the two Allen sources counted separately). Non-inferior at both levels in every pair except Cellina
  spike-free at study level (3 studies, upper +.0311, δ .03).
- **Positive controls** (`X_posctl`; run 2026-10-28 to 11-02, before registration): scVI, CellCharter, Novae,
  LightGBM, scANVI and Tangram pass; linear NCEM (the locked configuration) fails the planted-composition control;
  ED-NCEM passes.
- **Development effort** (`X_effort`): configurations evaluated on development units. Floor grids: F1 48, F2 225 (+ 173
  in reference and k sweeps), F3 144, F5 84. Learned models: 2–6 each.
- **F3 LightGBM with early stopping** (`X_f3_lgbm_es`; added after `results-v1.1`): up to 1,000 rounds, patience 30,
  validation on held-out FOVs or tiles; registered: fixed 100 rounds. Same 22 units.
  - All v3 feature groups: −.0211 (−.0234, −.0188).
  - The floor's features only: −.0237 (−.0265, −.0209).
  - Registered: −.0155. Median stopping round: 5–6.
- **F5 Tangram on normalised input** (`X_f5_tangram_norm`; added after `results-v1.1`): both sides normalised before
  mapping, everything else as registered (registered: raw counts).
  - Counts per 10,000: −.1044 (−.1290, −.0798).
  - log1p counts per 10,000: −.1094 (−.1391, −.0796).
  - Registered: −.1132. Predicted share immune non-myeloid / epithelial: .51 / .11 (raw), .49 / .10 (counts per
    10,000), .51 / .09 (log1p).
- **Specification curve** (`X_speccurve`; added after `results-v1.1`): every estimate of D already computed in the
  reports above under another choice of units, aggregation level, tuning, metric, floor, learned-model configuration or
  scale; no new run. Rows of one pair share units; no test across rows.
  - 146 specifications on the registered scale: D < 0 in 131; non-inferior in 135; upper bound < 0 in 117; lower bound
    > 0 in 3 (Cellina registered; Cellina, six CRC units; LightGBM − v2 floor on whole slides).
  - Not non-inferior: 11 specifications, each with 6 units or fewer; 7 in F4.
- **Common scale and meta-analysis over studies** (`X_meta`; added after `results-v1.1`): D as a percentage of the
  floor's above-chance score A; random-effects model with studies (data releases) as the random factor (REML, one
  within-study variance; 2–8 studies per pair).
  - Registered D, % of A: scVI −3.6, scANVI −6.6, CPA −12.0, Nicheformer −15.7, LightGBM −21.3, CellCharter −21.8,
    Novae −26.1, NCEM −48.9, Tangram −51.8, Cellina +1.1 (F4 official metric). Margins: 3.4–11.6.
  - Mean over studies non-inferior in 10 of 12 pairs; not non-inferior: CPA official, Cellina spike-free (3 studies
    each).
  - 95% prediction interval for a new study below the margin: scVI, LightGBM, Tangram; with unweighted study means:
    scVI, CPA spike-free.
- **Baseline ladder** (`X_ladder`; added after `results-v1.1`): naive baselines per family. New CPU runs on the
  confirmation crops: F1 log-normalised PCA 50, F2 GMM on expression PCA 50 without a neighbourhood term, F3 ridge on
  the neighbourhood type mix only (the three single-lineage F3 crops fail as in tier 1). Existing: F2 k-means on
  coordinates, F3 v2 floor, F4 rungs L1–L5 and the control mean, F5 logistic regression and InSituType without the
  spatial prior.
  - Floor − each naive baseline: lower bound > 0 in all five families.
  - Learned model − each naive baseline, lower bound > 0 against all: Cellina only. Lower bound ≤ 0 against at least
    one: scVI and Nicheformer vs log-normalised PCA (−.0021, −.0313); CellCharter and Novae vs k-means on coordinates;
    LightGBM vs v2 floor (−.0005); NCEM vs type-mix ridge (−.0027); scANVI vs InSituType without the prior (−.0100);
    CPA vs L3–L5; Tangram vs logistic regression and InSituType without the prior.

## Reproducing the analysis

Data are public (PROTOCOL §3) and not included. The per-crop outputs the analysis reads are in `results/p4/` (scores and
configurations as JSON / CSV; no labels, embeddings or counterfactuals). With the confirmation outputs under
`runs/p4` (the job list `registration/p4_jobs.txt`, run through `scripts/sched.sh`):

```
zsh scripts/p5_tier1_run.sh   # F1 rescore, F4 scores, tier 1 -> runs/p4/P5_corrected, interior re-score
zsh scripts/p5_tier2_run.sh   # tier 2 -> runs/p4/T2_corrected, secondary outputs, unit-size moderators
```

`P5_corrected.json` (the tier-1 report) has SHA-256 `6d2b2b18…`, recorded in the tier-2 report `T2_corrected.json`. Run
outputs carry the commit hash of the development repository in which they ran; the code as run is this repository at
`results-v1`.
