# Results: closed-form floors versus learned models in spatial omics

Confirmatory results of the protocol registered as DOI [10.5281/zenodo.23931264](https://doi.org/10.5281/zenodo.23931264)
(`registration/PROTOCOL.md`, tag `registration-v1`). Release `results-v1`, 2026-11-08, deposited on Zenodo as version 2
of the registration record, DOI [10.5281/zenodo.24044400](https://doi.org/10.5281/zenodo.24044400); version 1 is
unchanged. Every number below is read from
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
  above 0 (one-sided p = .035); its lead is inside the margin on both F4 metrics (equivalent).
- **Floor superiority** (upper bound < 0; no decision weight) holds on 10 of 12 pairs.
- **Units:** F3 lost 3 of 25 crops before scoring. Their label-free receiver identities form a single lineage, so no
  model has an other-lineage context (`results/DEVIATIONS.md`). Donors 7 and 170 are therefore missing, and donor 280
  is scored on sample 280N only.
- **First tier-1 run:** it left GSE1 222 out of F4 (it read 222's scores with the wrong held-out type), giving N = 9 / 4.
  The corrected report above is the registered analysis on every registered unit. The first run is kept as
  `results/reports/P5.{json,md}`, and its verdicts are the same (`results/DEVIATIONS.md`).
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
- **Leave-one-cancer-out:** the failures are in splits of 2–6 units. F4 non-inferiority for Cellina rests on the mouse
  sections, where D ≈ 0. On the six CRC slides Cellina leads by +.017 (official), close to the .03 margin. F5
  non-inferiority rests on the CRC units. On the three Bruker units (kidney, skin, breast) scANVI is level with the
  floor (D +.002, upper +.0525).
- **Nicheformer without the Yao animal** (in its pretraining material): D −.0388, upper −.0333.
- **Novae collapse rule** (effective rank < 5, or > 50% of cells in one domain): no mouse run collapsed (effective rank
  7.6–8.7, largest domain ≤ 39%). 5 of the 6 human zero-shot runs met the rule; human Novae is not a test pair.

## Tier 2: tuning transfer and capacity (diagnostic, no decision weight)

The transfer gap T = best configuration on the unit − locked configuration (≥ 0). The capacity gap C = best learned − best
floor. Models trained on the unit ran on the seeded subset (PROTOCOL §5). Selection on the scored unit makes the best
scores optimistic.

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

- **Development tuning did not handicap the learned models.** Their mean transfer gaps are at most .0034 outside F2,
  and .024 for CellCharter.
- **The F2 floor is the one under-tuned model.** The multi-scale secondary floor is best on 8 of 10 units. At best
  tuning the floor's lead grows in every F2 pair.
- **Fine-tuned Novae** scored below zero-shot on all four animals (balanced ARI .08–.20 against .29–.41).

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

- **F3:** LightGBM matches the v2 floor (other lineages only). The floor's F3 lead over LightGBM comes from the
  same-lineage rings added in v3, which carry tissue-scale gradients of the receiver's own state (PROTOCOL §2, F3).
  NCEM is below a ridge on its own inputs.
- **F5:** the floor assigns fewer cells to the immune non-myeloid class than scANVI (10% against 24%) and recovers fewer
  immune marker-positive cells (.26 against .45). Its AUROC lead comes with a more conservative immune call. Recall is
  undefined on crops where another class's 90th percentile is 0: it is defined on 19, 21 and 19 of the 21 F5 crops.
- **F2 geometry reference:** k-means on x, y alone scores .317 on average against the floor's .461. It is at or above
  the floor on GSE1 210 and 221.

## Exploratory analyses (run after the tier-1 report; no decision weight)

- **Interior-only re-score** (`results/reports/X_interior.md`): with cells within 100–500 µm of a crop's cut edge left
  unscored, D stays within .03 of the full-crop value for F2 and within .003 for F5.
  - CellCharter: −.100 (full crop), −.093 (100 µm), −.095 (250 µm), −.073 (500 µm).
  - Novae: −.115 to −.130.
  - scANVI: −.014 to −.015.
  - Tangram: −.111 to −.113.
- **Unit size** (`results/reports/X_moderators.md`): no association between D and the cells scored, except F5 scANVI
  (Spearman ρ = +.53, p = .017, no multiplicity correction). Unit size is confounded with platform there: the CosMx 1k
  samples are the small units and carry the floor's lead.
- **F5 at a matched call rate** (`results/reports/X_f5_matched.md`; added after `results-v1`). The
  held-out-marker AUROC can rise when a model calls fewer cells of a class. Here the floor calls, per class, as many
  cells as the comparator: its top cells by class posterior. Nothing is refitted.
  - scANVI − floor at scANVI's call rate: D −.0040 (one-sided 95% bounds −.0075, −.0005), against −.0143 registered.
    About 70% of the floor's lead over scANVI comes from its more conservative calling.
  - By class: immune non-myeloid +.025 (+.019, +.031; scANVI ranks these cells better), epithelial −.012, myeloid
    −.025.
  - By stratum: CosMx 1k +.002, CosMx WTx −.011.
  - Tangram at its own call rate: −.035 (−.048, −.022).
  - Ranking by class posterior raises the floor at its own call rate from .723 to .741, so matching does not handicap
    the floor. The reading is that scANVI and the floor are level on F5.
- **F1 with Euclidean k-NN** (`results/reports/X_f1_euclid.md`; added 2026-11-09). Every saved locked embedding was
  re-scored with Euclidean instead of cosine neighbours.
  - Crop means, cosine → Euclidean: scVI .154 → .111; floor .162 → .149.
  - D with Euclidean for both: scVI −.0376, Nicheformer −.0256; registered (cosine): −.0087 and −.0379.
  - The registered cosine metric is the more favourable choice for scVI.
- **Whole-slide arm** (`results/WHOLESLIDE_PLAN.md`; report `results/reports/X_wholeslide.md`): the six GSE1
  confirmation slides as whole slides (141k–741k cells), for the four pairs whose learned model trains on the unit. The
  plan was committed before any whole-slide output existed. All 68 jobs ran; none failed.
  - All four pairs are non-inferior on whole slides and no learned model is ahead after Holm.
  - D on whole slides vs crops of the same slides: scVI −.0076 vs −.0118; CellCharter −.155 vs −.140; LightGBM −.0071
    vs −.0126; scANVI −.0639 vs −.0015.
  - The F1 and F3 floor leads shrink at slide scale, as on the validation slide (PROTOCOL Appendix A).
  - LightGBM − v2 floor (no rings): +.0167 (one-sided 95% bounds +.0036, +.0297). On whole slides LightGBM beats the
    floor without the same-lineage rings.
- **F5 at slide scale** (`results/reports/X_f5_scale.md`; diagnostic added after the arm). scANVI's whole-slide loss
  comes from the ratio of query to reference cells, not from the model's capacity.
  - On GSE1 221, training for the crop's number of epochs on every cell changed nothing (.588).
  - Training on the reference plus 50,000 random cells of the slide, then labelling every cell, restores scANVI:
    D −.0106 (−.0253, +.0041) over the six slides, against −.0639 for the registered recipe and −.0015 on the crops.
  - The registered recipe (joint model of reference and the whole query) degrades when the query outnumbers the
    reference 6–11 to 1.
- **F4 by cancer and provenance** (`results/reports/X_f4_split.md`; subsets chosen after tier 1).
  - Cellina, official metric, CRC (6 slides): +.0165 (+.0017, +.0314). Ahead on 5 of 6 slides, above the .03 margin
    on 2 (110 +.044, run here; 120 +.031, the authors' run).
  - Mouse (4 animals): −.0010 (−.0051, +.0031).
  - CPA is below the floor in every subset.
- **Replication levels** (`results/reports/X_levels.md`): tier-1 D re-aggregated per slide (GSE303070's 11 donors sit
  on 5 slides) and per study. Every pair stays non-inferior at both levels except Cellina spike-free at study level
  (3 studies, upper +.0311 against .03). The two Allen sources (Zhuang lab, Yao et al.) count as separate studies.
- **Positive controls** (`results/reports/X_posctl.md`; run 2026-10-28 to 11-02, before registration). scVI,
  CellCharter, Novae, LightGBM, scANVI and Tangram pass. The locked NCEM configuration (linear NCEM) fails the
  planted-composition control; the ED-NCEM variant passes.
- **Development effort** (`results/reports/X_effort.md`): configurations evaluated on development units per model.
  Floor grids: F1 48, F2 225 (+ 173 in reference and k sweeps), F3 144, F5 84. Learned models: 2–6 each. The floors
  were tuned far more than the learned models before the lock; tier 2 (above) measures what tuning on the scored unit
  would add.
- **Not run in this release:** stress panels on the confirmation units, the protein shift-null on IF stains, and the F1
  depth-matched reference (PROTOCOL §6).

## Reproducing the analysis

Data are public (PROTOCOL §3) and not included. The per-crop outputs the analysis reads are in `results/p4/` (scores and
configurations as JSON / CSV; no labels, embeddings or counterfactuals). With the confirmation outputs under
`runs/p4` (the job list `registration/p4_jobs.txt`, run through `scripts/sched.sh`):

```
zsh scripts/p5_tier1_run.sh   # F1 rescore, F4 scores, tier 1 -> runs/p4/P5_corrected, interior re-score
zsh scripts/p5_tier2_run.sh   # tier 2 -> runs/p4/T2_corrected, secondary outputs, unit-size moderators
```

`P5_corrected.json` has SHA-256 `6d2b2b18…` and is the tier-1 report recorded in `T2_corrected.json`. The first tier-1
report, `P5.json` (`8c855932…`), is the one recorded in `T2.json`. The two tier-2 reports are otherwise identical. Run
outputs carry the commit hash of the development repository in which they ran; the code as run is this repository at
`results-v1`.
