# X_ceiling (exploratory)

F5: held-out-marker counts split by binomial thinning into halves A and B (5 draws x 2 directions). Ceiling = per
class, the cells with the highest half-A marker score called at the learned model's call rates, scored on half B.
Floor and learned labellings scored on half B. Means over units; headroom share = (learned - floor) / (ceiling at
the floor's call rates - floor). Above ceiling: mean over draws of (score - own-rate oracle) - 2 SE > 0.

| pair | units | floor | learned | ceiling | ceiling at floor rate | headroom share | crops learned > ceiling | crops floor > ceiling |
|---|---|---|---|---|---|---|---|---|
| F5|scanvi | 20 | 0.6868 | 0.6793 | 0.7117 | 0.7176 | -0.241 | 4/21 | 5/21 |
| F5|tangram | 20 | 0.6868 | 0.6011 | 0.6894 | 0.7176 | -2.778 | 0/21 | 5/21 |

Full-data scores equal the registered per-crop scores: True

F2: ceiling = out-of-fold logistic regression from the locked floor's features (ell 400 µm, λ .8) to the annotated
regions, 5 random cell folds, balanced ARI on the scored cells; quadrant folds (4 spatial blocks) reported beside it.
Floor and learned: balanced ARI of the saved partitions (seed mean). One crop per unit; means over units.

| pair | units | floor | learned | ceiling | ceiling, quadrant folds | headroom share | units learned > ceiling | units floor > ceiling |
|---|---|---|---|---|---|---|---|---|
| F2|cellcharter | 10 | 0.4606 | 0.3602 | 0.8853 | 0.3701 | -0.236 | 0/10 | 0/10 |
| F2|novae | 4 | 0.4392 | 0.3245 | 0.8753 | 0.2583 | -0.263 | 0/4 | 0/4 |

Floor equals the registered per-crop floor: True

F3: ceiling = per receiver, masked-gene counts split by binomial thinning (2 draws x 2 directions); per held-out
quadrant the full-data training fit (depth covariate, z-score, PCA top 10) applied to each half's test rows;
per PC sqrt of the Spearman-Brown corrected half-A / half-B Spearman rho; mean over PCs, folds, draws, receivers.
Floor and learned: registered per-crop unit scores. Crops -> units as registered; means over units.

| pair | units | floor | learned | ceiling | headroom share | units learned > ceiling | units floor > ceiling |
|---|---|---|---|---|---|---|---|
| F3|lightgbm | 22 | 0.0726 | 0.0572 | 0.6681 | -0.0260 | 0/22 | 0/22 |
| F3|ncem | 22 | 0.0726 | 0.0371 | 0.6681 | -0.0597 | 0/22 | 0/22 |

F1: oracle reference (not a bound) = neighbours from the expression PCA (50 PCs) of the held-out count half
itself, kNN-Poisson score as registered (count split seed 0). Floor and learned: registered rescored scores.
Share of reference gap = (learned - floor) / (reference - floor); n/a when reference - floor < .01.

| pair | units | floor | learned | reference | share of reference gap | units learned > reference | units floor > reference |
|---|---|---|---|---|---|---|---|
| F1|scvi | 25 | 0.1635 | 0.1548 | 0.1633 | n/a | 6/25 | 15/25 |
| F1|nicheformer | 25 | 0.1635 | 0.1256 | 0.1633 | n/a | 2/25 | 15/25 |

Checks before the round-1 units (results/p4/x_ceiling/checks): F5 on development unit gse1_232 (floor = v2
configuration), F2 on sim.domain_toy and on gse1_232, F3 on gse1_232 (development grid has no v3 floor; LightGBM
shown).

- f1_gse1_232: ceiling 0.0845, nicheformer 0.0121
- f2_gse1_232: ceiling 0.9910, quadrant folds 0.8009, floor 0.9292, cellcharter 0.4384
- f2_toy: ceiling 0.9880, quadrant folds 0.9732, floor 0.9831
- f3_gse1_232: ceiling 0.5473, lightgbm 0.0194
- f5_gse1_232: half B floor 0.5955, scanvi 0.6328, tangram 0.5962, oracle_floor 0.6455, oracle_scanvi 0.6597, oracle_tangram 0.6334
