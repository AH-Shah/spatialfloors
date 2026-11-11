# X_ladder: baseline ladder (exploratory)

Naive baselines per family. D = first − second per unit, one-sided 95% bounds; ahead = lower bound > 0.

## F1

| Rung | Role | Units | Mean score |
|---|---|---|---|
| log-normalised PCA 50 | naive | 25 | 0.1569 |
| floor (Pearson-residual PCA 50) | floor | 25 | 0.1635 |
| scVI | learned | 25 | 0.1548 |
| Nicheformer | learned | 25 | 0.1256 |

| Comparison | Units | D | Bounds | Ahead |
|---|---|---|---|---|
| floor (Pearson-residual PCA 50) − log-normalised PCA 50 | 25 | +0.0066 | +0.0036, +0.0096 | yes |
| scVI − log-normalised PCA 50 | 25 | -0.0021 | -0.0038, -0.0005 | no |
| Nicheformer − log-normalised PCA 50 | 25 | -0.0313 | -0.0350, -0.0277 | no |

## F2

| Rung | Role | Units | Mean score |
|---|---|---|---|
| k-means on coordinates | naive | 10 | 0.3169 |
| GMM on expression PCA 50 (no neighbourhood) | naive | 10 | 0.1482 |
| floor (expression + neighbourhood mean) | floor | 10 | 0.4606 |
| CellCharter | learned | 10 | 0.3602 |
| Novae zero-shot (mouse units) | learned | 4 | 0.3245 |

| Comparison | Units | D | Bounds | Ahead |
|---|---|---|---|---|
| floor (expression + neighbourhood mean) − k-means on coordinates | 10 | +0.1437 | +0.0147, +0.2727 | yes |
| CellCharter − k-means on coordinates | 10 | +0.0433 | -0.0564, +0.1430 | no |
| Novae zero-shot (mouse units) − k-means on coordinates | 4 | -0.0240 | -0.1136, +0.0657 | no |
| floor (expression + neighbourhood mean) − GMM on expression PCA 50 (no neighbourhood) | 10 | +0.3124 | +0.2222, +0.4027 | yes |
| CellCharter − GMM on expression PCA 50 (no neighbourhood) | 10 | +0.2121 | +0.1758, +0.2483 | yes |
| Novae zero-shot (mouse units) − GMM on expression PCA 50 (no neighbourhood) | 4 | +0.1334 | +0.0684, +0.1983 | yes |

## F3

| Rung | Role | Units | Mean score |
|---|---|---|---|
| ridge on the type mix only | naive | 22 | 0.0398 |
| floor without rings (v2) | naive | 22 | 0.0576 |
| floor (v3) | floor | 22 | 0.0726 |
| LightGBM | learned | 22 | 0.0572 |
| NCEM | learned | 22 | 0.0371 |

| Comparison | Units | D | Bounds | Ahead |
|---|---|---|---|---|
| floor (v3) − ridge on the type mix only | 22 | +0.0328 | +0.0281, +0.0375 | yes |
| LightGBM − ridge on the type mix only | 22 | +0.0173 | +0.0134, +0.0213 | yes |
| NCEM − ridge on the type mix only | 22 | -0.0027 | -0.0052, -0.0003 | no |
| floor (v3) − floor without rings (v2) | 22 | +0.0150 | +0.0120, +0.0179 | yes |
| LightGBM − floor without rings (v2) | 22 | -0.0005 | -0.0034, +0.0024 | no |
| NCEM − floor without rings (v2) | 22 | -0.0206 | -0.0246, -0.0166 | no |

## F4

| Rung | Role | Units | Mean score |
|---|---|---|---|
| control mean | naive | 10 | 0.4823 |
| L1 shared shift | naive | 10 | 0.4471 |
| L2 contamination | naive | 10 | 0.5688 |
| L3 niche ridge | naive | 10 | 0.7155 |
| L4 niche match | naive | 10 | 0.6986 |
| L5 pooled ridge | naive | 10 | 0.7810 |
| floor (L6 two-way pseudobulk) | floor | 10 | 0.8727 |
| Cellina | learned | 10 | 0.8823 |
| CPA | learned | 10 | 0.7683 |

| Comparison | Units | D | Bounds | Ahead |
|---|---|---|---|---|
| floor (L6 two-way pseudobulk) − control mean | 10 | +0.3905 | +0.3362, +0.4447 | yes |
| Cellina − control mean | 10 | +0.4000 | +0.3441, +0.4559 | yes |
| CPA − control mean | 10 | +0.2860 | +0.2427, +0.3293 | yes |
| floor (L6 two-way pseudobulk) − L1 shared shift | 10 | +0.4256 | +0.3668, +0.4844 | yes |
| Cellina − L1 shared shift | 10 | +0.4351 | +0.3730, +0.4972 | yes |
| CPA − L1 shared shift | 10 | +0.3211 | +0.2483, +0.3939 | yes |
| floor (L6 two-way pseudobulk) − L2 contamination | 10 | +0.3039 | +0.2387, +0.3692 | yes |
| Cellina − L2 contamination | 10 | +0.3135 | +0.2506, +0.3763 | yes |
| CPA − L2 contamination | 10 | +0.1995 | +0.1266, +0.2723 | yes |
| floor (L6 two-way pseudobulk) − L3 niche ridge | 10 | +0.1572 | +0.0922, +0.2221 | yes |
| Cellina − L3 niche ridge | 10 | +0.1667 | +0.1081, +0.2253 | yes |
| CPA − L3 niche ridge | 10 | +0.0527 | -0.0164, +0.1218 | no |
| floor (L6 two-way pseudobulk) − L4 niche match | 10 | +0.1741 | +0.0813, +0.2669 | yes |
| Cellina − L4 niche match | 10 | +0.1836 | +0.0958, +0.2714 | yes |
| CPA − L4 niche match | 10 | +0.0696 | -0.0329, +0.1721 | no |
| floor (L6 two-way pseudobulk) − L5 pooled ridge | 10 | +0.0917 | +0.0513, +0.1321 | yes |
| Cellina − L5 pooled ridge | 10 | +0.1013 | +0.0640, +0.1385 | yes |
| CPA − L5 pooled ridge | 10 | -0.0127 | -0.0409, +0.0155 | no |

## F5

| Rung | Role | Units | Mean score |
|---|---|---|---|
| logistic regression on the reference | naive | 20 | 0.6297 |
| InSituType without the spatial prior | naive | 20 | 0.7143 |
| floor (InSituType with the spatial prior) | floor | 20 | 0.7186 |
| scANVI | learned | 20 | 0.7042 |
| Tangram | learned | 20 | 0.6054 |

| Comparison | Units | D | Bounds | Ahead |
|---|---|---|---|---|
| floor (InSituType with the spatial prior) − logistic regression on the reference | 20 | +0.0888 | +0.0752, +0.1024 | yes |
| scANVI − logistic regression on the reference | 20 | +0.0745 | +0.0640, +0.0850 | yes |
| Tangram − logistic regression on the reference | 20 | -0.0244 | -0.0454, -0.0033 | no |
| floor (InSituType with the spatial prior) − InSituType without the spatial prior | 20 | +0.0043 | +0.0032, +0.0054 | yes |
| scANVI − InSituType without the spatial prior | 20 | -0.0100 | -0.0177, -0.0023 | no |
| Tangram − InSituType without the spatial prior | 20 | -0.1089 | -0.1381, -0.0797 | no |

## Summary

| Family | Floor ahead of every naive rung | Learned models ahead of every naive rung |
|---|---|---|
| F1 | yes | scVI: no; Nicheformer: no |
| F2 | yes | CellCharter: no; Novae zero-shot (mouse units): no |
| F3 | yes | LightGBM: no; NCEM: no |
| F4 | yes | Cellina: yes; CPA: no |
| F5 | yes | scANVI: no; Tangram: no |
