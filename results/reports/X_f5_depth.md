# F5 with depth held fixed (exploratory)

Depth-stratified = held-out-marker AUROC within each library-size quintile, mean over quintiles. Registered scores recomputed exactly on every crop.

| Pair | Score | Units | Floor | Learned | D (est) | one-sided 95% bounds | Non-inferior | Floor superior |
|---|---|---|---|---|---|---|---|---|
| F5 scanvi | registered | 20 | 0.719 | 0.704 | -0.0143 | -0.0225, -0.0062 | yes | yes |
| F5 scanvi | depth-stratified | 20 | 0.725 | 0.727 | +0.0015 | -0.0044, +0.0074 | yes | no |
| F5 tangram | registered | 20 | 0.719 | 0.605 | -0.1132 | -0.1431, -0.0833 | yes | yes |
| F5 tangram | depth-stratified | 20 | 0.725 | 0.573 | -0.1526 | -0.1798, -0.1253 | yes | yes |

Share of each class's calls in the deepest library-size quintile (mean over crops; 0.20 = depth-neutral):

| Model | epithelial | immune_nonmyeloid | myeloid |
|---|---|---|---|
| floors | 0.31 | 0.08 | 0.10 |
| scanvi | 0.34 | 0.08 | 0.10 |
| tangram | 0.70 | 0.09 | 0.21 |
