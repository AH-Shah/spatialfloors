# X_wholeslide: the six GSE1 confirmation slides as whole slides (exploratory; results/WHOLESLIDE_PLAN.md)

commit 672137e8d572abe79591b63012bfc0bbdfc5fd72-dirty

| Pair | Slides | D whole | One-sided 95% bounds | δ | Non-inferior | Learned superior (Holm) | D crop (same slides) | Δ = whole − crop (90% CI) | Floor lead change |
|---|---|---|---|---|---|---|---|---|---|
| F1|scvi | 6 | -0.0076 | -0.0123, -0.0030 | 0.0195 | True | False | -0.0118 | +0.0042 (+0.0001, +0.0083) | +36% |
| F2|cellcharter | 6 | -0.1553 | -0.2600, -0.0506 | 0.0500 | True | False | -0.1399 | -0.0154 (-0.0955, +0.0646) | -11% |
| F3|lightgbm | 6 | -0.0071 | -0.0132, -0.0009 | 0.0043 | True | False | -0.0126 | +0.0055 (-0.0009, +0.0120) | +44% |
| F5|scanvi | 6 | -0.0639 | -0.1069, -0.0208 | 0.0255 | True | False | -0.0015 | -0.0624 (-0.0980, -0.0267) | -4160% |
| F3 LightGBM - v2 floor | 6 | +0.0167 | +0.0036, +0.0297 | 0.0043 | – | – | – | – | – |
| F5 scANVI - floor at scANVI's call rate | 6 | +0.0051 | -0.0203, +0.0306 | 0.0255 | – | – | – | – | – |

Per slide (learned, floor): 110: F1 0.080/0.079, F2 0.421/0.540, F3 0.063/0.079, F3 LightGBM - v2 floor 0.063/0.044, F5 0.659/0.686, F5 scANVI - floor at scANVI's call rate 0.659/0.618; 120: F1 0.161/0.170, F2 0.324/0.340, F3 0.123/0.118, F3 LightGBM - v2 floor 0.123/0.082, F5 0.734/0.790, F5 scANVI - floor at scANVI's call rate 0.734/0.713; 210: F1 0.187/0.199, F2 0.235/0.283, F3 0.068/0.070, F3 LightGBM - v2 floor 0.068/0.054, F5 0.593/0.702, F5 scANVI - floor at scANVI's call rate 0.593/0.626; 221: F1 0.138/0.140, F2 0.207/0.341, F3 0.044/0.051, F3 LightGBM - v2 floor 0.044/0.035, F5 0.588/0.716, F5 scANVI - floor at scANVI's call rate 0.588/0.605; 222: F1 0.099/0.111, F2 0.459/0.795, F3 0.023/0.033, F3 LightGBM - v2 floor 0.023/0.030, F5 0.628/0.615, F5 scANVI - floor at scANVI's call rate 0.628/0.644; 242: F1 0.153/0.165, F2 0.474/0.754, F3 0.081/0.093, F3 LightGBM - v2 floor 0.081/0.057, F5 0.627/0.703, F5 scANVI - floor at scANVI's call rate 0.627/0.592
