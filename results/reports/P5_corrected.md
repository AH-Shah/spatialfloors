| pair | units | floor | learned | D (est) | upper | margin | non-inferior | equivalent | floor superior |
|---|---|---|---|---|---|---|---|---|---|
| F1|nicheformer | 25 | 0.163 | 0.126 | -0.0379 | -0.0327 | 0.0195 | True | False | True |
| F1|scvi | 25 | 0.163 | 0.155 | -0.0087 | -0.0062 | 0.0195 | True | True | True |
| F2|cellcharter | 10 | 0.461 | 0.360 | -0.1003 | -0.0193 | 0.0500 | True | False | True |
| F2|novae | 4 | 0.439 | 0.324 | -0.1147 | -0.0857 | 0.0500 | True | False | True |
| F3|lightgbm | 22 | 0.073 | 0.057 | -0.0155 | -0.0127 | 0.0043 | True | False | True |
| F3|ncem | 22 | 0.073 | 0.037 | -0.0355 | -0.0298 | 0.0043 | True | False | True |
| F4_pearson|cellina | 10 | 0.873 | 0.882 | +0.0095 | +0.0180 | 0.0300 | True | True | False |
| F4_pearson|cpa | 10 | 0.873 | 0.768 | -0.1045 | -0.0588 | 0.0300 | True | False | True |
| F4_pearson_clean|cellina | 5 | 0.887 | 0.896 | +0.0086 | +0.0268 | 0.0300 | True | True | False |
| F4_pearson_clean|cpa | 5 | 0.887 | 0.820 | -0.0677 | -0.0253 | 0.0300 | True | False | True |
| F5|scanvi | 20 | 0.719 | 0.704 | -0.0143 | -0.0062 | 0.0255 | True | True | True |
| F5|tangram | 20 | 0.719 | 0.605 | -0.1132 | -0.0833 | 0.0255 | True | False | True |

floor non-inferior on every pair: True

Secondary (no decision weight): F3 decomposition, D = first - second named

| comparison | units | D (est) | lower | upper |
|---|---|---|---|---|
| F3|ncem_vs_input_ridge | 22 | -0.0051 | -0.0066 | -0.0037 |
| F3|ncem_vs_v2_floor | 22 | -0.0206 | -0.0246 | -0.0166 |
| F3|v3_floor_vs_v2_floor | 22 | +0.0150 | +0.0120 | +0.0179 |
