# X_scale: floor vs learned across cell counts (exploratory)

commit 7385b73fd5275a30738dcd48e72e0f9c13da6df8

Windows by the registered crop rule at each cell count, gene panel and count split of the whole slide; D = learned - floor.

## F1|scvi

| cells | slides | mean D |
|---|---|---|
| 12500 | 6 | -0.0134 |
| 25000 | 6 | -0.0123 |
| 50000 | 6 | -0.0125 |
| 100000 | 6 | -0.0129 |
| 200000 | 5 | -0.0098 |
| 400000 | 4 | -0.0089 |
| whole | 6 | -0.0076 |

Slope of D per log10(cells): mean +0.0036, 95% CI +0.0001, +0.0072 (6 slides)

## F2|cellcharter

| cells | slides | mean D |
|---|---|---|
| 12500 | 6 | -0.0670 |
| 25000 | 6 | -0.1583 |
| 50000 | 6 | -0.1464 |
| 100000 | 6 | -0.1451 |
| 200000 | 5 | -0.1438 |
| 400000 | 4 | -0.0931 |
| whole | 6 | -0.1553 |

Slope of D per log10(cells): mean -0.0541, 95% CI -0.1163, +0.0081 (6 slides)

## F3|lightgbm

| cells | slides | mean D |
|---|---|---|
| 12500 | 6 | -0.0123 |
| 25000 | 6 | -0.0086 |
| 50000 | 6 | -0.0093 |
| 100000 | 6 | -0.0070 |
| 200000 | 5 | -0.0076 |
| 400000 | 4 | -0.0060 |
| whole | 6 | -0.0071 |

Slope of D per log10(cells): mean +0.0029, 95% CI -0.0041, +0.0099 (6 slides)

## F5|scanvi

| cells | slides | mean D |
|---|---|---|
| 12500 | 6 | -0.0093 |
| 25000 | 6 | -0.0083 |
| 50000 | 6 | -0.0012 |
| 100000 | 6 | -0.0220 |
| 200000 | 5 | -0.0591 |
| 400000 | 4 | -0.0876 |
| whole | 6 | -0.0639 |

Slope of D per log10(cells): mean -0.0415, 95% CI -0.0723, -0.0108 (6 slides)
