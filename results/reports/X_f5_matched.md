# X_f5_matched: F5 AUROC at the comparator's call rate (exploratory, after results-v1; no decision weight)

commit 0714d4f3e838ef7ef7cc2cecf4f14294e1fb8965; margin 0.025452985695741493. Floor top-n by class posterior, n = comparator's calls of that class.

| Comparison | Classes | Units | D | One-sided 95% bounds | Mean first | Mean second |
|---|---|---|---|---|---|---|
| scanvi − auroc_floor | mean | 20 | -0.0143 | -0.0225, -0.0062 | 0.709 | 0.723 |
| scanvi − auroc_floor_ownrate | mean | 20 | -0.0326 | -0.0423, -0.0228 | 0.709 | 0.741 |
| scanvi − auroc_floor_matched | mean | 20 | -0.0040 | -0.0075, -0.0005 | 0.709 | 0.712 |
| scanvi − auroc_floor | epithelial | 20 | -0.0013 | -0.0170, +0.0143 | 0.829 | 0.831 |
| scanvi − auroc_floor_ownrate | epithelial | 20 | -0.0125 | -0.0292, +0.0042 | 0.829 | 0.841 |
| scanvi − auroc_floor_matched | epithelial | 20 | -0.0115 | -0.0185, -0.0044 | 0.829 | 0.839 |
| scanvi − auroc_floor | immune_nonmyeloid | 20 | -0.0487 | -0.0690, -0.0285 | 0.579 | 0.626 |
| scanvi − auroc_floor_ownrate | immune_nonmyeloid | 20 | -0.0682 | -0.0933, -0.0432 | 0.579 | 0.647 |
| scanvi − auroc_floor_matched | immune_nonmyeloid | 20 | +0.0249 | +0.0190, +0.0309 | 0.579 | 0.554 |
| scanvi − auroc_floor | myeloid | 20 | +0.0071 | -0.0082, +0.0224 | 0.718 | 0.712 |
| scanvi − auroc_floor_ownrate | myeloid | 20 | -0.0170 | -0.0347, +0.0006 | 0.718 | 0.736 |
| scanvi − auroc_floor_matched | myeloid | 20 | -0.0254 | -0.0319, -0.0188 | 0.718 | 0.743 |
| tangram − auroc_floor | mean | 20 | -0.1132 | -0.1431, -0.0833 | 0.605 | 0.723 |
| tangram − auroc_floor_ownrate | mean | 20 | -0.1314 | -0.1634, -0.0995 | 0.605 | 0.741 |
| tangram − auroc_floor_matched | mean | 20 | -0.0347 | -0.0479, -0.0215 | 0.605 | 0.640 |
| tangram − auroc_floor | epithelial | 20 | -0.0916 | -0.1332, -0.0500 | 0.735 | 0.831 |
| tangram − auroc_floor_ownrate | epithelial | 20 | -0.1028 | -0.1421, -0.0635 | 0.735 | 0.841 |
| tangram − auroc_floor_matched | epithelial | 20 | -0.0433 | -0.0605, -0.0261 | 0.735 | 0.776 |
| tangram − auroc_floor | immune_nonmyeloid | 20 | -0.1279 | -0.1512, -0.1045 | 0.493 | 0.626 |
| tangram − auroc_floor_ownrate | immune_nonmyeloid | 20 | -0.1474 | -0.1759, -0.1188 | 0.493 | 0.647 |
| tangram − auroc_floor_matched | immune_nonmyeloid | 20 | -0.0082 | -0.0189, +0.0025 | 0.493 | 0.501 |
| tangram − auroc_floor | myeloid | 20 | -0.1201 | -0.1548, -0.0854 | 0.587 | 0.712 |
| tangram − auroc_floor_ownrate | myeloid | 20 | -0.1442 | -0.1823, -0.1061 | 0.587 | 0.736 |
| tangram − auroc_floor_matched | myeloid | 20 | -0.0526 | -0.0862, -0.0189 | 0.587 | 0.643 |

Calls summed over crops (comparator / floor): scanvi epithelial 385561 / 437554; scanvi immune_nonmyeloid 227755 / 86109; scanvi myeloid 99378 / 111379; tangram epithelial 93619 / 437554; tangram immune_nonmyeloid 438237 / 86109; tangram myeloid 141509 / 111379

Strata means of D (mean over classes): scanvi - auroc_floor [mean]: CosMx 1k -0.0258, CosMx WTx -0.0003; scanvi - auroc_floor_ownrate [mean]: CosMx 1k -0.0441, CosMx WTx -0.0185; scanvi - auroc_floor_matched [mean]: CosMx 1k +0.0018, CosMx WTx -0.0111; tangram - auroc_floor [mean]: CosMx 1k -0.1254, CosMx WTx -0.0983; tangram - auroc_floor_ownrate [mean]: CosMx 1k -0.1436, CosMx WTx -0.1165; tangram - auroc_floor_matched [mean]: CosMx 1k -0.0360, CosMx WTx -0.0331
