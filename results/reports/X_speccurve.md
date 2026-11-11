# X_speccurve: specification curve of the round-1 comparison (exploratory)

D = learned − floor under every choice already computed in the saved reports; one-sided 95% bounds; δ = the pair's locked margin. Rows of one pair share units and are not independent: descriptive only.

On the registered scale: 146 specifications; D < 0 in 131; floor non-inferior in 135; floor ahead (upper bound < 0) in 117; learned model ahead (lower bound > 0) in 3. 2 further specifications use another score scale (no margin).

| Pair | Specifications | Registered D | Min | Median | Max | δ | D < 0 | Non-inferior | Floor ahead | Learned ahead |
|---|---|---|---|---|---|---|---|---|---|---|
| F1|nicheformer | 16 | -0.0379 | -0.0731 | -0.0379 | -0.0117 | 0.0195 | 16 | 16 | 16 | 0 |
| F1|scvi | 16 | -0.0087 | -0.0510 | -0.0087 | -0.0023 | 0.0195 | 16 | 16 | 15 | 0 |
| F2|cellcharter | 11 | -0.1003 | -0.2464 | -0.1003 | -0.0411 | 0.0500 | 11 | 11 | 8 | 0 |
| F2|novae | 7 | -0.1147 | -0.1402 | -0.1296 | -0.1108 | 0.0500 | 7 | 7 | 7 | 0 |
| F3|lightgbm | 17 | -0.0155 | -0.0268 | -0.0155 | +0.0167 | 0.0043 | 16 | 16 | 15 | 1 |
| F3|ncem | 15 | -0.0355 | -0.0874 | -0.0355 | -0.0051 | 0.0043 | 15 | 15 | 15 | 0 |
| F4_pearson|cellina | 8 | +0.0095 | -0.0010 | +0.0102 | +0.0228 | 0.0300 | 1 | 5 | 0 | 2 |
| F4_pearson|cpa | 8 | -0.1045 | -0.1995 | -0.0941 | -0.0348 | 0.0300 | 8 | 7 | 5 | 0 |
| F4_pearson_clean|cellina | 5 | +0.0086 | -0.0012 | +0.0086 | +0.0232 | 0.0300 | 1 | 3 | 0 | 0 |
| F4_pearson_clean|cpa | 5 | -0.0677 | -0.0777 | -0.0677 | -0.0611 | 0.0300 | 5 | 4 | 4 | 0 |
| F5|scanvi | 19 | -0.0143 | -0.0639 | -0.0141 | +0.0051 | 0.0255 | 16 | 17 | 14 | 0 |
| F5|tangram | 19 | -0.1132 | -0.1526 | -0.1117 | -0.0347 | 0.0255 | 19 | 18 | 18 | 0 |

Bounds: report = as in the source report; symmetry = lower bound completed as 2 est − upper; computed here = the registered test on the saved per-unit D.

## F1|nicheformer

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| units | CosMx WTx only | P5_corrected | 9 | -0.0731 | -0.0858, -0.0604 | -3.75 | yes | computed here |
| aggregation | one value per slide | X_levels | 19 | -0.0430 | -0.0494, -0.0366 | -2.21 | yes | report |
| units | without mouse brain | P5_loco_crc | 21 | -0.0429 | -0.0491, -0.0367 | -2.20 | yes | symmetry |
| aggregation | one value per study | X_levels | 8 | -0.0427 | -0.0547, -0.0308 | -2.19 | yes | report |
| tuning | both sides at their best configuration (tier 2) | T2_corrected | 25 | -0.0421 | -0.0471, -0.0371 | -2.16 | yes | report |
| metric | Euclidean for the model, cosine for the floor | X_f1_euclid | 25 | -0.0389 | -0.0443, -0.0336 | -2.00 | yes | report |
| units | without units in the model's pretraining data | P5_corrected | – | -0.0388 | -0.0443, -0.0333 | -1.99 | yes | symmetry |
| registered | registered | P5_corrected | 25 | -0.0379 | -0.0432, -0.0327 | -1.95 | yes | report |
| units | without CRC | P5_loco_crc | 7 | -0.0378 | -0.0487, -0.0269 | -1.94 | yes | symmetry |
| metric | the model's better distance, cosine for the floor | X_f1_euclid | 25 | -0.0378 | -0.0430, -0.0325 | -1.94 | yes | report |
| units | without breast | P5_loco_crc | 24 | -0.0370 | -0.0424, -0.0317 | -1.90 | yes | symmetry |
| units | without skin SCC | P5_loco_crc | 24 | -0.0370 | -0.0424, -0.0316 | -1.90 | yes | symmetry |
| units | without RCC | P5_loco_crc | 24 | -0.0354 | -0.0405, -0.0303 | -1.82 | yes | symmetry |
| metric | Euclidean neighbours on both sides | X_f1_euclid | 25 | -0.0256 | -0.0310, -0.0202 | -1.31 | yes | report |
| units | CosMx 1k only | P5_corrected | 11 | -0.0219 | -0.0286, -0.0153 | -1.13 | yes | computed here |
| units | MERFISH mouse only | P5_corrected | 4 | -0.0117 | -0.0157, -0.0077 | -0.60 | yes | computed here |

## F1|scvi

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| metric | Euclidean for the model, cosine for the floor | X_f1_euclid | 25 | -0.0510 | -0.0566, -0.0454 | -2.62 | yes | report |
| metric | Euclidean neighbours on both sides | X_f1_euclid | 25 | -0.0376 | -0.0429, -0.0323 | -1.93 | yes | report |
| units | CosMx WTx only | P5_corrected | 9 | -0.0136 | -0.0177, -0.0096 | -0.70 | yes | computed here |
| tuning | both sides at their best configuration (tier 2) | T2_corrected | 14 | -0.0122 | -0.0162, -0.0082 | -0.63 | yes | report |
| units | without mouse brain | P5_loco_crc | 21 | -0.0100 | -0.0128, -0.0071 | -0.51 | yes | symmetry |
| aggregation | one value per study | X_levels | 8 | -0.0100 | -0.0127, -0.0072 | -0.51 | yes | report |
| aggregation | one value per slide | X_levels | 19 | -0.0093 | -0.0120, -0.0067 | -0.48 | yes | report |
| registered | registered | P5_corrected | 25 | -0.0087 | -0.0112, -0.0062 | -0.45 | yes | report |
| metric | the model's better distance, cosine for the floor | X_f1_euclid | 25 | -0.0087 | -0.0112, -0.0062 | -0.45 | yes | report |
| units | without CRC | P5_loco_crc | 7 | -0.0087 | -0.0113, -0.0061 | -0.45 | yes | symmetry |
| units | without skin SCC | P5_loco_crc | 24 | -0.0085 | -0.0111, -0.0059 | -0.44 | yes | symmetry |
| units | without breast | P5_loco_crc | 24 | -0.0083 | -0.0109, -0.0057 | -0.43 | yes | symmetry |
| units | without RCC | P5_loco_crc | 24 | -0.0083 | -0.0108, -0.0057 | -0.43 | yes | symmetry |
| scale | whole slides (6 GSE1) | X_wholeslide | 6 | -0.0076 | -0.0123, -0.0030 | -0.39 | yes | report |
| units | CosMx 1k only | P5_corrected | 11 | -0.0072 | -0.0118, -0.0026 | -0.37 | yes | computed here |
| units | MERFISH mouse only | P5_corrected | 4 | -0.0023 | -0.0068, +0.0023 | -0.12 | yes | computed here |

## F2|cellcharter

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| tuning | both sides at their best configuration (tier 2) | T2_corrected | 5 | -0.2464 | -0.4262, -0.0666 | -4.93 | yes | report |
| scale | whole slides (6 GSE1) | X_wholeslide | 6 | -0.1553 | -0.2600, -0.0506 | -3.11 | yes | report |
| floor | multi-scale floor | P5_secondary | 10 | -0.1521 | -0.2688, -0.0354 | -3.04 | yes | report |
| units | without mouse brain = CosMx WTx only | P5_loco_crc, P5_corrected | 6 | -0.1399 | -0.2831, +0.0033 | -2.80 | yes | symmetry |
| model | CellCharter on PCA instead of scVI | P5_secondary | 10 | -0.1149 | -0.1970, -0.0329 | -2.30 | yes | report |
| registered | registered | P5_corrected | 10 | -0.1003 | -0.1814, -0.0193 | -2.01 | yes | report |
| metric | scored 250 µm inside the crop edge | X_interior | 10 | -0.0949 | -0.1802, -0.0095 | -1.90 | yes | computed here |
| metric | scored 100 µm inside the crop edge | X_interior | 10 | -0.0926 | -0.1680, -0.0172 | -1.85 | yes | computed here |
| aggregation | one value per study | X_levels | 3 | -0.0756 | -0.1698, +0.0186 | -1.51 | yes | report |
| metric | scored 500 µm inside the crop edge | X_interior | 10 | -0.0727 | -0.1741, +0.0287 | -1.45 | yes | computed here |
| units | without CRC = MERFISH mouse only | P5_loco_crc, P5_corrected | 4 | -0.0411 | -0.0540, -0.0281 | -0.82 | yes | symmetry |

## F2|novae

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| tuning | both sides at their best configuration (tier 2) | T2_corrected | 4 | -0.1402 | -0.1724, -0.1080 | -2.80 | yes | report |
| floor | multi-scale floor | P5_secondary | 4 | -0.1402 | -0.1724, -0.1080 | -2.80 | yes | report |
| metric | scored 500 µm inside the crop edge | X_interior | 4 | -0.1298 | -0.1922, -0.0675 | -2.60 | yes | computed here |
| metric | scored 250 µm inside the crop edge | X_interior | 4 | -0.1296 | -0.1863, -0.0730 | -2.59 | yes | computed here |
| aggregation | one value per study | X_levels | 2 | -0.1218 | -0.2111, -0.0325 | -2.44 | yes | report |
| registered | registered | P5_corrected | 4 | -0.1147 | -0.1437, -0.0857 | -2.29 | yes | report |
| metric | scored 100 µm inside the crop edge | X_interior | 4 | -0.1108 | -0.1359, -0.0857 | -2.22 | yes | computed here |

## F3|lightgbm

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| units | MERFISH mouse only | P5_corrected | 4 | -0.0268 | -0.0411, -0.0126 | -6.20 | yes | computed here |
| model | early stopping, the floor's features | X_f3_lgbm_es | 22 | -0.0237 | -0.0265, -0.0209 | -5.46 | yes | report |
| units | without CRC | P5_loco_crc | 7 | -0.0217 | -0.0290, -0.0144 | -5.01 | yes | symmetry |
| model | early stopping, all features | X_f3_lgbm_es | 22 | -0.0211 | -0.0234, -0.0188 | -4.87 | yes | report |
| aggregation | one value per study | X_levels | 7 | -0.0188 | -0.0238, -0.0138 | -4.34 | yes | report |
| aggregation | one value per slide | X_levels | 18 | -0.0166 | -0.0196, -0.0136 | -3.83 | yes | report |
| tuning | both sides at their best configuration (tier 2) | T2_corrected | 22 | -0.0166 | -0.0194, -0.0137 | -3.82 | yes | report |
| units | without breast | P5_loco_crc | 21 | -0.0156 | -0.0185, -0.0128 | -3.61 | yes | symmetry |
| units | without skin SCC | P5_loco_crc | 21 | -0.0155 | -0.0184, -0.0126 | -3.58 | yes | symmetry |
| registered | registered | P5_corrected | 22 | -0.0155 | -0.0182, -0.0127 | -3.57 | yes | report |
| units | without RCC | P5_loco_crc | 21 | -0.0153 | -0.0182, -0.0125 | -3.54 | yes | symmetry |
| units | CosMx WTx only | P5_corrected | 9 | -0.0134 | -0.0159, -0.0108 | -3.08 | yes | computed here |
| units | without mouse brain | P5_loco_crc | 18 | -0.0129 | -0.0155, -0.0104 | -2.99 | yes | symmetry |
| units | CosMx 1k only | P5_corrected | 9 | -0.0125 | -0.0173, -0.0077 | -2.89 | yes | computed here |
| scale | whole slides (6 GSE1) | X_wholeslide | 6 | -0.0071 | -0.0132, -0.0009 | -1.64 | yes | report |
| floor | floor without ring features (v2) | P5_secondary | 22 | -0.0005 | -0.0034, +0.0024 | -0.11 | yes | report |
| scale | whole slides, floor without ring features (v2) | X_wholeslide | 6 | +0.0167 | +0.0036, +0.0297 | +3.84 | no | report |

## F3|ncem

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| units | MERFISH mouse only | P5_corrected | 4 | -0.0874 | -0.1144, -0.0604 | -20.17 | yes | computed here |
| units | without CRC | P5_loco_crc | 7 | -0.0606 | -0.0752, -0.0460 | -13.99 | yes | symmetry |
| aggregation | one value per study | X_levels | 7 | -0.0446 | -0.0566, -0.0326 | -10.30 | yes | report |
| model | encoder-decoder NCEM | P5_secondary | 22 | -0.0412 | -0.0466, -0.0359 | -9.52 | yes | report |
| aggregation | one value per slide | X_levels | 18 | -0.0376 | -0.0439, -0.0314 | -8.69 | yes | report |
| units | without breast | P5_loco_crc | 21 | -0.0365 | -0.0425, -0.0305 | -8.43 | yes | symmetry |
| units | without RCC | P5_loco_crc | 21 | -0.0362 | -0.0422, -0.0301 | -8.35 | yes | symmetry |
| registered | registered | P5_corrected | 22 | -0.0355 | -0.0413, -0.0298 | -8.21 | yes | report |
| units | without skin SCC | P5_loco_crc | 21 | -0.0355 | -0.0413, -0.0297 | -8.19 | yes | symmetry |
| tuning | both sides at their best configuration (tier 2) | T2_corrected | 11 | -0.0347 | -0.0429, -0.0264 | -8.00 | yes | report |
| units | CosMx 1k only | P5_corrected | 9 | -0.0279 | -0.0376, -0.0182 | -6.44 | yes | computed here |
| units | without mouse brain | P5_loco_crc | 18 | -0.0240 | -0.0297, -0.0184 | -5.55 | yes | symmetry |
| floor | floor without ring features (v2) | P5_secondary | 22 | -0.0206 | -0.0246, -0.0166 | -4.75 | yes | report |
| units | CosMx WTx only | P5_corrected | 9 | -0.0202 | -0.0272, -0.0132 | -4.66 | yes | computed here |
| floor | ridge on NCEM's own inputs | P5_corrected | 22 | -0.0051 | -0.0066, -0.0037 | -1.19 | yes | report |

## F4_pearson|cellina

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| units | without CRC = MERFISH mouse only = mouse (Allen) only | P5_loco_crc, P5_corrected, X_f4_split | 4 | -0.0010 | -0.0051, +0.0031 | -0.03 | yes | symmetry |
| aggregation | one value per study | X_levels | 3 | +0.0040 | -0.0150, +0.0229 | +0.13 | yes | report |
| units | prospective only | X_f4_split | 5 | +0.0081 | -0.0110, +0.0272 | +0.27 | yes | report |
| registered | registered | P5_corrected | 10 | +0.0095 | +0.0011, +0.0180 | +0.32 | yes | report |
| units | retrospective only | X_f4_split | 5 | +0.0109 | -0.0024, +0.0242 | +0.36 | yes | report |
| units | CRC retrospective only | X_f4_split | 4 | +0.0134 | -0.0040, +0.0308 | +0.45 | no | report |
| units | without mouse brain = CosMx WTx only = CRC (GSE1) only | P5_loco_crc, P5_corrected, X_f4_split | 6 | +0.0165 | +0.0017, +0.0314 | +0.55 | no | symmetry |
| units | CRC prospective only | X_f4_split | 2 | +0.0228 | -0.1075, +0.1532 | +0.76 | no | report |

## F4_pearson|cpa

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| units | CRC retrospective only | X_f4_split | 4 | -0.1995 | -0.2777, -0.1213 | -6.65 | yes | report |
| units | retrospective only | X_f4_split | 5 | -0.1636 | -0.2578, -0.0693 | -5.45 | yes | report |
| units | without mouse brain = CosMx WTx only = CRC (GSE1) only | P5_loco_crc, P5_corrected, X_f4_split | 6 | -0.1509 | -0.2293, -0.0725 | -5.03 | yes | symmetry |
| registered | registered | P5_corrected | 10 | -0.1045 | -0.1501, -0.0588 | -3.48 | yes | report |
| aggregation | one value per study | X_levels | 3 | -0.0837 | -0.1945, +0.0271 | -2.79 | yes | report |
| units | CRC prospective only | X_f4_split | 2 | -0.0536 | -0.3281, +0.2208 | -1.79 | no | report |
| units | prospective only | X_f4_split | 5 | -0.0454 | -0.0839, -0.0068 | -1.51 | yes | report |
| units | without CRC = MERFISH mouse only = mouse (Allen) only | P5_loco_crc, P5_corrected, X_f4_split | 4 | -0.0348 | -0.0709, +0.0013 | -1.16 | yes | symmetry |

## F4_pearson_clean|cellina

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| units | without CRC = MERFISH mouse only = mouse (Allen) only | P5_loco_crc, P5_corrected, X_f4_split | 3 | -0.0012 | -0.0066, +0.0042 | -0.04 | yes | symmetry |
| aggregation | one value per study | X_levels | 3 | +0.0066 | -0.0179, +0.0311 | +0.22 | no | report |
| registered | registered | P5_corrected | 5 | +0.0086 | -0.0096, +0.0268 | +0.29 | yes | report |
| units | prospective only | X_f4_split | 5 | +0.0086 | -0.0106, +0.0277 | +0.29 | yes | report |
| units | without mouse brain = CosMx WTx only = CRC (GSE1) only = CRC prospective only | P5_loco_crc, P5_corrected, X_f4_split | 2 | +0.0232 | -0.1089, +0.1554 | +0.77 | no | symmetry |

## F4_pearson_clean|cpa

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| units | without mouse brain = CosMx WTx only = CRC (GSE1) only = CRC prospective only | P5_loco_crc, P5_corrected, X_f4_split | 2 | -0.0777 | -0.3337, +0.1783 | -2.59 | no | symmetry |
| aggregation | one value per study | X_levels | 3 | -0.0699 | -0.0974, -0.0425 | -2.33 | yes | report |
| registered | registered | P5_corrected | 5 | -0.0677 | -0.1102, -0.0253 | -2.26 | yes | report |
| units | prospective only | X_f4_split | 5 | -0.0677 | -0.1021, -0.0334 | -2.26 | yes | report |
| units | without CRC = MERFISH mouse only = mouse (Allen) only | P5_loco_crc, P5_corrected, X_f4_split | 3 | -0.0611 | -0.1086, -0.0136 | -2.04 | yes | symmetry |

## F5|scanvi

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| scale | whole slides (6 GSE1) | X_wholeslide | 6 | -0.0639 | -0.1069, -0.0208 | -2.51 | yes | report |
| metric | floor thresholded at its own call rate | X_f5_matched | 20 | -0.0326 | -0.0423, -0.0228 | -1.28 | yes | report |
| units | CosMx 1k only | P5_corrected | 11 | -0.0258 | -0.0371, -0.0145 | -1.01 | yes | computed here |
| units | without RCC | P5_loco_crc | 19 | -0.0169 | -0.0248, -0.0090 | -0.66 | yes | symmetry |
| tuning | both sides at their best configuration (tier 2) | T2_corrected | 11 | -0.0166 | -0.0297, -0.0036 | -0.65 | yes | report |
| metric | scored 500 µm inside the crop edge | X_interior | 20 | -0.0152 | -0.0238, -0.0065 | -0.60 | yes | computed here |
| units | without skin SCC | P5_loco_crc | 19 | -0.0149 | -0.0235, -0.0062 | -0.58 | yes | symmetry |
| metric | scored 100 µm inside the crop edge | X_interior | 20 | -0.0143 | -0.0226, -0.0061 | -0.56 | yes | computed here |
| registered | registered | P5_corrected | 20 | -0.0143 | -0.0225, -0.0062 | -0.56 | yes | report |
| metric | scored 250 µm inside the crop edge | X_interior | 20 | -0.0141 | -0.0225, -0.0057 | -0.55 | yes | computed here |
| units | without breast | P5_loco_crc | 19 | -0.0138 | -0.0221, -0.0055 | -0.54 | yes | symmetry |
| floor | floor without the spatial prior | P5_secondary | 20 | -0.0100 | -0.0177, -0.0023 | -0.39 | yes | report |
| aggregation | one value per slide | X_levels | 14 | -0.0089 | -0.0176, -0.0003 | -0.35 | yes | report |
| aggregation | one value per study | X_levels | 5 | -0.0042 | -0.0298, +0.0215 | -0.16 | yes | report |
| metric | floor at the model's call rate | X_f5_matched | 20 | -0.0040 | -0.0075, -0.0005 | -0.16 | yes | report |
| units | CosMx WTx only | P5_corrected | 9 | -0.0003 | -0.0136, +0.0130 | -0.01 | yes | computed here |
| metric | AUROC within depth quintiles | X_f5_depth | 20 | +0.0015 | -0.0044, +0.0074 | +0.06 | yes | report |
| units | without CRC | P5_loco_crc | 3 | +0.0022 | -0.0481, +0.0525 | +0.09 | no | symmetry |
| scale | whole slides, floor at the model's call rate | X_wholeslide | 6 | +0.0051 | -0.0203, +0.0306 | +0.20 | no | report |
| metric | protein-gated macro-F1 (GSE1) | P5_secondary | 6 | -0.0529 | -0.0878, -0.0180 | other scale | – | report |

## F5|tangram

| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |
|---|---|---|---|---|---|---|---|---|
| metric | AUROC within depth quintiles | X_f5_depth | 20 | -0.1526 | -0.1798, -0.1253 | -5.99 | yes | report |
| metric | floor thresholded at its own call rate | X_f5_matched | 20 | -0.1314 | -0.1634, -0.0995 | -5.16 | yes | report |
| units | CosMx 1k only | P5_corrected | 11 | -0.1254 | -0.1740, -0.0768 | -4.93 | yes | computed here |
| tuning | both sides at their best configuration (tier 2) | T2_corrected | 11 | -0.1183 | -0.1575, -0.0791 | -4.65 | yes | report |
| units | without RCC | P5_loco_crc | 19 | -0.1181 | -0.1487, -0.0875 | -4.64 | yes | symmetry |
| units | without CRC | P5_loco_crc | 3 | -0.1158 | -0.2615, +0.0298 | -4.55 | no | symmetry |
| registered | registered | P5_corrected | 20 | -0.1132 | -0.1431, -0.0833 | -4.45 | yes | report |
| aggregation | one value per study | X_levels | 5 | -0.1125 | -0.1880, -0.0369 | -4.42 | yes | report |
| metric | scored 100 µm inside the crop edge | X_interior | 20 | -0.1118 | -0.1420, -0.0815 | -4.39 | yes | computed here |
| units | without breast | P5_loco_crc | 19 | -0.1117 | -0.1430, -0.0803 | -4.39 | yes | symmetry |
| metric | scored 500 µm inside the crop edge | X_interior | 20 | -0.1115 | -0.1425, -0.0805 | -4.38 | yes | computed here |
| metric | scored 250 µm inside the crop edge | X_interior | 20 | -0.1109 | -0.1417, -0.0800 | -4.36 | yes | computed here |
| units | without skin SCC | P5_loco_crc | 19 | -0.1094 | -0.1397, -0.0790 | -4.30 | yes | symmetry |
| model | input as log counts per 10k | X_f5_tangram_norm | 20 | -0.1094 | -0.1391, -0.0796 | -4.30 | yes | report |
| floor | floor without the spatial prior | P5_secondary | 20 | -0.1089 | -0.1381, -0.0797 | -4.28 | yes | report |
| aggregation | one value per slide | X_levels | 14 | -0.1078 | -0.1349, -0.0806 | -4.23 | yes | report |
| model | input as counts per 10k | X_f5_tangram_norm | 20 | -0.1044 | -0.1290, -0.0798 | -4.10 | yes | report |
| units | CosMx WTx only | P5_corrected | 9 | -0.0983 | -0.1349, -0.0617 | -3.86 | yes | computed here |
| metric | floor at the model's call rate | X_f5_matched | 20 | -0.0347 | -0.0479, -0.0215 | -1.36 | yes | report |
| metric | protein-gated macro-F1 (GSE1) | P5_secondary | 6 | -0.3971 | -0.4426, -0.3516 | other scale | – | report |
