Secondary outputs (PROTOCOL §2; no decision weight). D = first - second named; unit means with one-sided 95% bounds.

| comparison | units | mean first | mean second | D | lower | upper |
|---|---|---|---|---|---|---|
| F2 CellCharter - multi-scale floor | 10 | 0.360 | 0.512 | -0.1521 | -0.2688 | -0.0354 |
| F2 Novae - multi-scale floor (mouse) | 4 | 0.324 | 0.465 | -0.1402 | -0.1724 | -0.1080 |
| F2 multi-scale floor - locked floor | 10 | 0.512 | 0.461 | +0.0518 | -0.0379 | +0.1414 |
| F2 CellCharter on PCA - locked floor | 10 | 0.346 | 0.461 | -0.1149 | -0.1970 | -0.0329 |
| F2 CellCharter (scVI) - CellCharter on PCA | 10 | 0.360 | 0.346 | +0.0146 | +0.0002 | +0.0290 |
| F3 LightGBM - v2 floor | 22 | 0.057 | 0.058 | -0.0005 | -0.0034 | +0.0024 |
| F3 NCEM - v2 floor | 22 | 0.037 | 0.058 | -0.0206 | -0.0246 | -0.0166 |
| F3 ED-NCEM - floor | 22 | 0.031 | 0.073 | -0.0412 | -0.0466 | -0.0359 |
| F3 ED-NCEM - NCEM (linear interactions) | 22 | 0.031 | 0.037 | -0.0057 | -0.0076 | -0.0038 |
| F5 scANVI - InSituType without prior (v2 floor) | 20 | 0.704 | 0.714 | -0.0100 | -0.0177 | -0.0023 |
| F5 Tangram - InSituType without prior (v2 floor) | 20 | 0.605 | 0.714 | -0.1089 | -0.1381 | -0.0797 |
| F5 floor - InSituType without prior (spatial prior gain) | 20 | 0.719 | 0.714 | +0.0043 | +0.0032 | +0.0054 |
| F5 scANVI - logistic regression (logreg_C0.1_ref) | 20 | 0.704 | 0.630 | +0.0745 | +0.0640 | +0.0850 |
| F5 Tangram - logistic regression (logreg_C0.1_ref) | 20 | 0.605 | 0.630 | -0.0244 | -0.0454 | -0.0033 |
| F5 logistic regression (logreg_C0.1_ref) - floor | 20 | 0.630 | 0.719 | -0.0888 | -0.1024 | -0.0752 |
| F5 IF-gated macro-F1 (GSE1): scanvi - floor | 6 | 0.605 | 0.657 | -0.0529 | -0.0878 | -0.0180 |
| F5 IF-gated macro-F1 (GSE1): tangram - floor | 6 | 0.260 | 0.657 | -0.3971 | -0.4426 | -0.3516 |
| F5 IF-gated macro-F1 (GSE1): insitutype_no_prior - floor | 6 | 0.654 | 0.657 | -0.0033 | -0.0073 | +0.0007 |
| F5 IF-gated macro-F1 (GSE1): logreg - floor | 6 | 0.528 | 0.657 | -0.1292 | -0.1479 | -0.1105 |

F5 per model, mean over crops. Marker-positive recall: share of the cells in the top 10% of one class's held-out-marker score
(and of no other class's) that the model assigns to that class; defined on 19 / 21 / 19 crops (undefined where another class's 90th percentile is 0). Class share: predicted fraction.

| model | config | AUROC | IF-gated macro-F1 (GSE1) | recall epi / imm / mye | share epi / imm / mye / other |
|---|---|---|---|---|---|
| floor | insitutype_spatial | 0.723 | 0.657 | 0.84 / 0.26 / 0.47 | 0.51 / 0.10 / 0.13 / 0.26 |
| insitutype_no_prior | insitutype | 0.719 | 0.654 | 0.83 / 0.27 / 0.48 | 0.49 / 0.11 / 0.13 / 0.27 |
| logreg | logreg_C0.1_ref | 0.635 | 0.528 | 0.75 / 0.36 / 0.41 | 0.48 / 0.25 / 0.14 / 0.13 |
| scanvi | nl30_nb | 0.709 | 0.605 | 0.83 / 0.45 / 0.47 | 0.47 / 0.24 / 0.11 / 0.18 |
| tangram | uniform | 0.605 | 0.260 | 0.20 / 0.62 / 0.32 | 0.11 / 0.51 / 0.16 / 0.22 |

Novae zero-shot collapse rule (effective rank < 5 or > 50% of cells in one domain; test pair on mouse units only)

| crop | effective rank | largest domain share | collapsed |
|---|---|---|---|
| allen_animal1_C57BL6J-1.086 | 8.53 | 0.21 | False |
| allen_animal2_C57BL6J-2.039 | 8.72 | 0.17 | False |
| allen_animal3_C57BL6J-3.003 | 8.14 | 0.39 | False |
| allen_yao638850_C57BL6J-638850.36 | 7.59 | 0.24 | False |
| gse1_110 | 5.18 | 0.51 | True |
| gse1_120 | 5.20 | 0.49 | False |
| gse1_210 | 5.84 | 0.55 | True |
| gse1_221 | 4.46 | 0.74 | True |
| gse1_222 | 5.20 | 1.00 | True |
| gse1_242 | 5.18 | 0.59 | True |

F2 balanced ARI per crop

| crop | floor | floor_multiscale | cellcharter | cellcharter_pca | novae_zero_shot |
|---|---|---|---|---|---|
| allen_animal1_C57BL6J-1.086 | 0.400 | 0.448 | 0.367 | 0.368 | 0.296 |
| allen_animal2_C57BL6J-2.039 | 0.389 | 0.405 | 0.359 | 0.358 | 0.304 |
| allen_animal3_C57BL6J-3.003 | 0.542 | 0.552 | 0.490 | 0.460 | 0.408 |
| allen_yao638850_C57BL6J-638850.36 | 0.425 | 0.453 | 0.377 | 0.381 | 0.289 |
| gse1_110 | 0.502 | 0.248 | 0.332 | 0.343 | 0.231 |
| gse1_120 | 0.257 | 0.276 | 0.300 | 0.272 | 0.040 |
| gse1_210 | 0.451 | 0.690 | 0.361 | 0.376 | 0.066 |
| gse1_221 | 0.323 | 0.584 | 0.321 | 0.261 | 0.196 |
| gse1_222 | 0.662 | 0.818 | 0.214 | 0.178 | 0.000 |
| gse1_242 | 0.654 | 0.649 | 0.481 | 0.459 | 0.044 |
