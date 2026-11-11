# X_meta: common scale and random-effects meta-analysis over studies (exploratory)

All differences as % of the floor's above-chance score A (D = learned − floor; negative = floor ahead). δ = the locked margin on the same scale.

| Pair | Units | Floor | Chance | A | D (registered) | D % of A | One-sided 95% bounds, % | δ, % |
|---|---|---|---|---|---|---|---|---|
| F1|nicheformer | 25 | 0.163 | -0.078 | 0.241 | -0.0379 | -15.7 | -17.9, -13.6 | 8.1 |
| F1|scvi | 25 | 0.163 | -0.078 | 0.241 | -0.0087 | -3.6 | -4.6, -2.6 | 8.1 |
| F2|cellcharter | 10 | 0.461 | +0.000 | 0.461 | -0.1003 | -21.8 | -39.4, -4.2 | 10.9 |
| F2|novae | 4 | 0.439 | +0.000 | 0.439 | -0.1147 | -26.1 | -32.7, -19.5 | 11.4 |
| F3|lightgbm | 22 | 0.073 | +0.000 | 0.073 | -0.0155 | -21.3 | -25.1, -17.5 | 6.0 |
| F3|ncem | 22 | 0.073 | +0.000 | 0.073 | -0.0355 | -48.9 | -56.8, -41.1 | 6.0 |
| F4_pearson|cellina | 10 | 0.873 | +0.000 | 0.873 | +0.0095 | +1.1 | +0.1, +2.1 | 3.4 |
| F4_pearson|cpa | 10 | 0.873 | +0.000 | 0.873 | -0.1045 | -12.0 | -17.2, -6.7 | 3.4 |
| F4_pearson_clean|cellina | 5 | 0.887 | +0.000 | 0.887 | +0.0086 | +1.0 | -1.1, +3.0 | 3.4 |
| F4_pearson_clean|cpa | 5 | 0.887 | +0.000 | 0.887 | -0.0677 | -7.6 | -12.4, -2.9 | 3.4 |
| F5|scanvi | 20 | 0.719 | +0.500 | 0.219 | -0.0143 | -6.6 | -10.3, -2.8 | 11.6 |
| F5|tangram | 20 | 0.719 | +0.500 | 0.219 | -0.1132 | -51.8 | -65.5, -38.1 | 11.6 |

Random-effects model over studies (REML, one within-study variance). μ = mean difference over studies; τ = between-study SD; PI = 95% prediction interval for the true difference in a new study (three studies needed). Check = unweighted study means, mean and 95% prediction interval for a new study's mean. All in % of A.

| Pair | Studies | Units | μ | One-sided 95% bounds | τ | 95% PI | Non-inferior (μ) | PI below δ | Check: mean | Check: PI | Check: PI below δ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| F1|nicheformer | 8 | 25 | -17.6 | -27.1, -8.1 | 13.3 | -52.3, +17.1 | yes | no | -17.7 | -53.9, +18.4 | no |
| F1|scvi | 8 | 25 | -3.6 | -5.0, -2.3 | 0.4 | -5.7, -1.6 | yes | yes | -4.1 | -11.8, +3.5 | yes |
| F2|cellcharter | 3 | 10 | -21.8 | -49.8, +6.2 | 0.0 | -143.5, +100.0 | yes | no | -16.4 | -76.7, +43.8 | no |
| F2|novae | 2 | 4 | -26.1 | -43.8, -8.4 | 0.0 | – | yes | – | -27.7 | -98.6, +43.2 | no |
| F3|lightgbm | 7 | 22 | -24.4 | -32.8, -16.0 | 8.6 | -49.1, +0.3 | yes | yes | -25.9 | -60.9, +9.2 | no |
| F3|ncem | 7 | 22 | -60.9 | -97.6, -24.2 | 47.2 | -191.6, +69.8 | yes | no | -61.4 | -195.2, +72.4 | no |
| F4_pearson|cellina | 3 | 10 | +0.8 | -1.5, +3.2 | 0.9 | -14.6, +16.3 | yes | no | +0.5 | -5.9, +6.8 | no |
| F4_pearson|cpa | 3 | 10 | -10.2 | -25.5, +5.1 | 7.1 | -121.8, +101.3 | no | no | -9.6 | -47.0, +27.8 | no |
| F4_pearson_clean|cellina | 3 | 5 | +1.0 | -2.0, +3.9 | 0.2 | -12.2, +14.2 | no | no | +0.7 | -7.4, +8.9 | no |
| F4_pearson_clean|cpa | 3 | 5 | -7.6 | -12.9, -2.3 | 0.0 | -30.7, +15.5 | yes | no | -7.9 | -17.0, +1.2 | yes |
| F5|scanvi | 5 | 20 | -3.4 | -12.9, +6.2 | 7.6 | -31.4, +24.7 | yes | no | -1.9 | -35.8, +32.0 | no |
| F5|tangram | 5 | 20 | -51.8 | -68.4, -35.1 | 0.0 | -76.6, -26.9 | yes | yes | -51.5 | -138.5, +35.6 | no |

Study means of D (units):

- F1|nicheformer: Allen Yao -0.0163 (1); Allen Zhuang -0.0101 (3); Bruker breast multiomic -0.0589 (1); Bruker kidney -0.0983 (1); Bruker skin -0.0608 (1); GSE1 CRC (Crowell et al.) -0.0733 (6); GSE303070 -0.0219 (11); SPATCH -0.0023 (1)
- F1|scvi: Allen Yao +0.0021 (1); Allen Zhuang -0.0037 (3); Bruker breast multiomic -0.0185 (1); Bruker kidney -0.0192 (1); Bruker skin -0.0141 (1); GSE1 CRC (Crowell et al.) -0.0118 (6); GSE303070 -0.0072 (11); SPATCH -0.0071 (1)
- F2|cellcharter: Allen Yao -0.0483 (1); Allen Zhuang -0.0387 (3); GSE1 CRC (Crowell et al.) -0.1399 (6)
- F2|novae: Allen Yao -0.1359 (1); Allen Zhuang -0.1076 (3)
- F3|lightgbm: Allen Yao -0.0390 (1); Allen Zhuang -0.0228 (3); Bruker breast multiomic -0.0117 (1); Bruker kidney -0.0182 (1); Bruker skin -0.0148 (1); GSE1 CRC (Crowell et al.) -0.0126 (6); GSE303070 -0.0125 (9)
- F3|ncem: Allen Yao -0.1131 (1); Allen Zhuang -0.0788 (3); Bruker breast multiomic -0.0149 (1); Bruker kidney -0.0225 (1); Bruker skin -0.0374 (1); GSE1 CRC (Crowell et al.) -0.0178 (6); GSE303070 -0.0279 (9)
- F4_pearson|cellina: Allen Yao -0.0050 (1); Allen Zhuang +0.0004 (3); GSE1 CRC (Crowell et al.) +0.0165 (6)
- F4_pearson|cpa: Allen Yao -0.0807 (1); Allen Zhuang -0.0195 (3); GSE1 CRC (Crowell et al.) -0.1509 (6)
- F4_pearson_clean|cellina: Allen Yao -0.0034 (1); Allen Zhuang -0.0001 (2); GSE1 CRC (Crowell et al.) +0.0232 (2)
- F4_pearson_clean|cpa: Allen Yao -0.0809 (1); Allen Zhuang -0.0512 (2); GSE1 CRC (Crowell et al.) -0.0777 (2)
- F5|scanvi: Bruker breast multiomic -0.0242 (1); Bruker kidney +0.0345 (1); Bruker skin -0.0037 (1); GSE1 CRC (Crowell et al.) -0.0015 (6); GSE303070 -0.0258 (11)
- F5|tangram: Bruker breast multiomic -0.1425 (1); Bruker kidney -0.0192 (1); Bruker skin -0.1858 (1); GSE1 CRC (Crowell et al.) -0.0895 (6); GSE303070 -0.1254 (11)
