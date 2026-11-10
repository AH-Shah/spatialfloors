# Leave-one-cancer-out with SPATCH Xenium 5K counted as CRC

Source: `P5_corrected.json` (sha256 6d2b2b18339a), per-unit D of the registered tier-1 report. Registered labels reproduce every registered value; relabelled: `spatch_xenium5k` (SPATCH, tissue unverified -> CRC).

| Pair | Left out | Registered: units, est, upper, NI | Corrected: units, est, upper, NI |
|---|---|---|---|
| F1 nicheformer | CRC | 8, -0.0333, -0.0232, yes | 7, -0.0378, -0.0269, yes |
| F1 nicheformer | RCC | 24, -0.0354, -0.0303, yes | 24, -0.0354, -0.0303, yes |
| F1 nicheformer | SPATCH (tissue unverified) | 24, -0.0394, -0.0340, yes | - |
| F1 nicheformer | breast | 24, -0.0370, -0.0317, yes | 24, -0.0370, -0.0317, yes |
| F1 nicheformer | mouse brain | 21, -0.0429, -0.0367, yes | 21, -0.0429, -0.0367, yes |
| F1 nicheformer | skin SCC | 24, -0.0370, -0.0316, yes | 24, -0.0370, -0.0316, yes |
| F1 scvi | CRC | 8, -0.0085, -0.0060, yes | 7, -0.0087, -0.0061, yes |
| F1 scvi | RCC | 24, -0.0083, -0.0057, yes | 24, -0.0083, -0.0057, yes |
| F1 scvi | SPATCH (tissue unverified) | 24, -0.0088, -0.0063, yes | - |
| F1 scvi | breast | 24, -0.0083, -0.0057, yes | 24, -0.0083, -0.0057, yes |
| F1 scvi | mouse brain | 21, -0.0100, -0.0071, yes | 21, -0.0100, -0.0071, yes |
| F1 scvi | skin SCC | 24, -0.0085, -0.0059, yes | 24, -0.0085, -0.0059, yes |

Pairs whose splits do not change (SPATCH not among their units): F2 cellcharter, F3 lightgbm, F3 ncem, F4_pearson cellina, F4_pearson cpa, F4_pearson_clean cellina, F4_pearson_clean cpa, F5 scanvi, F5 tangram.

| Pair | Non-inferior in every split: registered | corrected |
|---|---|---|
| F1 nicheformer | yes | yes |
| F1 scvi | yes | yes |
| F2 cellcharter | yes | yes |
| F3 lightgbm | yes | yes |
| F3 ncem | yes | yes |
| F4_pearson cellina | no | no |
| F4_pearson cpa | yes | yes |
| F4_pearson_clean cellina | no | no |
| F4_pearson_clean cpa | no | no |
| F5 scanvi | no | no |
| F5 tangram | no | no |
