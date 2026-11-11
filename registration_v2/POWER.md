# Replication round: predicted power (from round-1 pilot data)

sd = residual sd of the per-unit D in round 1; N = planned round-2 units; δ = locked margin. Non-inferiority: one-sided 95% upper bound of mean D < δ. Superiority: Holm first step (α / 7). In brackets: sd × 1.5.

| Pair | Round-1 mean D | sd | δ | N (df) | NI power at round-1 mean | at D = 0 | at D = +δ/2 | Superiority power at D = +δ | D for 80% superiority power |
|---|---|---|---|---|---|---|---|---|---|
| F1|scvi | -0.0087 | 0.0072 | 0.0195 | 15 (13) | 1.00 (1.00) | 1.00 (1.00) | 1.00 (0.95) | 1.00 (1.00) | +0.0069 (+0.0104) |
| F1|nicheformer | -0.0379 | 0.0153 | 0.0195 | 15 (13) | 1.00 (1.00) | 1.00 (0.93) | 0.76 (0.47) | 0.97 (0.68) | +0.0147 (+0.0221) |
| F2|cellcharter | -0.1003 | 0.1378 | 0.0500 | 5 (4) | 0.64 (0.39) | 0.17 (0.12) | 0.10 (0.08) | 0.03 (0.02) | +0.3292 (+0.4938) |
| F3|lightgbm | -0.0155 | 0.0074 | 0.0043 | 15 (13) | 1.00 (1.00) | 0.69 (0.41) | 0.28 (0.18) | 0.33 (0.13) | +0.0072 (+0.0108) |
| F3|ncem | -0.0355 | 0.0155 | 0.0043 | 15 (13) | 1.00 (1.00) | 0.27 (0.17) | 0.13 (0.10) | 0.07 (0.04) | +0.0149 (+0.0224) |
| F5|scanvi | -0.0143 | 0.0211 | 0.0255 | 13 (11) | 1.00 (1.00) | 0.99 (0.86) | 0.65 (0.39) | 0.90 (0.52) | +0.0223 (+0.0335) |
| F5|tangram | -0.1132 | 0.0771 | 0.0255 | 13 (11) | 1.00 (0.99) | 0.30 (0.18) | 0.14 (0.10) | 0.08 (0.04) | +0.0818 (+0.1227) |

Joint non-inferiority power (all seven pairs): at the round-1 means 0.64 (0.38); at D = 0 0.01 (0.00); at D = +δ/2 0.00 (0.00).
