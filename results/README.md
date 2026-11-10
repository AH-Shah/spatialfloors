# Results files

Summary: [`../RESULTS.md`](../RESULTS.md). Deviations and additions after registration: [`DEVIATIONS.md`](DEVIATIONS.md).

| Path | Content |
|---|---|
| `reports/P5_corrected.{json,md}` | Tier 1, the registered analysis (`scripts/p5_analyze.py`) |
| `reports/T2_corrected.{json,md}` | Tier 2 (`scripts/p5_tier2.py`) |
| `reports/P5_secondary.{json,md}` | Secondary outputs (`scripts/p5_secondary.py`) |
| `reports/X_*.{json,md}` | Exploratory analyses (`scripts/x_*.py`); whole-slide arm: `X_wholeslide`, `X_f5_scale` |
| `WHOLESLIDE_PLAN.md` | Plan of the whole-slide arm, committed before any whole-slide output |
| `data_sha256.tsv` | Size and SHA-256 of every input file (`scripts/data_hashes.sh`) |
| `p4/f1_*.json`, `p4/f1_rescored.csv` | F1 runs; the primary score is the rescored table |
| `p4/f2_*.json`, `p4/f2v3/` | F2 runs (the balanced ARI is computed from saved partitions) and the multi-scale secondary floor |
| `p4/f3_*.json`, `p4/v2/` | F3 runs (per receiver and fold) and the v2 floor |
| `p4/f4/` | F4 scores per task and seed (Cellina, CPA, L6) and L6 per task |
| `p4/f5_*.json` | F5 runs (held-out-marker AUROC per class, class shares, IF-gated macro-F1) |
| `p4/t2/` | Tier-2 runs |
