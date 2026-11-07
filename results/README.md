# Results files

The summary of the results is in [`../RESULTS.md`](../RESULTS.md). Events, crash fixes and deviations after registration
are in [`EXECUTION_LOG.md`](EXECUTION_LOG.md).

| Path | Content |
|---|---|
| `reports/P5_corrected.{json,md}` | Tier 1, the registered analysis on every registered unit (`scripts/p5_analyze.py`) |
| `reports/P5.{json,md}` | First tier-1 run (2026-11-06): GSE1 222 missing from F4; kept for the record |
| `reports/T2_corrected.{json,md}` | Tier 2 (`scripts/p5_tier2.py`), recording `P5_corrected.json` |
| `reports/T2.{json,md}` | Tier 2 as first run, recording `P5.json`; otherwise identical to `T2_corrected` |
| `reports/P5_secondary.{json,md}` | Secondary outputs (`scripts/p5_secondary.py`), including per-crop F2 balanced ARI |
| `reports/X_interior.{json,md}`, `reports/X_moderators.{json,md}` | Exploratory: interior-only re-score, unit size as a moderator |
| `p4/f1_<crop>_<model>_<config>.json`, `p4/f1_rescored.csv` | F1 runs; the primary score is the rescored table |
| `p4/f2_<crop>_<model>.json` | F2 runs (plain ARI and run metadata; the balanced ARI is computed from saved partitions) |
| `p4/f3_<crop>_<model>.json`, `p4/v2/` | F3 runs (per receiver and per fold), and the v2 floor |
| `p4/f5_<crop>_<model>.json` | F5 runs (held-out-marker AUROC per class, class shares, IF-gated macro-F1) |
| `p4/f2v3/` | F2 multi-scale secondary floor |
| `p4/f4/f4_<slide or section>.csv`, `p4/f4/l6/` | F4 scores per task and seed (Cellina, CPA, L6), and L6 per task on every F4 slide / section |
| `p4/t2/<family>_<model>_<config>/`, `p4/t2/f1/`, `p4/t2/f5_floor_prior/` | Tier-2 runs |

`p4/` mirrors the run directory `runs/p4` without partitions, embeddings, counterfactuals and logs. Every output records
the commit of the development repository it ran from (`spatialfloors_commit`). The code as run is this repository at tag
`results-v1`.
