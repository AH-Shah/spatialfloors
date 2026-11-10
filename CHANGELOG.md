# Changes since `results-v1` (2026-11-08)

The registered verdicts and every tier-1 and tier-2 number are unchanged; the tier-1 and tier-2 reports recompute
byte-identically from the committed outputs. Each change has a row in `results/EXECUTION_LOG.md`.

## Exploratory analyses (no decision weight)
- F5 held-out-marker AUROC at the comparator's call rate (`X_f5_matched`); F1 with Euclidean k-NN (`X_f1_euclid`).
- Whole-slide arm on the six GSE1 confirmation slides (plan `results/WHOLESLIDE_PLAN.md`, committed before any
  whole-slide output; report `X_wholeslide`), and the F5 slide-scale diagnostic (`X_f5_scale`).
- F4 by cancer and by provenance (`X_f4_split`); replication at donor, slide and study level (`X_levels`); positive
  controls of the learned models in one table (`X_posctl`); configurations explored in development (`X_effort`).

## Analysis code
- `scripts/p5_analyze.py` checks the registered unit count of every pair (`SF_EXPECT`,
  `manifests/expected_units_tier1.json`); output unchanged.

## Reproducibility and documentation
- One data root (`SF_DATA`, `src/spatialfloors/paths.py`); no hard-coded paths.
- `DATA.md`: sources of every dataset, the names used in the text, and the SHA-256 of every input file
  (`results/data_sha256.tsv`).
- `floors` and `r` environments locked for Linux x86-64 too; macOS arm64 lock unchanged.
- SPEC: the locked v3 floors at the top; F3 far-ring docstring corrected to 150–400 µm.
- InSituType is called with `pixi run --frozen`; `test_genes` skips without the HGNC file.
- `AI_DISCLOSURE.md`.
