# Registration package

Confirmatory protocol and locked inputs for the `spatialfloors` benchmark (closed-form floors vs learned spatial-omics
models), fixed before any confirmation metric was computed (except the F4 retrospective comparison, `PROTOCOL.md` §3).
DOI 10.5281/zenodo.23931264; code https://github.com/AH-Shah/spatialfloors, tag `registration-v1`.

| Path | Content |
|---|---|
| `PROTOCOL.md` | Aim, hypotheses, families, data roles, analysis, tier-2 diagnostic, deviation policy, compute, development history |
| `configs/` | Locked configurations (`locked_v3.json`, NCEM test configuration `locked_v3_ncem.json`; earlier locks v1, v2 for the record), tier-2 design and unit subset (`tier2.json`), power (`power_v3.json`), model versions (`models.json`), F4 sections (`allen.json`), GSE303070 donor map, IF gate / QC, F5 held-out markers (`panels_index.json`) and new-cancer references (`refs_*.json`) |
| `units/units.tsv` | Every crop: unit, role, families, platform stratum, cells, SHA-256 of its cell-ID manifest |
| `units/crops/`, `units/panels/` | Cell-ID manifests per crop (with the crop index) and gene panels per unit; in the git repository the cell-ID manifests are `manifests/crops/*.csv.gz` (same bytes) |
| `p4_jobs.txt` | Sealed confirmation job list: tier 1 (gates G7, P4LOCK) and tier 2 (gate T2) |
| `spatialfloors-<commit>.tar.gz` | Code archive of the release commit (`git archive`; run outputs and internal notes excluded) |
| `SHA256SUMS` | SHA-256 of every file above |

Verify: `shasum -a 256 -c SHA256SUMS` (in this folder). The archive can be regenerated from the repository with
`git archive --format=tar.gz --prefix=spatialfloors-<commit>/ <commit>`; the copies in `configs/` and `units/` are
byte-identical to `manifests/` in the archive.

Release checks (development and validation data only): `pytest` (34 passed, 1 slow test skipped), `spatialfloors check`
(L6 bit-identical to the reference run), and one rerun per family with the release code compared with the development
outputs (`scripts/verify_dev.sh`, `scripts/verify_compare.py`: 16 of 16 identical); tier-1 and tier-2 analysis dry run on
development outputs (`scripts/verify_fixture.sh`).

Built with `scripts/registration_build.py copy` and `seal <commit>` (environment `floors`, `pixi.lock`).
