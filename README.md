# spatialfloors

Pre-registered confirmatory benchmark of closed-form floors against published learned models in spatial omics, in five
task families: cell identity (F1), spatial domains (F2), spatial context → receiver state (F3), held-out cell-type
counterfactuals (F4) and label transfer (F5).

- Protocol: [`registration/PROTOCOL.md`](registration/PROTOCOL.md); registration DOI
  [10.5281/zenodo.23931264](https://doi.org/10.5281/zenodo.23931264), tag `registration-v1`.
- Results: [`RESULTS.md`](RESULTS.md); release `results-v1.1`, DOI [10.5281/zenodo.24109807](https://doi.org/10.5281/zenodo.24109807) (first release
  `results-v1`: DOI [10.5281/zenodo.24044400](https://doi.org/10.5281/zenodo.24044400)). The floors were non-inferior to the learned model in all 12 registered pairs, and no learned model was
  superior after Holm correction.
- Deviations and additions after registration: [`results/DEVIATIONS.md`](results/DEVIATIONS.md).
- Data sources: [`DATA.md`](DATA.md).

## Layout

| Path | Content |
|---|---|
| `src/spatialfloors/` | floors (`floors/`), metrics (`metrics/`), learned-model wrappers (`models/`), data loading (`data/`), statistics (`stats.py`) |
| `scripts/` | data preparation (`p1_panels.py`, `p2_cache_crop.py`, `v_bruker_*`, `ref_prep.py`, `f4_crop.py`, `f4_merfish_prep.py`); model runs (`f1_fit.py`, `p3_f2.py`, `p3_f3.py`, `p3_f3_ncem_auth.py`, `p3_f5.py`, `f4_run.py`, `f4_l6.py`); confirmation job list (`p4_build.py`) and job pool (`sched.sh`); tier 1 (`p5_tier1_run.sh`, `p5_analyze.py`, `p5_loco_crc.py`), tier 2 (`p5_tier2_run.sh`, `p5_tier2.py`), secondary outputs (`p5_secondary.py`); positive controls (`p3_posctl_*`); exploratory analyses (`x_*`, whole-slide arm `xw_*`) |
| `manifests/` | locked configurations, crop cell-ID lists, gene panels |
| `registration/` | the registration package (protocol, configurations, units, sealed job list, `SHA256SUMS`) |
| `results/` | reports and per-crop outputs ([index](results/README.md)) |
| `tests/` | unit tests |

Code used only during development (landscape sweeps, lock derivation, validation gates) is at tags `registration-v1`
and `results-v1`.

## Environments and tests

Environments are defined in `pixi.toml` / `pixi.lock`: `floors` for the floors and the analysis, `r` for InSituType,
and one environment per learned model. `floors` and `r` are locked for macOS arm64 and Linux x86-64, the learned-model
environments for macOS arm64 (they ran on Apple MPS). On macOS arm64 the floors and the analysis reproduce the committed
outputs exactly; on Linux, expect agreement within floating-point tolerance.

```
pixi run --frozen -e r r-setup        # once: Bioconductor data package and InSituType (commit f3ef0dd)
pixi run --frozen -e floors pytest
```

Data are public and not included. Put them under one root and set `SF_DATA` to it (default: `data/` in the
repository); [`DATA.md`](DATA.md) lists every source and the SHA-256 of every input file. Run outputs record the commit
of the repository they ran from (`spatialfloors_commit`).

## Citation and licence

See `CITATION.cff`. Code: MIT (`LICENSE`). Registration package on Zenodo: CC BY 4.0. Use of AI tools:
[`AI_DISCLOSURE.md`](AI_DISCLOSURE.md).
