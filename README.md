# spatialfloors

Pre-registered confirmatory benchmark of closed-form floors against published learned models in spatial omics, in five
task families: cell identity (F1), spatial domains (F2), spatial context → receiver state (F3), held-out cell-type
counterfactuals (F4) and label transfer (F5).

- Protocol: [`registration/PROTOCOL.md`](registration/PROTOCOL.md). Registration: DOI
  [10.5281/zenodo.23931264](https://doi.org/10.5281/zenodo.23931264), tag `registration-v1`.
- Results: [`RESULTS.md`](RESULTS.md), tag `results-v1`. The floors were non-inferior to the learned model in all 12
  registered pairs, and no learned model was superior after Holm correction. Execution events and deviations:
  [`results/EXECUTION_LOG.md`](results/EXECUTION_LOG.md).
- Metric definitions and their tests: [`SPEC.md`](SPEC.md).

## Layout

| Path | Content |
|---|---|
| `src/spatialfloors/` | floors (`floors/`), metrics (`metrics/`), learned-model wrappers (`models/`), stress panels (`stress/`), statistics (`stats.py`) |
| `scripts/` | runners per phase: `p0_*`–`p3_*` development, `v_*` / `v3_*` validation and the v3 lock, `p4_build.py` (confirmation job list), `sched.sh` (job pool), `p5_tier1_run.sh` / `p5_analyze.py` (tier 1), `p5_tier2_run.sh` / `p5_tier2.py` (tier 2), `p5_secondary.py` (secondary outputs), `x_*.py` (exploratory), `registration_build.py` (package) |
| `manifests/` | locked configurations, crop cell-ID lists, gene panels, references |
| `registration/` | the registration package (protocol, configuration and unit copies, sealed job list, `SHA256SUMS`) |
| `results/` | tier-1, tier-2, secondary and exploratory reports, per-crop run outputs, execution log ([index](results/README.md)) |
| `tests/` | unit tests (`pytest`) |

## Environments and checks

Environments are defined in `pixi.toml` / `pixi.lock` (`floors` for the floors and the analysis; one environment per
learned model). On macOS arm64:

```
pixi run -e floors pytest
pixi run -e floors spatialfloors check
```

At `results-v1`, `pytest` gives 33 passed and 1 skipped (the slow `spatialfloors check` test, which needs the GSE1 data).

Data are public (sources in the protocol, §3) and not included. Commit hashes recorded inside the lock files and run
stamps refer to the development repository. The registered state is this repository at tag `registration-v1`, and the
code as run for the results is at `results-v1`.

## Citation and licence

See `CITATION.cff`. Code: MIT (`LICENSE`). Registration package on Zenodo: CC BY 4.0.
