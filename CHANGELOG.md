# Changes since `results-v1`

`results-v1.1` (2026-11-10, DOI 10.5281/zenodo.24109807).

The registered verdicts and every tier-1 and tier-2 number are unchanged; both reports recompute byte-identically.

- Unit label corrected: SPATCH Xenium 5K (colon adenocarcinoma) counts as CRC in leave-one-cancer-out
  (`scripts/p5_loco_crc.py` -> `results/reports/P5_loco_crc.md`); the registered code had labelled it "tissue
  unverified". Only F1's split without CRC changes; every split stays non-inferior.
- Exploratory analyses and the whole-slide arm (`RESULTS.md`, `results/DEVIATIONS.md`).
- The tier-1 analysis checks the registered unit count of every pair.
- One data root (`SF_DATA`); `DATA.md` with every source and the SHA-256 of every input file.
- `floors` and `r` environments locked for Linux too; `r-setup` installs InSituType at the commit used.
- Development-only code removed (it remains at tags `registration-v1` and `results-v1`); documentation and formatting
  of the remaining code cleaned, with the code itself unchanged.
