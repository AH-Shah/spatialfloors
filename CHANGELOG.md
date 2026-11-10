# Changes since `results-v1`

The registered verdicts and every tier-1 and tier-2 number are unchanged; both reports recompute byte-identically.

- Exploratory analyses and the whole-slide arm (`RESULTS.md`, `results/DEVIATIONS.md`).
- The tier-1 analysis checks the registered unit count of every pair.
- One data root (`SF_DATA`); `DATA.md` with every source and the SHA-256 of every input file.
- `floors` and `r` environments locked for Linux too; `r-setup` installs InSituType at the commit used.
- Development-only code removed (it remains at tags `registration-v1` and `results-v1`); documentation and formatting
  of the remaining code cleaned, with the code itself unchanged.
