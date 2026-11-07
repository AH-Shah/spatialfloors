# After the tier-1 report: tier-2 diagnostic (PROTOCOL §5) -> runs/p4/T2_corrected.{json,md}, secondary outputs (PROTOCOL §2)
# -> runs/p4/P5_secondary.{json,md}, and the exploratory unit-size moderators -> runs/p4/X_moderators.{json,md}.
# Usage: zsh scripts/p5_tier2_run.sh (after scripts/p5_tier1_run.sh and every tier-2 job)
set -e; cd ${0:A:h}/..; export PYTHONPATH=src PIXI_FROZEN=true OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4
CROPS=($(awk -F'\t' '$3 == "confirmation" {print $1}' registration/units/units.tsv))
echo "[$(date +%T)] start $(git rev-parse --short HEAD) crops ${#CROPS}"
SF_T1NAME=P5_corrected.json SF_T2NAME=T2_corrected pixi run -e floors python scripts/p5_tier2.py runs/p4 $CROPS; echo "[$(date +%T)] p5_tier2 done"
pixi run -e floors python scripts/p5_secondary.py runs/p4 $CROPS; echo "[$(date +%T)] p5_secondary done"
pixi run -e floors python scripts/x_moderators.py runs/p4 P5_corrected; echo "[$(date +%T)] x_moderators done"
echo "[$(date +%T)] TIER2_ANALYSIS_DONE"
