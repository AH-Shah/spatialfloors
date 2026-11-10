# Tier-1 analysis on the sealed confirmation outputs: F1 rescore, F4 scores per slide / section, p5_analyze (confirmation crops + F4
# retrospective) -> runs/p4/P5_corrected.{json,md}, then the exploratory interior re-score. Usage: zsh scripts/p5_tier1_run.sh
# GSE1 222 holds out TVA (PROTOCOL §2, F4), so its scores are read with SF_HOLDOUT=TVA.
set -e; cd ${0:A:h}/..; export PYTHONPATH=src PIXI_FROZEN=true OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4
CROPS=($(awk -F'\t' '$3 == "confirmation" {print $1}' registration/units/units.tsv))
F4=(110 222 C57BL6J-1.069 C57BL6J-1.084 C57BL6J-1.086 C57BL6J-3.003 C57BL6J-3.011 C57BL6J-3.012 C57BL6J-638850.36 C57BL6J-638850.38 C57BL6J-638850.40)
RETRO=(retro:120 retro:210 retro:221 retro:242 retro:C57BL6J-2.036 retro:C57BL6J-2.039 retro:C57BL6J-2.041)
echo "[$(date +%T)] start $(git rev-parse --short HEAD) crops ${#CROPS}"
SF_OUT=runs/p4 pixi run -e floors python scripts/p3_f1_rescore.py; echo "[$(date +%T)] f1 rescore done"
for s in $F4; do SF_F4_RUN=runs/p4/f4 SF_BLIND=1 SF_HOLDOUT=$([ $s = 222 ] && echo TVA || echo CRC) pixi run -e floors python scripts/f4_score.py $s; echo "[$(date +%T)] f4_score $s done"; done
SF_EXPECT=manifests/expected_units_tier1.json SF_P5NAME=P5_corrected pixi run -e floors python scripts/p5_analyze.py runs/p4 $CROPS $RETRO; echo "[$(date +%T)] p5_analyze done"
pixi run -e floors python scripts/x_interior.py runs/p4 $CROPS; echo "[$(date +%T)] x_interior done"
echo "[$(date +%T)] TIER1_ANALYSIS_DONE"
