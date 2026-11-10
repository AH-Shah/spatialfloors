#!/bin/zsh
# B_ratio (scANVI on ref + 50k random query cells, predict every cell) on the other five whole slides (221 done).
# 2026-11-10 user order; exploratory. Log: runs/xw/diag/x_f5_scale_all.log (ends X_F5_SCALE_ALL_DONE).
cd /Users/Hemato/project_3/antig/spatialfloors
export PIXI_FROZEN=true OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 SF_VARIANTS=B_ratio
for s in 222 242 210 110 120; do
  /Users/crocodilus/.pixi/bin/pixi run -e scvi python scripts/x_f5_scale.py $s 2>&1 | grep -E 'gse1_|Traceback|Error' || echo "[$(date +%T)] FAILED $s"
done
echo "[$(date +%T)] X_F5_SCALE_ALL_DONE"
