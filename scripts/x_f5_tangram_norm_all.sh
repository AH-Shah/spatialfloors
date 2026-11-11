# Exploratory Tangram on normalised input, every F5 confirmation crop (one GPU job at a time), then the report.
# Usage (cwd = the runs root): zsh <repo>/scripts/x_f5_tangram_norm_all.sh
set -u; P=${0:A:h}/..; export PYTHONPATH=$P/src PIXI_FROZEN=true OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PIXI=${PIXI:-pixi}; OUT=runs/p4/x_f5tn
CROPS=($(awk -F'\t' '$3 == "confirmation" && $4 ~ /F5/ {print $1}' $P/registration/units/units.tsv))
echo "[$(date +%T)] start, ${#CROPS} crops"
for u in $CROPS; do
  [ -f $OUT/f5tn_$u.json ] && continue
  $PIXI run --manifest-path $P/pixi.toml -e scvi python $P/scripts/x_f5_tangram_norm.py unit $u $OUT || echo "[$(date +%T)] FAILED $u"
  echo "[$(date +%T)] done $u"
done
$PIXI run --manifest-path $P/pixi.toml -e floors python $P/scripts/x_f5_tangram_norm.py report $OUT runs/p4 && echo "X_F5TN_DONE"
