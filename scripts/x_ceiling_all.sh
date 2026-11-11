#!/bin/zsh
# Exploratory ceilings on the round-1 confirmation crops of one family, from the round-1 runs root.
# Usage: zsh scripts/x_ceiling_all.sh <round-1 runs root> <family> -> results/p4/x_ceiling/<family>_<crop>.json, log beside
set -eu
P=${0:A:h}/..; R=$1; FAM=$2; OUT=$P/results/p4/x_ceiling; mkdir -p $OUT
CROPS=($(awk -F'\t' -v f=${FAM:u} '$3 == "confirmation" && $4 ~ f {print $1}' $P/registration/units/units.tsv))
cd $R; LAB=runs/p4/labels; [[ $FAM = f[13] ]] && LAB=runs/p4
for c in $CROPS; do
  PIXI_FROZEN=true PYTHONPATH=$P/src OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
    pixi run --manifest-path $P/pixi.toml -e floors python $P/scripts/x_ceiling.py $FAM $LAB $c $OUT 2>&1 | grep -v -i warn
done
echo "X_CEILING_${FAM:u}_DONE"
