# Exploratory baseline ladder: the naive baselines that need a new CPU run on the round-1 confirmation crops, then the
# report. F1 log-normalised PCA 50, F2 GMM on expression PCA 50, F3 ridge on the type mix only.
# Usage (cwd = the runs root): zsh <repo>/scripts/x_ladder_all.sh
set -u; P=${0:A:h}/..; export PYTHONPATH=$P/src PIXI_FROZEN=true OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PIXI=${PIXI:-pixi}; OUT=runs/p4/x_ladder; RUN=($PIXI run --manifest-path $P/pixi.toml -e floors python)
crops() { awk -F'\t' -v f=$1 '$3 == "confirmation" && $4 ~ f {print $1}' $P/registration/units/units.tsv }
echo "[$(date +%T)] start"
for u in $(crops F1); do
  [ -f $OUT/f1_${u}_expr_pca_pc50.json ] || $RUN $P/scripts/f1_fit.py $u expr_pca pc50 --out $OUT || echo "[$(date +%T)] FAILED f1 $u"
done
echo "[$(date +%T)] F1 done"
for u in $(crops F2); do
  [ -f $OUT/f2_${u}_expr_gmm.json ] || $RUN $P/scripts/x_ladder.py f2 $u $OUT || echo "[$(date +%T)] FAILED f2 $u"
done
echo "[$(date +%T)] F2 done"
for u in $(crops F3); do
  [ -f $OUT/f3_${u}_ridge.json ] || SF_OUT=$OUT SF_CONFIG=T:ell25:pc10:a1000 $RUN $P/scripts/p3_f3.py $u ridge || echo "[$(date +%T)] FAILED f3 $u"
done
echo "[$(date +%T)] F3 done"
$RUN $P/scripts/x_ladder.py report $OUT runs/p4 runs/p4 && echo "X_LADDER_ALL_DONE"
