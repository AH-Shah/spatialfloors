#!/bin/zsh
# Release check on development units (no confirmation data): re-run one P4 job per family with the release code into
# runs/verify and compare with the development outputs (scripts/verify_compare.py). F4 L6 = `spatialfloors check`.
# Usage: zsh scripts/verify_dev.sh > runs/verify/verify.log 2>&1
cd /Users/Hemato/project_3/antig/spatialfloors; export PATH=$HOME/.pixi/bin:$PATH PIXI_FROZEN=true OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
V=runs/verify; mkdir -p $V; G=gse303070_148191A_278N; H=gse303070_148191A_118N; A=allen_animal4_C57BL6J-4.003
run() { echo "== $*"; "$@" || echo "FAILED: $*"; }
for u in $G $H; do run pixi run -e floors python scripts/f1_fit.py $u expr_embed pearson_none_pc50 --out $V; done
run pixi run -e floors python scripts/f1_fit.py $G expr_embed pearson_z_pc100 --out $V
run env SF_OUT=$V pixi run -e floors python scripts/p3_f1_rescore.py
for u in gse1_232 $A; do run env SF_OUT=$V SF_CONFIG=ell400_lam0.8 pixi run -e floors python scripts/p3_f2.py $u floor; done
run env SF_OUT=$V/f2v3 SF_CONFIG=v3:ms50_200_800+c200 pixi run -e floors python scripts/p3_f2.py gse1_232 floor
for c in v3 v3_a100; do run env SF_OUT=$V/f3_$c SF_CONFIG=$c SF_F3_FEATS=v3 pixi run -e floors python scripts/p3_f3.py $G ridge; done
run env SF_OUT=$V/v2 SF_CONFIG='T+ctx:ell25:pc10:a1000' pixi run -e floors python scripts/p3_f3.py $G ridge
run env SF_OUT=$V SF_CONFIG=leaves15_lr0.05 SF_F3_FEATS=v3 pixi run -e scvi python scripts/p3_f3.py $G lightgbm
run env SF_OUT=$V SF_CONFIG=insitutype_spatial pixi run -e floors python scripts/p3_f5.py $G floors
run env SF_OUT=$V/f5_gse1 pixi run -e floors python scripts/p3_f5.py gse1_232 floors
run env SF_OUT=$V/t2_prior SF_T1=$V SF_BETAS=1,.5,2 pixi run -e floors python scripts/t2_f5_prior.py $G
echo "VERIFY_DONE $(date +%T)"
