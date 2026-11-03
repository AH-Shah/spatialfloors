#!/bin/zsh
# Dry run of the confirmatory analysis on development units: arrange development outputs (runs/p3, p3v2, p3v3, v3) and
# the release re-runs (scripts/verify_dev.sh) in the P4 layout under runs/verify, then run tier 1 (p5_analyze.py, and the
# pre-release version for comparison) and tier 2 (p5_tier2.py, every unit treated as subset). No confirmation data.
# Usage: zsh scripts/verify_fixture.sh [old_p5.py] > runs/verify/fixture.log 2>&1
cd /Users/Hemato/project_3/antig/spatialfloors; export PATH=$HOME/.pixi/bin:$PATH PIXI_FROZEN=true OMP_NUM_THREADS=4
R=runs/verify; G=gse303070_148191A_278N; H=gse303070_148191A_118N; A=allen_animal4_C57BL6J-4.003; P=$PWD
py() { pixi run -e floors python "$@"; }
mkdir -p $R/labels $R/t2/f1 $R/t2/f2_floor_dev/labels $R/t2/f3_lightgbm_dev
py -c "
import pandas as pd
u = ['$G', '$H']; a, b = pd.read_csv('runs/p3/f1_rescored.csv'), pd.read_csv('runs/p3v2/f1_rescored.csv')
b[b.unit.isin(u) & b.model.isin(['scvi', 'nicheformer'])].to_csv('$R/f1_rescored_learned.csv', index=False)
pd.concat([a[a.unit.isin(u) & (a.model.isin(['scvi']) | (a.config == 'pearson_z_pc100'))], b[b.unit.isin(u) & b.model.isin(['scvi', 'nicheformer'])]]).to_csv('$R/t2/f1/f1_rescored.csv', index=False)"
for m in cellcharter novae; do for u in gse1_232 $A; do [ -e runs/p3v2/labels/f2_${u}_$m.npz ] && ln -sf $P/runs/p3v2/labels/f2_${u}_$m.npz $R/labels/; done; done
for u in gse1_232 $A; do ln -sf $P/runs/p3/labels/f2_${u}_floor.npz $R/t2/f2_floor_dev/labels/; done
mkdir -p $R/v2; ln -sf $P/runs/p3v2/f3_${G}_ridge.json $R/v2/; ln -sf $P/$R/f3_v3/f3_${G}_ridge.json $R/; ln -sf $P/runs/p3v3/f3_${G}_ncemA.json $R/; ln -sf $P/runs/v3/lgbm/f3_${G}_lightgbm.json $R/t2/f3_lightgbm_dev/
for c in v3_a300 v3_a3000 v3_a10000; do env SF_OUT=$R/t2/f3_floor_$c SF_CONFIG=$c SF_F3_FEATS=v3 pixi run -e floors python scripts/p3_f3.py $G ridge; done
ln -sf $P/$R/f5_gse1/f5_gse1_232_floors.json $R/; for u in $G gse1_232; do for m in scanvi tangram; do ln -sf $P/runs/p3v2/f5_${u}_$m.json $R/; done; done
env SF_OUT=$R/t2/f5_floor_prior SF_T1=$R/f5_gse1 pixi run -e floors python scripts/t2_f5_prior.py gse1_232
U=($G $H gse1_232 $A)
py scripts/p5_analyze.py $R $U
[ -n "$1" ] && env SF_LOCKED=manifests/locked_v3.json SF_P5NAME=P5_prerelease pixi run -e floors python $1 $R $U
env SF_T2_SUBSET=all pixi run -e floors python scripts/p5_tier2.py $R $U
echo "FIXTURE_DONE $(date +%T)"
