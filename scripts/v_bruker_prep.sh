#!/bin/zsh
# Bruker confirmation / validation slide -> crop manifest (metadata only) -> crop h5ad -> panel -> cache
# Usage: zsh scripts/v_bruker_prep.sh <unit> <channel map>
set -e; cd ${0:a:h}/..
pixi run -e floors python scripts/v_crops.py $1
pixi run -e floors python scripts/v_bruker_h5ad.py $1 $2
pixi run -e floors python scripts/p1_panels.py $1
pixi run -e floors python scripts/p2_cache_crop.py $1
