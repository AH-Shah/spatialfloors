#!/bin/zsh
# Download the first-ranked CELLxGENE dataset per new tissue (registration_v2/refs/selection.json, scripts/v2_ref_prep.py
# select) to data/raw/refs/<tissue>_<dataset id>/<file>.h5ad; size checked against the index, sha256 to the log.
# Usage: zsh scripts/v2_ref_fetch.sh <tissue ...>   (log: data/logs/v2_ref_fetch.log; last line V2_REF_FETCH_DONE <status>)
set -u; cd ${0:A:h}/..
LOG=data/logs/v2_ref_fetch.log; st=0
for t in "$@"; do
  read id url bytes < <(python3 -I -c "import json,sys;x=json.load(open('registration_v2/refs/selection.json'))['tissues'][sys.argv[1]][0];print(x['dataset_id'],x['h5ad'],x['h5ad_bytes'])" $t)
  d=data/raw/refs/${t}_${id}; f=$d/${url:t}; mkdir -p $d
  echo "[$(date +%T)] $t $id $url" >> $LOG
  caffeinate -ims curl -sS -f -L -C - -o $f $url >> $LOG 2>&1 || { echo "[$(date +%T)] $t curl failed" >> $LOG; st=1; continue; }
  n=$(stat -f %z $f); [ $n = $bytes ] || { echo "[$(date +%T)] $t size $n != $bytes" >> $LOG; st=1; continue; }
  echo "[$(date +%T)] $t ok $n bytes sha256 $(shasum -a 256 $f | cut -d' ' -f1)" >> $LOG
done
echo "V2_REF_FETCH_DONE $st" >> $LOG
