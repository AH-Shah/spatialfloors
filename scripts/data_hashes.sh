#!/bin/zsh
# SHA-256 of every input file under the data root that the code reads (DATA.md; 2026-11-10, external review R2).
# Usage: zsh scripts/data_hashes.sh <data root> -> runs/data_sha256.tsv (path relative to the root, bytes, sha256)
cd ${0:A:h}/..; R=${1:-data}; O=runs/data_sha256.tsv; : > $O.tmp
for d in raw/GSE1 raw/GSE0 raw/GSE303070 raw/spatch_xenium_5k raw/allen_zhuang raw/allen_yao_638850 raw/bruker_wtx raw/refs \
         raw/hgnc_2026-10-06 raw/ensembl_orthologs raw/novae_p5 derived; do
  find -L $R/$d -type f ! -name '*.log' ! -name '*DONE' ! -name '.DS_Store' | sort | while read -r f; do
    printf "%s\t%s\t%s\n" "${f#$R/}" "$(stat -L -f %z $f)" "$(shasum -a 256 $f | cut -d' ' -f1)" >> $O.tmp
  done
done
for f in $R/pelka.h5ad $R/processed/p5/labeled_h5ad/GSE303070_*.h5ad(N); do
  printf "%s\t%s\t%s\n" "${f#$R/}" "$(stat -L -f %z $f)" "$(shasum -a 256 $f | cut -d' ' -f1)" >> $O.tmp
done
mv $O.tmp $O; echo "[$(date +%T)] DATA_HASHES_DONE $(wc -l < $O) files"
