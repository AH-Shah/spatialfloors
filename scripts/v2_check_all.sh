# Replication round: check every fetched file (scripts/v2_fetch_all.sh) and write registration_v2/data_sha256.tsv.
# Zenodo zip members: size and CRC-32 against the zip directory. Plain files: size against Content-Length.
# NSCLC tarballs: gzip -t. Then one line per file: sha256, bytes, path under $SF_DATA/raw, source URL.
# Usage: zsh scripts/v2_check_all.sh   (prints V2_CHECK_DONE <number of BAD or FAILED lines>)
setopt pipefail; set -u; cd ${0:A:h}/..; RAW=${SF_DATA:-data}/raw; PY="python3 -I scripts/v2_fetch.py"
KEEP='(cell_feature_matrix\.h5|cells\.parquet|gene_panel\.json|experiment\.xenium|metrics_summary\.csv|annot|patholog|\.geojson|[Rr]eadme)'
OUT=registration_v2/data_sha256.tsv; TMP=$(mktemp); LOG=$(mktemp)
row() { print -r -- "$(shasum -a 256 $RAW/$1 | cut -d' ' -f1)\t$(stat -f %z $RAW/$1)\t$1\t$2" >> $TMP }
Z=https://zenodo.org/api/records/18000256/files
for t in Breast Colorectal Lung Lymphoma Ovarian; do
  $=PY verify $Z/xenium5K_$t.zip/content $RAW/cervilla2026_xenium5k "$KEEP" | tee -a $LOG || echo "FAILED cervilla $t" | tee -a $LOG
  for f in $RAW/cervilla2026_xenium5k/xenium5K_$t/*(.N); do row ${f#$RAW/} "$Z/xenium5K_$t.zip (member ${f:t})"; done
done
X=https://cf.10xgenomics.com/samples/xenium/3.0.0
for n in Xenium_Prime_Human_Skin_FFPE Xenium_Prime_Human_Prostate_FFPE Xenium_Prime_Ovarian_Cancer_FFPE_XRrun Xenium_Prime_Breast_Cancer_FFPE Xenium_Prime_Human_Lung_Cancer_FFPE; do
  for s in cell_feature_matrix.h5 cells.parquet gene_panel.json he_imagealignment.csv he_image.ome.tif; do
    $=PY size $X/$n/${n}_$s $RAW/tenx_xenium5k/$n/$s | tee -a $LOG || echo "FAILED 10x $n $s" | tee -a $LOG
    row tenx_xenium5k/$n/$s $X/$n/${n}_$s
  done
done
B=https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed
for s in Lung5_Rep1 Lung5_Rep2 Lung5_Rep3 Lung6 Lung9_Rep1 Lung9_Rep2 Lung12 Lung13; do
  $=PY size "$B/$s/$s+SMI+Flat+data.tar.gz" $RAW/cosmx_nsclc/$s.tar.gz | tee -a $LOG || echo "FAILED nsclc $s" | tee -a $LOG
  gzip -t $RAW/cosmx_nsclc/$s.tar.gz && echo "OK gzip $s" | tee -a $LOG || echo "BAD gzip $s" | tee -a $LOG
  row cosmx_nsclc/$s.tar.gz "$B/$s/$s+SMI+Flat+data.tar.gz"
done
{ print "sha256\tbytes\tpath\tsource"; cat $TMP } > $OUT
echo "V2_CHECK_DONE $(grep -c '^BAD\|^FAILED' $LOG) bad, $(grep -c '^OK' $LOG) ok, $(($(wc -l < $OUT) - 1)) files in $OUT"
