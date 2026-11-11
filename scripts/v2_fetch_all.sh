# Replication round: fetch every selected dataset (registration_v2/SELECTION.md) into $SF_DATA/raw.
# Count matrices, cell tables, gene panels and annotation files only. Usage: zsh scripts/v2_fetch_all.sh
set -u; cd ${0:A:h}/..; RAW=${SF_DATA:-data}/raw; PY="python3 -I scripts/v2_fetch.py"
KEEP='(cell_feature_matrix\.h5|cells\.parquet|gene_panel\.json|experiment\.xenium|metrics_summary\.csv|annot|patholog|\.geojson|[Rr]eadme)'
for t in Breast Colorectal Lung Lymphoma Ovarian; do
  $=PY get https://zenodo.org/api/records/18000256/files/xenium5K_$t.zip/content $RAW/cervilla2026_xenium5k "$KEEP" || echo "FAILED cervilla $t"
done
X=https://cf.10xgenomics.com/samples/xenium/3.0.0
for n in Xenium_Prime_Human_Skin_FFPE Xenium_Prime_Human_Prostate_FFPE Xenium_Prime_Ovarian_Cancer_FFPE_XRrun Xenium_Prime_Breast_Cancer_FFPE Xenium_Prime_Human_Lung_Cancer_FFPE; do
  for s in cell_feature_matrix.h5 cells.parquet gene_panel.json; do
    $=PY file $X/$n/${n}_$s $RAW/tenx_xenium5k/$n/$s || echo "FAILED 10x $n $s"
  done
done
B=https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed
for s in Lung5_Rep1 Lung5_Rep2 Lung5_Rep3 Lung6 Lung9_Rep1 Lung9_Rep2 Lung12 Lung13; do
  $=PY file "$B/$s/$s+SMI+Flat+data.tar.gz" $RAW/cosmx_nsclc/$s.tar.gz || echo "FAILED nsclc $s"
done
echo "V2_FETCH_DONE"
