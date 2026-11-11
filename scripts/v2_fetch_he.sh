# Replication round: H&E image and its alignment to the Xenium coordinates for the five 10x sections (pathologist
# annotation, registration_v2/ANNOTATION_PROTOCOL.md). Resumes partial files; checks each size against Content-Length.
# Usage: zsh scripts/v2_fetch_he.sh   (prints V2_HE_DONE <number of BAD or FAILED lines>)
set -u; cd ${0:A:h}/..; RAW=${SF_DATA:-data}/raw; PY="python3 -I scripts/v2_fetch.py"; BAD=0
X=https://cf.10xgenomics.com/samples/xenium/3.0.0
for n in Xenium_Prime_Human_Skin_FFPE Xenium_Prime_Human_Prostate_FFPE Xenium_Prime_Ovarian_Cancer_FFPE_XRrun Xenium_Prime_Breast_Cancer_FFPE Xenium_Prime_Human_Lung_Cancer_FFPE; do
  for s in he_imagealignment.csv he_image.ome.tif; do
    for try in 1 2 3; do
      $=PY file $X/$n/${n}_$s $RAW/tenx_xenium5k/$n/$s && $=PY size $X/$n/${n}_$s $RAW/tenx_xenium5k/$n/$s && break
      echo "retry $try $n $s"
    done || { echo "FAILED $n $s"; BAD=$((BAD + 1)) }
  done
done
echo "V2_HE_DONE $BAD"
