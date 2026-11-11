# Replication round: pathologist region annotation (F2 ground truth)

Draft, 2026-11-11. Part of registration v2. The labels it produces are deposited, with their sha256, before any F2
score is computed on these sections.

## Purpose

F2 (spatial domains) is scored against regions drawn by people who see only the tissue. Two pathologists annotate
each section independently; their agreement is the human ceiling for F2, and their consensus is the reference the
floor and the learned models are scored against (SELECTION.md, rule point 7).

## Sections

The five 10x Genomics Xenium Prime 5K FFPE sections with a post-run H&E of the same section: skin (melanoma),
prostate adenocarcinoma, ovarian carcinoma, breast carcinoma, lung carcinoma. Files per section:
`<name>_he_image.ome.tif` and `<name>_he_imagealignment.csv` (3 × 3 affine from H&E pixels to Xenium µm), fetched by
`scripts/v2_fetch_he.sh`, checksums in `data_sha256.tsv`. One pilot section first (prostate, the smallest tumour
image after skin) to time the work and fix the protocol; the protocol is frozen after the pilot and before the rest.

## Blinding

- Annotators see the H&E image only: no gene expression, no cell segmentation, no cluster, domain or cell-type
  output of any method, and not each other's annotation.
- They know the organ and that the section is a primary tumour; they are not told the study's hypotheses.
- The images are given to them as copies under neutral names (`section_A` … `section_E`); the key stays with the
  authors.

## Classes

One shared set, plus at most two tissue-specific classes agreed in the pilot (listed in the frozen version). A class
is drawn only where it is present.

| Class | Definition |
|---|---|
| Tumour | Invasive or in situ malignant cells, including tumour cells within stroma when they dominate the area |
| Tumour stroma | Desmoplastic or reactive stroma inside or at the edge of the tumour, without dominant tumour cells |
| Normal tissue | Non-neoplastic parenchyma of the organ (glands, epidermis, alveoli, ovarian stroma) |
| Immune aggregate | Dense lymphoid infiltrate or tertiary lymphoid structure, at least ~200 µm across |
| Necrosis | Coagulative or comedo necrosis |
| Other | Fat, large vessels, muscle, nerve, mucin, haemorrhage |
| Exclude | Folds, tears, out-of-focus or artefact areas (cells inside are not scored) |

## Procedure

- Software: QuPath (version recorded), brush and polygon tools at 10–20× equivalent magnification.
- Annotate the whole section; areas left unannotated are not scored. Minimum region size ~100 µm across, except
  immune aggregates as above.
- Export: one GeoJSON per section and annotator (QuPath "Export objects as GeoJSON", pixel coordinates of the H&E
  image), with the class names above. Record the time spent.

## From polygons to cells

`scripts/v2_annot_to_cells.py` (to write): H&E pixel polygons → Xenium µm with the alignment matrix; each cell takes
the class of the polygon that contains its centroid; cells in no polygon or in "Exclude" are unlabelled.

## Agreement and consensus

- Per section: Cohen's κ over cells labelled by both annotators, and the adjusted Rand index between the two label
  vectors (the F2 metric applied to the second annotator against the first: the human ceiling).
- Consensus: cells where both agree keep the class; disagreeing cells are unlabelled (not scored). If the consensus
  leaves fewer than two regions with ≥ 2% of the labelled cells, the section does not enter F2 (rule point 7).
- Crop: the round-1 region rule (`data.crops.choose_crop` with regions) on the consensus labels; that crop is then the
  section's crop for every family.

## Deposit

GeoJSONs, per-cell consensus labels, agreement table and their sha256 go into `registration_v2/annotations/` and are
part of the sealed registration before any F2 job is opened.
