# Data

All inputs are public. They are not in this repository. The code reads them under one data root, `SF_DATA` (default:
`data/` in the repository, a directory or a symlink; `src/spatialfloors/paths.py`). SHA-256 and size of every file
read: [`results/data_sha256.tsv`](results/data_sha256.tsv).

The registration, the code and every output use short IDs (first column). They are kept unchanged, because the sealed
registration (`registration/`, DOI 10.5281/zenodo.23931264) and the crop manifests' SHA-256 values refer to them. Text,
tables and figures meant for readers use the names in the second column.

## Spatial datasets

| ID in code and registration | Name used in the text | Role | Source | Files read (under `raw/`) |
|---|---|---|---|---|
| `gse1_<sample>` (110, 120, 210, 221, 222, 242; 231, 232) | Crowell CRC `<sample>` | confirmation (110–242); validation (231); development (232) | Crowell et al., CosMx whole-transcriptome colon sections from normal mucosa through tubulovillous adenoma to carcinoma; Zenodo [10.5281/zenodo.15574384](https://doi.org/10.5281/zenodo.15574384) (CC-BY-4.0), code github.com/HelenaLC/WTx-TVA | `GSE1/cosmx_wt_<sample>.h5ad` (the record's `<sample>.h5ad`) |
| `gse1_<sample>w` | Crowell CRC `<sample>`, whole slide | exploratory whole-slide arm; validation (231w) | as above | as above |
| `gse0_S0` | Bruker colon WTx | development | Bruker Spatial Biology, CosMx Human Whole Transcriptome Colon Dataset (sigmoid adenocarcinoma, FFPE): [dataset page](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-human-whole-transcriptome-colon-dataset/); files `https://objects.liquidweb.services/smi-public/wtx_manuscript/colon_discovery/` | `GSE0/S0_exprMat_file`, `S0_metadata_file`, `S0_fov_positions_file` |
| `Donor<n>` (crops `gse303070_<slide>_<sample>`) | GSE303070 donor `<n>` | confirmation (11 donors); development (donors 46, 118, 278) | GEO [GSE303070](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303070): "A spatial atlas of colorectal cancer reveals the influence of stromal niches on tumour differentiation" (CosMx 1k) | `GSE303070/GSE303070_All_cells_metadata.csv`; `processed/p5/labeled_h5ad/GSE303070_CosMx_Slide_*.h5ad` (made from `GSE303070/Processed`) |
| `spatch_xenium5k` | SPATCH Xenium 5K (colon adenocarcinoma) | confirmation (F1 only; CRC in leave-one-cancer-out) | SPATCH benchmark (Xenium 5K section with adjacent CODEX), Nature Communications 16 (2025), doi [10.1038/s41467-025-64292-3](https://doi.org/10.1038/s41467-025-64292-3); SPATCH web server http://spatch.pku-genomics.org/ (raw data also NGDC GSA-Human HRA011129). Downloaded 2025-05-27 (file dates). The study has colon, liver and ovarian tumours; this file is the colon sample (CDX1 / CDX2 in 27% / 21% of cells, ovarian and hepatocyte markers absent) | `spatch_xenium_5k/spatch_xenium_5k.h5ad` |
| `allen_animal<k>` (k = 1–4; F4 retrospective `allen2`) | Zhuang MERFISH mouse `<k>` | confirmation (1–3); development (4) | Zhang et al. 2023, Nature, doi [10.1038/s41586-023-06808-9](https://doi.org/10.1038/s41586-023-06808-9) (Zhuang lab whole mouse brain MERFISH; Allen Brain Cell Atlas Zhuang-ABCA-1 to 4), as curated by CZ CELLxGENE Discover, collection `0cca8620-8dee-45d0-aef5-23f032a5cf09`. Dataset versions (`https://datasets.cellxgene.cziscience.com/<id>.h5ad`): animal 1 `aeaa6e77-7f6e-4d8c-b217-4c1e777cc0be`, animal 2 `93c3bb97-ea05-4ee0-a760-a1508cd04612`, animal 3 `0cd6f9fe-4054-480e-a4b7-ccce8d67c4cb`, animal 4 `4b5c682f-4b5d-4dcb-81c1-606fef2b601e` (read from each file's `uns/citation`) | `allen_zhuang/WB_MERFISH_animal<k>_*.h5ad` |
| `allen_yao638850` | Yao MERFISH mouse 638850 | confirmation | Yao et al. 2023, Allen Brain Cell Atlas dataset MERFISH-C57BL6J-638850: `https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/` `expression_matrices/MERFISH-C57BL6J-638850/20230830/C57BL6J-638850-raw.h5ad`, `metadata/MERFISH-C57BL6J-638850-CCF/20231215/views/cell_metadata_with_parcellation_annotation.csv`, `metadata/MERFISH-C57BL6J-638850/20241115/gene.csv` | `allen_yao_638850/` (these three files) |
| `bruker_kidney`, `bruker_skin` | Bruker WTx kidney RCC, skin SCC | confirmation | Bruker Spatial Biology CosMx whole-transcriptome FFPE datasets ([list](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/)); files `https://objects.liquidweb.services/smi-public/wtx_manuscript/kidney_discovery/` and `.../wtx_manuscript/skin_scc/` (local `md5sums.txt` identical to the bucket's) | `bruker_wtx/{kidney_rcc,skin_scc}` |
| `bruker_mobreast` | Bruker breast (multiomic) | confirmation | Bruker CosMx Human Multiomic Breast Dataset ([page](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-human-multiomic-breast-dataset/)); files `https://objects.liquidweb.services/smi-public/multiomic_breast/Flatfiles_RNA.zip`, `Flatfiles_Protein.zip` | `bruker_wtx/multiomic_breast/flat` |
| `bruker_breast` | Bruker WTx breast IDC | validation | Bruker CosMx Human Whole Transcriptome Breast Dataset ([page](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-human-whole-transcriptome-breast-dataset/)); files `https://objects.liquidweb.services/smi-public/wtx_manuscript/breast_discovery/` | `bruker_wtx/breast_discovery` |
| – | Novae tutorial slide | positive control (Novae) | Xenium human colon cancer P5 with the Novae authors' domains, Hugging Face dataset `prism-oncology/novae` | `novae_p5/` |

## Single-cell references (F5)

| Reference | Used for | Source | File |
|---|---|---|---|
| Pelka et al. 2021, CRC atlas | CRC units (Crowell CRC, GSE303070, Bruker colon) | GEO [GSE178341](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178341) | `pelka.h5ad` |
| Wu et al. 2021, breast | Bruker breast | GEO [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) | `raw/refs/breast_wu2021/` |
| Krishna et al. 2021, ccRCC | Bruker kidney | CZ CELLxGENE, `https://datasets.cellxgene.cziscience.com/6602fa5a-d663-4fe0-8dc3-90191c3a013b.h5ad` | `raw/refs/ccrcc_krishna2021/` |
| Ji et al. 2020, cSCC | Bruker skin | GEO [GSE144240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144240) | `raw/refs/cscc_ji2020/` |

## Annotation files

| File | Source |
|---|---|
| `raw/hgnc_2026-10-06/hgnc_complete_set_2026-10-06.txt` | HGNC complete set, downloaded 2026-10-06 (genenames.org) |
| `raw/ensembl_orthologs/mouse_human_orthologs.tsv` | Ensembl release 116 BioMart export (`https://www.ensembl.org/biomart/martservice`): mouse "Gene stable ID" → "Human gene stable ID" (`src/spatialfloors/models/nicheformer.py`) |

## External code (`external/`, not in this repository)

| Directory | Repository and commit |
|---|---|
| `external/cellina` | github.com/PMBio/cellina @ dd860d0 |
| `external/cellina-reproducibility` | github.com/PMBio/cellina-reproducibility @ 6bb060b; the F4 retrospective units read its `results/loo_summary_*_DEG_50_camready.csv` |
| `external/cpa` | github.com/theislab/cpa @ fbd7c02 |
| `external/nicheformer` | github.com/theislab/nicheformer @ 485cadb (tokenization notebooks, technology means) |
| `external/hf_nicheformer` | Hugging Face model `theislab/Nicheformer` (weights; SHA-256 in `manifests/models.json`) |
| `external/ncem_env`, `external/ncem_benchmarks_src` | NCEM 0.1.5 (`pip install ncem==0.1.5`, TensorFlow 2.15) and github.com/theislab/ncem_benchmarks |
| R environment (pixi env `r`) | InSituType 2.0, github.com/Nanostring-Biostats/InSituType @ f3ef0dd |

