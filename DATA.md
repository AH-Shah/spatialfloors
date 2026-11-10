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
| `gse0_S0` | Bruker colon WTx | development | Bruker Spatial Biology, CosMx Human Whole Transcriptome Colon Dataset (sigmoid adenocarcinoma, FFPE): [dataset page](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-human-whole-transcriptome-colon-dataset/) | `GSE0/S0_exprMat_file`, `S0_metadata_file`, `S0_fov_positions_file` |
| `Donor<n>` (crops `gse303070_<slide>_<sample>`) | GSE303070 donor `<n>` | confirmation (11 donors); development (donors 46, 118, 278) | GEO [GSE303070](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303070): "A spatial atlas of colorectal cancer reveals the influence of stromal niches on tumour differentiation" (CosMx 1k) | `GSE303070/GSE303070_All_cells_metadata.csv`; `processed/p5/labeled_h5ad/GSE303070_CosMx_Slide_*.h5ad` (made from `GSE303070/Processed`) |
| `spatch_xenium5k` | SPATCH Xenium 5K | confirmation (F1 only) | SPATCH benchmark, Nature Communications 16 (2025), doi [10.1038/s41467-025-64292-3](https://doi.org/10.1038/s41467-025-64292-3); data from the SPATCH web server. TODO: ask user for the exact download URL and the tissue | `spatch_xenium_5k/spatch_xenium_5k.h5ad` |
| `allen_animal<k>` (k = 1–4; F4 retrospective `allen2`) | Zhuang MERFISH mouse `<k>` | confirmation (1–3); development (4) | Zhang et al. 2023 (Zhuang lab), whole mouse brain MERFISH, Allen Brain Cell Atlas datasets Zhuang-ABCA-1 to 4 ([ABC atlas](https://alleninstitute.github.io/abc_atlas_access/descriptions/Zhuang-ABCA-2.html)). TODO: ask user for the exact file URLs of `WB_MERFISH_animal<k>_*.h5ad` | `allen_zhuang/WB_MERFISH_animal{1,3,4}_*.h5ad`; animal 2 under `derived/` |
| `allen_yao638850` | Yao MERFISH mouse 638850 | confirmation | Yao et al. 2023, Allen Brain Cell Atlas dataset MERFISH-C57BL6J-638850 | `allen_yao_638850/C57BL6J-638850-raw.h5ad`, `cell_metadata_with_parcellation_annotation.csv`, `gene.csv` |
| `bruker_mobreast`, `bruker_skin`, `bruker_kidney` | Bruker WTx breast (multiomic), skin SCC, kidney RCC | confirmation | Bruker Spatial Biology CosMx whole-transcriptome FFPE datasets ([list](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/)). TODO: ask user for the exact dataset pages | `bruker_wtx/{multiomic_breast,skin_scc,kidney_rcc}` |
| `bruker_breast` | Bruker WTx breast IDC | validation | as above (breast discovery dataset). TODO: exact page | `bruker_wtx/breast_discovery` |
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
| `raw/ensembl_orthologs/mouse_human_orthologs.tsv` | Ensembl BioMart mouse–human orthologs. TODO: ask user for the Ensembl release |

## Derived inputs (`derived/`)

Made by the investigators' earlier project, before this study. Each is listed with its SHA-256 in `results/data_sha256.tsv`.

| File | What it is |
|---|---|
| `WB_MERFISH_animal2_coronal.h5ad` | Zhuang MERFISH mouse 2, the same format as the other `WB_MERFISH_*` files |
| `gse303070_fov_meta.csv` | GSE303070 FOV-to-donor table, used to cut the donor crops |
| `perturb_opt/` | L6 floor reference outputs on Crowell CRC 231 and 232 (F4 validation check and a test fixture) |

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

Vendored inside this repository with their source hashes: `src/spatialfloors/floors/l6.py` (the investigators' earlier
L6 baseline, verbatim) and `src/spatialfloors/floors/validate.py` (verbatim from the investigators' earlier library).
