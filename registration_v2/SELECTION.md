# Replication round: dataset selection

Draft, 2026-11-10. Written before any file of the datasets below was opened. It becomes part of registration v2, which
is deposited before any score is computed on these units. The floors, learned-model configurations, margins, crop rule
and analysis of registration v1 (DOI 10.5281/zenodo.23931264) are used unchanged.

## Rule

A dataset enters if it meets every point. Only metadata were used to apply the rule: tissue, platform, panel size,
cell count, licence and the list of files.

1. Imaging-based spatial transcriptomics with a per-cell count matrix and cell coordinates from the data producer.
2. Human primary tumour tissue. (Mouse brain enters only for F2, point 7.)
3. Panel of at least 900 genes, the smallest human panel of v1 (CosMx 1k). Smaller panels were already set aside
   during the development of v1.
4. Public without an access agreement or registration, under a licence that allows reuse (CC BY, or the producer's
   open data release).
5. The study and its specimens had no role in v1: development, validation, confirmation or an exploratory analysis.
6. At least 10,000 cells in a section after the registered gene and cell filters.
7. F2 needs region labels drawn by a pathologist or anatomist, not derived from expression clustering, with at least
   two regions that each hold 2% of the labelled cells. Datasets without such labels do not enter F2.

Units and families:

- Unit = patient. Sections or replicates of one patient are averaged into one unit, as in v1.
- A specimen profiled on several platforms enters once, on the eligible platform with the fewest v1 confirmation units
  (Xenium 5K, then CosMx 1k).
- F1 and F3: every unit. F3 units whose label-free identities form a single lineage are not scored (v1, DEVIATIONS).
- F5: carcinomas only, because the classes are epithelial, immune non-myeloid and myeloid. Tissues with a v1 reference
  reuse it (CRC, breast); a new tissue gets one CELLxGENE reference by the v1 rule (`scripts/ref_prep.py`). A tissue
  whose reference fails that rule's QC does not enter F5.
- F4: not replicated.

## Search

2026-11-10: web searches for public Xenium, CosMx and MERFISH / MERSCOPE tumour data (vendor dataset pages, Zenodo,
GEO, platform-comparison studies), and the datasets already on disk. Sources read: the 10x Genomics, Bruker Spatial
Biology and Vizgen dataset pages, Zenodo records 18000256, 17986017 and 8327576, and the CORA record 10.34810/DATA3045.

## Selected

| Source | Platform | Units | Access |
|---|---|---|---|
| Cervilla et al. 2026, Genome Biology 27:22 (one FFPE primary tumour per cancer type) | Xenium Prime 5K | breast, colorectal, lung, lymphoma, ovarian (5 patients) | Zenodo 10.5281/zenodo.18000256, CC BY 4.0 |
| CosMx NSCLC FFPE dataset (He et al. 2022, Nat Biotechnol) | CosMx 960-plex | Lung 5 (3 sections), 6, 9 (2 sections), 12, 13 (5 patients) | Bruker Spatial Biology / NanoString open data release |
| 10x Genomics Xenium Prime 5K FFPE datasets | Xenium Prime 5K | skin melanoma, prostate adenocarcinoma, ovarian cancer, breast cancer, lung cancer; cervical cancer if its files are found (one donor each) | 10x Genomics datasets, CC BY 4.0 |

15 or 16 units from three sources. F5-eligible (carcinoma): breast 2, ovarian 2, lung 7, colorectal 1, prostate 1,
cervix 0 or 1.

## Examined and not selected

| Dataset | Reason (rule point) |
|---|---|
| Cervilla et al., bladder | Profiled on Xenium Multi-Tissue (377 genes) only; CosMx 1k slides to be checked for a bladder region (3) |
| Cervilla et al., Xenium Multi-Tissue and CosMx data of the five selected tumours | Same specimens; one platform per specimen |
| Vizgen MERSCOPE FFPE Human Immuno-oncology (8 tumour types, 16 datasets) | 500-gene panel (3); registration form (4) |
| 10x Xenium Prime 5K reactive lymph node | Not a tumour (2) |
| Tissue-microarray platform comparisons (GEO GSE308145–GSE308148) | Cores far below 10,000 cells (6) |
| SPATCH liver and ovarian sections | Study used in v1 (5) |
| STARmap PLUS mouse CNS atlas (Zenodo 8327576) | Region labels come from expression clustering (7); CC BY-NC-ND |
| GSE280314, GSE308699, GSE293124 (on disk) | Set aside during v1 development: panels of 422 and 541 genes, and metastases without a matching reference |

F2 has no eligible dataset so far. The pathology annotations of the 10x ovarian dataset are checked against point 7
once its files are listed; if they pass, that unit enters F2.

## Stop rule

If fewer than 8 units or fewer than 3 sources remain after the metric-blind preparation (file checks, panels, crops,
references), the round stops and that is reported. The rule above is not relaxed.
