# spatialfloors SPEC (locked at tag `p0-lock`)

Metric, split and statistics definitions; each item names the function that implements it and the test that pins it.
The locked configurations, units and decision rules of the confirmatory study are registered in
`registration/PROTOCOL.md` and `manifests/locked_v3.json`; where this file names an earlier floor configuration, the
protocol governs. Section references (PLAN §n) point to the development record summarised in the protocol's appendix.

**Locked floors of the confirmatory study (v3; added 2026-11-10). These supersede the floor configurations described in
the family sections below, which record earlier development stages.**
- F1: analytic Pearson residuals on X_a, centred (not z-scored), PCA 50 (`pearson_none_pc50`).
- F2: full-covariance GMM on [PCA 50, mean over the 200 nearest neighbours], Gaussian ℓ 400 µm, λ .8.
- F3: ridge (α 1000) on the v2 features (human: type mix + bleed-safe context of other lineages, ℓ 25 µm, 10
  components; mouse: context only, ℓ 100 µm, 50 components) plus two same-lineage rings of the receiver's own genes:
  50–150 µm among the 600 nearest cells and 150–400 µm among the 1,500 nearest (median reach 220–330 µm on CosMx).
  The note in `locked_v3.json` says 150–250 µm; that file is hashed in the tier-2 manifest and is left unchanged.
- F4: L6 two-way (`floors.l6.twoway`), unchanged.
- F5: InSituType ML with a spatial prior on its log-likelihoods (15 nearest neighbours, 5 iterations;
  `floors.f5.spatial_prior`).
- GAT (F3) was dropped as a test pair before confirmation (PLAN §18, 2026-11-02); `models.gat` is kept because the
  development landscape and the F3 positive controls used it.

## Units and statistics — `spatialfloors.stats`
- Unit = patient / donor / animal. Per unit, paired D = metric(learned) − metric(floor), one primary metric per family.
- Model: D ~ 1 + S_c, where S_c is the stratum indicator(s) mean-centred over units, so the intercept equals the grand
  mean of D (`test_stats_intercept_is_mean`). df = N − 2 with two strata, N − 1 with one stratum (no S column)
  (`test_stats_df_rule`, `test_stats_one_stratum`).
- Floor non-inferior: one-sided 95% upper bound of the intercept < margin. Learned superior: one-sided 95% lower bound
  > 0, Holm over the 5 families. Exact t; Wilcoxon signed-rank as sensitivity; no bootstrap; moderators descriptive.
- Tests are per (family, learned model, primary metric) pair (10; 12 with F4's co-primary metrics): floor
  non-inferiority = intersection-union over all pairs; learned superiority = Holm over all pairs (PLAN §11).
- Margin caps: F1 .03, F2 .05, F3 .03, F4 .03, F5 .05. Locked margin = min(cap, 10% of the median development-unit
  above-chance floor score) (`stats.margin_from_dev`, review 2026-10-28; F4 keeps .03), written to manifests/locked.json.
- D is averaged per unit before testing (`stats.aggregate_units`); strata = platform (`stats.stratum`: CosMx WTx, CosMx
  1k, Xenium 5K, MERFISH mouse); df = N − number of strata present (`test_aggregate_units_and_margin`).

## F1 identity — `metrics.f1` (`count_split`, `knn_poisson_score`)
- Floor (v2, 2026-10-29, user-approved; PLAN §14, HYPOTHESES MP12): `floors.simple.expr_embed` = analytic Pearson
  residuals (theta 100, clip sqrt(n)) on X_a panel genes, z-scored (clip 10), PCA 100 (`p3_f1.py <u> expr_embed pearson_z_pc100`).
  v1 (`expr_pca`, log1p-CP10k, PCs {30, 50}) kept as a reported reference.
- X = X_a + X_b, X_a ~ Binomial(X, .5) elementwise, seed 0 (numpy default_rng(0)).
- HVG (scanpy seurat_v3, 2000) and every fit (embeddings, scalers, PCA) on X_a only.
- k = 15 cosine neighbours per embedding, self excluded.
- L_{a,i}, L_{b,i}: total counts of cell i over all filtered genes (data.genes) in X_a, X_b; the split is drawn on the
  full crop, then columns are restricted (scripts/f1_fit.py; P3 runs re-scored by scripts/p3_f1_rescore.py).
- mu_hat_{b,i,g} = L_{b,i} · ( (1−λ) · (1/k) Σ_{j∈N_i} X_{a,j,g}/L_{a,j} + λ · pbar_{a,g} ), λ = .5 (fixed before p0-lock:
  the PLAN's λ = 1/(k+1) gave negative scores at ~32 X_a counts/cell on gse1_232, Expr PCA −.155 vs chance −.466; λ = .5: .074 vs −.096),
  pbar_{a,g} = Σ_i X_{a,i,g} / Σ_i L_{a,i}. Cells with L_{a,j} = 0 contribute a zero profile.
- Dev_Poisson(x, mu) = 2 [x log(x/mu) − (x − mu)], with 0·log 0 = 0.
- Score = 1 − Σ Dev(X_b, mu_hat) / Σ Dev(X_b, L_b · pbar_a), sums over HVGs and cells.
- Genes with pbar_{a,g} = 0 are dropped (never the case for HVGs chosen on X_a).
- Chance level (`chance_level`, rows of E permuted) is reported next to every score. It is NOT 0: a k-cell neighbour
  mean is a noisier rate than pbar, so a random embedding scores < 0 (synthetic: −.05 to −.10). Paired D between
  embeddings on the same cells / k cancels it. (P0 finding; PLAN §10 "shuffled = chance" read as "= permutation chance".)
- Tests: perfect embedding (neighbours = true group) ≫ chance and rises with depth (→ 1 only as λ → 0); random embeddings agree
  with the permutation chance level and are < 0; finite with all-zero genes/cells; self never in N_i.
- Assumption: under NB the two halves of one cell correlate (shared rate); the predictor uses only other cells, so
  this does not leak.
- Nicheformer embedding (`bench.smoke._nf_masked`): mean over non-pad gene tokens, each cell cut to its own length
  (128-token buckets). Non-pad outputs equal the full 1,500-token forward (max |diff| 0.0 on 4 cells); the authors'
  get_embeddings also averages ~1,300 padding positions. P2 dev (gse1_232, 3,000 cells): F1 .006 masked vs −.001
  authors' pooling (chance −.09), median cosine .76; 16 vs 3.2 cells/s. Full crops, no cell subsample.
  Inputs per platform (`models.nicheformer.inputs`): CosMx (Pelka symbol → Ensembl, cosmx mean), Xenium (file Ensembl ids,
  xenium mean, token "Xenium"), MERFISH mouse (Ensembl 116 orthologs → human, merfish mean, species mouse); vocabulary
  share: Allen 95%, SPATCH 99.5%.

## F2 domains — `metrics.f2` (`ari_at_true_k`)
- ARI vs annotated regions at k = number of regions present in the crop. CellCharter: its own GMM
  (cc.tl.Cluster); Novae: assign_domains; floor: full-covariance GMM (`gmm_labels`, v2 2026-10-29; k-means split
  gse1_232's larger elongated region, PLAN §12), grid ell {100, 200} um x lam {.8, 1}, 50 PCs (v3 2026-10-29: centre of the
  PLAN §14 landscape plateau, gmm_full ell 100-400 x lam .8-1 x 30-50 PCs within 10% of the best balanced ARI on both dev units;
  locked config = plateau centre ell 200, lam .8, 50 PCs). 10 clustering seeds (0..9), ARI averaged. Reported with balanced
  ARI (N2 T.ARI_bal) and per-region recall (user decision 2026-10-29, row 359: the Allen small-region loss is shared by
  every method at true k).
- Scored cells (`scored_regions`): the crop's present regions (≥ 2% of the slide's labelled cells), minus non-region
  labels ('n/a', 'unassigned', …); k = number of present regions; every cell is clustered.
- Tests: perfect = 1, random labels ≈ 0.
- Floor (decision 2026-10-28, row 347; `floors.simple.domain_features` defaults): [Expr PCA, all-type neighbour
  mean of it (Gaussian ell µm, 200 NN, same FOV, self excluded)], each block centred and divided by its total SD,
  block weights √(1−λ) / √λ; grid ell {25, 100 µm} × λ {.5, .8}, 30 PCs (review 2026-10-28, PLAN §11).
- Novae: {zero_shot, fine_tune}; novae-human-0 b8c0a5d / novae-mouse-0 3a23be9; spatial_neighbors radius 80 µm;
  assign_domains(n_domains = k). CellCharter: scVI once (seed 0), GMM seeds 0–9. The PLAN's bleed-safe (other-lineage) context is blind to
  composition-only domains (sim.domain_toy ARI ~0 vs .98) and stays a diagnostic.

## F3 niche → receiver state — `metrics.f3` (`quadrant_blocks`, `edge_buffer`, `block_cv`)
- 4 contiguous quadrants (split at the median x and median y); each is held out in turn.
- Graphs / neighbourhoods are built within a block (`restrict_graph_to_blocks`).
- Cells closer than 2 × receptive field to a block edge (the median lines) are not scored.
- Bleed-safe residualisation (floors.bleedsafe: local_foreign_mask, own_states, ctx_bleedsafe); the receiver-state PCA
  (top 10 PCs) and every scaler are fitted on training quadrants only — asserted by index in `block_cv`.
- Score = mean held-out Spearman over the receiver's top-10 state PCs, averaged over the 4 folds.
- Shared design (`floors.f3.prepare`, P2): receivers = Leiden ids on Expr50 of the panel genes with >= 1,000 cells;
  lineages and local foreign mask computed once per crop (label-free); receptive field 50 µm for every model (buffer
  100 µm). Target = receiver's masked genes renormalised on them, FOV-centred; depth (log detected masked genes)
  removed with train-only coefficients, then z-scored on train. Unit score = mean over receivers.
- Floor = RidgeCV (alphas 10^-2..10^4) on [kernel id mix T, bleed-safe context]; kernel Gaussian ell µm over 200 NN,
  same FOV, same quadrant, other lineages; context SVD and gene scale fitted outside the held-out quadrant.
- GAT (`models.gat`, P3): message passing over the floor kernel's top-15 neighbours (same FOV, same quadrant, other
  lineages; no self loops) on [one-hot id, the receiver id's context SVD], zeroed for the receiver's lineage; MLP head
  on [GAT output, floor features] plus a linear map of the floor features (contains its floor). Grid: attention
  {GATConv, GATv2Conv} × layers {1, 2}, hidden 128 (decision 2026-10-28). Kept despite failing the quadratic control
  (row 349, decision 2026-10-28): its control profile (linear .18, product .10, quadratic .03, bleed ~0) is reported.
- Secondary F3 readouts (reported, not tested): PC1 Spearman and top-3 mean (dev: ~4× and ~2× the 10-PC mean's range).
- Toys (tests/test_f3_floor.py): no coupling or bleed only → |score| < .02; planted linear coupling → PC1 ρ > .15
  (unit score only ~.025: one coupled axis of 10 PCs); nonlinear (sender state × third-type presence) → ridge ~0.

## F4 held-out-type counterfactual — `floors.l6`, `metrics.f4`
- Cellina / CPA at the authors' configs (external/cellina-reproducibility@6bb060b, cellina@dd860d0 = v1.2.0,
  cpa@fbd7c02); adata_crc_config (control REF, held-out CRC); 222 uses adata_crc_holdout_config (control REF,
  held-out TVA).
- Slides: authors' whole-slide protocol, except GSE1 110 (695k cells): its ~50k crop for Cellina, CPA and L6 alike
  (P0 decision 2026-10-28, memory).
- Floor L6 = `floors.l6.twoway` (tools/87, verbatim), locked.
- Metrics: official pearson over the top-50 |LFC| ground-truth DEGs (Cellina: log2((m+1e-6)/(c+1e-6))) and spike-free
  pearson (pseudocount .1: log2((m+.1)/(c+.1)), its own top 50) — co-primary. Implemented in `floors.l6.score`.
- `spatialfloors check` reproduces N2/trials/perturb/opt/crc_232.csv and crc_232_lfc.npz (row 323) bit-identically.
- Reporting: retrospective N = 5 (GSE1 120, 210, 221, 242 + Allen animal 2), prospective N = 5 (GSE1 110, 222 +
  Allen animals 1, 3, 4 or the substitute in manifests/allen.json), combined N = 10.

## F5 label transfer — `metrics.f5` (`if_gate`, `macro_f1_gated`)
- Classes: epithelial / immune non-myeloid / myeloid from conservative IF gates (confident tails only; coverage
  reported). GSE1 myeloid = CD45+ ∧ CD68_CK8_18+ ∧ PanCK−. Thresholds and the validation on GSE1 232
  (held-out CD68 / CD163 / C1QC counts) in manifests/if_gate_gse1_232.json; failure → 2-class gate, reported.
- scANVI: joint scVI on a stratified Pelka subsample (≤ 80k cells, ≤ 4k per ClusterMidway) + the un-split crop,
  batch_key = "modality", then SCANVI.from_scvi_model.
- Primary (decision 2026-10-28, row 348): `heldout_auroc` on every F5 unit — per class, AUROC of log1p CP10k of its
  summed held-out markers (on the platform; never inputs; library size over filtered genes) for cells predicted that
  class vs every other cell, "other" included (review 2026-10-28); < 50 cells on either side → .5; mean over 3 classes.
- scANVI / Tangram reference = Pelka raw counts (.raw.X, `floors.f5.reference(counts=True)`), stratified subsample.
- Primary F5 floor = InSituType ML (one config); logistic regression (4 configs) is a secondary floor (PLAN §11). Per-class predicted counts
  reported. Secondary: macro-F1 on IF-gated cells, GSE1 2-class gate only.
- Floors (`floors.f5`, P2): inputs = panel − all F5 held-out markers ∩ Pelka. InSituType ML (profiles = Pelka mean
  linear expression per ClusterMidway, neg = NegPrb mean or .05); logistic regression on the stratified Pelka subsample,
  C {.1, 1} × scaling {own, ref}. ClusterMidway → class: Epi, EpiT → epithelial; Plasma, TCD4, TCD8, B, Tgd, TZBTB16,
  NK, ILC → immune non-myeloid; Macro, Mono, DC, Granulo → myeloid; rest → other (a miss on gated cells).
- GSE303070 IF QC (manifests/if_qc_gse303070.json, model-free vs the authors' labels): no unit has all three channels
  at AUROC >= .75; the 3-class gate is not used (decision 2026-10-28, row 348).

## Gene panel — `data.genes` (P1)
- Remove MT-, RPL/RPS, negative / system probes, unannotated lncRNA / antisense / pseudogenes; seurat_v3 HVG 2000 (on
  X_a where count-split applies); targeted panels keep all filtered genes. No Pelka markers, no protein data in inputs.
- Operational rules (P1, `data.genes`, tests/test_genes.py): annotation = HGNC complete set 2026-10-06 (sha256 8bf6f686…);
  previous / alias symbols resolve when unique. Removed: probes (NEGATIVE*, SystemControl*, Blank*, FalseCode*, …), MT-,
  HGNC groups "L / S ribosomal proteins" (cytoplasmic RPL/RPS; MRPL/MRPS and RPS6K* kept), HGNC pseudogene locus types,
  lncRNAs without a functional name (LINC#, -AS/-DT/-IT/-OT), and symbols absent from HGNC that are clone / Ensembl / LOC
  ids. Mouse (targeted): probes, mt-, Rpl/Rps. Whole-transcriptome = > 10,000 features (GSE1, GSE0) -> HVG 2000 on X_a
  (count split seed 0 of the crop); every other panel keeps all filtered genes. F4 keeps the authors' own HVG step.
- F5 held-out markers (secondary AUROC, WTx; removed from F5 inputs): epithelial EPCAM KRT8 KRT18 KRT19 CDH1;
  immune non-myeloid CD3E CD2 MS4A1 CD79A MZB1; myeloid CD68 CD163 C1QC CD14 LYZ. Panels: manifests/panels/.

## Crops — `data.crops`
- One contiguous ~50k-cell crop per slide; every present annotated region ≥ 20% of the crop's cells, else the crop
  maximising the minimum region share (recorded). Chosen from coordinates + region labels only, before any model runs.
  Manifests (cell IDs + window) in manifests/crops/.

## Reproducibility
- Every output carries the spatialfloors commit hash (`spatialfloors.provenance.stamp`).
- CPU floors are bit-identical across runs in the `floors` env; GPU (MPS) tolerance is documented per model.
