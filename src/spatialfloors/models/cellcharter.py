"""F2 learned model: CellCharter: scVI on counts -> Delaunay + remove_long_links ->
aggregate_neighbors(n_layers) -> cc.tl.Cluster (GMM) at the true k. No harmonising PCA, no flavor="seurat" HVG step.
Grid: n_layers {1, 3} x scVI n_latent {10, 30}. fit(rep="pca"): the published PCA-input variant."""

import numpy as np
import scipy.sparse as sp


def fit(X, xy_um, k, n_latent=10, n_layers=3, accelerator="mps", seed=0, max_epochs=None, gmm_seeds=None, rep="scvi"):
    """X counts (cells x panel genes). scVI once (embedding seed `seed`), then the GMM per clustering seed (PROTOCOL §2:
    10
    seeds; gmm_seeds=None -> [seed]). Returns (labels, info), labels = list of arrays when gmm_seeds is given."""
    import anndata as ad, scvi, squidpy as sq, cellcharter as cc

    a = ad.AnnData(sp.csr_matrix(X, dtype=np.float32))
    a.obsm["spatial"] = np.asarray(xy_um, float)
    scvi.settings.seed = seed
    m = None
    if rep == "pca":
        from spatialfloors.floors.simple import expr_pca

        a.obsm["X_scVI"] = expr_pca(X, 50, seed)  # key kept so the CellCharter steps below are shared
    else:
        scvi.model.SCVI.setup_anndata(a)
        m = scvi.model.SCVI(a, n_latent=n_latent)
        m.train(max_epochs=max_epochs, accelerator=accelerator, devices=1, enable_progress_bar=False)
        a.obsm["X_scVI"] = m.get_latent_representation()
    sq.gr.spatial_neighbors(a, coord_type="generic", delaunay=True, spatial_key="spatial")
    cc.gr.remove_long_links(a)
    cc.gr.aggregate_neighbors(a, n_layers=n_layers, use_rep="X_scVI", out_key="X_cellcharter")
    labs = []
    for s in gmm_seeds if gmm_seeds is not None else [seed]:
        g = cc.tl.Cluster(
            n_clusters=k, random_state=s, trainer_params=dict(accelerator="cpu", enable_progress_bar=False)
        )
        g.fit(a, use_rep="X_cellcharter")
        labs.append(np.asarray(g.predict(a, use_rep="X_cellcharter")))
    info = dict(
        epochs=len(m.history["elbo_train"]) if m is not None else 0, dim=int(a.obsm["X_cellcharter"].shape[1]), rep=rep
    )
    return (labs if gmm_seeds is not None else labs[0]), info
