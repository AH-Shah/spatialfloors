"""F5 learned models (PLAN §5): scANVI (joint scVI on reference + un-split query, batch_key="modality", then
SCANVI.from_scvi_model; query unlabelled) and Tangram mode="clusters" (reference ClusterMidway mean profiles mapped to
query cells; label = argmax of the cell's mapping probability). Inputs: integer counts on the shared genes (reference =
floors.f5.reference(..., counts=True), Pelka .raw.X; query = the un-split crop).
Grids: scANVI inherits scVI's (n_latent, likelihood); Tangram density prior {uniform, rna_count_based}."""
import numpy as np
import scipy.sparse as sp


def scanvi(Xr, yr, Xq, n_latent=10, likelihood="nb", accelerator="mps", seed=0, max_epochs=None, scanvi_epochs=None):
    import anndata as ad, scvi, pandas as pd
    a = ad.AnnData(sp.vstack([sp.csr_matrix(Xr), sp.csr_matrix(Xq)]).tocsr().astype(np.float32))
    a.obs["modality"] = ["ref"] * Xr.shape[0] + ["query"] * Xq.shape[0]; a.obs["label"] = pd.Categorical(list(map(str, yr)) + ["Unknown"] * Xq.shape[0])
    scvi.settings.seed = seed; scvi.model.SCVI.setup_anndata(a, batch_key="modality", labels_key="label")
    m = scvi.model.SCVI(a, n_latent=n_latent, gene_likelihood=likelihood); m.train(max_epochs=max_epochs, accelerator=accelerator, devices=1, enable_progress_bar=False)
    s = scvi.model.SCANVI.from_scvi_model(m, unlabeled_category="Unknown"); s.train(max_epochs=scanvi_epochs, accelerator=accelerator, devices=1, enable_progress_bar=False)
    return np.asarray(s.predict(a[Xr.shape[0]:]))


def tangram(Xr, yr, Xq, density="rna_count_based", device="mps", seed=0, num_epochs=1000):
    import anndata as ad, pandas as pd, tangram as tg
    genes = [f"g{i}" for i in range(Xr.shape[1])]
    r = ad.AnnData(sp.csr_matrix(Xr, dtype=np.float32), var=pd.DataFrame(index=genes)); r.obs["label"] = pd.Categorical(list(map(str, yr)))
    q = ad.AnnData(sp.csr_matrix(Xq, dtype=np.float32), var=pd.DataFrame(index=genes))
    tg.pp_adatas(r, q, genes=genes)
    m = tg.map_cells_to_space(r, q, mode="clusters", cluster_label="label", density_prior=density, num_epochs=num_epochs, device=device, random_state=seed, verbose=False)
    P = m.X; names = np.asarray(m.obs["label"].astype(str)); return names[np.asarray(P).argmax(0)]
