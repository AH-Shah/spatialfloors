"""Load one crop (manifests/crops/<name>.csv.gz) as AnnData: X = raw counts (csr), var_names = gene symbols,
obs = x_um, y_um, region (or NA), IF means where measured. Usage: load_crop("gse1_110")."""
import json
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

SF = Path(__file__).resolve().parents[3]; MAN = SF / "manifests/crops"
RAW = Path("/Users/Hemato/project_3/data/raw"); PROC = Path("/Users/Hemato/project_3/data/processed/p5/labeled_h5ad")
COSMX_UM = 0.12028  # um per CosMx px
ALLEN = {"animal1": RAW / "allen_zhuang/WB_MERFISH_animal1_coronal.h5ad", "animal2": Path("/Users/Hemato/project_3/antig/N2/trials/perturb/WB_MERFISH_animal2_coronal.h5ad"),
         "animal3": RAW / "allen_zhuang/WB_MERFISH_animal3_sagittal.h5ad", "animal4": RAW / "allen_zhuang/WB_MERFISH_animal4_sagittal.h5ad"}  # 4 = mouse dev unit (review 2026-10-28)
YAO = RAW / "allen_yao_638850"
BRUKER = {"bruker_breast": RAW / "bruker_wtx/breast_discovery.h5ad", "bruker_mobreast": RAW / "bruker_wtx/multiomic_breast.h5ad",
          "bruker_skin": RAW / "bruker_wtx/skin_scc.h5ad", "bruker_kidney": RAW / "bruker_wtx/kidney_rcc.h5ad"}  # crop h5ads (scripts/v_bruker_h5ad.py)


BRUKER_SRC = {"bruker_breast": RAW / "bruker_wtx/breast_discovery", "bruker_mobreast": RAW / "bruker_wtx/multiomic_breast/flat",
              "bruker_skin": RAW / "bruker_wtx/skin_scc", "bruker_kidney": RAW / "bruker_wtx/kidney_rcc"}


def bruker_obs(u):
    """Cells of a Bruker WTx slide from its metadata flat file (index = cell id c_<slide>_<fov>_<cell_ID>; x_px, y_px global
    CosMx px; fov; cell_ID). Used to choose the crop before any expression is read (scripts/v_crops.py)."""
    mf = next(BRUKER_SRC[u].glob("*metadata*.csv*")); m = pd.read_csv(mf, usecols=lambda c: c in ("fov", "cell_ID", "cell", "cell_id", "CenterX_global_px", "CenterY_global_px"))
    idx = (m["cell"] if "cell" in m else m["cell_id"]).astype(str).values
    return pd.DataFrame(dict(x_px=m.CenterX_global_px.values.astype(float), y_px=m.CenterY_global_px.values.astype(float), fov=m.fov.astype(str).values, cell_ID=m.cell_ID.values), index=idx)


def index():
    return {k: v for k, v in json.loads((MAN / "index.json").read_text()).items() if isinstance(v, dict)}


def _take(path, ids, key=None, raw=False):
    """Rows of a backed h5ad by id (obs_names or obs[key]), in manifest order; counts as csr float32."""
    a = ad.read_h5ad(path, backed="r"); names = pd.Index(a.obs[key].astype(str).values if key else a.obs_names.astype(str))
    pos = names.get_indexer(pd.Index(ids)); assert (pos >= 0).all(), f"{path}: {(pos < 0).sum()} manifest ids missing"
    o = np.argsort(pos); rows = pos[o]; X = (a.raw.X if raw else a.X)[rows]; X = sp.csr_matrix(X, dtype=np.float32)
    inv = np.empty_like(o); inv[o] = np.arange(len(o)); X = X[inv]; obs = a.obs.iloc[pos].copy(); obsm = {k: np.asarray(a.obsm[k])[rows][inv] for k in ("spatial", "X_spatial_coords") if k in a.obsm}
    var = (a.raw.var if raw else a.var).copy(); a.file.close(); return X, obs, var, obsm


def load_crop(name):
    rec = index()[name]; ids = pd.read_csv(MAN / f"{name}.csv.gz").cell_id.astype(str).values
    if name.startswith("gse1_"):
        X, obs, var, _ = _take(RAW / f"GSE1/cosmx_wt_{name[5:].rstrip('w')}.h5ad", ids); genes = var.index.astype(str)  # <sid>w = whole slide (crop check)
        xy = obs[["CenterX_global_px", "CenterY_global_px"]].values * COSMX_UM; region = obs["typ"].astype(str).str.extract(r"(REF|TVA|CRC)", expand=False).values
        ifc = {c: obs[f"Mean.{c}"].values for c in ("PanCK", "CD45", "CD68_CK8_18") if f"Mean.{c}" in obs}
    elif name == "gse0_S0":
        X, obs, var, _ = _take(PROC / "GSE0_S0.h5ad", ids, key="cell_id"); genes = var.index.astype(str)
        xy = obs[["CenterX_global_px", "CenterY_global_px"]].values * COSMX_UM; region = np.full(len(ids), None); ifc = {}
    elif name.startswith("gse303070_"):
        slide = "_".join(name.split("_")[1:-1])  # gse303070_<slide>_<sample>, slide may be Tumour_A
        f = PROC / f"GSE303070_CosMx_Slide_{slide.replace('Tumour', 'Tumor')}.h5ad"
        key = [f"{b.split('_')[-1]}_{b.split('_')[-2]}" for b in ids]  # Barcode <slide>_<cell_ID>_<fov> -> h5ad <fov>_<cell_ID>
        X, obs, var, _ = _take(f, key); genes = var.index.astype(str)
        M = pd.read_csv(RAW / "GSE303070/GSE303070_All_cells_metadata.csv", index_col=False, usecols=["Barcode", "nCount_RNA", "CenterX_global_px", "CenterY_global_px", "Mean.PanCK", "Mean.CD45", "Mean.CD68", "Mean.CD3"]).set_index("Barcode").loc[ids]
        tot = np.asarray(X.sum(1)).ravel()  # metadata nCount_RNA = total over all h5ad features (probes included)
        assert np.mean(np.abs(tot - M.nCount_RNA.values) <= .5) > .99, f"{name}: h5ad counts do not match the metadata nCount_RNA"
        xy = M[["CenterX_global_px", "CenterY_global_px"]].values * COSMX_UM; region = np.full(len(ids), None)
        ifc = {c: M[f"Mean.{c}"].values for c in ("PanCK", "CD45", "CD68", "CD3")}
    elif name in BRUKER:
        X, obs, var, _ = _take(BRUKER[name], ids); genes = var.index.astype(str)
        xy = obs[["x_px", "y_px"]].values * COSMX_UM; region = np.full(len(ids), None)
        ifc = {c: obs[f"Mean.{c}"].values for c in ("PanCK", "CD45", "CD68", "CD68_CK8_18", "CD3") if f"Mean.{c}" in obs}
    elif name == "spatch_xenium5k":
        X, obs, var, obsm = _take(RAW / "spatch_xenium_5k/spatch_xenium_5k.h5ad", ids); genes = var.index.astype(str)
        xy = obsm["spatial"]; region = np.full(len(ids), None); ifc = {}  # Xenium coordinates are in um
    elif name.startswith("allen_animal") or name.startswith("allen2_"):
        an = "animal2" if name.startswith("allen2_") else name.split("_")[1]
        X, obs, var, obsm = _take(ALLEN[an], ids, raw=True); genes = var["feature_name"].astype(str).values
        xy = obsm["X_spatial_coords"][:, :2]; region = obs["major_brain_region"].astype(str).values; ifc = {}  # already um (span ~5,000 um per section)
    elif name.startswith("allen_yao638850_"):
        X, obs, var, _ = _take(YAO / "C57BL6J-638850-raw.h5ad", ids); genes = var["gene_symbol"].astype(str).values if "gene_symbol" in var else var.index.astype(str)
        m = pd.read_csv(YAO / "cell_metadata_with_parcellation_annotation.csv", usecols=["cell_label", "x_section", "y_section", "parcellation_division", "parcellation_category"]).set_index("cell_label").loc[ids]
        xy = m[["x_section", "y_section"]].values * 1000
        region = np.select([m.parcellation_division == "TH", m.parcellation_division == "Isocortex", m.parcellation_category == "fiber tracts"], ["Thalamus", "Isocortex", "Fiber_tracts"], default=m.parcellation_division.astype(str)); ifc = {}
    else:
        raise KeyError(name)
    assert np.allclose(X.data, np.round(X.data)), f"{name}: X is not integer counts"
    fov = obs["fov"].astype(str).values if "fov" in obs and name.startswith(("gse1_", "gse0_", "gse303070_", "bruker_")) else np.full(len(ids), "0")
    o = pd.DataFrame(dict(x_um=xy[:, 0].astype(float), y_um=xy[:, 1].astype(float), region=region, fov=fov), index=ids)
    for c, v in ifc.items(): o[f"IF_{c}"] = np.asarray(v, float)
    a = ad.AnnData(X, obs=o, var=pd.DataFrame(index=pd.Index(genes, name="symbol")))
    a.var_names_make_unique(); a.uns["unit"] = dict(name=name, **{k: v for k, v in rec.items() if isinstance(v, (str, int, float, bool)) or v is None}); return a
