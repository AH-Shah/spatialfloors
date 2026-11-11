"""Replication round (registration v2): load one section of the selected datasets (registration_v2/SELECTION.md) in the
round-1 cache layout, so the round-1 fit and score scripts run unchanged on it.
AnnData: X = raw counts (csr float32) over the genes and the negative-control features (the gene filter,
data.genes.filter_reason, removes the controls; F5 reads them as negatives); var index = symbol, var.feature_type,
var.gene_id (Ensembl id where the producer gives one); obs = x_um, y_um, region (None until annotations exist), fov,
IF_<marker> means where measured; uns["unit"] = section, unit, source, platform.
  Xenium Prime 5K (10x output bundle): cell_feature_matrix.h5 (features x cells, CSC), cells.parquet (centroids, µm).
    One FOV per section ("0"), as for the round-1 Xenium unit.
  CosMx NSCLC (flat files from scripts/v2_nsclc_extract.py): <S>_exprMat_file.csv (rows = fov, cell_ID; cell_ID 0 =
    transcripts outside cells, dropped), <S>_metadata_file.csv (global centroids in pixels, 0.18 µm per pixel per the
    dataset README, IF means PanCK, CD45, CD3). Cell id = <section>_<fov>_<cell_ID>.
obs(section) reads coordinates only, so a crop can be chosen before any expression is read."""

from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

from spatialfloors.paths import RAW

NSCLC_UM_PER_PX = 0.18  # SMI Data File ReadMe v1 (NanoString, 2021-10-04): pixel edge 180 nm
XENIUM_KEEP = ("Gene Expression", "Negative Control Probe", "Negative Control Codeword")
NSCLC_IF = ("PanCK", "CD45", "CD3")
_CV = RAW / "cervilla2026_xenium5k"
_TX = RAW / "tenx_xenium5k"
_NS = RAW / "cosmx_nsclc"
# section -> (platform, directory, unit = patient, source, tissue, cancer)
SECTIONS = {
    "cervilla_breast": ("xenium", _CV / "xenium5K_Breast", "cervilla_breast", "Cervilla 2026", "breast", "carcinoma"),
    "cervilla_crc": ("xenium", _CV / "xenium5K_Colorectal", "cervilla_crc", "Cervilla 2026", "colon", "carcinoma"),
    "cervilla_lung": ("xenium", _CV / "xenium5K_Lung", "cervilla_lung", "Cervilla 2026", "lung", "carcinoma"),
    "cervilla_lymphoma": (
        "xenium",
        _CV / "xenium5K_Lymphoma",
        "cervilla_lymphoma",
        "Cervilla 2026",
        "lymph node",
        "lymphoma",
    ),  # fmt: skip
    "cervilla_ovarian": ("xenium", _CV / "xenium5K_Ovarian", "cervilla_ovarian", "Cervilla 2026", "ovary", "carcinoma"),
    "tenx_skin": ("xenium", _TX / "Xenium_Prime_Human_Skin_FFPE", "tenx_skin", "10x Genomics", "skin", "melanoma"),
    "tenx_prostate": (
        "xenium",
        _TX / "Xenium_Prime_Human_Prostate_FFPE",
        "tenx_prostate",
        "10x Genomics",
        "prostate",
        "carcinoma",
    ),  # fmt: skip
    "tenx_ovarian": (
        "xenium",
        _TX / "Xenium_Prime_Ovarian_Cancer_FFPE_XRrun",
        "tenx_ovarian",
        "10x Genomics",
        "ovary",
        "carcinoma",
    ),  # fmt: skip
    "tenx_breast": (
        "xenium",
        _TX / "Xenium_Prime_Breast_Cancer_FFPE",
        "tenx_breast",
        "10x Genomics",
        "breast",
        "carcinoma",
    ),  # fmt: skip
    "tenx_lung": (
        "xenium",
        _TX / "Xenium_Prime_Human_Lung_Cancer_FFPE",
        "tenx_lung",
        "10x Genomics",
        "lung",
        "carcinoma",
    ),  # fmt: skip
} | {
    f"nsclc_{s.lower()}": (
        "nsclc",
        _NS / s,
        f"nsclc_{s.split('_')[0].lower()}",
        "He 2022 (NanoString NSCLC)",
        "lung",
        "carcinoma",
    )  # fmt: skip
    for s in ("Lung5_Rep1", "Lung5_Rep2", "Lung5_Rep3", "Lung6", "Lung9_Rep1", "Lung9_Rep2", "Lung12", "Lung13")
}
PLATFORM = dict(xenium="Xenium Prime 5K", nsclc="CosMx 960-plex")


def _nsclc_files(d):
    return {k: next(d.glob(f"*_{k}_file.csv")) for k in ("exprMat", "metadata")}


def obs(section):
    """Cells of a section: index = cell id; x_um, y_um, fov. Coordinates only (no expression is read)."""
    kind, d = SECTIONS[section][:2]
    if kind == "xenium":
        c = pd.read_parquet(d / "cells.parquet", columns=["cell_id", "x_centroid", "y_centroid"])
        return pd.DataFrame(
            dict(x_um=c.x_centroid.values.astype(float), y_um=c.y_centroid.values.astype(float), fov="0"),
            index=pd.Index(c.cell_id.astype(str).values, name="cell_id"),
        )
    m = pd.read_csv(_nsclc_files(d)["metadata"], usecols=["fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"])
    return pd.DataFrame(
        dict(
            x_um=m.CenterX_global_px.values * NSCLC_UM_PER_PX,
            y_um=m.CenterY_global_px.values * NSCLC_UM_PER_PX,
            fov=m.fov.astype(str).values,
        ),
        index=pd.Index([f"{section}_{f}_{c}" for f, c in zip(m.fov, m.cell_ID)], name="cell_id"),
    )


def _xenium(d):
    import h5py

    with h5py.File(d / "cell_feature_matrix.h5", "r") as f:
        m = f["matrix"]
        nf, nc = (int(x) for x in m["shape"][:])
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=(nf, nc)).T.tocsr()
        cells = [x.decode() for x in m["barcodes"][:]]
        var = pd.DataFrame(
            dict(
                feature_type=[x.decode() for x in m["features/feature_type"][:]],
                gene_id=[x.decode() for x in m["features/id"][:]],
            ),
            index=pd.Index([x.decode() for x in m["features/name"][:]], name="symbol"),
        )
    keep = np.flatnonzero(var.feature_type.isin(XENIUM_KEEP).values)
    return X[:, keep], pd.Index(cells), var.iloc[keep].copy(), {}


def _nsclc(section, d):
    f = _nsclc_files(d)
    e = pd.read_csv(f["exprMat"])
    e = e[e.cell_ID != 0]
    ids = pd.Index([f"{section}_{a}_{b}" for a, b in zip(e.fov, e.cell_ID)])
    feat = [c for c in e.columns if c not in ("fov", "cell_ID")]
    X = sp.csr_matrix(e[feat].values)
    var = pd.DataFrame(
        dict(feature_type=["Negative Control Probe" if c.startswith("NegPrb") else "Gene Expression" for c in feat]),
        index=pd.Index(feat, name="symbol"),
    )
    var["gene_id"] = ""
    m = pd.read_csv(f["metadata"], usecols=["fov", "cell_ID"] + [f"Mean.{c}" for c in NSCLC_IF])
    m.index = pd.Index([f"{section}_{a}_{b}" for a, b in zip(m.fov, m.cell_ID)])
    return X, ids, var, {c: m[f"Mean.{c}"] for c in NSCLC_IF}


def load(section, ids=None):
    """One section (or the cells `ids` of it, in that order) in the round-1 cache layout."""
    kind, d, unit, source, tissue, cancer = SECTIONS[section]
    X, cells, var, ifc = _xenium(d) if kind == "xenium" else _nsclc(section, d)
    o = obs(section)
    ids = o.index if ids is None else pd.Index(ids).astype(str)
    pos = cells.get_indexer(ids)
    assert (pos >= 0).all(), f"{section}: {(pos < 0).sum()} cells without counts"
    assert o.index.get_indexer(ids).min() >= 0, f"{section}: cells without coordinates"
    X = sp.csr_matrix(X[pos], dtype=np.float32)
    assert np.allclose(X.data, np.round(X.data)), f"{section}: X is not integer counts"
    ob = o.loc[ids, ["x_um", "y_um"]].copy()
    ob["region"] = None
    ob["fov"] = o.loc[ids, "fov"].values
    for c, v in ifc.items():
        ob[f"IF_{c}"] = np.asarray(v.reindex(ids).values, float)
    a = ad.AnnData(X, obs=ob, var=var)
    a.var_names_make_unique()
    a.uns["unit"] = dict(name=section, unit=unit, source=source, platform=PLATFORM[kind], tissue=tissue, cancer=cancer)
    return a
