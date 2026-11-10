"""Whole-slide arm (results/WHOLESLIDE_PLAN.md, exploratory): lean cache, panel and manifest of one whole GSE1 slide,
streamed from the on-disk CSC raw file (data.units.load_crop of a 741k-cell slide needs ~23 GB; this stays near 6 GB).
Same definitions as the registered pipeline:
- cells: every cell of the slide, in on-disk order (= load_crop of a whole-slide manifest, as gse1_231w);
- X_a = Binomial(X, .5) with numpy default_rng(0) over the slide's non-zeros in CSR (row-major) order
  (= metrics.f1.count_split of that load_crop matrix); X_b = X - X_a;
- panel = data.genes.panel on X_a (filter, seurat_v3 HVG 2,000; scanpy 1.11.5 code path as floors/l6_lean.py);
- obs as load_crop (x_um, y_um, region from typ, fov, IF means).
Cache (runs/p2/cache/gse1_<sid>w.h5ad): X = counts on panel + F5 held-out markers on the platform + negative probes
+ SF_REST_FILTERED (counts of every other filtered gene, so library sizes over filtered genes are unchanged; the column
is never a feature: every model reads the panel). runs/p2/cache/gse1_<sid>w_split.npz: X_a on the same columns.
--check (sid 231): rebuild gse1_231w this way under the name gse1_231wl and compare with runs/p2/cache/gse1_231w.h5ad
(count_split restricted to the same columns, rest column) and manifests/panels/gse1_231w.json; exit 1 on any mismatch.
Usage (env floors): python scripts/xw_prepare.py <sid> [--check]"""

import json, sys, time
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
from anndata.io import read_elem
from anndata.utils import make_index_unique
from scanpy.preprocessing._highly_variable_genes import _sum_and_sum_squares_clipped
from scanpy.preprocessing._utils import _get_mean_var
from spatialfloors.data.genes import filter_reason, N_HVG, HGNC, HGNC_SHA
from spatialfloors.data.units import RAW, COSMX_UM
from spatialfloors.provenance import stamp

sid = sys.argv[1]
CHECK = "--check" in sys.argv
u = f"gse1_{sid}w" + ("l" if CHECK else "")
GCH = 500
BCH = 50_000_000
F5_HELDOUT = json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"]
REST = "SF_REST_FILTERED"
assert filter_reason([REST])[0] == "", "the rest column must count as a filtered gene"
t0 = time.time()
log = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)
f = h5py.File(RAW / f"GSE1/cosmx_wt_{sid}.h5ad", "r")
Xh = f["X"]
assert Xh.attrs["encoding-type"] == "csc_matrix"
n, G = map(int, Xh.attrs["shape"])
ip = Xh["indptr"][:].astype(np.int64)
obs = read_elem(f["obs"])
var = read_elem(f["var"])
genes = make_index_unique(pd.Index(var.index.astype(str))).values.astype(
    str
)  # = AnnData.var_names_make_unique in load_crop
reason = filter_reason(genes, "human")
filt = reason == ""
log(u, n, "cells", G, "genes", ip[-1], "non-zeros")


def chunks():
    """(g0, g1, rows int64, vals float64) per window of GCH on-disk columns."""
    for g0 in range(0, G, GCH):
        g1 = min(G, g0 + GCH)
        lo, hi = ip[g0], ip[g1]
        yield g0, g1, Xh["indices"][lo:hi].astype(np.int64), Xh["data"][lo:hi]


def positions(rows, cnt, rowptr):
    """CSR position of each element of a column window (elements in column-major order); advances cnt."""
    o = np.argsort(rows, kind="stable")
    rs = rows[o]
    st = np.r_[0, np.flatnonzero(np.diff(rs)) + 1]
    rk = np.empty(len(rs), np.int64)
    rk[o] = np.arange(len(rs)) - np.repeat(st, np.diff(np.r_[st, len(rs)]))
    pos = rowptr[rows] + cnt[rows] + rk
    cnt += np.bincount(rows, minlength=n)
    return pos


# pass A: non-zeros per row -> CSR row pointer
rnz = np.zeros(n, np.int64)
for _, _, r, _ in chunks():
    rnz += np.bincount(r, minlength=n)
rowptr = np.r_[0, np.cumsum(rnz)]
nnz = int(rowptr[-1])
assert nnz == ip[-1]
log("pass A")
# pass B: values in CSR order (uint16), then the binomial split in that order
V = np.zeros(nnz, np.uint16)
cnt = np.zeros(n, np.int64)
for _, _, r, v in chunks():
    assert np.all(v == np.round(v)) and v.max(initial=0) < 65535 and v.min(initial=0) >= 0
    V[positions(r, cnt, rowptr)] = v
log("pass B")
rng = np.random.default_rng(0)
A = np.empty(nnz, np.uint16)
for s in range(0, nnz, BCH):
    A[s : s + BCH] = rng.binomial(V[s : s + BCH].astype(np.int64), 0.5)
del V
log("split")
# pass C1: per-gene mean / var of X_a (filtered genes, float32 as panel), library sizes over filtered genes
fi = np.flatnonzero(filt)
col = -np.ones(G, np.int64)
col[fi] = np.arange(len(fi))
mean = np.zeros(len(fi))
varr = np.zeros(len(fi))
totX = np.zeros(n)
totA = np.zeros(n)
cnt[:] = 0
blocks = {}


def window(g0, g1, r, v):
    a = A[positions(r, cnt, rowptr)].astype(np.float32)
    ptr = ip[g0 : g1 + 1] - ip[g0]
    return sp.csc_matrix((v.astype(np.float32), r, ptr), shape=(n, g1 - g0)), sp.csc_matrix(
        (a, r, ptr), shape=(n, g1 - g0)
    )


for g0, g1, r, v in chunks():
    Xw, Aw = window(g0, g1, r, v)
    k = np.flatnonzero(filt[g0:g1])
    if not len(k):
        continue
    Bw = sp.csr_matrix(Aw[:, k])
    m, vv = _get_mean_var(Bw)
    mean[col[g0 + k]], varr[col[g0 + k]] = m, vv
    totX += np.asarray(Xw[:, k].sum(1)).ravel()
    totA += np.asarray(Aw[:, k].sum(1)).ravel()
log("pass C1")
from skmisc.loess import loess  # scanpy _highly_variable_genes_seurat_v3, one batch

nc = varr > 0
est = np.zeros(len(fi))
mdl = loess(np.log10(mean[nc]), np.log10(varr[nc]), span=0.3, degree=2)
mdl.fit()
est[nc] = mdl.outputs.fitted_values
reg_std = np.sqrt(10**est)
clip_val = reg_std * np.sqrt(n) + mean
sq, sm = np.zeros(len(fi)), np.zeros(len(fi))
cnt[:] = 0
for g0, g1, r, v in chunks():  # pass C2: clipped sums
    _, Aw = window(g0, g1, r, v)
    k = np.flatnonzero(filt[g0:g1])
    if not len(k):
        continue
    Bw = sp.csr_matrix(Aw[:, k])
    a_, b_ = _sum_and_sum_squares_clipped(
        Bw.indices, Bw.data, n_cols=Bw.shape[1], clip_val=clip_val[col[g0 + k]], nnz=Bw.nnz
    )
    sq[col[g0 + k]], sm[col[g0 + k]] = a_, b_
ngv = ((1 / ((n - 1) * np.square(reg_std))) * (n * np.square(mean) + sq - 2 * sm * mean)).reshape(1, -1)
rk = np.argsort(np.argsort(-ngv, axis=1), axis=1).astype(np.float32)
nb_ = np.sum((rk < N_HVG).astype(int), axis=0)
rk[rk >= N_HVG] = np.nan
med = np.ma.median(np.ma.masked_invalid(rk), axis=0).filled(np.nan)
top = (
    pd.DataFrame(dict(r=med, b=nb_)).sort_values(["r", "b"], ascending=[True, False], na_position="last").index[:N_HVG]
)
hv = np.zeros(len(fi), bool)
hv[np.asarray(top)] = True
panel_idx = fi[hv]
log("pass C2: panel", len(panel_idx))
# pass D: cache columns
mk = sorted({g for v in F5_HELDOUT.values() for g in v} & set(genes))
probes = [i for i, g in enumerate(genes) if g.startswith(("NegPrb", "Negative"))]
K = np.unique(np.r_[panel_idx, [int(np.flatnonzero(genes == g)[0]) for g in mk], probes]).astype(np.int64)
inK = np.zeros(G, bool)
inK[K] = True
cX, cA = [], []
kX = np.zeros(n)
kA = np.zeros(n)
cnt[:] = 0
for g0, g1, r, v in chunks():
    Xw, Aw = window(g0, g1, r, v)
    k = np.flatnonzero(inK[g0:g1])
    if not len(k):
        continue
    cX.append(Xw[:, k])
    cA.append(Aw[:, k])
    kf = k[filt[g0 + k]]
    kX += np.asarray(Xw[:, kf].sum(1)).ravel()
    kA += np.asarray(Aw[:, kf].sum(1)).ravel()
restX, restA = totX - kX, totA - kA
assert restX.min() >= 0 and restA.min() >= 0
del A
X = sp.hstack(cX + [sp.csc_matrix(restX[:, None].astype(np.float32))], format="csc", dtype=np.float32)
del cX
X = X.tocsr()
Xa = sp.hstack(cA + [sp.csc_matrix(restA[:, None].astype(np.float32))], format="csc", dtype=np.float32)
del cA
Xa = Xa.tocsr()
vn = np.r_[genes[K], [REST]]
log("pass D: cache", X.shape, X.nnz, "non-zeros")
reg = obs["typ"].astype(str).str.extract(r"(REF|TVA|CRC)", expand=False)
o = pd.DataFrame(
    dict(
        x_um=obs["CenterX_global_px"].values.astype(float) * COSMX_UM,
        y_um=obs["CenterY_global_px"].values.astype(float) * COSMX_UM,
        region=reg.values,
        fov=obs["fov"].astype(str).values,
    ),
    index=obs.index.astype(str),
)
for c in ("PanCK", "CD45", "CD68_CK8_18"):
    if f"Mean.{c}" in obs:
        o[f"IF_{c}"] = obs[f"Mean.{c}"].values.astype(float)
f.close()
o["region"] = o.region.astype("category")
lab = reg.dropna()
shares = lab.value_counts(normalize=True)
rec = dict(
    n_slide=n,
    n_crop=n,
    whole_slide=True,
    slide_shares=shares.round(4).to_dict(),
    present=sorted(shares.index[shares >= 0.02]),
    dataset="GSE1",
    role="exploratory_wholeslide",
    region_key="typ",
    note="whole-slide arm (results/WHOLESLIDE_PLAN.md); lean cache scripts/xw_prepare.py",
)
pan = dict(
    mode="wtx_hvg",
    genes=genes[panel_idx].tolist(),
    removed={k_: int(v_) for k_, v_ in pd.Series(reason[reason != ""]).value_counts().items()},
    n_input=G,
    f5_heldout_in_panel={c: sorted(set(genes[panel_idx]) & set(m_)) for c, m_ in F5_HELDOUT.items()},
    f5_heldout_on_platform={c: sorted(set(genes) & set(m_)) for c, m_ in F5_HELDOUT.items()},
    unit=u,
    species="human",
    n_cells=n,
    n_panel=len(panel_idx),
    annotation=dict(hgnc=str(HGNC), sha256=HGNC_SHA),
    counts_per_cell_Xa=round(float(totA.sum() / n), 1),
    lean_cache=dict(columns=len(vn), rest=REST),
)
if CHECK:  # compare with the registered whole-slide path on gse1_231w (cache built by load_crop, panel by p1_panels)
    from spatialfloors.metrics.f1 import count_split

    ref = ad.read_h5ad("runs/p2/cache/gse1_231w.h5ad")
    assert (ref.obs_names.values == o.index.values).all(), "cell order"
    RA, _ = count_split(ref.X, seed=0)
    rv = ref.var_names.values.astype(str)
    assert (rv == genes).all(), "gene names"
    ok = {
        "panel": pan["genes"] == json.loads(Path("manifests/panels/gse1_231w.json").read_text())["genes"],
        "X": (sp.csr_matrix(ref.X[:, K]) != X[:, :-1]).nnz == 0,
        "Xa": (sp.csr_matrix(RA[:, K]) != Xa[:, :-1]).nnz == 0,
    }
    rf = filt & ~inK
    ok["rest_X"] = np.array_equal(np.asarray(ref.X[:, rf].sum(1)).ravel(), restX)
    ok["rest_Xa"] = np.array_equal(np.asarray(RA[:, rf].sum(1)).ravel(), restA)
    ok["obs"] = (
        all(
            np.allclose(ref.obs[c].values.astype(float), o[c].values.astype(float), equal_nan=True)
            for c in ("x_um", "y_um", "IF_PanCK", "IF_CD68_CK8_18")
        )
        and (ref.obs.fov.astype(str).values == o.fov.values).all()
        and (ref.obs.region.astype(str).values == o.region.astype(str).values).all()
    )
    log("CHECK", ok)
    del ref, RA
    if not all(ok.values()):
        sys.exit(1)
a = ad.AnnData(X, obs=o, var=pd.DataFrame(index=pd.Index(vn, name="symbol")))
a.uns["unit"] = dict(name=u, **{k_: v_ for k_, v_ in rec.items() if isinstance(v_, (str, int, float, bool))})
a.write_h5ad(f"runs/p2/cache/{u}.h5ad", compression="gzip")
sp.save_npz(f"runs/p2/cache/{u}_split.npz", Xa)
pd.Series(o.index.values, name="cell_id").to_csv(f"manifests/crops/{u}.csv.gz", index=False)
for p_, key, val in (
    (Path("manifests/crops/index.json"), u, rec),
    (
        Path("manifests/panels/index.json"),
        u,
        {k_: pan[k_] for k_ in ("mode", "n_input", "n_panel", "removed", "n_cells", "counts_per_cell_Xa")},
    ),
):
    d_ = json.loads(p_.read_text())
    d_[u] = stamp(val)
    p_.write_text(json.dumps(d_, indent=1, default=float))
Path(f"manifests/panels/{u}.json").write_text(json.dumps(stamp(pan), indent=1))
log("written", u, "present", rec["present"], "XW_PREPARE_DONE")
