"""Memory-lean twin of l6.load_crc for whole GSE1 slides (execution-only, EXECUTION_LOG 2026-11-06): the same obs filter,
sc.pp.filter_cells / filter_genes (min_counts 3) and seurat_v3 HVG (2000, scanpy 1.11.5 code path) computed over gene
chunks of the on-disk csc X, then only the HVG columns of the kept cells are loaded. l6.load_crc on slide 120 (741k cells,
969M non-zeros, float64) peaks at ~24 GB; this peaks at a few GB. Identity with l6.load_crc: scripts/l6_lean_check.py."""
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
from anndata.io import read_elem
from scanpy.preprocessing._highly_variable_genes import _sum_and_sum_squares_clipped
from scanpy.preprocessing._utils import _get_mean_var
from . import l6, _labels_to_coarse as lm

CHUNK = 1000  # genes per chunk


def _chunks(X, rows, genes):
    """(offset into `genes`, csr float64 block rows x genes-in-window) over windows of CHUNK consecutive on-disk columns;
    `genes` and `rows` sorted column / row ids."""
    ip = X["indptr"][:]; g0 = 0
    for w in range(0, X.attrs["shape"][1], CHUNK):
        g = genes[(genes >= w) & (genes < w + CHUNK)]
        if not len(g): continue
        lo, hi = ip[g[0]], ip[g[-1] + 1]; ptr = ip[g[0]:g[-1] + 2] - lo
        B = sp.csc_matrix((X["data"][lo:hi], X["indices"][lo:hi], ptr), shape=(X.attrs["shape"][0], g[-1] - g[0] + 1))[:, g - g[0]]
        yield g0, sp.csr_matrix(B[rows], dtype=np.float64); g0 += len(g)


def load_crc_lean(sid, raw=None, ntop=2000, span=.3):
    with h5py.File(Path(raw or l6.RAW) / f"cosmx_wt_{sid}.h5ad", "r") as f:
        obs = read_elem(f["obs"]); X = f["X"]; assert X.attrs["encoding-type"] == "csc_matrix"; n, G = X.attrs["shape"]
        ct = obs["ist"].astype(str).map(lm.LABEL_TO_COARSE); dom = obs["typ"].astype(str).str.extract(r"(REF|TVA|CRC)", expand=False)
        keep = np.flatnonzero((ct.notna() & dom.notna()).values); allg = np.arange(G)
        lib = np.zeros(len(keep))
        for _, B in _chunks(X, keep, allg): lib += np.asarray(B.sum(1)).ravel()
        cells = keep[lib >= 3]  # filter_cells(min_counts=3)
        gsum = np.zeros(G)
        for g0, B in _chunks(X, cells, allg): gsum[g0:g0 + B.shape[1]] = np.asarray(B.sum(0)).ravel()
        genes = allg[gsum >= 3]  # filter_genes(min_counts=3)
        mean, var = np.zeros(len(genes)), np.zeros(len(genes))
        for g0, B in _chunks(X, cells, genes): m, v = _get_mean_var(B); mean[g0:g0 + len(m)], var[g0:g0 + len(m)] = m, v
        from skmisc.loess import loess  # seurat_v3, one batch (scanpy _highly_variable_genes_seurat_v3)
        not_const = var > 0; est = np.zeros(len(genes)); model = loess(np.log10(mean[not_const]), np.log10(var[not_const]), span=span, degree=2)
        model.fit(); est[not_const] = model.outputs.fitted_values; reg_std = np.sqrt(10 ** est); N = len(cells); clip_val = reg_std * np.sqrt(N) + mean
        sq, s = np.zeros(len(genes)), np.zeros(len(genes))
        for g0, B in _chunks(X, cells, genes):
            a, b = _sum_and_sum_squares_clipped(B.indices, B.data, n_cols=B.shape[1], clip_val=clip_val[g0:g0 + B.shape[1]], nnz=B.nnz); sq[g0:g0 + len(a)], s[g0:g0 + len(b)] = a, b
        ngv = ((1 / ((N - 1) * np.square(reg_std))) * (N * np.square(mean) + sq - 2 * s * mean)).reshape(1, -1)
        rk = np.argsort(np.argsort(-ngv, axis=1), axis=1).astype(np.float32); nb = np.sum((rk < ntop).astype(int), axis=0)
        rk[rk >= ntop] = np.nan; med = np.ma.median(np.ma.masked_invalid(rk), axis=0).filled(np.nan)
        df = pd.DataFrame(dict(highly_variable_rank=med, highly_variable_nbatches=nb), index=pd.RangeIndex(len(genes)))
        top = df.sort_values(["highly_variable_rank", "highly_variable_nbatches"], ascending=[True, False], na_position="last").index[:ntop]
        hv = np.zeros(len(genes), bool); hv[np.asarray(top)] = True; hvg = genes[hv]
        C = sp.hstack([B for _, B in _chunks(X, cells, hvg)], format="csr", dtype=np.float64)
    xy = obs[["CenterX_global_px", "CenterY_global_px"]].values.astype(float)[cells]
    l6.log(f"crc_{sid} loaded (lean)", (len(keep), G))
    return C, xy, ct.values[cells].astype(str), dom.values[cells].astype(str)
