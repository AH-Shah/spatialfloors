"""Bruker CosMx WTx flat files (exprMat + metadata csv.gz) -> counts h5ad of the unit's CROP (data.units.BRUKER); the crop is
chosen first from the metadata coordinates (scripts/v_crops.py), so the 400-650k-cell slides never sit in memory (audit
2026-10-30: the whole-slide matrix hit 16 GB). obs index = metadata `cell` (c_<slide>_<fov>_<cell_ID>); x_px / y_px =
CenterX/Y_global_px; fov; IF means renamed by the dataset's channel map; negative / system probes kept (data.genes filters
them; F5 uses Negative<n> as InSituType's background). Channel-map QC: each IF mean vs its marker gene's CP10k (Spearman).
Usage (env floors): python scripts/v_bruker_h5ad.py <unit> <channel map, e.g. G=PanCK,R=CD45,Y=CD68,B=Membrane>"""
import sys
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.stats import spearmanr
from spatialfloors.data.units import BRUKER, BRUKER_SRC

u, cmap = sys.argv[1], dict(kv.split("=") for kv in sys.argv[2].split(","))
d = BRUKER_SRC[u]; ef = next(d.glob("*exprMat*.csv*")); mf = next(d.glob("*metadata*.csv*"))
ids = pd.read_csv(f"manifests/crops/{u}.csv.gz").cell_id.astype(str).values
m = pd.read_csv(mf); m.index = (m["cell"] if "cell" in m else m["cell_id"]).astype(str).values; m = m.loc[ids]
want = pd.Index(m.fov.astype(str) + "_" + m.cell_ID.astype(str)); blocks, keys, genes = [], [], None
for ch in pd.read_csv(ef, chunksize=4000, dtype=np.int32):
    genes = [c for c in ch.columns if c not in ("fov", "cell_ID")] if genes is None else genes
    k = ch.fov.astype(str) + "_" + ch.cell_ID.astype(str); ch = ch[k.isin(want).values]
    if len(ch): keys.append(ch.fov.astype(str) + "_" + ch.cell_ID.astype(str)); blocks.append(sp.csr_matrix(ch[genes].values.astype(np.float32)))
X = sp.vstack(blocks).tocsr(); key = pd.Index(pd.concat(keys).values); pos = key.get_indexer(want); assert (pos >= 0).all(), "crop cells missing in exprMat"
X = X[pos]
obs = pd.DataFrame(dict(x_px=m.CenterX_global_px.values.astype(float), y_px=m.CenterY_global_px.values.astype(float), fov=m.fov.astype(str).values, nCount_RNA=m.nCount_RNA.values), index=ids)
for c, name in cmap.items():
    if f"Mean.{c}" in m: obs[f"Mean.{name}"] = m[f"Mean.{c}"].values.astype(float)
a = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="symbol"))); a.var_names_make_unique()
tot = np.asarray(X.sum(1)).ravel(); print("counts vs metadata nCount_RNA: share within 1 count", float(np.mean(np.abs(tot - m.nCount_RNA.values) <= 1)), flush=True)
import json
from pathlib import Path
VAL = json.loads(Path("manifests/roles_v2.json").read_text())["validation"]  # confirmation crops get no data check here (their QC is the smoke window, v_smoke.py)
cp = sp.csr_matrix(sp.diags(1e4 / np.maximum(tot, 1)) @ X); gi = {g: i for i, g in enumerate(genes)}
for name, gs in () if u not in VAL else (("PanCK", ["KRT8", "KRT18", "KRT19", "EPCAM"]), ("CD45", ["PTPRC"]), ("CD68", ["CD68"])):
    if f"Mean.{name}" in obs: print("IF check", name, {g: round(float(spearmanr(obs[f"Mean.{name}"], cp[:, gi[g]].toarray().ravel())[0]), 3) for g in gs if g in gi}, flush=True)
a.write_h5ad(BRUKER[u], compression="gzip"); print(u, a.shape, "IF", [c for c in obs if c.startswith("Mean.")])
