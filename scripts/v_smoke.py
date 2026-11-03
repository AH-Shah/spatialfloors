"""Metric-blind smoke window for a new-cancer confirmation slide (PLAN §17 G5): ~5k contiguous cells OUTSIDE the unit's crop
(densest square among the remaining cells), read from the flat files. Checks, no model and no metric: the IF channel map
(each IF mean vs its marker gene's CP10k, Spearman), detection of the F5 held-out markers of the unit's reference
(share of cells with >= 1 count), and the reference's class coverage (manifests/refs/<tissue>.json).
Usage (env floors): python scripts/v_smoke.py <unit> <channel map> -> runs/p2/smoke_<unit>.json"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.stats import spearmanr
from spatialfloors.data.crops import choose_crop
from spatialfloors.data.units import BRUKER_SRC, bruker_obs
from spatialfloors.floors.f5 import UNIT_REF
from spatialfloors.provenance import stamp

u, cmap = sys.argv[1], dict(kv.split("=") for kv in sys.argv[2].split(","))
o = bruker_obs(u); crop = set(pd.read_csv(f"manifests/crops/{u}.csv.gz").cell_id.astype(str)); rest = o[~o.index.isin(crop)]
idx, rec = choose_crop(rest[["x_px", "y_px"]].values, None, n_target=5000); w = rest.iloc[idx]
d = BRUKER_SRC[u]; ef = next(d.glob("*exprMat*.csv*")); mf = next(d.glob("*metadata*.csv*"))
want = pd.Index(w.fov.astype(str) + "_" + w.cell_ID.astype(str)); blocks, keys, genes = [], [], None
for ch in pd.read_csv(ef, chunksize=4000, dtype=np.int32):
    genes = [c for c in ch.columns if c not in ("fov", "cell_ID")] if genes is None else genes
    k = ch.fov.astype(str) + "_" + ch.cell_ID.astype(str); ch = ch[k.isin(want).values]
    if len(ch): keys.append(ch.fov.astype(str) + "_" + ch.cell_ID.astype(str)); blocks.append(sp.csr_matrix(ch[genes].values.astype(np.float32)))
X = sp.vstack(blocks).tocsr()[pd.Index(pd.concat(keys).values).get_indexer(want)]
m = pd.read_csv(mf); m.index = (m["cell"] if "cell" in m else m["cell_id"]).astype(str).values; m = m.loc[w.index]
tot = np.maximum(np.asarray(X.sum(1)).ravel(), 1); cp = sp.csr_matrix(sp.diags(1e4 / tot) @ X); gi = {g: i for i, g in enumerate(genes)}
out = dict(unit=u, n_cells=int(X.shape[0]), outside_crop=True, counts_per_cell=float(tot.mean()), IF={}, markers={}, reference=None)
for c, name in cmap.items():
    if f"Mean.{c}" in m and name in ("PanCK", "CD45", "CD68"):
        gs = {"PanCK": ["KRT8", "KRT18", "KRT19", "EPCAM"], "CD45": ["PTPRC"], "CD68": ["CD68"]}[name]
        out["IF"][name] = {g: float(spearmanr(m[f"Mean.{c}"].values, cp[:, gi[g]].toarray().ravel())[0]) for g in gs if g in gi}
t = UNIT_REF.get(u); r = json.loads(Path(f"manifests/refs/{t}.json").read_text()) if t and Path(f"manifests/refs/{t}.json").exists() else None
if r:
    out["reference"] = dict(tissue=t, class_counts=r["class_counts"], qc_pass=r["qc_pass"])
    out["markers"] = {c: {g: float((X[:, gi[g]].toarray() > 0).mean()) for g in gs if g in gi} for c, gs in r["heldout_markers"].items()}
Path(f"runs/p2/smoke_{u}.json").write_text(json.dumps(stamp(out), indent=1, default=float)); print(json.dumps(out, default=float)[:1500])
