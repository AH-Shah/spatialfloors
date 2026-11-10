"""Re-score every saved F1 embedding with the protocol's library sizes: L = totals over the unit's FILTERED genes
(data.genes
filter_reason == ""), count split seed 0 of the full crop (as fitted) then columns restricted.  Embeddings are
unchanged.
Usage (env floors): python scripts/p3_f1_rescore.py [unit ...] -> runs/p3/f1_rescored.csv (all units) or
runs/p3/f1_rescored_validation.csv."""

import sys
import glob, json, re
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd
from spatialfloors.data.genes import filter_reason
from spatialfloors.metrics.f1 import count_split, knn_poisson_score, chance_level
from spatialfloors.provenance import commit

rows, cache = [], {}
import os

IN = os.environ.get("SF_OUT", "runs/p3")  # P4 / P5: SF_OUT=runs/p4
for f in sorted(glob.glob(f"{IN}/f1_*.json")):
    r = json.loads(Path(f).read_text())
    if re.search(r"_n\d+\.json$", f) or "scores" not in r:
        continue  # profiles
    u = r["unit"]
    if sys.argv[1:] and u not in sys.argv[1:]:
        continue
    if u not in cache:
        a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
        keep = np.flatnonzero(filter_reason(a.var_names, "human" if not u.startswith("allen") else "mouse") == "")
        Xa, Xb = count_split(a.X, seed=0)
        Xa, Xb = Xa[:, keep], Xb[:, keep]
        vn = a.var_names[keep]
        g = vn.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"])
        assert (g >= 0).all()
        cache = {u: (Xa, Xb, g)}
    Xa, Xb, g = cache[u]
    z = np.load(f.replace(".json", "_emb.npz"))
    keys = (
        {"final": "-1_masked", "penultimate": "-2_masked"}
        if r["model"] == "nicheformer"
        else {k: k for k in r["scores"]}
    )
    for k, zk in keys.items():
        E = z[zk]
        rows.append(
            dict(
                unit=u,
                model=r["model"],
                config=k,
                score=float(knn_poisson_score(E, Xa, Xb, g)),
                chance=chance_level(E, Xa, Xb, g, n_perm=1),
                score_allL=r["scores"][k]["score"],
            )
        )
        print(rows[-1], flush=True)
pd.DataFrame(rows).assign(spatialfloors_commit=commit()).to_csv(
    f"{IN}/f1_rescored_validation.csv" if sys.argv[1:] else f"{IN}/f1_rescored.csv", index=False
)
