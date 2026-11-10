"""Exploratory.
F1 re-score of the saved locked embeddings (floor pearson_none_pc50, scVI nl50/nl100, Nicheformer final layer) with
Euclidean instead of cosine k-NN; everything else as metrics.f1 (count split seed 0, filtered-gene library sizes, panel
genes, k 15, lambda .5), as scripts/p3_f1_rescore.py. D = learned - floor (both Euclidean), and learned (Euclidean) -
floor (cosine, registered). Unit means and one-sided 95% bounds as tier 1.
Usage (env floors): python scripts/x_f1_euclid.py runs/p4 <confirmation crops> -> runs/p4/X_f1_euclid.{json,md}"""

import sys, json
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd
from sklearn.neighbors import NearestNeighbors
import spatialfloors.metrics.f1 as M
from spatialfloors import stats
from spatialfloors.data.genes import filter_reason
from spatialfloors.provenance import commit


def knn_euclid(E, k=M.K):
    nb = NearestNeighbors(n_neighbors=k + 1).fit(E).kneighbors(E, return_distance=False)
    self_ = nb == np.arange(len(E))[:, None]
    drop = np.where(self_.any(1), np.argmax(self_, 1), k)
    keep = np.ones_like(nb, bool)
    keep[np.arange(len(E)), drop] = False
    return nb[keep].reshape(len(E), k)


R = Path(sys.argv[1])
U = sys.argv[2:]
L = json.loads(Path("manifests/locked_v3.json").read_text())["locked"]
donors = {
    f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"]
    for d in json.loads(Path("manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]
}
unit_of = lambda c: donors.get(c, c)
reg = pd.read_csv(R / "f1_rescored.csv")
COS = M.knn_cosine
rows = []
for i, u in enumerate(U):
    sp_ = "mouse" if u.startswith("allen") else "human"
    emb = {
        "floor": (f"expr_embed_{L['F1']['expr_embed'][sp_]['config']}", L["F1"]["expr_embed"][sp_]["config"]),
        "scvi": (f"scvi_{L['F1']['scvi'][sp_]['config']}", L["F1"]["scvi"][sp_]["config"]),
        "nicheformer": ("nicheformer_layers", "-1_masked"),
    }
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
    keep = np.flatnonzero(filter_reason(a.var_names, sp_) == "")
    Xa, Xb = M.count_split(a.X, seed=0)
    Xa, Xb = Xa[:, keep], Xb[:, keep]
    vn = a.var_names[keep]
    del a
    g = vn.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"])
    assert (g >= 0).all()
    for m, (stem, key) in emb.items():
        E = np.load(R / f"f1_{u}_{stem}_emb.npz")[key]
        cfgname = {"floor": emb["floor"][1], "scvi": emb["scvi"][1], "nicheformer": "final"}[m]
        r0 = reg[(reg.unit == u) & (reg.config == cfgname)]
        assert len(r0) == 1, (u, m)
        s_cos = float(r0.score.iloc[0])
        if i < 2:  # sanity: cosine re-score equals the registered rescored value
            M.knn_cosine = COS
            assert abs(M.knn_poisson_score(E, Xa, Xb, g) - s_cos) < 1e-9, f"{u} {m}: cosine re-score differs"
        M.knn_cosine = knn_euclid
        s_euc = float(M.knn_poisson_score(E, Xa, Xb, g))
        M.knn_cosine = COS
        rows.append(dict(crop=u, model=m, cosine=s_cos, euclid=s_euc))
        print(rows[-1], flush=True)

df = pd.DataFrame(rows)
W = df.pivot(index="crop", columns="model", values=["cosine", "euclid"])
margin = L["F1"]["margin"]
out = dict(
    commit=commit(),
    margin=margin,
    per_crop=rows,
    tests={},
    means={f"{k}_{m}": float(W[k][m].mean()) for k in ("cosine", "euclid") for m in ("floor", "scvi", "nicheformer")},
)
md = [
    "# X_f1_euclid: F1 with Euclidean k-NN (exploratory, after results-v1; no decision weight)",
    "",
    f"commit {commit()}; margin {margin}.",
    "",
    "Crop means: " + "; ".join(f"{k} {v:.4f}" for k, v in out["means"].items()),
    "",
    "| Comparison | Units | D | One-sided 95% bounds | Non-inferior | Floor superior |",
    "|---|---|---|---|---|---|",
]
for m in ("scvi", "nicheformer"):
    for name, a_, b_ in (
        (f"{m} euclid − floor euclid", W["euclid"][m], W["euclid"]["floor"]),
        (f"{m} euclid − floor cosine", W["euclid"][m], W["cosine"]["floor"]),
        (f"{m} cosine − floor cosine (registered)", W["cosine"][m], W["cosine"]["floor"]),
        (f"{m} best − floor cosine", np.maximum(W["euclid"][m], W["cosine"][m]), W["cosine"]["floor"]),
    ):
        units, D, strata = stats.aggregate_units(dict(zip(W.index, a_ - b_)), {c: unit_of(c) for c in W.index})
        t = stats.family_test(D, strata, margin=margin)
        out["tests"][name] = dict(
            n_units=len(D),
            est=t["est"],
            lower=t["lower"],
            upper=t["upper"],
            noninferior=t["noninferior"],
            floor_superior=t["floor_superior"],
            n_learned_ahead=int((D > 0).sum()),
        )
        md.append(
            f"| {name} | {len(D)} | {t['est']:+.4f} | {t['lower']:+.4f}, {t['upper']:+.4f} | {t['noninferior']} | {t['floor_superior']} |"
        )
(R / "X_f1_euclid.json").write_text(json.dumps(out, indent=1))
(R / "X_f1_euclid.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
