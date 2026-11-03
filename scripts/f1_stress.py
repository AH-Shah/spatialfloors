"""F1 development stress panel (PLAN §13): does an embedding win the kNN-Poisson score by encoding identity, or by
something the score happens to reward? Development crops only, saved P3 embeddings (runs/p3/f1_<unit>_*_emb.npz).
  agreement   score - chance at k {5, 15, 50} (chance = rows permuted, same k)
  references  depth only (1-D log library size; N2 row R2: library size alone beat a model on an imputation metric),
              spatial smoothing (neighbour mean of Expr PCA 50, 25 um; bench REVIEW #22: smoothing inflates scores),
              random Gaussian embedding
  degeneracy  depth leak (CV R2 of log library size from the embedding); share of the score's deviance gain carried by
              the top 10 genes; share of genes improved (a few dominant genes = the score is a marker detector)
  detail      Allen: kNN (k 15) macro recall over subclasses >= 50 cells, recall on the rarest quartile of those types,
              k-means AMI at k = #types (N2 rows R3 / R4: ARI was EpiT compaction); GSE1: neighbour-predicted IF stain vs
              the cell's own IF (orthogonal truth, bench REVIEW #7), Spearman per stain
Usage (env floors, OMP_NUM_THREADS=4): python scripts/f1_stress.py <unit> [--refs: references only, merged into the existing json]
  -> runs/p3/f1_stress_<unit>.json"""
import glob, json, re, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.linear_model import RidgeCV
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.model_selection import cross_val_predict
from spatialfloors.data.genes import filter_reason
from spatialfloors.data.units import ALLEN
from spatialfloors.floors.simple import expr_pca, neighbour_mean
from spatialfloors.metrics.f1 import count_split, knn_poisson_score, chance_level, knn_cosine
from spatialfloors.provenance import stamp

u = sys.argv[1]; a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); sp_ = "mouse" if u.startswith("allen") else "human"
keep = np.flatnonzero(filter_reason(a.var_names, sp_) == ""); Xa, Xb = count_split(a.X, seed=0); Xa, Xb = Xa[:, keep], Xb[:, keep]
g = a.var_names[keep].get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"]); assert (g >= 0).all()
La = np.asarray(Xa.sum(1)).ravel(); xy = a.obs[["x_um", "y_um"]].values.astype(float); log = lambda *s: print(*s, flush=True)
emb = {}
for f in sorted(glob.glob(f"runs/p3/f1_{u}_*_emb.npz")):
    if re.search(r"_n\d+_emb", f): continue
    z = np.load(f); m = Path(f).name[len(f"f1_{u}_"):-len("_emb.npz")]
    for key in z.files:
        if key.endswith("_authors"): continue
        emb[f"{m.split('_')[0]}:{key}"] = z[key]
E50 = expr_pca(Xa[:, g], 50); rng = np.random.default_rng(0)
zl = (np.log1p(La) - np.log1p(La).mean()) / np.log1p(La).std()
ref = {"ref:depth_only": np.c_[zl, np.ones(len(La))],  # cosine kNN on (z, 1) = neighbours of similar depth
       "ref:smooth_expr50_25um": neighbour_mean(xy, a.obs.fov.values, E50, 25.)[0],
       "ref:random": rng.normal(size=(len(La), 30))}
labels = None
if u.startswith("allen"):
    st = ad.read_h5ad(ALLEN[u.split("_")[1]], backed="r").obs.loc[a.obs_names, "subclass_transfer"].astype(str).values
    t, c = np.unique(st, return_counts=True); ok = t[c >= 50]; labels = np.where(np.isin(st, ok), st, None)
    rare = set(ok[np.argsort(c[c >= 50])[: max(1, len(ok) // 4)]])
ifc = [c for c in a.obs.columns if c.startswith("IF_")]
OUTF = Path(f"runs/p3/f1_stress_{u}.json")
if "--refs" in sys.argv: emb = {}
rep = json.loads(OUTF.read_text()) if "--refs" in sys.argv else dict(unit=u, n_cells=int(a.n_obs), embeddings={})
for name, E in {**emb, **ref}.items():
    t = time.time(); E = np.asarray(E, float); r = {}
    for k in (5, 15, 50):
        r[f"score_k{k}"] = float(knn_poisson_score(E, Xa, Xb, g, k=k)); r[f"chance_k{k}"] = chance_level(E, Xa, Xb, g, k=k, n_perm=1); r[f"excess_k{k}"] = r[f"score_k{k}"] - r[f"chance_k{k}"]
    _, num, den, gg = knn_poisson_score(E, Xa, Xb, g, per_gene=True); gain = den - num
    r["top10_gene_share"] = float(np.sort(gain)[::-1][:10].sum() / max(gain.sum(), 1e-12)) if gain.sum() > 0 else float("nan"); r["genes_improved"] = float((gain > 0).mean())
    r["depth_leak_r2"] = float(1 - np.mean((cross_val_predict(RidgeCV(), E, np.log1p(La), cv=5) - np.log1p(La)) ** 2) / np.var(np.log1p(La)))
    nb = knn_cosine(E, 15)
    if labels is not None:
        lab = labels; ok_ = lab != None; vote = np.array([max(set(v), key=list(v).count) if len(v) else None for v in (lab[n][lab[n] != None] for n in nb)], dtype=object)
        rec = {c: float(np.mean(vote[lab == c] == c)) for c in np.unique(lab[ok_].astype(str))}
        r["knn_macro_recall"] = float(np.mean(list(rec.values()))); r["knn_rare_recall"] = float(np.mean([rec[c] for c in rare]))
        r["kmeans_ami_types"] = float(adjusted_mutual_info_score(lab[ok_].astype(str), KMeans(len(rec), n_init=3, random_state=0).fit_predict(E)[ok_]))
    for c in ifc:
        v = a.obs[c].values.astype(float); r[f"{c}_nbr_spearman"] = float(spearmanr(v[nb].mean(1), v)[0])
    r["s"] = round(time.time() - t, 1); rep["embeddings"][name] = r; log(name, {k: round(v, 3) for k, v in r.items()})
OUTF.write_text(json.dumps(stamp(rep), indent=1, default=float))
