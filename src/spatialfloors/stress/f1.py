"""F1 stress panel (PLAN §13) as a reusable context: built once per (unit, perturbation), profiles any embedding.
  agreement   score - chance at k {5, 15, 50} (chance = rows permuted; the same for every embedding, computed once)
  degeneracy  share of the deviance gain carried by the top 10 genes, share of genes improved; depth leak (CV R2 of
              log library size from the embedding)
  detail      Allen: kNN (k 15) macro recall over subclasses >= 50 cells, rarest-quartile recall, k-means AMI at k = #types;
              IF units: neighbour-mean IF stain vs the cell's own IF (orthogonal truth), Spearman per stain
Floor variants (floor_embed): normalisation x gene scaling x n_pc, the knobs of the expr_pca floor."""
import json
from pathlib import Path
import anndata as ad, numpy as np, scipy.sparse as sp
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.model_selection import cross_val_predict
from spatialfloors.data.genes import filter_reason
from spatialfloors.data.units import ALLEN
from spatialfloors.floors.simple import expr_embed
from spatialfloors.metrics.f1 import count_split, knn_poisson_score, chance_level, knn_cosine

KS = (5, 15, 50)


def load(u):
    """Unit inputs on the scored gene set: raw counts (filtered genes), panel column indices, xy, fov, labels (Allen), IF."""
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); sp_ = "mouse" if u.startswith("allen") else "human"
    keep = np.flatnonzero(filter_reason(a.var_names, sp_) == ""); X = sp.csr_matrix(a.X[:, keep])
    g = a.var_names[keep].get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"]); assert (g >= 0).all()
    labels = None
    if u.startswith("allen"):
        st = ad.read_h5ad(ALLEN[u.split("_")[1]], backed="r").obs.loc[a.obs_names, "subclass_transfer"].astype(str).values
        t, c = np.unique(st, return_counts=True); labels = np.where(np.isin(st, t[c >= 50]), st, None)
    IF = {c: a.obs[c].values.astype(float) for c in a.obs.columns if c.startswith("IF_")}
    return dict(unit=u, X=X, g=g, xy=a.obs[["x_um", "y_um"]].values.astype(float), fov=a.obs.fov.astype(str).values, labels=labels, IF=IF)


floor_embed = expr_embed  # the F1 floor family lives in floors.simple


class Panel:
    def __init__(self, d, split_seed=0, cells=None, genes=None, X=None):
        """d from load(); X: perturbed counts (rows = cells, all filtered genes); cells / genes: kept indices."""
        c = np.arange(d["X"].shape[0]) if cells is None else np.asarray(cells); X = d["X"][c] if X is None else X
        self.g = d["g"] if genes is None else d["g"][np.asarray(genes)]
        self.Xa, self.Xb = count_split(X, seed=split_seed); self.La = np.asarray(self.Xa.sum(1)).ravel()
        self.labels = None if d["labels"] is None else d["labels"][c]; self.IF = {k: v[c] for k, v in d["IF"].items()}
        self.chance = {k: chance_level(np.random.default_rng(0).normal(size=(len(c), 10)), self.Xa, self.Xb, self.g, k=k, n_perm=1) for k in KS}
        if self.labels is not None:
            ok = self.labels[self.labels != None]; t, n = np.unique(ok.astype(str), return_counts=True)
            self.rare = set(t[np.argsort(n)[: max(1, len(t) // 4)]])

    def __call__(self, E, ks=KS):
        E = np.asarray(E, float); r = {}
        for k in ks:
            if k == 15:
                s, num, den, _ = knn_poisson_score(E, self.Xa, self.Xb, self.g, k=15, per_gene=True); gain = den - num
                r["top10_gene_share"] = float(np.sort(gain)[::-1][:10].sum() / gain.sum()) if gain.sum() > 0 else float("nan"); r["genes_improved"] = float((gain > 0).mean())
            else: s = knn_poisson_score(E, self.Xa, self.Xb, self.g, k=k)
            r[f"score_k{k}"] = float(s); r[f"excess_k{k}"] = float(s) - self.chance[k]
        ll = np.log1p(self.La); r["depth_leak_r2"] = float(1 - np.mean((cross_val_predict(RidgeCV(), E, ll, cv=5) - ll) ** 2) / np.var(ll))
        nb = knn_cosine(E, 15)
        if self.labels is not None:
            lab = self.labels; ok = lab != None
            vote = np.array([max(set(v), key=list(v).count) if len(v) else None for v in (lab[n][lab[n] != None] for n in nb)], dtype=object)
            rec = {c: float(np.mean(vote[lab == c] == c)) for c in np.unique(lab[ok].astype(str))}
            r["knn_macro_recall"] = float(np.mean(list(rec.values()))); r["knn_rare_recall"] = float(np.mean([rec[c] for c in self.rare if c in rec]))
            r["kmeans_ami_types"] = float(adjusted_mutual_info_score(lab[ok].astype(str), KMeans(len(rec), n_init=3, random_state=0).fit_predict(E)[ok]))
        for c, v in self.IF.items(): r[f"{c}_nbr_spearman"] = float(spearmanr(v[nb].mean(1), v)[0])
        return r
