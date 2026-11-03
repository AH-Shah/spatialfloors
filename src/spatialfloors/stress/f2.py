"""F2 stress panel (PLAN §13) as a reusable context: built once per (unit, cell subset), scores any partition.
  agreement   ARI, balanced ARI (N2 T.ARI_bal), AMI, homogeneity, completeness, matched macro-F1, min region recall
  geometry    ARI on boundary cells (other region < 30 um) vs interior cells (> 150 um); ARI vs the region map rotated
              180 deg about the crop centre (blob-geometry null); excess over a label-free spatial-blob partition (k-means on xy)
  degeneracy  cluster-size entropy / log k, largest cluster share, seed stability (mean pairwise ARI)
  detail      the same features clustered at a fine k: AMI vs fine labels (Allen CCF regions >= 1% of cells) and vs cell
              types (Allen subclass; GSE1 Leiden ids), nesting = homogeneity of fine clusters w.r.t. regions"""
import json
from pathlib import Path
import anndata as ad, numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from sklearn import metrics as M
from sklearn.cluster import KMeans
from spatialfloors.data.units import index, ALLEN
from spatialfloors.metrics.f2 import scored_regions, ari_bal


def load(u):
    """Unit inputs: panel counts, xy, fov, region strings, present regions, fine labels (or None), cell types."""
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); p = json.loads(Path(f"manifests/panels/{u}.json").read_text())
    X = a.X[:, a.var_names.get_indexer(p["genes"])]; xy = a.obs[["x_um", "y_um"]].values.astype(float); fov = a.obs.fov.astype(str).values
    if u.startswith("allen"):
        o = ad.read_h5ad(ALLEN[u.split("_")[1]], backed="r").obs.loc[a.obs_names, ["ccf_region_name", "subclass_transfer"]].astype(str)
        fine, ctype = o.ccf_region_name.values, o.subclass_transfer.values
        keep = [r for r, c in zip(*np.unique(fine, return_counts=True)) if c >= .01 * len(fine) and r not in ("nan", "unassigned")]
        fine = np.where(np.isin(fine, keep), fine, None)
    else:
        from spatialfloors.floors import bleedsafe as bs
        from spatialfloors.floors.simple import lognorm
        fine, ctype = None, bs.leiden_ids(bs.expr50(lognorm(X), 0), seed=0).astype(str)
    return dict(unit=u, X=X, xy=xy, fov=fov, region=a.obs.region.astype(str).values, present=index()[u]["present"], fine=fine, ctype=ctype)


def hungarian(yt, yp):
    C = M.cluster.contingency_matrix(yt, yp); r, c = linear_sum_assignment(-C); f1, rec = [], []
    for i, j in zip(r, c):
        tp = C[i, j]; rec.append(tp / C[i].sum()); f1.append(2 * tp / (C[i].sum() + C[:, j].sum()))
    rec += [0.] * (C.shape[0] - len(r)); f1 += [0.] * (C.shape[0] - len(r)); return float(np.mean(f1)), float(np.min(rec))


class Panel:
    def __init__(self, d, cells=None):
        """d from load(); cells: kept cell indices (subsample / FOV hold-out), None = all."""
        c = np.arange(len(d["xy"])) if cells is None else np.asarray(cells)
        self.xy, self.ctype = d["xy"][c], np.asarray(d["ctype"])[c]; self.fine = None if d["fine"] is None else np.asarray(d["fine"], dtype=object)[c]
        self.m, self.y = scored_regions(d["region"][c], d["present"]); self.k = len(np.unique(self.y))
        self.kf = len({x for x in self.fine if x is not None}) if self.fine is not None else 4 * self.k
        ys = np.asarray(d["region"][c], dtype=object); ys[~self.m] = None; idx = np.flatnonzero(self.m); dmin = np.full(len(c), np.inf)
        for r in np.unique(self.y):
            t = cKDTree(self.xy[idx[self.y != r]]); mr = idx[self.y == r]; dmin[mr] = t.query(self.xy[mr])[0]
        self.ys, self.bnd, self.inn = ys, self.m & (dmin < 30), self.m & (dmin > 150)
        c0 = (self.xy.min(0) + self.xy.max(0)) / 2; self.rot = ys[idx[cKDTree(self.xy[idx]).query(2 * c0 - self.xy)[1]]].astype(str)
        self.rot_ok = self.m & np.isin(self.rot, np.unique(self.y))
        self.blob = KMeans(self.k, n_init=10, random_state=0).fit_predict((self.xy - self.xy.mean(0)) / self.xy.std(0)).astype(str)
        self.blob_ari = float(M.adjusted_rand_score(self.y, self.blob[self.m]))

    def one(self, lab):
        lab = np.asarray(lab).astype(str); s = lab[self.m]; f1, rmin = hungarian(self.y, s); sz = np.unique(lab, return_counts=True)[1] / len(lab)
        return dict(ari=M.adjusted_rand_score(self.y, s), ami=M.adjusted_mutual_info_score(self.y, s), homogeneity=M.homogeneity_score(self.y, s),
                    completeness=M.completeness_score(self.y, s), matched_f1=f1, min_region_recall=rmin,
                    ari_boundary=M.adjusted_rand_score(self.ys[self.bnd].astype(str), lab[self.bnd]),
                    ari_interior=M.adjusted_rand_score(self.ys[self.inn].astype(str), lab[self.inn]),
                    ari_rotated_null=M.adjusted_rand_score(self.rot[self.rot_ok], lab[self.rot_ok]), ari_vs_blob=M.adjusted_rand_score(self.blob, lab),
                    size_entropy=float(-(sz * np.log(sz)).sum() / np.log(max(len(sz), 2))), largest_share=float(sz.max()), ari_bal=self.ari_bal(s))

    def ari_bal(self, s):
        return ari_bal(self.y, s)

    def __call__(self, L):
        """L: seeds x cells. Seed-averaged panel + seed sd of ARI, seed stability, excess over the blob partition."""
        out = [self.one(lab) for lab in L]; r = {key: float(np.mean([o[key] for o in out])) for key in out[0]}
        r["ari_sd"] = float(np.std([o["ari"] for o in out]))
        r["seed_stability"] = float(np.mean([M.adjusted_rand_score(L[i], L[j]) for i in range(len(L)) for j in range(i + 1, len(L))])) if len(L) > 1 else float("nan")
        r["excess_over_blob"] = r["ari"] - self.blob_ari; return r

    def detail(self, F, clus, seed=0):
        """Fine-k partition of the same features: AMI vs fine regions / cell types, nesting in the coarse regions."""
        lab = clus(F, self.kf, seed).astype(str)
        d = dict(k_fine=self.kf, ami_celltype=float(M.adjusted_mutual_info_score(self.ctype, lab)), nesting_homogeneity=float(M.homogeneity_score(lab[self.m], self.y)))
        if self.fine is not None: f = self.fine != None; d["ami_fine_regions"] = float(M.adjusted_mutual_info_score(self.fine[f].astype(str), lab[f]))
        return d
