"""Exploratory: ceilings per family and unit, from saved round-1 outputs, and headroom share
= (learned - floor) / (ceiling - floor).
F5: the held-out-marker counts of every cell are split by binomial thinning (p = .5) into halves A and B (5 draws, both
directions). Oracle at a labelling's call rates: per class c, the n_c cells with the highest half-A score of c's
markers are called c (n_c = cells the labelling calls c); its held-out-marker AUROC is computed on half B (metric
.heldout_auroc rule per class: AUROC of the half-B score, called vs not called; < 50 cells on either side -> .5).
Floor and learned labellings are re-scored on half B. Oracles at every labelling's call rates are saved
(scripts/x_ceiling_report.py: headroom uses the oracle at the floor's call rates; the own-rate oracle is the check).
Usage (env floors, round-1 runs root as working directory):
  python scripts/x_ceiling.py f5 <labels dir> <unit> <out dir> [--cfg name=file_stem:key ...]
  default --cfg: floor=floors:insitutype_spatial scanvi=scanvi:nl30_nb tangram=tangram:uniform (locked)
  -> <out dir>/f5_<unit>.json"""

import argparse, json
from pathlib import Path
import anndata as ad, numpy as np, scipy.sparse as sp
from sklearn.metrics import roc_auc_score
from spatialfloors.floors import f5
from spatialfloors.metrics.f5 import CLASSES, MIN_PRED, heldout_auroc, unit_inputs
from spatialfloors.provenance import stamp

DRAWS = 5
LOCKED = ("floor=floors:insitutype_spatial", "scanvi=scanvi:nl30_nb", "tangram=tangram:uniform")


def class_scores(X, var_names, markers, lib):
    """n x classes log1p CP10k of each class's summed held-out markers (library = full row sum)."""
    vn = list(map(str, var_names))
    S = np.zeros((X.shape[0], len(CLASSES)))
    for k, c in enumerate(CLASSES):
        gi = [vn.index(g) for g in markers[c] if g in vn]
        S[:, k] = np.log1p(np.asarray(X[:, gi].sum(1)).ravel() / lib * 1e4) if gi else 0.0
    return S


def auroc_called(called, s):
    pos = called.sum()
    if min(pos, len(called) - pos) < MIN_PRED:
        return 0.5
    return float(roc_auc_score(called, s))


def f5_unit(u, labdir, cfg, out):
    MK = f5.for_unit(u, json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"])
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
    X, vn, _, _ = unit_inputs(u, a)
    X = sp.csr_matrix(X)
    vl = list(map(str, vn))
    gi = sorted({vl.index(g) for c in CLASSES for g in MK[c] if g in vl})
    M = X[:, gi].toarray().astype(np.int64)
    lib = np.maximum(np.asarray(X.sum(1)).ravel(), 1)
    lib_rest = lib - M.sum(1)
    lab = {
        n: f5.coarse(
            np.load(Path(f"{stem}.npz") if "/" in stem else Path(labdir) / f"f5_{u}_{stem}.npz", allow_pickle=True)[key]
        )
        for n, (stem, key) in ((n, v.split(":")) for n, v in (x.split("=") for x in cfg))
    }
    full = {n: heldout_auroc(c, X, vn, MK)[0] for n, c in lab.items()}
    rng = np.random.default_rng(0)
    rec = {n: [] for n in lab} | {f"oracle_{n}": [] for n in lab}
    sub = lambda Mh: sp.csr_matrix((Mh.ravel(), (np.repeat(np.arange(len(Mh)), len(gi)), np.tile(np.arange(len(gi)), len(Mh)))), shape=(len(Mh), len(gi)))
    for d in range(DRAWS):
        A = rng.binomial(M, 0.5)
        for HA, HB in ((A, M - A), (M - A, A)):
            libA, libB = np.maximum(lib_rest + HA.sum(1), 1), np.maximum(lib_rest + HB.sum(1), 1)
            vg = np.array(vl)[gi]
            SA, SB = class_scores(sub(HA), vg, MK, libA), class_scores(sub(HB), vg, MK, libB)
            for n, c in lab.items():
                rec[n].append(float(np.mean([auroc_called(c == k, SB[:, j]) for j, k in enumerate(CLASSES)])))
                orc = []
                for j, k in enumerate(CLASSES):
                    called = np.zeros(len(c), bool)
                    called[np.argsort(-SA[:, j], kind="stable")[: int((c == k).sum())]] = True
                    orc.append(auroc_called(called, SB[:, j]))
                rec[f"oracle_{n}"].append(float(np.mean(orc)))
    r = dict(
        unit=u,
        n_cells=int(X.shape[0]),
        markers_on_platform={c: [g for g in MK[c] if g in vl] for c in CLASSES},
        cfg=list(cfg),
        full=full,
        half_b_mean={k: float(np.mean(v)) for k, v in rec.items()},
        half_b_sd={k: float(np.std(v)) for k, v in rec.items()},
        half_b_draws=rec,
        call_share={n: {k: float((c == k).mean()) for k in (*CLASSES, "other")} for n, c in lab.items()},
    )
    h = r["half_b_mean"]
    r["headroom_share"] = {
        n: (h[n] - h["floor"]) / (h[f"oracle_{n}"] - h["floor"]) if h[f"oracle_{n}"] != h["floor"] else None
        for n in lab
        if n != "floor"
    }
    Path(out).mkdir(parents=True, exist_ok=True)
    (Path(out) / f"f5_{u}.json").write_text(json.dumps(stamp(r), indent=1, default=float))
    print(u, {k: round(v, 4) for k, v in full.items()}, {k: round(v, 4) for k, v in h.items()}, flush=True)


def oof_logreg(F, y, blocks):
    """Out-of-fold multinomial logistic regression (standardised features), one fold per spatial block."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    pred = np.empty(len(y), dtype=object)
    for b in np.unique(blocks):
        te = blocks == b
        clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
        pred[te] = clf.fit(F[~te], y[~te]).predict(F[te])
    return pred.astype(str)


def f2_unit(u, labdir, cfg, out):
    """F2 ceiling: out-of-fold logistic regression from the locked floor's features (floors.simple.domain_features,
    ell 400 um, lambda .8) to the annotated regions, scored with balanced ARI on the scored cells; folds = 5 random
    cell folds (ceiling) and 4 spatial quadrants (ceiling_quadrant_folds, reported; the random-fold version was chosen
    for the headroom after the quadrant version scored below the floor on 7 of 10 crops). Floor and learned: balanced ARI of their saved partitions (seed mean). unit "toy": sim.domain_toy,
    floor = GMM on the same features (10 seeds)."""
    from spatialfloors.floors.simple import domain_features
    from spatialfloors.metrics.f2 import SEEDS, ari_bal, ari_bal_saved, gmm_labels, scored_regions
    from spatialfloors.metrics.f3 import quadrant_blocks
    from sklearn.metrics import adjusted_rand_score

    if u == "toy":
        from spatialfloors import sim

        X, xy, fov, dom = sim.domain_toy(seed=0)
        m, regl = np.ones(len(dom), bool), np.asarray(dom).astype(str)
    else:
        from spatialfloors.data.units import index

        a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
        g = a.var_names.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"])
        X, xy, fov = a.X[:, g], a.obs[["x_um", "y_um"]].values, a.obs.fov.values
        m, regl = scored_regions(a.obs.region.astype(str).values, index()[u]["present"])
    F, _ = domain_features(X, xy, fov, 50, 400.0, lam=0.8)
    pred = oof_logreg(F[m], regl, np.random.default_rng(0).integers(0, 5, int(m.sum())))
    pq = oof_logreg(F[m], regl, quadrant_blocks(xy[m])[0])
    r = dict(unit=u, n_scored=int(m.sum()), k=int(len(np.unique(regl))), ceiling=ari_bal(regl, pred),
             ceiling_ari=float(adjusted_rand_score(regl, pred)), ceiling_quadrant_folds=ari_bal(regl, pq), scores={})
    if u == "toy":
        r["scores"]["floor"] = float(np.mean([ari_bal(regl, gmm_labels(F, r["k"], s).astype(str)) for s in SEEDS]))
    else:
        for x in cfg:
            n, v = x.split("=")
            stem, key = v.split(":")
            p = Path(f"{stem}.npz") if "/" in stem else Path(labdir) / f"f2_{u}_{stem}.npz"
            if p.exists() and key in np.load(p, allow_pickle=True).files:
                r["scores"][n] = ari_bal_saved(p, key, m, regl)
    Path(out).mkdir(parents=True, exist_ok=True)
    (Path(out) / f"f2_{u}.json").write_text(json.dumps(stamp(r), indent=1, default=float))
    print(u, round(r["ceiling"], 4), {k: round(v, 4) for k, v in r["scores"].items()}, flush=True)


def f3_unit(u, runs, out, draws=2):
    """F3 ceiling: per receiver of the floor's design (floors.f3.prepare on the frozen ids), the counts of its masked
    genes are split by binomial thinning into halves A and B (2 draws x 2 directions). Per held-out quadrant (block_cv
    layout, edge buffer), the full-data training fit (depth covariate, z-score, PCA top 10) is applied to the test rows
    of each half; per PC, Spearman rho between half-A and half-B scores, Spearman-Brown r = 2 rho / (1 + rho), ceiling
    = sqrt(max(r, 0)) (the largest correlation a predictor can have with the observed state PC); mean over PCs, folds,
    draws, then receivers. Floor and learned: the registered per-crop unit scores (<runs>/f3_<crop>_<model>.json)."""
    from scipy.stats import spearmanr
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from spatialfloors.floors import bleedsafe as bs, f3
    from spatialfloors.metrics.f3 import NPC, edge_buffer

    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
    g = a.var_names.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"])
    Xc = sp.csr_matrix(a.X[:, g])
    D = f3.prepare(Xc, a.obs[["x_um", "y_um"]].values, a.obs.fov.values, ids=np.load(f"runs/p2/f3_ids/{u}.npy"))
    rng = np.random.default_rng(0)

    def half_target(C, fovr):
        y = sp.csr_matrix(C.astype(np.float32))
        y.data = np.log1p(y.data)
        x = bs.renorm_log(y)
        return bs.group_center(x.toarray(), fovr), bs.group_center(np.log1p(np.diff(x.indptr)).astype(float)[:, None], fovr)

    per = {}
    for r in f3.receivers(D):
        rows, Y, dep = f3.target(D, r)
        C = Xc[rows][:, np.flatnonzero(D["mask"][r])].toarray().astype(np.int64)
        fovr = D["fov"][rows]
        blocks, buf = D["blocks"][rows], edge_buffer(D["xy"][rows], D["split"], f3.RF_UM)
        vals = []
        for _ in range(draws):
            Ah = rng.binomial(C, 0.5)
            for HA, HB in ((Ah, C - Ah), (C - Ah, Ah)):
                (YA, dA), (YB, dB) = half_target(HA, fovr), half_target(HB, fovr)
                for q in range(4):
                    tr, te = np.flatnonzero(blocks != q), np.flatnonzero((blocks == q) & ~buf)
                    if len(te) < 20 or len(tr) < 20:
                        continue
                    Cm = np.c_[np.ones(len(Y)), dep]
                    beta = np.linalg.lstsq(Cm[tr], Y[tr], rcond=None)[0]
                    sy = StandardScaler().fit((Y - Cm @ beta)[tr])
                    pca = PCA(min(NPC, Y.shape[1]), random_state=0).fit(sy.transform((Y - Cm @ beta)[tr]))
                    PA = pca.transform(sy.transform(YA[te] - np.c_[np.ones(len(te)), dA[te]] @ beta))
                    PB = pca.transform(sy.transform(YB[te] - np.c_[np.ones(len(te)), dB[te]] @ beta))
                    rho = np.array([spearmanr(PA[:, j], PB[:, j])[0] for j in range(PA.shape[1])])
                    vals.append(float(np.nanmean(np.sqrt(np.clip(2 * rho / (1 + rho), 0, None)))))
        per[int(r)] = float(np.mean(vals)) if vals else float("nan")
    L = json.loads(Path("manifests/locked_v3.json").read_text())["locked"]["F3"]
    NC = json.loads(Path("manifests/locked_v3_ncem.json").read_text())["locked"]
    spc = "mouse" if u.startswith("allen") else "human"
    sc = {}
    for n, stem, c in (("floor", "ridge", L["ridge"][spc]["config"]), ("lightgbm", "lightgbm", L["lightgbm"][spc]["config"]),
                       ("ncem", "ncemA", NC[spc]["config"])):
        p = Path(runs) / f"f3_{u}_{stem}.json"
        gr = json.loads(p.read_text())["grid"] if p.exists() else {}
        if c in gr:
            sc[n] = gr[c]["unit"]
    v = [x for x in per.values() if np.isfinite(x)]
    r = dict(unit=u, receivers=per, ceiling=float(np.mean(v)) if v else None, scores=sc)
    Path(out).mkdir(parents=True, exist_ok=True)
    (Path(out) / f"f3_{u}.json").write_text(json.dumps(stamp(r), indent=1, default=float))
    print(u, r["ceiling"], sc, flush=True)


def f1_unit(u, runs, out):
    """F1 oracle reference (not a bound): neighbours from the expression PCA (floors.simple.expr_pca, 50 PCs) of the
    held-out half X_b itself, scored with metrics.f1.knn_poisson_score exactly as the rescored embeddings (count split
    seed 0, filtered genes, panel HVGs). Floor and learned: the registered rescored scores (<runs>/f1_rescored.csv)."""
    import pandas as pd
    from spatialfloors.data.genes import filter_reason
    from spatialfloors.floors.simple import expr_pca
    from spatialfloors.metrics.f1 import count_split, knn_poisson_score

    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
    spc = "mouse" if u.startswith("allen") else "human"
    keep = np.flatnonzero(filter_reason(a.var_names, spc) == "")
    Xa, Xb = count_split(a.X, seed=0)
    Xa, Xb = Xa[:, keep], Xb[:, keep]
    g = a.var_names[keep].get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"])
    ref = float(knn_poisson_score(expr_pca(Xb[:, g]), Xa, Xb, g))
    L = json.loads(Path("manifests/locked_v3.json").read_text())["locked"]["F1"]
    d = pd.read_csv(Path(runs) / "f1_rescored.csv")
    sc = {}
    for n, m in (("floor", "expr_embed"), ("scvi", "scvi"), ("nicheformer", "nicheformer")):
        s = d[(d.unit == u) & (d.model == m) & (d.config == L[m][spc]["config"])].score
        if len(s):
            sc[n] = float(s.iloc[0])
    r = dict(unit=u, ceiling=ref, scores=sc)
    Path(out).mkdir(parents=True, exist_ok=True)
    (Path(out) / f"f1_{u}.json").write_text(json.dumps(stamp(r), indent=1, default=float))
    print(u, round(ref, 4), {k: round(v, 4) for k, v in sc.items()}, flush=True)


F2_LOCKED = {
    "human": ("floor=floor:ell400_lam0.8", "cellcharter=cellcharter:nl10_L3"),
    "mouse": ("floor=floor:ell400_lam0.8", "cellcharter=cellcharter:nl30_L3", "novae=novae:zero_shot"),
}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("fam", choices=("f5", "f2", "f3", "f1"))
    ap.add_argument("labdir")
    ap.add_argument("unit")
    ap.add_argument("out")
    ap.add_argument("--cfg", nargs="*")
    A = ap.parse_args()
    if A.fam == "f5":
        f5_unit(A.unit, A.labdir, A.cfg or list(LOCKED), A.out)
    elif A.fam == "f3":  # labdir = the runs dir holding f3_<crop>_<model>.json
        f3_unit(A.unit, A.labdir, A.out)
    elif A.fam == "f1":  # labdir = the runs dir holding f1_rescored.csv
        f1_unit(A.unit, A.labdir, A.out)
    else:
        f2_unit(A.unit, A.labdir, A.cfg or list(F2_LOCKED["mouse" if A.unit.startswith("allen") else "human"]), A.out)
