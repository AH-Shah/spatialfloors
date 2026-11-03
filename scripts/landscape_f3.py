"""F3 floor characterisation landscape (PLAN §14). Development units only. One JSON line per landscape point
(runs/landscape/f3_<unit>.jsonl, resumable), each with the stress.f3 panel averaged over receivers.
  grid     ell {10, 25, 50, 100} um x context PCs {10, 30, 50, 100} x features {T, ctx, T+ctx} x ridge alpha {CV (P3), .1,
           10, 1000}; T has no PCs. Context SVD is fitted once per (receiver, fold) at 100 components on training quadrants
           and sliced (leading components of a randomized SVD; P3 fits k components directly)
  null     features shuffled within FOV x quadrant at the anchors -> above_null
  seeds    prepare seed 1, 2 (Leiden ids, context SVD) at the anchors
  perturb  thin50, thin25, sub50, fov75, genes75 x seeds {0, 1} at the anchors; receiver ids fixed from the unperturbed
           data (re-deriving them changes the receiver set, e.g. 278N 8 -> 5 receivers at thin25, so unit means stop being comparable)
  protein  IF units: stress.f3.Protein shift-null check at the anchors and for LightGBM
  learned  P3 GAT / LightGBM unit scores (runs/p3/f3_<unit>_<model>.json); LightGBM (P3-best config) on the panel, its
           shuffled null and thin50 / sub50 / genes75 (seed 0). GAT under perturbations: scripts/landscape_learned_f3.py (MPS)
Usage (pixi env scvi, OMP_NUM_THREADS=4): python scripts/landscape_f3.py <unit> [grid null seeds perturb protein learned ...]"""
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge, RidgeCV
import scipy.sparse as sp
from spatialfloors.floors import f3, bleedsafe as bs
from spatialfloors.metrics import f3 as m3
from spatialfloors.provenance import stamp
from spatialfloors.stress import perturb
from spatialfloors.stress.f3 import panel, summary, receiver_panel, Protein

u = sys.argv[1]; parts = sys.argv[2:] or ["grid", "null", "seeds", "perturb", "protein", "learned"]
a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); g = a.var_names.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"])
X0, XY, FOV = sp.csr_matrix(a.X[:, g]), a.obs[["x_um", "y_um"]].values.astype(float), a.obs.fov.astype(str).values
IF = [c for c in a.obs.columns if c.startswith("IF_")]
OUT = Path("runs/landscape"); fo = OUT / f"f3_{u}.jsonl"; log = lambda *s: print(*s, flush=True)
done = {json.loads(l)["id"] for l in fo.read_text().splitlines()} if fo.exists() else set()
ELL, NPC, FEATS, ALPHA = (10., 25., 50., 100.), (10, 30, 50, 100), ("T", "ctx", "T+ctx"), ("cv", .1, 10., 1000.)
ANCHORS = [("T+ctx", 25., 50, "cv"), ("T+ctx", 10., 30, "cv"), ("ctx", 25., 50, "cv"), ("T", 25., 0, "cv")]  # P3 floor grid corners + feature ablations
key = lambda ft, ell, npc, al: f"{ft}|ell{ell:g}|pc{npc}|a{al}"


def write(rec):
    with fo.open("a") as f: f.write(json.dumps(stamp(rec), default=float) + "\n")
    done.add(rec["id"]); log(rec["id"], {k: round(v, 4) for k, v in rec.items() if k in ("top10", "pc1", "gene_median", "worst_top10", "above_null", "prot_pred_sig", "prot_true_sig")}, rec.get("s"))


def fitter(al):
    if al == "cv": return f3.ridge_fit_predict
    return lambda Xtr, Ytr, Xte: Ridge(alpha=al).fit(Xtr, Ytr).predict(Xte)


def ctx_cache(D, r, rows, kers, seed=0, kmax=100):
    """{ell: {q: kernel mean of the bleed-safe context (rows x kmax)}}: SVD fitted per fold on the training quadrants (as floors.f3.context_fn)."""
    gg = np.flatnonzero(~D["mask"][r]); k = min(kmax, len(gg) - 1); out = {e: {} for e in kers}
    if k < 1: return {e: {q: np.zeros((len(rows), kmax), np.float32) for q in range(4)} for e in kers}  # crash fix: < 2 foreign genes (thin25) -> no context
    xg = bs.renorm_log(D["x"][:, gg])
    for q in range(4):
        fit = np.flatnonzero(D["blocks"] != q); fit = np.sort(np.random.default_rng(seed).choice(fit, min(len(fit), 100000), replace=False))
        m1 = np.asarray(xg[fit].mean(0)).ravel(); sd = np.sqrt(np.maximum(np.asarray(xg[fit].multiply(xg[fit]).mean(0)).ravel() - m1 ** 2, 1e-6))
        xs = sp.csr_matrix(xg @ sp.diags(1 / sd)); E = bs.group_center(TruncatedSVD(k, random_state=seed).fit(xs[fit]).transform(xs), D["fov"]).astype(np.float32)
        for e, (nb, w, _) in kers.items(): out[e][q] = bs.kernel_mean(nb, w, E, rows)
    return out


def score_points(D, pts, shuffle_seed=None, full_out=False):
    """pts: list of (feats, ell, npc, alpha). Returns {pt: {receiver: panel}} (and {pt: score_floor-style out} if full_out)."""
    ells = sorted({p[1] for p in pts}); kers = {}
    for e in ells: nb, w = f3.kernel(D, e); kers[e] = (nb, w, f3.id_mix(D, nb, w))
    res, outs = {p: {} for p in pts}, {p: {} for p in pts}
    for r in f3.receivers(D):
        rows, Y, dep = f3.target(D, r); C = ctx_cache(D, r, rows, kers, seed=D.get("seed", 0)) if any(p[0] != "T" for p in pts) else None
        for p in pts:
            ft, e, npc, al = p; T = kers[e][2][rows]
            Xf = {"T": lambda q, tr: T, "ctx": lambda q, tr: C[e][q][:, :npc], "T+ctx": lambda q, tr: np.c_[T, C[e][q][:, :npc]]}[ft]
            if shuffle_seed is not None: Xf = f3._shuffled(Xf, D["fov"][rows], D["blocks"][rows], shuffle_seed)
            det = {"full": True}
            try: s, fo_ = m3.block_cv(Xf, Y, D["xy"][rows], fitter(al), f3.RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep, detail=det)
            except ValueError as e:  # crash fix 2026-10-29: no non-constant feature column (e.g. thin25 merges every id into one lineage -> T empty)
                log("degenerate", p, r, str(e)[:80]); s, fo_, det = float("nan"), [], {}
            res[p][int(r)] = receiver_panel(s, det)
            if full_out: outs[p][int(r)] = (s, fo_, None, det)
    return (res, outs) if full_out else res


def prep(X, xy, fov, seed=0, ids=None):
    D = f3.prepare(X, xy, fov, ids=ids, seed=seed); D["seed"] = seed; return D


D0 = prep(X0, XY, FOV)
if "grid" in parts:
    pts = [(ft, e, 0 if ft == "T" else n, al) for ft in FEATS for e in ELL for n in ((0,) if ft == "T" else NPC) for al in ALPHA]
    pts = [p for p in pts if f"grid|{key(*p)}" not in done]
    if pts:
        t = time.time(); R = score_points(D0, pts); s = round((time.time() - t) / len(pts), 1)
        for p in pts: write(dict(id=f"grid|{key(*p)}", unit=u, part="grid", features=p[0], ell=p[1], n_pc=p[2], alpha=str(p[3]), perturb="none", **summary(R[p]), per_receiver=R[p], s=s))
NULL = {}
if "null" in parts or "learned" in parts:
    nul = [p for p in ANCHORS if f"null|{key(*p)}" not in done]
    if nul:
        R = score_points(D0, nul, shuffle_seed=0)
        for p in nul: write(dict(id=f"null|{key(*p)}", unit=u, part="null", features=p[0], ell=p[1], n_pc=p[2], alpha=str(p[3]), perturb="none", **summary(R[p]), per_receiver=R[p]))
    for l in fo.read_text().splitlines():
        r = json.loads(l)
        if r["part"] == "null": NULL[key(r["features"], r["ell"], r["n_pc"], r["alpha"] if r["alpha"] == "cv" else float(r["alpha"]))] = {int(k): v for k, v in r["per_receiver"].items()}
    with fo.open() as f: G = [json.loads(l) for l in f]
    for p in ANCHORS:
        pid = f"anchor|{key(*p)}"
        if pid in done or key(*p) not in NULL: continue
        gl = next((r for r in G if r["id"] == f"grid|{key(*p)}"), None)
        if gl: per = {int(k): v for k, v in gl["per_receiver"].items()}; write(dict(id=pid, unit=u, part="anchor", features=p[0], ell=p[1], n_pc=p[2], alpha=str(p[3]), perturb="none", **summary(per, NULL[key(*p)])))
if "seeds" in parts:
    for sd in (1, 2):
        pts = [p for p in ANCHORS if f"seed|{key(*p)}|s{sd}" not in done]
        if not pts: continue
        D = prep(X0, XY, FOV, sd); R = score_points(D, pts)
        for p in pts: write(dict(id=f"seed|{key(*p)}|s{sd}", unit=u, part="split", features=p[0], ell=p[1], n_pc=p[2], alpha=str(p[3]), split_seed=sd, perturb="none", **summary(R[p])))
if "perturb" in parts:
    for name in perturb.NAMES[1:]:
        for ps in (0, 1):
            pts = [p for p in ANCHORS if f"perturb|{key(*p)}|{name}|ps{ps}" not in done]
            if not pts: continue
            Xp, c, gi = perturb.apply(name, X0, FOV, ps); D = prep(Xp, XY[c], FOV[c], ids=D0["ids"][c]); R = score_points(D, pts)
            for p in pts: write(dict(id=f"perturb|{key(*p)}|{name}|ps{ps}", unit=u, part="perturb", features=p[0], ell=p[1], n_pc=p[2], alpha=str(p[3]), perturb=name, perturb_seed=ps, **summary(R[p])))


def lgbm_fp():
    cfg = max(json.loads(Path(f"runs/p3/f3_{u}_lightgbm.json").read_text())["grid"].items(), key=lambda kv: kv[1]["unit"])[0]
    lv, lr = int(cfg.split("_")[0][6:]), float(cfg.split("_lr")[1])

    def fp(Xtr, Ytr, Xte):
        from lightgbm import LGBMRegressor
        return np.column_stack([LGBMRegressor(num_leaves=lv, learning_rate=lr, n_estimators=100, n_jobs=4, random_state=0, verbose=-1).fit(Xtr, Ytr[:, j]).predict(Xte) for j in range(Ytr.shape[1])])
    return cfg, fp


PROT = None
if "protein" in parts and IF:
    lib = np.asarray(X0.sum(1)).ravel(); PROT = Protein(D0, a.obs[IF].values.astype(float), lib)
    for p in ANCHORS:
        pid = f"protein|{key(*p)}"
        if pid in done: continue
        _, O = score_points(D0, [p], full_out=True); pr = PROT(O[p])
        write(dict(id=pid, unit=u, part="protein", features=p[0], ell=p[1], n_pc=p[2], alpha=str(p[3]), perturb="none", **{k: v for k, v in pr.items() if k != "per_receiver"}, per_receiver=pr["per_receiver"]))
if "learned" in parts:
    for m in ("gat", "lightgbm"):
        f = Path(f"runs/p3/f3_{u}_{m}.json")
        if not f.exists(): continue
        for c, v in json.loads(f.read_text())["grid"].items():
            pid = f"learned|{m}|{c}|p3"
            if pid not in done: write(dict(id=pid, unit=u, part="learned", model=m, config=c, perturb="none", source="p3", top10=v["unit"], worst_top10=min(x["score"] for x in v["per_receiver"].values())))
    cfg, fp = lgbm_fp()
    for name, ps, shuf in (("none", 0, None), ("none", 0, 0), ("thin50", 0, None), ("sub50", 0, None), ("genes75", 0, None)):
        pid = f"learned|lightgbm|{cfg}|{name}{'|shuffled' if shuf is not None else ''}"
        if pid in done: continue
        t = time.time(); Xp, c, gi = perturb.apply(name, X0, FOV, ps); D = D0 if name == "none" else prep(Xp, XY[c], FOV[c], ids=D0["ids"][c])
        O = f3.score_floor(D, ell_um=25., n_pc=50, fit_predict=fp, full=True, shuffle_seed=shuf); per = panel(O)
        rec = dict(id=pid, unit=u, part="learned", model="lightgbm", config=cfg, perturb=name if shuf is None else "shuffled", perturb_seed=ps, **summary(per), per_receiver=per)
        if name == "none" and shuf is None and IF:
            if PROT is None: PROT = Protein(D0, a.obs[IF].values.astype(float), np.asarray(X0.sum(1)).ravel())
            pr = PROT(O); rec.update({k: v for k, v in pr.items() if k != "per_receiver"})
        rec["s"] = round(time.time() - t, 1); write(rec)
