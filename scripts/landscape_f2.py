"""F2 floor characterisation landscape (PLAN §14). Development units only. One JSON line per landscape point
(runs/landscape/f2_<unit>.jsonl, resumable: points already written are skipped), each with the full §13 panel.
  grid     ell {10, 25, 50, 100, 200, 400} um x lam {0, .2, .5, .8, 1} x n_pc {10, 30, 50} x clusterer {kmeans, gmm_full,
           gmm_diag}, 5 clustering seeds; lam 0 (expression only) once per n_pc / clusterer; fine-k detail at n_pc 30
  k        true k + {-2, -1, +1, +2} at the anchor points (how the partition degrades when k is wrong)
  perturb  thin50, thin25, sub50, fov75, genes75 x 3 perturbation seeds x 3 clustering seeds at the anchor points
  learned  P3 partitions of the learned models (runs/p3*/labels) on the same panel; their perturbation runs, when present
           (runs/landscape/labels/f2_<unit>_<model>_<perturb>_s<seed>.npz), on the matching perturbed panel
Usage (env cellcharter, OMP_NUM_THREADS=4): python scripts/landscape_f2.py <unit> [grid|k|perturb|learned ...]"""
import glob, json, sys, time
from pathlib import Path
import numpy as np
from sklearn.mixture import GaussianMixture
from spatialfloors.floors.simple import expr_pca, neighbour_mean
from spatialfloors.metrics.f2 import kmeans_labels, gmm_labels
from spatialfloors.provenance import stamp
from spatialfloors.stress import perturb
from spatialfloors.stress.f2 import load, Panel

u = sys.argv[1]; parts = sys.argv[2:] or ["grid", "k", "perturb", "learned"]; d0 = load(u); P0 = Panel(d0)
OUT = Path("runs/landscape"); OUT.mkdir(parents=True, exist_ok=True); fo = OUT / f"f2_{u}.jsonl"
done = {json.loads(l)["id"] for l in fo.read_text().splitlines()} if fo.exists() else set(); log = lambda *s: print(*s, flush=True)
ELL, LAM, NPC, SEEDS = (10., 25., 50., 100., 200., 400.), (0., .2, .5, .8, 1.), (10, 30, 50), range(5)
CLUS = {"kmeans": kmeans_labels, "gmm_full": gmm_labels,
        "gmm_diag": lambda E, k, s: GaussianMixture(k, covariance_type="diag", reg_covar=1e-4, random_state=s).fit(E).predict(E)}
ANCHORS = [("kmeans", 25., .5, 30), ("gmm_full", 100., .5, 30), ("gmm_full", 100., .8, 30), ("gmm_full", 200., .8, 30), ("gmm_full", 200., .5, 30)]


def zb(Mx): Mx = Mx - Mx.mean(0); return Mx / max(float(np.sqrt((Mx ** 2).sum(1).mean())), 1e-6)


class Feats:
    """Floor features for one (perturbed) input, cached per n_pc and (n_pc, ell)."""
    def __init__(self, X, xy, fov): self.X, self.xy, self.fov, self.E, self.C = X, xy, fov, {}, {}

    def __call__(self, npc, ell, lam):
        if npc not in self.E: self.E[npc] = zb(expr_pca(self.X, npc))
        if lam == 0: return self.E[npc]
        if (npc, ell) not in self.C: self.C[npc, ell] = zb(neighbour_mean(self.xy, self.fov, self.E[npc], ell)[0])
        return self.C[npc, ell] if lam == 1 else np.hstack([np.sqrt(1 - lam) * self.E[npc], np.sqrt(lam) * self.C[npc, ell]]).astype(np.float32)


def write(rec):
    with fo.open("a") as f: f.write(json.dumps(stamp(rec), default=float) + "\n")
    done.add(rec["id"]); log(rec["id"], {k: round(rec[k], 3) for k in ("ari", "ami", "matched_f1", "min_region_recall", "excess_over_blob", "seed_stability") if k in rec}, rec.get("s"))


def point(pid, F, P, clus, k, seeds, detail=False, **meta):
    if pid in done: return
    t = time.time(); L = np.stack([CLUS[clus](F, k, s) for s in seeds]); r = dict(id=pid, unit=u, clusterer=clus, k=k, **meta, **P(L))
    if detail: r.update(P.detail(F, CLUS[clus]))
    r["s"] = round(time.time() - t, 1); write(r)


F0 = Feats(d0["X"], d0["xy"], d0["fov"])
if "grid" in parts:
    for npc in NPC:
        for clus in CLUS:
            for ell in ELL:
                for lam in LAM:
                    if lam == 0 and ell != ELL[0]: continue
                    point(f"grid|{clus}|ell{ell:g}|lam{lam:g}|pc{npc}", F0(npc, ell, lam), P0, clus, P0.k, SEEDS, detail=npc == 30,
                          part="grid", ell=ell if lam else None, lam=lam, n_pc=npc, perturb="none")
if "k" in parts:
    for clus, ell, lam, npc in ANCHORS:
        for dk in (-2, -1, 1, 2):
            if P0.k + dk >= 2: point(f"k|{clus}|ell{ell:g}|lam{lam:g}|pc{npc}|dk{dk:+d}", F0(npc, ell, lam), P0, clus, P0.k + dk, SEEDS, part="k", ell=ell, lam=lam, n_pc=npc, dk=dk, perturb="none")
if "perturb" in parts:
    for name in perturb.NAMES[1:]:
        for ps in range(3):
            X, c, g = perturb.apply(name, d0["X"], d0["fov"], ps); P = P0 if len(c) == len(d0["xy"]) else Panel(d0, c); Fp = Feats(X, d0["xy"][c], d0["fov"][c])
            for clus, ell, lam, npc in ANCHORS:
                point(f"perturb|{clus}|ell{ell:g}|lam{lam:g}|pc{npc}|{name}|ps{ps}", Fp(npc, ell, lam), P, clus, P.k, range(3), part="perturb", ell=ell, lam=lam, n_pc=npc, perturb=name, perturb_seed=ps)
if "learned" in parts:
    for f in sorted(set(glob.glob(f"runs/p3*/labels/f2_{u}_*.npz") + glob.glob(f"runs/p3/relabel/labels/f2_{u}_*.npz"))):
        if "floor" in f: continue
        z = np.load(f, allow_pickle=True); model = Path(f).stem.split("_")[-1]
        for c in z.files:
            pid = f"learned|{model}|{c}|none"
            if pid not in done: t = time.time(); write(dict(id=pid, unit=u, part="learned", model=model, config=c, perturb="none", **P0(z[c]), s=round(time.time() - t, 1)))
    for f in sorted(glob.glob(f"runs/landscape/labels/f2_{u}_*.npz")):  # learned perturbation / seed runs
        z = np.load(f, allow_pickle=True); model, name, ps = Path(f).stem[len(f"f2_{u}_"):].rsplit("_", 2); ps = int(ps[1:])
        c = perturb.apply(name, d0["X"], d0["fov"], ps)[1] if name != "none" else np.arange(len(d0["xy"])); P = P0 if len(c) == len(d0["xy"]) else Panel(d0, c)
        pid = f"learned|{model}|{name}|ps{ps}"
        if pid not in done: t = time.time(); write(dict(id=pid, unit=u, part="learned", model=model, perturb=name, perturb_seed=ps, **P(z["labels"]), s=round(time.time() - t, 1)))
