"""F5 floor characterisation landscape (PLAN §14). Development units only. One JSON line per landscape point
(runs/landscape/f5_<unit>.jsonl, resumable), each with the stress.f5 panel.
  ist      InSituType ML: reference profile {all cells (P3), stratified subsample seeds 0-2, coarse (fine profiles averaged
           within each coarse class)} x background {NegPrb mean (P3), constant .05, 2 x NegPrb} x posterior threshold
           {0, .5, .8, .95} (below it: unassigned -> "other"; call rate vs AUROC, row O2)
  logreg   C {.01, .1, 1, 10} x scaling {own, ref} x reference subsample seed {0, 1, 2}
  perturb  thin50, thin25, sub50, genes75 x seeds {0, 1} at the anchors
  learned  P3 scANVI / Tangram calls (runs/p3*/labels/f5_<unit>_*.npz); perturbation runs when present
           (runs/landscape/labels/f5_<unit>_<model>_<perturb>_s<seed>.npz, key labels)
Usage (env floors): python scripts/landscape_f5.py <unit> [learned ist logreg perturb ...]"""
import glob, json, sys, time
from pathlib import Path
import numpy as np
from spatialfloors.floors import f5
from spatialfloors.provenance import stamp
from spatialfloors.stress.f5 import load, inputs, Panel

u = sys.argv[1]; parts = sys.argv[2:] or ["learned", "ist", "logreg", "perturb"]; d0 = load(u)
OUT = Path("runs/landscape"); fo = OUT / f"f5_{u}.jsonl"; OUT.mkdir(parents=True, exist_ok=True); WD = OUT / "ist" / u
done = {json.loads(l)["id"] for l in fo.read_text().splitlines()} if fo.exists() else set(); log = lambda *s: print(*s, flush=True)
THR, CS = (0., .5, .8, .95), (.01, .1, 1., 10.)
PROF_ALL = f5.reference(d0["genes"], rows="all"); REF = {}


def write(rec):
    with fo.open("a") as f: f.write(json.dumps(stamp(rec), default=float) + "\n")
    done.add(rec["id"]); log(rec["id"], {k: round(v, 3) for k, v in rec.items() if k in ("heldout_auroc", "worst_class", "depth_strat_mean", "excess_over_threshold", "gated_macro_f1", "share_other")}, rec.get("s"))


def ref(genes, s):
    if (len(genes), s) not in REF: REF[len(genes), s] = f5.reference(genes, seed=s)
    return REF[len(genes), s]


def profile(kind, genes):
    """(genes x labels profile, names) on `genes` for a profile kind."""
    gi = [d0["genes"].index(g) for g in genes]; P, names = np.asarray(PROF_ALL[0])[gi], list(PROF_ALL[1])
    if kind == "all": return P, names
    if kind.startswith("strat"):
        X, y = ref(genes, int(kind[-1])); nm = np.unique(y); return np.column_stack([np.asarray(X[y == n].mean(0)).ravel() for n in nm]), list(nm)
    grp = {}
    for i, n in enumerate(names): grp.setdefault(f5.COARSE.get(n, n), []).append(i)
    return np.column_stack([P[:, ix].mean(1) for ix in grp.values()]), [names[ix[0]] for ix in grp.values()]  # named by a member: coarse() maps it


def ist_points(tag, Xq, genes, neg, P_, profiles=("all", "strat_s0", "strat_s1", "strat_s2", "coarse"), negs=("negprb", "const05", "negprb_x2"), thr=THR, **meta):
    for kind in profiles:
        for nk in negs:
            ids = [f"ist|{kind}|{nk}|t{t:g}|{tag}" for t in thr]
            if all(i in done for i in ids): continue
            if nk != "const05" and neg is None: continue
            t0 = time.time(); prof, names = profile(kind, genes)
            nv = None if nk == "const05" else (neg if nk == "negprb" else 2 * neg)
            lab, prob = f5.insitutype(Xq, prof, names, nv, workdir=WD / f"{kind}_{nk}_{tag.replace('|', '_')}"); s = round(time.time() - t0, 1)
            for t, pid in zip(thr, ids):
                if pid not in done: write(dict(id=pid, unit=u, part=meta.get("part", "ist"), method="insitutype", profile=kind, neg=nk, threshold=t, **meta.get("m", {}),
                                               **P_(np.where(prob >= t, lab, "unassigned")), s=s))


def lr_points(tag, Xq, genes, P_, cs=CS, scs=("own", "ref"), seeds=(0, 1, 2), **meta):
    for s in seeds:
        Xr, yr = ref(genes, s)
        for C in cs:
            for sc in scs:
                pid = f"logreg|C{C:g}|{sc}|s{s}|{tag}"
                if pid in done: continue
                t0 = time.time(); write(dict(id=pid, unit=u, part=meta.get("part", "logreg"), method="logreg", C=C, scaling=sc, ref_seed=s, **meta.get("m", {}),
                                             **P_(f5.logreg(Xr, yr, Xq, C, sc)), s=round(time.time() - t0, 1)))


Xq, genes, Xf, vf, neg, gate, _ = inputs(d0); P0 = Panel(Xq, genes, Xf, vf, gate, *profile("all", genes))
if "learned" in parts:
    for f in sorted(set(glob.glob(f"runs/p3*/labels/f5_{u}_*.npz") + glob.glob(f"runs/p3/relabel/labels/f5_{u}_*.npz"))):
        if "floors" in f: continue
        z = np.load(f, allow_pickle=True); model = Path(f).stem.split("_")[-1]
        for c in z.files:
            pid = f"learned|{model}|{c}|none"
            if pid not in done: t0 = time.time(); write(dict(id=pid, unit=u, part="learned", model=model, config=c, perturb="none", **P0(z[c]), s=round(time.time() - t0, 1)))
    for f in sorted(glob.glob(f"runs/landscape/labels/f5_{u}_*.npz")):
        model, name, ps = Path(f).stem[len(f"f5_{u}_"):].rsplit("_", 2); ps = int(ps[1:]); pid = f"learned|{model}|{name}|ps{ps}"
        if pid in done: continue
        q, g, xf, v, _, gt, _ = inputs(d0, name, ps); P = Panel(q, g, xf, v, gt, *profile("all", g))
        write(dict(id=pid, unit=u, part="learned", model=model, perturb=name, perturb_seed=ps, **P(np.load(f, allow_pickle=True)["labels"])))
if "ist" in parts: ist_points("none", Xq, genes, neg, P0, m=dict(perturb="none"))
if "logreg" in parts: lr_points("none", Xq, genes, P0, m=dict(perturb="none"))
if "perturb" in parts:
    for name in ("thin50", "thin25", "sub50", "genes75"):
        for ps in (0, 1):
            q, g, xf, v, nq, gt, _ = inputs(d0, name, ps); P = Panel(q, g, xf, v, gt, *profile("all", g)); tag = f"{name}|ps{ps}"; m = dict(perturb=name, perturb_seed=ps)
            ist_points(tag, q, g, nq, P, profiles=("all",), negs=("negprb",) if nq is not None else ("const05",), thr=(0., .8), part="perturb", m=m)
            lr_points(tag, q, g, P, cs=(.1, 1.), scs=("ref",), seeds=(0,), part="perturb", m=m)
