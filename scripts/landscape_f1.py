"""F1 floor characterisation landscape (PLAN §14). Development units only. One JSON line per landscape point
(runs/landscape/f1_<unit>.jsonl, resumable), each with the full stress.f1 panel.
  grid     norm {logcp10k, log1p, sqrtcp10k, pearson} x scale {z, none} x n_pc {5, 10, 20, 30, 50, 100}; count-split seed 0
  refs     depth only, spatial smoothing (neighbour mean of the P3 floor, 25 um), random (the score a trivial embedding gets)
  split    count-split seeds 1..4 at the anchors (the score's own randomness)
  perturb  thin50, thin25, sub50, fov75, genes75 (perturbation seed 0) at the anchors + refs
  learned  P3 embeddings (runs/p3/f1_<unit>_*_emb.npz, trained on split seed 0) on the base panel; their perturbation runs,
           when present (runs/landscape/emb/f1_<unit>_<model>_<perturb>_s<seed>.npz), on the matching panel
Usage (env floors, OMP_NUM_THREADS=4): python scripts/landscape_f1.py <unit> [learned refs grid split perturb ...]"""
import glob, json, re, sys, time
from pathlib import Path
import numpy as np
from spatialfloors.floors.simple import neighbour_mean
from spatialfloors.provenance import stamp
from spatialfloors.stress import perturb
from spatialfloors.stress.f1 import load, floor_embed, Panel

u = sys.argv[1]; parts = sys.argv[2:] or ["learned", "refs", "grid", "split", "perturb"]; d0 = load(u)
OUT = Path("runs/landscape"); fo = OUT / f"f1_{u}.jsonl"; OUT.mkdir(parents=True, exist_ok=True)
done = {json.loads(l)["id"] for l in fo.read_text().splitlines()} if fo.exists() else set(); log = lambda *s: print(*s, flush=True)
NORM, SCALE, NPC = ("logcp10k", "log1p", "sqrtcp10k", "pearson"), ("z", "none"), (5, 10, 20, 30, 50, 100)
ANCHORS = [("logcp10k", "z", 50), ("logcp10k", "z", 30), ("pearson", "z", 50), ("sqrtcp10k", "z", 50)]


def write(rec):
    with fo.open("a") as f: f.write(json.dumps(stamp(rec), default=float) + "\n")
    done.add(rec["id"]); log(rec["id"], {k: round(v, 3) for k, v in rec.items() if isinstance(v, float) and k.startswith(("excess_k15", "genes_imp", "depth", "knn", "IF_PanCK"))}, rec.get("s"))


def point(pid, P, fn, **meta):
    if pid in done: return
    t = time.time(); E = fn(); write(dict(id=pid, unit=u, **meta, **P(E), s=round(time.time() - t, 1)))


def refs(P, Xa, xy, fov, tag, meta):
    zl = np.log1p(P.La); zl = (zl - zl.mean()) / zl.std()
    point(f"ref|depth_only|{tag}", P, lambda: np.c_[zl, np.ones(len(zl))], part="ref", ref="depth_only", **meta)
    point(f"ref|smooth25|{tag}", P, lambda: neighbour_mean(xy, fov, floor_embed(Xa[:, P.g], 50), 25.)[0], part="ref", ref="smooth25", **meta)
    point(f"ref|random|{tag}", P, lambda: np.random.default_rng(0).normal(size=(len(zl), 30)), part="ref", ref="random", **meta)


P0 = None
def base():
    global P0
    if P0 is None: P0 = Panel(d0, split_seed=0)
    return P0


if "learned" in parts:
    for f in sorted(glob.glob(f"runs/p3/f1_{u}_*_emb.npz")):
        if re.search(r"_n\d+_emb", f): continue
        m = Path(f).name[len(f"f1_{u}_"):-len("_emb.npz")]; z = np.load(f)
        for key in z.files:
            if key.endswith("_authors") or m.startswith("expr"): continue
            point(f"learned|{m}|{key}|none|split0", base(), lambda: z[key], part="learned", model=m.split("_")[0], config=key, perturb="none", split_seed=0)
    for f in sorted(glob.glob(f"runs/landscape/emb/f1_{u}_*.npz")):
        model, name, ps = Path(f).stem[len(f"f1_{u}_"):].rsplit("_", 2); ps = int(ps[1:]); pid = f"learned|{model}|{name}|ps{ps}"
        if pid in done: continue
        X, c, gg = perturb.apply(name, d0["X"], d0["fov"], ps) if name != "none" else (d0["X"], None, None)
        if name.startswith("genes"): X, gg = d0["X"], np.flatnonzero(np.isin(np.arange(len(d0["g"])), perturb.subsample(len(d0["g"]), .75, ps)))
        P = Panel(d0, 0, cells=c if name.startswith(("sub", "fov")) else None, genes=gg if name.startswith("genes") else None, X=X)
        point(pid, P, lambda: np.load(f)["emb"], part="learned", model=model, perturb=name, perturb_seed=ps, split_seed=0)
if "refs" in parts: refs(base(), base().Xa, d0["xy"], d0["fov"], "none|split0", dict(perturb="none", split_seed=0))
if "grid" in parts:
    for npc in NPC:
        for norm in NORM:
            for sc in SCALE:
                point(f"grid|{norm}|{sc}|pc{npc}|none|split0", base(), lambda: floor_embed(base().Xa[:, base().g], npc, norm, sc), part="grid", norm=norm, scale=sc, n_pc=npc, perturb="none", split_seed=0)
if "split" in parts:
    for ss in range(1, 5):
        P = Panel(d0, split_seed=ss)
        for norm, sc, npc in ANCHORS: point(f"split|{norm}|{sc}|pc{npc}|none|split{ss}", P, lambda: floor_embed(P.Xa[:, P.g], npc, norm, sc), part="split", norm=norm, scale=sc, n_pc=npc, perturb="none", split_seed=ss)
        refs(P, P.Xa, d0["xy"], d0["fov"], f"none|split{ss}", dict(perturb="none", split_seed=ss))
if "perturb" in parts:
    for name in perturb.NAMES[1:]:
        if name.startswith("genes"): c, gg, X = None, perturb.subsample(len(d0["g"]), .75, 0), None  # panel genes dropped (scored and embedded)
        else: X, c, _ = perturb.apply(name, d0["X"], d0["fov"], 0); gg = None; c = c if name.startswith(("sub", "fov")) else None
        P = Panel(d0, 0, cells=c, genes=gg, X=X); xy, fov = (d0["xy"], d0["fov"]) if c is None else (d0["xy"][c], d0["fov"][c])
        for norm, sc, npc in ANCHORS: point(f"perturb|{norm}|{sc}|pc{npc}|{name}|ps0", P, lambda: floor_embed(P.Xa[:, P.g], npc, norm, sc), part="perturb", norm=norm, scale=sc, n_pc=npc, perturb=name, perturb_seed=0, split_seed=0)
        refs(P, P.Xa, xy, fov, f"{name}|ps0", dict(perturb=name, perturb_seed=0, split_seed=0))
