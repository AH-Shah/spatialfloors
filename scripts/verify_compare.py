"""Compare the release re-runs (runs/verify, scripts/verify_dev.sh) with the development outputs: grid entries (scores,
per receiver, shares) must be equal; timing ('s', 'fit_s') and provenance stamps are ignored. Saved partitions / labels /
embeddings must be bit-identical. Usage: python scripts/verify_compare.py -> stdout, exit 1 on any difference"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

V = Path("runs/verify"); G, H, A = "gse303070_148191A_278N", "gse303070_148191A_118N", "allen_animal4_C57BL6J-4.003"
SKIP = {"s", "fit_s", "spatialfloors_commit", "spatialfloors_dirty", "created", "host", "python", "train", "info"}


def strip(x):
    if isinstance(x, dict): return {k: strip(v) for k, v in x.items() if k not in SKIP}
    if isinstance(x, list): return [strip(v) for v in x]
    return x


def grid(new, ref, keys=None):
    a, b = json.loads(Path(new).read_text())["grid"], json.loads(Path(ref).read_text())["grid"]
    ks = keys or sorted(set(a) & set(b)); return all(strip(a[k]) == strip(b[k]) for k in ks), ks


def arrays(new, ref, keys=None):
    a, b = np.load(new, allow_pickle=True), np.load(ref, allow_pickle=True); ks = keys or sorted(set(a.files) & set(b.files))
    return all(np.array_equal(a[k], b[k]) for k in ks), ks


checks = [
 ("F1 floor 278N json", lambda: (strip(json.loads((V / f"f1_{G}_expr_embed_pearson_none_pc50.json").read_text())["scores"]) == strip(json.loads(Path(f"runs/p3v2/f1_{G}_expr_embed_pearson_none_pc50.json").read_text())["scores"]), ["scores"])),
 ("F1 floor 278N embedding", lambda: arrays(V / f"f1_{G}_expr_embed_pearson_none_pc50_emb.npz", f"runs/p3v2/f1_{G}_expr_embed_pearson_none_pc50_emb.npz")),
 ("F1 floor 118N embedding", lambda: arrays(V / f"f1_{H}_expr_embed_pearson_none_pc50_emb.npz", f"runs/p3v2/f1_{H}_expr_embed_pearson_none_pc50_emb.npz")),
 ("F1 v1 floor 278N embedding", lambda: arrays(V / f"f1_{G}_expr_embed_pearson_z_pc100_emb.npz", f"runs/p3/f1_{G}_expr_embed_pearson_z_pc100_emb.npz")),
 ("F1 rescored primary", lambda: (lambda n, r: (bool(np.array_equal(n.score.values, r.score.values)) and len(n) == 2, ["score"]))(
     *[d[(d.model == "expr_embed") & (d.config == "pearson_none_pc50") & d.unit.isin([G, H])].sort_values("unit") for d in (pd.read_csv(V / "f1_rescored.csv"), pd.read_csv("runs/p3v2/f1_rescored.csv"))])),
 ("F2 floor gse1_232", lambda: grid(V / "f2_gse1_232_floor.json", "runs/p3v2/f2_gse1_232_floor.json")),
 ("F2 floor gse1_232 partitions", lambda: arrays(V / "labels/f2_gse1_232_floor.npz", "runs/p3v2/labels/f2_gse1_232_floor.npz")),
 ("F2 floor allen 4.003 partitions", lambda: arrays(V / f"labels/f2_{A}_floor.npz", f"runs/p3v2/labels/f2_{A}_floor.npz")),
 ("F2 v3 secondary gse1_232", lambda: grid(V / "f2v3/f2_gse1_232_floor.json", "runs/v3/p3f2test/f2_gse1_232_floor.json")),
 ("F3 floor v3 278N", lambda: grid(V / f"f3_v3/f3_{G}_ridge.json", f"runs/v3/p3f3test/f3_{G}_ridge.json")),
 ("F3 floor v2 278N", lambda: grid(V / f"v2/f3_{G}_ridge.json", f"runs/p3v2/f3_{G}_ridge.json")),
 ("F3 LightGBM v3 278N", lambda: grid(V / f"f3_{G}_lightgbm.json", f"runs/v3/lgbm/f3_{G}_lightgbm.json", ["leaves15_lr0.05"])),
 ("F5 floors 278N", lambda: grid(V / f"f5_{G}_floors.json", f"runs/v3/p3f5test/f5_{G}_floors.json")),
 ("F5 floors 278N labels", lambda: arrays(V / f"labels/f5_{G}_floors.npz", f"runs/v3/p3f5test/labels/f5_{G}_floors.npz")),
 ("F5 floors gse1_232 (gate path)", lambda: grid(V / "f5_gse1/f5_gse1_232_floors.json", "runs/p3v2/f5_gse1_232_floors.json")),
 ("F5 tier-2 prior beta 1 = locked floor", lambda: (strip(json.loads((V / f"t2_prior/f5_{G}_prior.json").read_text())["grid"]["prior_b1"]) == strip(json.loads((V / f"f5_{G}_floors.json").read_text())["grid"]["insitutype_spatial"]), ["prior_b1"])),
]
bad = 0
for name, fn in checks:
    try: ok, ks = fn()
    except Exception as e: ok, ks = False, [f"{type(e).__name__}: {e}"]
    bad += not ok; print(("SAME " if ok else "DIFF ") + name, ks if not ok or len(ks) < 6 else f"{len(ks)} keys")
print("all identical" if not bad else f"{bad} differ"); sys.exit(1 if bad else 0)
