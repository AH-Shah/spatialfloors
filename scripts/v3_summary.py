"""PLAN §19 acceptance tables for the v3 floor candidates (development units; runs/v3). Usage: python3 scripts/v3_summary.py [f1 f2 f3 f4 f5]"""
import glob, json, sys
from pathlib import Path
import numpy as np

FAM = sys.argv[1:] or ["f1", "f2", "f3", "f4", "f5"]; sh = lambda u: u.replace("gse303070_148191A_", "").replace("allen_animal4_C57BL6J-", "a")
J = lambda p: json.loads(Path(p).read_text())
HUMAN = ["118N", "278N", "278T", "46T", "gse0_S0", "gse1_232"]
if "f3" in FAM:
    print("## F3 unit mean (T / lock = v2 / v3 / v3_null), v3 > lock receivers")
    R = {sh(J(f)["unit"]): J(f) for f in sorted(glob.glob("runs/v3/f3_*.json"))}
    for u, r in R.items(): m = r["unit_mean"]; print(f"  {u:9} " + " ".join(f"{k} {v:.4f}" for k, v in m.items()) + f" | {r['n_v3_gt_lock']}/{r['n_receivers']}")
    H = [u for u in HUMAN if u in R]
    if H:
        d = np.array([R[u]["unit_mean"]["v3"] - R[u]["unit_mean"]["lock"] for u in H]); dn = np.array([R[u]["unit_mean"]["v3_null"] - R[u]["unit_mean"]["lock"] for u in H])
        print(f"  human ({len(H)}): lock {np.mean([R[u]['unit_mean']['lock'] for u in H]):.4f} v3 {np.mean([R[u]['unit_mean']['v3'] for u in H]):.4f} (+{d.mean():.4f}); "
              f"v3 >= lock on {int((d >= 0).sum())}/{len(H)}; null - lock mean {dn.mean():+.4f}")
if "f5" in FAM:
    print("## F5 held-out AUROC | marker-positive recall e/i/m | gated F1 ; acceptance vs v2")
    for f in sorted(glob.glob("runs/v3/f5_*.json")):
        r = J(f); s = r["summaries"]; v2 = s["v2"]; print(f"  {sh(r['unit'])} (v2 rerun {r['v2_reproduced']:.3f}, Python replica {r.get('python_replica_agreement', float('nan')):.3f})")
        for k, x in s.items():
            dr = min(x["marker_pos_recall"][c] - v2["marker_pos_recall"][c] for c in v2["marker_pos_recall"])
            print(f"    {k:22} {x['heldout_auroc']:.4f} ({x['heldout_auroc'] - v2['heldout_auroc']:+.4f}) | " + "/".join(f"{x['marker_pos_recall'][c]:.2f}" for c in x["marker_pos_recall"])
                  + f" (worst {dr:+.3f}) | {x.get('gated_macro_f1', float('nan')):.3f}")
if "f4" in FAM:
    import pandas as pd
    print("## F4 crc mean over tasks: official / spike-free pearson; worst task vs L6 (official)")
    for f in sorted(glob.glob("runs/v3/f4/*_v3.csv")):
        d = pd.read_csv(f); t = d.pivot_table(index="model_name", columns="holdout_celltype", values="pearson"); c = d.pivot_table(index="model_name", columns="holdout_celltype", values="pearson_clean")
        print(" ", Path(f).stem)
        for m in t.index: print(f"    {m:22} {t.loc[m].mean():.3f} / {c.loc[m].mean():.3f} | worst {float((t.loc[m] - t.loc['L6_twoway_pseudobulk']).min()):+.3f}")
if "f2" in FAM:
    print("## F2 balanced ARI per unit (no Potts / + Potts), mean")
    for d_ in ("runs/v3", "runs/v3/f2r2"):
        R = {sh(J(f)["unit"]): J(f)["partitions"] for f in sorted(glob.glob(f"{d_}/f2_*.json"))}
        if not R: continue
        names = sorted({k for r in R.values() for k in r}); print(" ", d_, list(R))
        for k in names: v = [R[u][k]["ari_bal"] for u in R if k in R[u]]; print(f"    {k:30} " + " ".join(f"{x:.3f}" for x in v) + f" | mean {np.mean(v):.3f}")
if "f1" in FAM:
    print("## F1 excess_k15 | rarest-quartile recall | within-type dhom (v3 candidates, runs/probe/f1_<u>_v3.json)")
    for f in sorted(glob.glob("runs/probe/f1_*_v3.json")):
        r = J(f); print(" ", sh(r["unit"]))
        for k, x in r["embeddings"].items(): print(f"    {k:22} {x['excess_k15']:.4f} | {x['recall_by_freq_quartile'][0]:.3f} | {x['depth_homophily_within_type']:.3f}")
