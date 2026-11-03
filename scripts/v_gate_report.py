"""Validation gate report (PLAN §17, G1-G7) for the validation units -> runs/landscape/VALIDATION_GATE.md + .json.
G1 jobs done and results finite; G2 floor landscape: the locked floor vs the best eligible config on each validation unit
(eligible = the set the lock chose from: F1 / F2 / F3 grid, F5 InSituType at threshold 0 + logistic regression), within
EPS .02 = on the plateau, guards at the locked config vs their development range; if a floor leaves the plateau, the §14
rule over development + validation units (max of the worst-unit share of the best) names the lock-v2 candidate; G4 CPU
determinism; G5 reference QC and IF channel checks; G6 power per pair from development + validation D; G7 crop check
(whole slide vs crop D, |difference| <= margin). Usage: python scripts/v_gate_report.py <validation unit ...>"""
import glob, json, re, sys
from pathlib import Path
import numpy as np
from scipy import stats as st

V = sys.argv[1:]; EPS = .02; L = json.loads(Path("manifests/locked.json").read_text())["locked"]
DEV = ["gse1_232", "gse0_S0", "gse303070_148191A_118N", "gse303070_148191A_278N", "gse303070_148191A_278T", "gse303070_148191A_46T", "allen_animal4_C57BL6J-4.003"]
FAM = dict(f1=dict(primary="excess_k15", locked="grid|pearson|z|pc100|none|split0", guards=("excess_k5", "excess_k50", "genes_improved"), elig=lambda r: r["part"] == "grid"),
           f2=dict(primary="ari_bal", locked="grid|gmm_full|ell200|lam0.8|pc50", guards=("min_region_recall", "ari_boundary", "excess_over_blob", "seed_stability"), elig=lambda r: r["part"] == "grid"),
           f3=dict(primary="top10", locked="grid|T+ctx|ell25|pc50|acv", guards=("pc1", "gene_median", "worst_top10"), elig=lambda r: r["part"] == "grid"),
           f5=dict(primary="heldout_auroc", locked=None, guards=("worst_class", "depth_strat_mean", "excess_over_threshold"),
                   elig=lambda r: (r["part"] == "ist" and "|t0|" in r["id"] and r["id"].endswith("|none")) or (r["part"] == "logreg" and r["id"].endswith("|none"))))
rep, md = dict(units=V, G2={}), ["# Validation gate (PLAN §17)", ""]


def rows(f, u):
    p = Path(f"runs/landscape/{f}_{u}.jsonl")
    return {r["id"]: r for r in map(json.loads, open(p)) if FAM[f]["elig"](r) and r.get(FAM[f]["primary"]) == r.get(FAM[f]["primary"])} if p.exists() else {}


def locked_id(f, R):
    if FAM[f]["locked"]: return FAM[f]["locked"]
    return next((i for i in ("ist|all|negprb|t0|none", "ist|all|const05|t0|none") if i in R), None)  # P3: NegPrb means when the panel has them


md += ["## G2 floor landscape", "| family | unit | locked | best eligible | share | on plateau | guards below dev range |", "|---|---|---|---|---|---|---|"]
for f, c in FAM.items():
    allu = {u: rows(f, u) for u in DEV + V}; allu = {u: R for u, R in allu.items() if R}
    for u in V:
        R = allu.get(u)
        if not R: md.append(f"| {f} | {u} | (no landscape) | | | | |"); continue
        lid = locked_id(f, R); p = c["primary"]; b = max(R.values(), key=lambda r: r[p]); lv = R[lid][p] if lid in R else float("nan")
        dev_rng = {g: min((allu[d][locked_id(f, allu[d])].get(g, np.nan) for d in DEV if d in allu and locked_id(f, allu[d]) in allu[d]), default=np.nan) for g in c["guards"]}
        low = [g for g in c["guards"] if R.get(lid, {}).get(g, np.nan) < dev_rng[g] - 1e-9]
        ok = bool(b[p] - lv <= EPS); rep["G2"][f"{f}|{u}"] = dict(locked=lid, locked_value=lv, best=b["id"], best_value=b[p], on_plateau=ok, guards_below_dev=low)
        md.append(f"| {f} | {u} | {lv:.3f} | {b[p]:.3f} ({b['id']}) | {lv / b[p]:.3f} | {ok} | {', '.join(low) or '-'} |")
    us = list(allu); ids = set.intersection(*(set(allu[u]) for u in us)) if us else set(); best = {u: max(r[c['primary']] for r in allu[u].values()) for u in us}
    if ids and any(not v["on_plateau"] for k, v in rep["G2"].items() if k.startswith(f + "|")):
        key = lambda i: min(allu[u][i][c["primary"]] / best[u] for u in us); cand = max(ids, key=key)
        rep["G2"][f"{f}|lock_v2_candidate"] = dict(config=cand, worst_share=key(cand)); md.append(f"| {f} | lock-v2 candidate (§14 rule, dev + validation) | {cand} | worst share {key(cand):.3f} | | | |")
md += ["", "## G4 determinism"] + [f"- {Path(x).stem}: {'pass' if Path(x).suffix == '.done' else 'FAIL'}" for x in sorted(glob.glob("runs/sched/G4*_*.done") + glob.glob("runs/sched/G4*_*.fail"))]
md += ["", "## G5 references and IF"]
for f in sorted(glob.glob("manifests/refs/*.json")):
    r = json.loads(Path(f).read_text()); md.append(f"- {r['tissue']}: {r['source']}; classes {r['class_counts']}; QC {'pass' if r['qc_pass'] else 'FAIL'}; markers {r['heldout_markers']}")
for f in sorted(glob.glob("runs/p2/v_bruker_*_convert.log")):
    md += [f"- {Path(f).stem}: " + l.strip() for l in open(f) if l.startswith("IF check")]
# G6 / G7 from P5-style per-unit D (runs/p3: development + validation)
p5 = Path("runs/p3/P5_gate.json")
if p5.exists():
    P = json.loads(p5.read_text())["pairs"]; pr = Path("runs/v_retro/P5_retro.json")  # F4 retrospective D (N2 L6 + camera-ready) for G6
    P.update(json.loads(pr.read_text())["pairs"] if pr.exists() else {}); NC = dict(F1=22 + 3, F2=10, F3=21 + 3, F5=17 + 3, F4=10)  # confirmation units per family (GSE303070 = 11 donors; +3 Bruker for F1 / F3 / F5; F4 retro + prospective)
    md += ["", "## G6 power (development + validation D; confirmation N)", "| pair | units | mean D | sd | margin | N conf | power |", "|---|---|---|---|---|---|---|"]
    for k, v in P.items():
        f = k.split("|")[0].split("_")[0]; D = np.array([d for u, d in v["per_unit"].items() if not u.endswith("w")]); n = NC.get(f)
        if n is None or len(D) < 3: continue
        sd = D.std(ddof=1); se = sd / np.sqrt(n); df = n - 3; tq = st.t.ppf(.95, df); pw = float(st.nct.cdf(-tq, df, (D.mean() - v["margin"]) / se)) if se > 0 else float("nan")
        md.append(f"| {k} | {len(D)} | {D.mean():+.4f} | {sd:.4f} | {v['margin']:.4f} | {n} | {pw:.2f} |"); rep.setdefault("G6", {})[k] = dict(mean=float(D.mean()), sd=float(sd), power=pw)
    md += ["", "## G7 crop check (gse1_231 whole slide vs crop)", "| pair | D crop | D whole | |diff| | margin | pass |", "|---|---|---|---|---|---|"]
    for k, v in P.items():
        a, b = v["per_unit"].get("gse1_231"), v["per_unit"].get("gse1_231w")
        if a is None or b is None: continue
        ok = abs(b - a) <= v["margin"]; rep.setdefault("G7", {})[k] = dict(crop=a, whole=b, passed=bool(ok))
        md.append(f"| {k} | {a:+.4f} | {b:+.4f} | {abs(b - a):.4f} | {v['margin']:.4f} | {ok} |")
# G7 for F4 (gse1_231): our 50k crop (runs/f4/f4_231.csv, seed mean) vs the whole slide (N2 L6 opt/crc_231.csv + Cellina camera-ready)
import pandas as pd
c4 = Path("runs/f4/f4_231.csv"); w6 = Path("/Users/Hemato/project_3/antig/N2/trials/perturb/opt/crc_231.csv"); cam = Path("external/cellina-reproducibility/results/loo_summary_crc_DEG_50_camready.csv")
if c4.exists() and w6.exists():
    d = pd.read_csv(c4); W = pd.read_csv(w6); W = W[W.model_name == "L6_twoway_pseudobulk"].pearson.mean(); C = pd.read_csv(cam); C = C[C.sid == "crc_231"]
    md += ["", "## G7 F4 crop check (gse1_231, official pearson)", "| model | crop: learned - L6 | whole: learned - L6 | |diff| | margin | pass |", "|---|---|---|---|---|---|"]
    for m in ("cellina", "cpa"):
        dc = d[d.model == m].groupby("holdout").pearson.mean().mean() - d[d.model == "L6"].pearson.mean(); dw = C[C.model_name == m].pearson.mean() - W
        ok = abs(dc - dw) <= L["F4"]["margin"]; rep.setdefault("G7", {})[f"F4_pearson|{m}"] = dict(crop=float(dc), whole=float(dw), passed=bool(ok))
        md.append(f"| {m} | {dc:+.4f} | {dw:+.4f} | {abs(dc - dw):.4f} | {L['F4']['margin']:.3f} | {ok} |")
    l6c = d[d.model == "L6"].pearson.mean(); md.append(f"| L6 alone (crop vs whole, same code) | {l6c:.4f} | {W:.4f} | {abs(l6c - W):.4f} | | |")
    md.append("Learned whole-slide values are the authors' camera-ready runs; our P0 MPS re-runs of Cellina on crc_232 sat .03-.11 below camera-ready per type, so this row mixes crop and hardware / seed effects.")
# Audit: CellCharter configs on the F2 primary (balanced ARI; configs were chosen on ARI before row 359), development units
cc = {}
for f in glob.glob("runs/p3/relabel/labels/f2_*_cellcharter.npz") + glob.glob("runs/p3_ccbal*/labels/f2_*_cellcharter.npz"):
    uu = Path(f).name[3:-len("_cellcharter.npz")]
    if uu not in DEV: continue
    import anndata as ad
    from spatialfloors.metrics.f2 import scored_regions
    from spatialfloors.stress.f2 import Panel as F2P
    a_ = ad.read_h5ad(f"runs/p2/cache/{uu}.h5ad", backed="r"); mk, rg = scored_regions(a_.obs.region.astype(str).values, json.loads(Path("manifests/crops/index.json").read_text())[uu]["present"]); a_.file.close()
    P_ = F2P.__new__(F2P); P_.y = rg; z = np.load(f, allow_pickle=True)
    for c in z.files: Lb = np.asarray(z[c]); Lb = Lb[None] if Lb.ndim == 1 else Lb; cc.setdefault(uu, {})[c] = float(np.mean([P_.ari_bal(l[mk].astype(str)) for l in Lb]))
if cc:
    md += ["", "## Audit: CellCharter config on balanced ARI (development)"] + [f"- {uu}: " + ", ".join(f"{c} {v:.3f}" for c, v in sorted(d_.items(), key=lambda kv: -kv[1])) + f" (locked {L['F2']['cellcharter']['mouse' if uu.startswith('allen') else 'human']['config']})" for uu, d_ in cc.items()]
    rep["cellcharter_ari_bal"] = cc
Path("runs/landscape/VALIDATION_GATE.json").write_text(json.dumps(rep, indent=1, default=float)); Path("runs/landscape/VALIDATION_GATE.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
