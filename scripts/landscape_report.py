"""PLAN §14 landscape report: runs/landscape/<fam>_<unit>.jsonl -> runs/landscape/LANDSCAPE_<FAM>.md (markdown tables).
Per unit: (a) profile of anchors, references and learned models; (b) shape: per-knob marginals of the primary metric over the
grid and the plateau (grid points within EPS of the best primary whose every guard metric is >= its grid 75th percentile);
(c) stability: spread of each metric over seeds / perturbations at each anchor (and for the learned models);
(d) metric validity: Spearman over grid points of the primary with every other panel metric, and the grid points where
the primary is in its top 10% while a guard metric is in its bottom half (primary rewarding a loss of detail).
Usage: python scripts/landscape_report.py f1|f2|f3|f5 [units ...]"""
import glob, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

fam = sys.argv[1]; EPS = .02
CFG = {"f2": dict(primary=__import__("os").environ.get("SF_PRIMARY", "ari"), knobs=("clusterer", "ell", "lam", "n_pc"),
                  metrics=("ari", "ari_bal", "ami", "matched_f1", "min_region_recall", "ari_boundary", "ari_interior", "ari_rotated_null", "excess_over_blob",
                           "size_entropy", "largest_share", "seed_stability", "ari_sd", "ami_fine_regions", "ami_celltype", "nesting_homogeneity"),
                  guards=("min_region_recall", "ari_boundary", "excess_over_blob", "ami_fine_regions", "seed_stability"),
                  anchor_key=lambda r: f"{r['clusterer']}|ell{r['ell']:g}|lam{r['lam']:g}|pc{r['n_pc']}"),
       "f1": dict(primary="excess_k15", knobs=("norm", "scale", "n_pc"),
                  metrics=("excess_k5", "excess_k15", "excess_k50", "top10_gene_share", "genes_improved", "depth_leak_r2", "knn_macro_recall",
                           "knn_rare_recall", "kmeans_ami_types", "IF_PanCK_nbr_spearman", "IF_CD45_nbr_spearman", "IF_CD68_CK8_18_nbr_spearman"),
                  guards=("excess_k5", "excess_k50", "genes_improved", "knn_rare_recall", "IF_PanCK_nbr_spearman"),
                  anchor_key=lambda r: f"{r['norm']}|{r['scale']}|pc{r['n_pc']}"),
       "f5": dict(primary="heldout_auroc", knobs=("method", "profile", "neg", "threshold", "C", "scaling"),
                  metrics=("heldout_auroc", "worst_class", "gated_macro_f1", "depth_strat_mean", "excess_over_threshold", "top_depth_share",
                           "share_epithelial", "share_myeloid", "share_immune_nonmyeloid", "share_other"),
                  guards=("worst_class", "gated_macro_f1", "depth_strat_mean", "excess_over_threshold"),
                  anchor_key=lambda r: (f"ist|{r['profile']}|{r['neg']}|t{r['threshold']:g}" if r.get("method") == "insitutype" else f"logreg|C{r['C']:g}|{r['scaling']}|s{r['ref_seed']}")),
       "f3": dict(primary="top10", knobs=("features", "ell", "n_pc", "alpha"),
                  metrics=("top10", "pc1", "top3", "gene_median", "r2", "worst_top10", "above_null", "prot_true_sig", "prot_pred_sig", "prot_pred_r"),
                  guards=("pc1", "gene_median", "worst_top10"),
                  anchor_key=lambda r: f"{r['features']}|ell{r['ell']:g}|pc{r['n_pc']}|a{r['alpha']}")}[fam]
units = sys.argv[2:] or [Path(f).stem[len(fam) + 1:] for f in sorted(glob.glob(f"runs/landscape/{fam}_*.jsonl"))]
f3 = lambda v: "-" if v is None or v != v else f"{v:.3f}"
out = [f"# {fam.upper()} floor landscape (PLAN §14)", f"primary {CFG['primary']}; guards {', '.join(CFG['guards'])}; plateau eps {EPS}", ""]
for u in units:
    R = [json.loads(l) for l in Path(f"runs/landscape/{fam}_{u}.jsonl").read_text().splitlines()]
    if fam == "f3":  # GAT perturbation runs (landscape_learned.py) + anchor lines carry above_null / protein onto the grid rows
        R += [dict(json.loads(Path(f).read_text()), part="learned", id=Path(f).stem) for f in sorted(glob.glob(f"runs/landscape/learned_f3/{u}_*.json"))]
        extra = {r["id"].split("|", 1)[1]: r for r in R if r.get("part") in ("anchor", "protein")}
        for r in R:
            if r.get("part") == "grid":
                for kk in ("anchor", "protein"):
                    e = next((x for x in R if x.get("part") == kk and x["id"] == f"{kk}|{r['id'].split('|', 1)[1]}"), None)
                    if e: r.update({m: e[m] for m in ("above_null", "prot_true_sig", "prot_pred_sig", "prot_pred_r") if m in e})
    if not R: continue
    mets = [m for m in CFG["metrics"] if any(m in r for r in R)]; G = [r for r in R if r.get("part") in ("grid", "ist", "logreg")]; pm = CFG["primary"]
    out += [f"## {u}  ({len(R)} points, {len(G)} grid)", ""]
    # (a) profile
    rows = [(f"learned {r['model']}:{r.get('config', '')}{' (P3)' if r.get('source') == 'p3' else ''}", r) for r in R if r.get("part") == "learned" and r.get("perturb") in ("none", "shuffled") and "ps" not in r["id"] and "_s1" not in r["id"]]
    rows += [(f"ref {r['ref']}", r) for r in R if r.get("part") == "ref" and r.get("perturb") == "none" and r.get("split_seed", 0) == 0]
    if G:
        best = max(G, key=lambda r: r[pm]); rows.insert(0, (f"floor best-{pm} {CFG['anchor_key'](best)}", best))
    out += ["### (a) profile", "| row | " + " | ".join(mets) + " |", "|" + "---|" * (len(mets) + 1)] + [f"| {n} | " + " | ".join(f3(r.get(m)) for m in mets) + " |" for n, r in rows] + [""]
    if not G: continue
    # (b) shape
    out += ["### (b) shape: mean / max of the primary per knob value (over the rest of the grid)"]
    for kn in CFG["knobs"]:
        by = defaultdict(list)
        for r in G: by[str(r.get(kn))].append(r[pm])
        out.append(f"- {kn}: " + "; ".join(f"{v} {np.mean(x):.3f}/{np.max(x):.3f}" for v, x in sorted(by.items(), key=lambda kv: kv[0])))
    gb = {g: float(np.quantile([r[g] for r in G if g in r and r[g] == r[g]], .75)) for g in CFG["guards"] if any(g in r and r[g] == r[g] for r in G)}  # guard floor: grid 75th pct
    plat = [r for r in G if r[pm] >= best[pm] - EPS and all(r.get(g, np.nan) >= gb[g] for g in gb if r.get(g) == r.get(g))]
    near = [r for r in G if r[pm] >= best[pm] - EPS]
    out += [f"- within {EPS} of best {pm} ({best[pm]:.3f}): {len(near)} / {len(G)} points; also every guard >= its grid 75th pct: {len(plat)}"]
    out += ["  - " + ", ".join(sorted(CFG["anchor_key"](r) for r in plat))[:1500] if plat else "  - none: no point is near-best on the primary and every guard at once", ""]
    # (c) stability
    out += ["### (c) stability over seeds / perturbations (mean, range of the primary; worst guard drop vs unperturbed)"]
    groups = defaultdict(list)
    for r in R:
        if r.get("part") in ("perturb", "split", "k"): groups[(r["part"], CFG["anchor_key"](r) if r.get("part") != "learned" else r["model"])].append(r)
        if r.get("part") == "learned" and r.get("perturb") not in (None, "shuffled") and r.get("source") != "p3": groups[("learned", r["model"])].append(r)
    for (p, key), rs in sorted(groups.items()):
        v = [r[pm] for r in rs if pm in r]
        lab = {"perturb": lambda r: r["perturb"], "split": lambda r: f"split{r['split_seed']}", "k": lambda r: f"dk{r['dk']:+d}", "learned": lambda r: f"{r['perturb']}{r.get('perturb_seed', '')}"}[p]
        out.append(f"- {p} {key}: {pm} mean {np.mean(v):.3f} range [{min(v):.3f}, {max(v):.3f}] | " + ", ".join(f"{lab(r)} {r[pm]:.3f}" for r in rs if pm in r))
    out += [""]
    # (d) metric validity
    out += ["### (d) metric validity over the grid: Spearman of the primary with each panel metric"]
    out.append("- " + "; ".join(f"{m} {spearmanr([r[pm] for r in G if m in r and r[m] == r[m]], [r[m] for r in G if m in r and r[m] == r[m]])[0]:+.2f}"
                                for m in mets if m != pm and sum(m in r and r[m] == r[m] for r in G) > 5))
    q = np.quantile([r[pm] for r in G], .9); flags = []
    for g in gb:
        med = np.median([r[g] for r in G if g in r and r[g] == r[g]])
        flags += [f"{CFG['anchor_key'](r)} ({g} {r[g]:.3f} < median {med:.3f})" for r in G if r[pm] >= q and g in r and r[g] == r[g] and r[g] < med]
    out += [f"- top-10% {pm} points with a guard in its bottom half: {len(flags)}"] + [f"  - {x}" for x in flags[:15]] + [""]
Path(f"runs/landscape/LANDSCAPE_{fam.upper()}.md").write_text("\n".join(out)); print("\n".join(out[:3]), f"-> runs/landscape/LANDSCAPE_{fam.upper()}.md")
