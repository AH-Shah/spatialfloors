"""Collect the positive-control results into one table (exploratory report). Reads runs/p3/posctl_*.json only.
Usage: python scripts/x_posctl.py -> runs/p4/X_posctl.{json,md}"""

import json
from pathlib import Path

P = Path("runs/p3")
j = lambda n: json.loads((P / f"posctl_{n}.json").read_text())
f1, f2, f3, f5, nov, nc, na = j("f1"), j("f2"), j("f3"), j("f5"), j("novae"), j("ncem"), j("ncem_auth_r46")
rows = [
    (
        "F1",
        "scVI",
        "planted types (5 types, 200 genes)",
        "F1 ≥ floor − .05 and ≥ chance + .20",
        f"scVI {f1['scvi']:.3f}, floor {f1['floor']:.3f}, chance {f1['chance']:.3f}",
        f1["pass"],
    ),
    (
        "F1",
        "Nicheformer",
        "vocabulary gate (cannot read synthetic genes)",
        "≥ 80% panel in vocabulary; F1 ≥ chance + .05 on every dev crop",
        "see runs/p3/f1_*",
        None,
    ),
    (
        "F2",
        "CellCharter",
        "2 domains differing only in neighbourhood",
        "ARI ≥ .80 (k = 2, seeds 0–2)",
        f"ARI {f2['mean']:.3f}",
        f2["pass"],
    ),
    (
        "F2",
        "Novae",
        "authors' tutorial slide (Xenium colon P5)",
        "ARI vs published domains ≥ .50; eff. rank ≥ 5; largest domain ≤ 50%",
        f"ARI {nov['ari_vs_published']:.3f}, rank {nov['effective_rank']:.1f}, largest {nov['largest_domain_share']:.2f}",
        nov["pass"],
    ),
    (
        "F3",
        "LightGBM",
        "planted niche toys (quadratic, couple, bleed)",
        "quadratic ≥ ridge + .08; couple ≥ .15; bleed null",
        str(f3["checks"]["lightgbm"]),
        f3["checks"]["lightgbm"]["pass"],
    ),
]
for k, v in na["checks"].items():
    rows.append(
        (
            "F3",
            f"NCEM authors' setup {k}" + (" (locked: linA)" if k == "linA" else ""),
            "planted toys (compo, couple, bleed)",
            "compo / couple recovered; bleed null",
            str(v),
            v["pass"],
        )
    )
for k, v in nc["checks"].items():
    rows.append(("F3", f"NCEM {k} (earlier port, superseded)", "planted niche toys", "as LightGBM", str(v), v["pass"]))
rows += [
    (
        "F5",
        "scANVI",
        "planted-type toy (query depth ×0.3, platform factor)",
        "coarse macro-F1 ≥ .90",
        f"{f5['scanvi']['macro_f1']:.3f}",
        f5["scanvi"]["passed"],
    ),
    ("F5", "Tangram", "same toy", "coarse macro-F1 ≥ .90", f"{f5['tangram']['macro_f1']:.3f}", f5["tangram"]["passed"]),
]
md = [
    "# X_posctl: positive controls of the learned models (run before confirmation; rules fixed before each run)",
    "",
    "| Family | Model | Control | Pass rule | Observed | Pass |",
    "|---|---|---|---|---|---|",
]
md += [
    f"| {a} | {b} | {c} | {d} | {e} | {'yes' if p else ('no' if p is False else 'see note')} |"
    for a, b, c, d, e, p in rows
]
md += [
    "",
    "F4 (Cellina, CPA): the control is the reproduction of the authors' published GSE1 results (runs/p0_cellina; PROTOCOL §3).",
]
Path("runs/p4/X_posctl.json").write_text(
    json.dumps([dict(zip(("family", "model", "control", "rule", "observed", "pass"), r)) for r in rows], indent=1)
)
Path("runs/p4/X_posctl.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
