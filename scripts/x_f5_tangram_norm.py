"""Exploratory: does the registered Tangram result depend on its raw-count input?
The registered run maps the reference to the query on raw integer counts (cluster mode, uniform density prior). Here
both sides are normalised first, everything else as registered (models.transfer.tangram, same genes, same scoring):
  cp10k     counts per 10,000
  logcp10k  log1p of counts per 10,000
Usage (env scvi, cwd = the runs root):
  python scripts/x_f5_tangram_norm.py unit <crop> <out dir>      -> <out dir>/f5tn_<crop>.json (+ labels)
  python scripts/x_f5_tangram_norm.py report <out dir> <p4 dir>  -> <out dir>/X_f5_tangram_norm.{json,md}"""

import json, sys, time
from pathlib import Path
import numpy as np

PUB = Path(__file__).resolve().parents[1]
L5 = json.loads((PUB / "manifests/locked_v3.json").read_text())["locked"]["F5"]
VARIANTS = ("cp10k", "logcp10k")

if sys.argv[1] == "unit":
    import anndata as ad, scipy.sparse as sp, torch
    from spatialfloors.floors import f5
    from spatialfloors.metrics.f5 import score_classes, unit_inputs
    from spatialfloors.models import transfer
    from spatialfloors.provenance import stamp

    u, OUT = sys.argv[2], Path(sys.argv[3])
    (OUT / "labels").mkdir(parents=True, exist_ok=True)
    MK = f5.for_unit(u, json.loads((PUB / "manifests/panels/index.json").read_text())["f5_heldout_markers"])
    held = {x for v in MK.values() for x in v}
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
    pg = set(f5.pelka_genes())
    genes = [
        x for x in json.loads((PUB / f"manifests/panels/{u}.json").read_text())["genes"] if x not in held and x in pg
    ]
    Xq = sp.csr_matrix(a[:, genes].X)
    Xf, vf, gate, cls = unit_inputs(u, a)
    Xr, yr = f5.reference(genes, counts=True)

    def norm(X, log):
        X = sp.csr_matrix(X, dtype=np.float64)
        X = sp.diags(1e4 / np.maximum(np.asarray(X.sum(1)).ravel(), 1)) @ X
        return (X.log1p() if log else X).tocsr()

    rep, LAB = (
        dict(unit=u, n_cells=int(a.n_obs), n_genes=len(genes), density=L5["tangram"]["human"]["config"], grid={}),
        {},
    )
    for v in VARIANTS:
        t = time.time()
        lab = transfer.tangram(norm(Xr, v == "logcp10k"), yr, norm(Xq, v == "logcp10k"), density=rep["density"])
        LAB[v] = np.asarray(lab).astype(str)
        rep["grid"][v] = dict(score_classes(f5.coarse(lab), Xf, vf, MK, gate, cls), s=round(time.time() - t, 1))
        torch.mps.empty_cache()
        print(
            u,
            v,
            round(rep["grid"][v]["heldout_auroc"], 3),
            rep["grid"][v]["pred_share"],
            rep["grid"][v]["s"],
            flush=True,
        )
    np.savez_compressed(OUT / "labels" / f"f5tn_{u}.npz", **LAB)
    (OUT / f"f5tn_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
    print("F5TN_UNIT_DONE", u, flush=True)
else:
    import pandas as pd
    from spatialfloors import stats

    OUT, P4 = Path(sys.argv[2]), Path(sys.argv[3])
    donors = {
        f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"]
        for d in json.loads((PUB / "manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]
    }
    cfg = lambda m: L5[m]["human"]["config"]
    rows = []
    for f in sorted(OUT.glob("f5tn_*.json")):
        d = json.loads(f.read_text())
        u = d["unit"]
        g = lambda m: json.loads((P4 / f"f5_{u}_{m}.json").read_text())["grid"][cfg(m)]
        rows.append(
            dict(
                crop=u,
                floor=g("floors")["heldout_auroc"],
                raw=g("tangram")["heldout_auroc"],
                share_raw=g("tangram")["pred_share"],
                **{v: d["grid"][v]["heldout_auroc"] for v in VARIANTS},
                **{f"share_{v}": d["grid"][v]["pred_share"] for v in VARIANTS},
            )
        )
    df = pd.DataFrame(rows)
    unit_of = {c: donors.get(c, c) for c in df.crop}
    out = dict(n_crops=len(df), pairs={})
    md = [
        "# F5 Tangram on normalised input (exploratory)",
        "",
        "D = Tangram − floor per unit (held-out-marker AUROC); registered test (platform strata, margin of F5).",
        "",
        "| Tangram input | Units | Mean AUROC | D | One-sided 95% bounds | Non-inferior | Predicted share, epi / imm / mye / other |",
        "|---|---|---|---|---|---|---|",
    ]
    for key, name in (
        ("raw", "raw counts (registered)"),
        ("cp10k", "counts per 10,000"),
        ("logcp10k", "log1p counts per 10,000"),
    ):
        units, D, strata = stats.aggregate_units(dict(zip(df.crop, df[key] - df.floor)), unit_of)
        t = stats.family_test(D, strata, margin=L5["margin"])
        sh = pd.DataFrame(list(df[f"share_{key}"])).mean()
        out["pairs"][key] = {k: t[k] for k in ("est", "lower", "upper", "noninferior")} | dict(
            n_units=len(D),
            mean_auroc=float(stats.aggregate_units(dict(zip(df.crop, df[key])), unit_of)[1].mean()),
            pred_share=sh.to_dict(),
            per_unit=dict(zip(units, map(float, D))),
        )
        md.append(
            f"| {name} | {len(D)} | {out['pairs'][key]['mean_auroc']:.3f} | {t['est']:+.4f} | {t['lower']:+.4f}, {t['upper']:+.4f} | "
            f"{'yes' if t['noninferior'] else 'no'} | "
            + " / ".join(f"{sh[k]:.2f}" for k in ("epithelial", "immune_nonmyeloid", "myeloid", "other"))
            + " |"
        )
    (OUT / "X_f5_tangram_norm.json").write_text(json.dumps(out, indent=1, default=float))
    (OUT / "X_f5_tangram_norm.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
