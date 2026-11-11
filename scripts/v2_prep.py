"""Replication round, metric-blind preparation of one section: crop, cache and gene panel by the round-1 rules.
Crop: data.crops.choose_crop on coordinates only (densest square of 50,000 cells without region labels; sections of up
to 55,000 cells whole). Cache: data.units_v2.load on the crop's cells, written as <root>/runs/p2/cache/<section>.h5ad
(round-1 layout). Panel: data.genes.panel on X_a of the count split seed 0 (Xenium 5K and CosMx 960-plex are targeted
panels: every filtered gene), as scripts/p1_panels.py. No model is fitted and no score is computed.
Sections that may enter F2 (the five 10x sections with H&E) wait for the pathologists' regions, which set their crop.
Usage (env floors): python scripts/v2_prep.py <root> <section ...>
  -> <root>/manifests/crops/<section>.csv.gz + index.json, <root>/manifests/panels/<section>.json + index.json,
     <root>/runs/p2/cache/<section>.h5ad"""

import json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from spatialfloors.data import units_v2
from spatialfloors.data.crops import choose_crop
from spatialfloors.data.genes import HGNC, HGNC_SHA, panel
from spatialfloors.metrics.f1 import count_split
from spatialfloors.provenance import stamp

F5_HELDOUT = {  # scripts/p1_panels.py
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "immune_nonmyeloid": ["CD3E", "CD2", "MS4A1", "CD79A", "MZB1"],
    "myeloid": ["CD68", "CD163", "C1QC", "CD14", "LYZ"],
}
ROOT = Path(sys.argv[1])
CR, PN, CA = ROOT / "manifests/crops", ROOT / "manifests/panels", ROOT / "runs/p2/cache"
for d in (CR, PN, CA):
    d.mkdir(parents=True, exist_ok=True)
J = lambda p: json.loads(p.read_text()) if p.exists() else {}
cidx, pidx = J(CR / "index.json"), J(PN / "index.json")
for s in sys.argv[2:]:
    t = time.time()
    o = units_v2.obs(s)
    sel, rec = choose_crop(o[["x_um", "y_um"]].values)
    ids = o.index[sel]
    pd.Series(np.asarray(ids), name="cell_id").to_csv(CR / f"{s}.csv.gz", index=False)
    a = units_v2.load(s, ids)
    a.obs["region"] = a.obs.region.astype("category")
    a.write_h5ad(CA / f"{s}.h5ad", compression="gzip")
    Xa, _ = count_split(a.X, seed=0)
    p = panel(Xa, a.var_names.values, "human")
    p.pop("idx")
    gs = set(p["genes"])
    p.update(
        unit=s,
        species="human",
        n_cells=int(a.n_obs),
        n_panel=len(p["genes"]),
        f5_heldout_in_panel={c: sorted(gs & set(m)) for c, m in F5_HELDOUT.items()},
        f5_heldout_on_platform={c: sorted(set(a.var_names) & set(m)) for c, m in F5_HELDOUT.items()},
        annotation=dict(hgnc=str(HGNC), sha256=HGNC_SHA),
    )
    (PN / f"{s}.json").write_text(json.dumps(stamp(p), indent=1))
    cidx[s] = stamp(rec | a.uns["unit"] | dict(n_cells=int(a.n_obs)))
    pidx[s] = stamp({k: p[k] for k in ("mode", "n_input", "n_panel", "removed", "n_cells")})
    print(s, rec["n_slide"], "->", rec["n_crop"], "whole" if rec["whole_slide"] else f"half-side {rec['half_side']:.0f} um",
          "panel", p["n_panel"], f"{time.time() - t:.0f}s", flush=True)  # fmt: skip
pidx["f5_heldout_markers"] = F5_HELDOUT
(CR / "index.json").write_text(json.dumps(cidx, indent=1, default=float))
(PN / "index.json").write_text(json.dumps(pidx, indent=1, default=float))
print("V2_PREP_DONE", flush=True)
