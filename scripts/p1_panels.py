"""Gene panel per unit (PROTOCOL §3) -> manifests/panels/<unit>.json + index.json.
X_a = count split seed 0 (metrics.f1.count_split) of the crop; filter (data.genes); WTx panels: seurat_v3 HVG 2000 on
X_a;
targeted panels: all filtered genes. F4 is not covered (Cellina / CPA / L6 use the authors' own HVG step, PROTOCOL §2).
F5 held-out markers (secondary AUROC, WTx) are removed from F5 inputs; their overlap with each panel is recorded.
Usage (env floors): python scripts/p1_panels.py [unit ...]"""

import json, sys, time
from pathlib import Path
from spatialfloors.data.genes import panel, HGNC, HGNC_SHA
from spatialfloors.data.units import index, load_crop, YAO
from spatialfloors.metrics.f1 import count_split
from spatialfloors.provenance import stamp

OUT = Path("manifests/panels")
OUT.mkdir(parents=True, exist_ok=True)
F5_HELDOUT = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "immune_nonmyeloid": ["CD3E", "CD2", "MS4A1", "CD79A", "MZB1"],
    "myeloid": ["CD68", "CD163", "C1QC", "CD14", "LYZ"],
}
units = sys.argv[1:] or [u for u in index() if not (u.startswith("allen_yao") and not (YAO / "RAW_DONE").exists())]
summ = json.loads((OUT / "index.json").read_text()) if (OUT / "index.json").exists() else {}
for u in units:
    t = time.time()
    a = load_crop(u)
    species = "mouse" if u.startswith("allen") else "human"
    Xa, _ = count_split(a.X, seed=0)
    p = panel(Xa, a.var_names.values, species)
    p.pop("idx")
    gs = set(p["genes"])
    p["f5_heldout_in_panel"] = {c: sorted(gs & set(m)) for c, m in F5_HELDOUT.items()} if species == "human" else None
    p["f5_heldout_on_platform"] = (
        {c: sorted(set(a.var_names) & set(m)) for c, m in F5_HELDOUT.items()} if species == "human" else None
    )
    p.update(
        unit=u,
        species=species,
        n_cells=int(a.n_obs),
        n_panel=len(p["genes"]),
        annotation=dict(hgnc=str(HGNC), sha256=HGNC_SHA) if species == "human" else "symbol rules",
        counts_per_cell_Xa=round(float(Xa.sum() / a.n_obs), 1),
    )
    (OUT / f"{u}.json").write_text(json.dumps(stamp(p), indent=1))
    summ[u] = stamp(
        {k: p[k] for k in ("mode", "n_input", "n_panel", "removed", "n_cells", "counts_per_cell_Xa")}
    )  # per-unit stamp
    print(u, summ[u], f"{time.time() - t:.0f}s", flush=True)
summ["f5_heldout_markers"] = F5_HELDOUT
(OUT / "index.json").write_text(
    json.dumps(summ if "spatialfloors_commit" in summ else stamp(summ), indent=1)
)  # first stamp kept
