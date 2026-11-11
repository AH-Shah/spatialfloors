"""Replication round: file inventory of every fetched section (no expression value is summarised, no model, no score).
Per section: cells, features by type (genes and control features), table columns, coordinate unit and extent.
Xenium Prime 5K (10x output bundle): cell_feature_matrix.h5 (features x cells, feature_type per feature),
cells.parquet (centroids in µm), gene_panel.json (panel identity).
CosMx NSCLC (flat files, scripts/v2_nsclc_extract.py): exprMat (one row per fov and cell_ID; cell_ID 0 = transcripts
outside cells; NegPrb columns = negative probes), metadata (centroids in pixels; 0.18 µm per pixel, dataset README).
Usage (env floors): python -I scripts/v2_inventory.py <raw dir> registration_v2
  -> registration_v2/inventory.json, registration_v2/INVENTORY.md"""

import collections, json, sys
from pathlib import Path
import h5py, pandas as pd, pyarrow.parquet as pq

RAW, OUT = Path(sys.argv[1]), Path(sys.argv[2])
NSCLC_UM_PER_PX = 0.18  # SMI Data File ReadMe v1 (NanoString, 2021-10-04): "pixel edge length is 180nm"
XENIUM = [("cervilla2026_xenium5k", d.name) for d in sorted((RAW / "cervilla2026_xenium5k").glob("xenium5K_*"))] + [
    ("tenx_xenium5k", d.name) for d in sorted((RAW / "tenx_xenium5k").glob("Xenium_*"))
]
inv = {}
for src, name in XENIUM:
    d = RAW / src / name
    with h5py.File(d / "cell_feature_matrix.h5", "r") as f:
        m = f["matrix"]
        n_feat, n_cells = (int(x) for x in m["shape"][:])
        types = collections.Counter(x.decode() for x in m["features/feature_type"][:])
        assert sum(types.values()) == n_feat and len(m["barcodes"]) == n_cells
    p = pq.ParquetFile(d / "cells.parquet")
    c = p.read(columns=["x_centroid", "y_centroid", "segmentation_method"]).to_pandas()
    assert len(c) == n_cells, (name, len(c), n_cells)
    panel = json.loads((d / "gene_panel.json").read_text())["payload"]["panel"]
    base = panel.get("type", {}).get("data", {}).get("base_panel_identity", {})
    inv[f"{src}/{name}"] = dict(
        platform="Xenium Prime 5K",
        cells=n_cells,
        features=dict(types),
        genes=types["Gene Expression"],
        control_features=n_feat - types["Gene Expression"],
        panel=panel["identity"]["name"],
        panel_design=panel["identity"].get("design_id"),
        base_panel=base.get("name"),
        cell_columns=[x.name for x in p.schema_arrow],
        coordinate_unit="µm (x_centroid, y_centroid)",
        extent_mm=[round(float(c[k].max() - c[k].min()) / 1000, 2) for k in ("x_centroid", "y_centroid")],
        segmentation=dict(c.segmentation_method.value_counts()),
    )
for d in sorted((RAW / "cosmx_nsclc").iterdir()):
    if not d.is_dir():
        continue
    s = d.name
    cols = pd.read_csv(d / f"{s}_exprMat_file.csv", nrows=0).columns.tolist()
    ids = pd.read_csv(d / f"{s}_exprMat_file.csv", usecols=["fov", "cell_ID"])
    meta = pd.read_csv(d / f"{s}_metadata_file.csv")
    fov = pd.read_csv(d / f"{s}_fov_positions_file.csv")
    feat = [x for x in cols if x not in ("fov", "cell_ID")]
    neg = [x for x in feat if x.startswith("NegPrb")]
    cells = ids[ids.cell_ID != 0]
    assert not cells.duplicated().any() and len(cells) == len(meta) == len(meta[["fov", "cell_ID"]].drop_duplicates())
    inv[f"cosmx_nsclc/{s}"] = dict(
        platform="CosMx SMI 960-plex (NSCLC FFPE release)",
        cells=int(len(meta)),
        expr_rows=int(len(ids)),
        expr_rows_outside_cells=int((ids.cell_ID == 0).sum()),
        features={"gene": len(feat) - len(neg), "NegPrb": len(neg)},
        genes=len(feat) - len(neg),
        control_features=len(neg),
        fovs=int(meta.fov.nunique()),
        fov_positions_rows=int(len(fov)),
        fov_position_columns=fov.columns.tolist(),
        cell_columns=meta.columns.tolist(),
        coordinate_unit=f"pixels (CenterX_global_px, CenterY_global_px); {NSCLC_UM_PER_PX} µm per pixel (README)",
        extent_mm=[
            round(float(meta[k].max() - meta[k].min()) * NSCLC_UM_PER_PX / 1000, 2)
            for k in ("CenterX_global_px", "CenterY_global_px")
        ],
    )
(OUT / "inventory.json").write_text(json.dumps(inv, indent=1, default=int))
md = [
    "# Replication round: file inventory",
    "",
    "Counts of cells, features and columns as fetched; no expression value was summarised and no model was run "
    "(`scripts/v2_inventory.py`). Checksums: `data_sha256.tsv`.",
    "",
    "| Section | Platform | Cells | Genes | Control features | Extent (mm) | Coordinates |",
    "|---|---|---|---|---|---|---|",
]
for k, v in inv.items():
    md.append(
        f"| {k} | {v['platform']} | {v['cells']:,} | {v['genes']:,} | {v['control_features']:,} | "
        f"{v['extent_mm'][0]} × {v['extent_mm'][1]} | {v['coordinate_unit']} |"
    )
md += ["", "Feature types, panels and columns:", ""]
for k, v in inv.items():
    extra = f"panel {v['panel']}" + (f" (base {v['base_panel']})" if v.get("base_panel") else "") if "panel" in v else (
        f"{v['fovs']} FOVs; {v['expr_rows_outside_cells']} expression rows with cell_ID 0"
    )  # fmt: skip
    md.append(f"- {k}: " + ", ".join(f"{a} {b}" for a, b in v["features"].items()) + f"; {extra}; cell table columns: " + ", ".join(v["cell_columns"]))  # fmt: skip
(OUT / "INVENTORY.md").write_text("\n".join(md) + "\n")
print("\n".join(md[:40]))
print("V2_INVENTORY_DONE", len(inv))
