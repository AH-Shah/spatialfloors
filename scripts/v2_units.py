"""Replication round: the unit table (registration_v2/units/units.tsv), one row per section, from data.units_v2 and
the crop manifests of scripts/v2_prep.py. Families: F1 and F3 every unit; F5 carcinomas (SELECTION.md); F2 only after
the pathologists' regions (status "awaiting regions" for the five 10x sections with H&E, whose crop the regions set).
Unit = patient: sections of one patient share the unit and are averaged (stats.aggregate_units). Stratum = platform.
Usage (env floors): python scripts/v2_units.py <root>   (root = the round-2 runs root)"""

import hashlib, json, sys
from pathlib import Path
import pandas as pd
from spatialfloors import stats
from spatialfloors.data.units_v2 import PLATFORM, SECTIONS

PUB = Path(__file__).resolve().parents[1]
ROOT = Path(sys.argv[1])
idx = json.loads((ROOT / "manifests/crops/index.json").read_text())
rows = []
for s, (kind, _, unit, source, tissue, cancer) in SECTIONS.items():
    m = ROOT / f"manifests/crops/{s}.csv.gz"
    fam = ["F1", "F3"] + (["F5"] if cancer == "carcinoma" else [])
    rows.append(
        dict(
            crop=s,
            unit=unit,
            role="confirmation_v2",
            families=",".join(fam),
            stratum=stats.stratum(s),
            platform=PLATFORM[kind],
            source=source,
            tissue=tissue,
            cancer=cancer,
            n_cells=idx[s]["n_crop"] if s in idx else "",
            status="prepared" if m.exists() else "awaiting regions (F2 candidate; crop set by the regions)",
            crop_manifest_sha256=hashlib.sha256(m.read_bytes()).hexdigest() if m.exists() else "",
        )
    )
U = pd.DataFrame(rows)
out = PUB / "registration_v2/units"
out.mkdir(parents=True, exist_ok=True)
U.to_csv(out / "units.tsv", sep="\t", index=False)
print(U.drop(columns="crop_manifest_sha256").to_string(index=False))
print("units", U.unit.nunique(), "sections", len(U), "F5 units", U[U.families.str.contains("F5")].unit.nunique())
