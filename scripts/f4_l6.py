"""F4 floor L6 (floors.l6.run_slide, tools/87 verbatim). GSE1 slide <sid>: control REF, held-out CRC (SF_HOLDOUT=TVA for
222), the 5 Cellina types, niche bandwidth 100 um (as `spatialfloors check` on crc_232); reads the F4 crop
<run>/raw/cosmx_wt_<sid>.h5ad when it exists (f4_crop.py), else the whole slide. MERFISH section <label>: <run>/raw/<label>.h5ad
(f4_merfish_prep.py), control Thalamus, held-out Isocortex + Fiber_tracts, l6.M_TYPES, 100 um (tools/87 animal 2 call).
Usage (env floors): SF_F4_RUN=<run> python scripts/f4_l6.py <sid|section> -> <run>/l6/<crc_sid|section>.csv + _lfc.npz"""
import os, sys
from pathlib import Path
from spatialfloors.floors import l6, l6_lean

sid = sys.argv[1]; R = Path(os.environ.get("SF_F4_RUN", "runs/f4")); out = R / "l6"; out.mkdir(parents=True, exist_ok=True)
if sid.isdigit():
    crop = (R / f"raw/cosmx_wt_{sid}.h5ad").exists()  # whole slides: memory-lean loader, identical output (scripts/l6_lean_check.py; EXECUTION_LOG 2026-11-06)
    l6.run_slide(f"crc_{sid}", l6.load_crc(sid, str(R / "raw")) if crop else l6_lean.load_crc_lean(sid), 100 / 0.12028,
                 ["REF"], [os.environ.get("SF_HOLDOUT", "CRC")], l6.C_TYPES, out)
else:
    l6.M_SECTIONS = [sid]  # load_merfish iterates this module list (animal 2 sections in tools/87)
    for s, data in l6.load_merfish(R / f"raw/{sid}.h5ad"): l6.run_slide(s, data, 100.0, ["Thalamus"], ["Isocortex", "Fiber_tracts"], l6.M_TYPES, out)
