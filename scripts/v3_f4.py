"""F4 floor v3 candidates (PLAN §19, floors.l6.run_twoway_v3) on one GSE1 slide: control REF, held-out CRC (SF_HOLDOUT=TVA
for 222), the 5 Cellina types; the F4 crop (<run>/raw/cosmx_wt_<sid>.h5ad) when it exists, else the whole slide, as
scripts/f4_l6.py. Development: crc_232 (choice); validation: 231 crop / whole slide (run once with the chosen floor).
Usage (env floors): [SF_F4_RUN=runs/f4] [SF_F4_WHOLE=1] python scripts/v3_f4.py <sid> -> runs/v3/f4/crc_<sid>_{whole,crop}_v3.csv"""
import os, sys
from pathlib import Path
from spatialfloors.floors import l6

sid = sys.argv[1]; R = Path(os.environ.get("SF_F4_RUN", "runs/f4")); out = Path("runs/v3/f4"); out.mkdir(parents=True, exist_ok=True)
whole = os.environ.get("SF_F4_WHOLE") == "1" or not (R / f"raw/cosmx_wt_{sid}.h5ad").exists()
data = l6.load_crc(sid, None if whole else str(R / "raw"))
l6.run_twoway_v3(f"crc_{sid}" + ("_whole" if whole else "_crop"), data, 1 / 0.12028, ["REF"], [os.environ.get("SF_HOLDOUT", "CRC")], l6.C_TYPES, out)
print("V3F4_DONE", sid, flush=True)
