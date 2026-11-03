"""F4 v3 diagnostic: the similarity weights floors.l6.type_similarity gives each held-out type on a GSE1 slide (why L9
equals L6 on crc_232). Usage (env floors): python scripts/v3_f4_sim.py <sid>"""
import sys
import numpy as np
from spatialfloors.floors import l6

C, xy, ct, dom = l6.load_crc(sys.argv[1]); lib = np.asarray(C.sum(1)).ravel(); gs = l6.group_stats(C, ct, dom, lib)
for h in l6.C_TYPES:
    w = l6.type_similarity(gs, h, ["CRC"], ["REF"]); print(h, {f"{k[0]}|{k[1]}": round(float(v), 3) for k, v in zip(gs[0], w) if k[1] == "REF"})
