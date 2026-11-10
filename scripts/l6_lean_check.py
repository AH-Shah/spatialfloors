"""Structure-only identity check of l6_lean.load_crc_lean against l6.load_crc on GSE1 slides (raw data only; nothing
fitted
or scored). Usage (env floors, PYTHONPATH=src): python scripts/l6_lean_check.py <sid ...>"""

import resource, sys, time
import numpy as np
from spatialfloors.floors import l6, l6_lean

for sid in sys.argv[1:]:
    t = time.time()
    L = l6_lean.load_crc_lean(sid)
    tl = time.time() - t
    pk = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30
    O = l6.load_crc(sid)
    A, B = O[0], L[0]
    same = (
        A.shape == B.shape
        and np.array_equal(A.indptr, B.indptr)
        and np.array_equal(A.indices, B.indices)
        and np.array_equal(A.data, B.data)
        and np.array_equal(O[1], L[1])
        and np.array_equal(O[2], L[2])
        and np.array_equal(O[3], L[3])
    )
    print(
        f"{sid} shape {A.shape} lean {B.shape} identical {same} lean {tl:.0f}s peak_rss_after_lean {pk:.1f}G",
        flush=True,
    )
