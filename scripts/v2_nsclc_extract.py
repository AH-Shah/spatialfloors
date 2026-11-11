"""Replication round: take the three flat files of each CosMx NSCLC section out of its tarball (expression matrix,
cell metadata, FOV positions). Images and the transcript table stay in the archive. Members are written by base name
into <raw>/cosmx_nsclc/<section>/, so no path stored in the archive is followed.
Usage: python -I scripts/v2_nsclc_extract.py <raw dir>     (prints V2_NSCLC_EXTRACT_DONE)"""

import re, sys, tarfile
from pathlib import Path

RAW = Path(sys.argv[1]) / "cosmx_nsclc"
KEEP = re.compile(r"^[A-Za-z0-9_]+_(exprMat|metadata|fov_positions)_file\.csv$")
for tgz in sorted(RAW.glob("*.tar.gz")):
    out = RAW / tgz.name.removesuffix(".tar.gz")
    out.mkdir(exist_ok=True)
    with tarfile.open(tgz, "r|gz") as t:
        for m in t:
            name = m.name.rsplit("/", 1)[-1]
            if not (m.isfile() and KEEP.match(name)):
                continue
            dst = out / name
            if not (dst.exists() and dst.stat().st_size == m.size):
                with t.extractfile(m) as src, open(dst, "wb") as f:
                    while chunk := src.read(8 << 20):
                        f.write(chunk)
            print("EXTRACTED", dst, m.size, flush=True)
print("V2_NSCLC_EXTRACT_DONE", flush=True)
