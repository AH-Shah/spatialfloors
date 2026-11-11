"""Replication round (registration v2): fetch the count matrices and cell tables of the selected public datasets.
Zip archives on Zenodo are read through HTTP range requests, so only the members that match a pattern are transferred
(the count matrix, the cell table, the gene panel and any annotation file; never images or transcript tables).
Usage:
  python scripts/v2_fetch.py list <zip url>                        -> member names and sizes
  python scripts/v2_fetch.py get <zip url> <out dir> <regex>       -> extract matching members into <out dir>
  python scripts/v2_fetch.py file <url> <out path>                 -> plain download (resumes a partial file)
  python scripts/v2_fetch.py size <url> <path>                     -> local size against the server's Content-Length
  python scripts/v2_fetch.py verify <zip url> <out dir> <regex>    -> local size and CRC-32 against the zip directory
The two checks print one OK or BAD line per file and exit 1 on any BAD."""

import io, re, sys, urllib.request, zipfile, zlib
from pathlib import Path

UA = {"User-Agent": "python-urllib"}


class RangeFile(io.RawIOBase):
    """Read-only seekable file over HTTP range requests, with one read-ahead block."""

    def __init__(self, url, block=4 << 20):
        self.url, self.block, self.pos, self.buf, self.buf_at = url, block, 0, b"", -1
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD", headers=UA), timeout=120) as r:
            self.size = int(r.headers["Content-Length"])
            self.url = r.url  # after redirects

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def _fetch(self, a, n):
        req = urllib.request.Request(self.url, headers=UA | {"Range": f"bytes={a}-{min(a + n, self.size) - 1}"})
        for attempt in range(5):
            try:
                with urllib.request.urlopen(req, timeout=300) as r:
                    assert r.status == 206, f"no range support: HTTP {r.status}"
                    return r.read()
            except Exception as e:  # transient network errors: retry
                if attempt == 4:
                    raise
                print("retry", attempt + 1, e, flush=True)

    def read(self, n=-1):
        n = self.size - self.pos if n < 0 else min(n, self.size - self.pos)
        if n <= 0:
            return b""
        if not (self.buf_at <= self.pos and self.pos + n <= self.buf_at + len(self.buf)):
            self.buf_at, self.buf = self.pos, self._fetch(self.pos, max(n, self.block))
        out = self.buf[self.pos - self.buf_at : self.pos - self.buf_at + n]
        self.pos += len(out)
        return out

    def readinto(self, b):
        d = self.read(len(b))
        b[: len(d)] = d
        return len(d)


mode = sys.argv[1]
if mode == "file":
    url, out = sys.argv[2], Path(sys.argv[3])
    out.parent.mkdir(parents=True, exist_ok=True)
    have = out.stat().st_size if out.exists() else 0
    req = urllib.request.Request(url, headers=UA | ({"Range": f"bytes={have}-"} if have else {}))
    try:
        with urllib.request.urlopen(req, timeout=300) as r, open(out, "ab" if have and r.status == 206 else "wb") as f:
            while chunk := r.read(8 << 20):
                f.write(chunk)
    except urllib.error.HTTPError as e:
        if e.code != 416:  # 416 = already complete
            raise
    print("FETCHED", out, out.stat().st_size, flush=True)
elif mode == "size":
    with urllib.request.urlopen(urllib.request.Request(sys.argv[2], method="HEAD", headers=UA), timeout=120) as r:
        want = int(r.headers["Content-Length"])
    have = Path(sys.argv[3]).stat().st_size
    print("OK" if have == want else "BAD", sys.argv[3], have, want, flush=True)
    sys.exit(have != want)
else:
    z = zipfile.ZipFile(io.BufferedReader(RangeFile(sys.argv[2]), buffer_size=1 << 20))
    if mode == "list":
        for i in z.infolist():
            print(f"{i.file_size / 1e6:10.1f} MB  {i.filename}")
    elif mode == "verify":
        out, pat, bad = Path(sys.argv[3]), re.compile(sys.argv[4]), 0
        for i in z.infolist():
            if pat.search(i.filename) and not i.is_dir():
                crc, size = 0, 0
                with open(out / i.filename, "rb") as f:
                    while chunk := f.read(8 << 20):
                        crc, size = zlib.crc32(chunk, crc), size + len(chunk)
                ok = size == i.file_size and crc == i.CRC
                bad += not ok
                print(
                    "OK" if ok else "BAD", out / i.filename, size, i.file_size, f"{crc:08x}", f"{i.CRC:08x}", flush=True
                )
        sys.exit(bad > 0)
    else:
        out, pat = Path(sys.argv[3]), re.compile(sys.argv[4])
        for i in z.infolist():
            if pat.search(i.filename) and not i.is_dir():
                dst = out / i.filename
                dst.parent.mkdir(parents=True, exist_ok=True)
                if not (dst.exists() and dst.stat().st_size == i.file_size):
                    with z.open(i) as src, open(dst, "wb") as f:
                        while chunk := src.read(8 << 20):
                            f.write(chunk)
                print("FETCHED", dst, i.file_size, flush=True)
