"""spatialfloors CLI. `check`: re-run the vendored tools/87 ladder on crc_232 and assert bit-identity with row 323."""
import argparse, json, sys, tempfile, time
from pathlib import Path
import numpy as np, pandas as pd

from spatialfloors.paths import DERIVED
ROW323 = DERIVED / "perturb_opt"


def check(args):
    from .floors import l6
    from .provenance import stamp
    out = Path(args.out or tempfile.mkdtemp(prefix="sf_check_")); out.mkdir(parents=True, exist_ok=True); t0 = time.time()
    l6.run_slide("crc_232", l6.load_crc("232", args.raw), 100 / 0.12028, ["REF"], ["CRC"], l6.C_TYPES, out)
    ref, new = pd.read_csv(ROW323 / "crc_232.csv"), pd.read_csv(out / "crc_232.csv")
    csv_ok = ref.shape == new.shape and bool((ref.columns == new.columns).all()) and ref.equals(new)
    zr, zn = np.load(ROW323 / "crc_232_lfc.npz"), np.load(out / "crc_232_lfc.npz")
    keys_ok = sorted(zr.files) == sorted(zn.files); bad = [k for k in zr.files if k not in zn.files or not np.array_equal(zr[k], zn[k], equal_nan=True)]
    num = new.select_dtypes("number").to_numpy() - ref.select_dtypes("number").to_numpy() if ref.shape == new.shape else None
    rep = stamp(dict(csv_bit_identical=csv_ok, npz_keys_match=keys_ok, npz_mismatched_arrays=len(bad), npz_arrays=len(zr.files),
                     csv_max_abs_diff=float(np.nanmax(np.abs(num))) if num is not None else None, seconds=round(time.time() - t0, 1), out=str(out)))
    print(json.dumps(rep, indent=1)); (out / "CHECK.json").write_text(json.dumps(rep, indent=1))
    return 0 if csv_ok and keys_ok and not bad else 1


def main(argv=None):
    p = argparse.ArgumentParser("spatialfloors"); sp_ = p.add_subparsers(dest="cmd", required=True)
    c = sp_.add_parser("check"); c.add_argument("--raw", default=None); c.add_argument("--out", default=None); c.set_defaults(fn=check)
    a = p.parse_args(argv); return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
