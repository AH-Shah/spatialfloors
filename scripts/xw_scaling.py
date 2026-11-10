"""Whole-slide arm, feasibility review (no fitting): peak memory and wall time of every learned model and floor
against the number of cells, from the scheduler's per-job footprint logs (runs/sched/<id>.log.mem: 'HH:MM:SS GB'
every 5 s; mem_guard 'peak_footprint' line), then a linear fit (time and memory ~ a + b * cells) and extrapolation to
the six GSE1 whole slides. Durations include contention with parallel jobs (an upper bound on the dedicated time).
Usage: python3 scripts/xw_scaling.py -> runs/xw/XW_SCALING.md"""

import csv, glob, json, re
from pathlib import Path
import numpy as np

idx = {k: v for k, v in json.loads(Path("manifests/crops/index.json").read_text()).items() if isinstance(v, dict)}
ncell = {k: v.get("n_crop") for k, v in idx.items()}
for r in csv.DictReader(open("registration/units/units.tsv"), delimiter="\t"):
    ncell[r["crop"]] = int(r["n_cells"])
ncell.update(
    {
        "gse1_231w": 298151,
        "gse1_231wl": 298151,
        "gse1_120w": 741155,
        "gse1_110w": 694553,
        "gse1_210w": 617498,
        "gse1_242w": 474471,
        "gse1_221w": 370074,
        "gse1_222w": 140542,
    }
)
MODEL = {
    "1s": "scVI (F1)",
    "1n": "Nicheformer (F1)",
    "1f": "F1 floor",
    "2c": "CellCharter (F2)",
    "2n": "Novae (F2)",
    "2f": "F2 floor",
    "3l": "LightGBM (F3)",
    "3n": "NCEM (F3)",
    "3r": "F3 floor v3",
    "3g": "GAT (F3)",
    "5s": "scANVI (F5)",
    "5t": "Tangram (F5)",
    "5f": "F5 floors",
}
W231 = {
    "w1s_231": "1s",
    "w2c_231": "2c",
    "w5s_231": "5s",
    "w3l_231": "3l",
    "w3r_231": "3r",
    "w1f_231": "1f",
    "w2f_231": "2f",
    "w5f_231": "5f",
    "w3g_231": "3g",
}


def parse(f):
    t, m = [], []
    for line in open(f):
        p = line.split()
        if len(p) >= 2 and re.match(r"\d\d:\d\d:\d\d$", p[0]):
            try:
                m.append(float(p[1]))
            except ValueError:
                continue
            h, mi, s = map(int, p[0].split(":"))
            t.append(h * 3600 + mi * 60 + s)
    if len(t) < 2:
        return None
    t = np.array(t)
    t = t + 86400 * np.cumsum(np.r_[0, np.diff(t) < -3600])
    d = np.diff(t)  # past midnight
    return d[(d > 0) & (d <= 30)].sum() / 60, max(
        m
    )  # active minutes: 5 s samples only (gaps = machine asleep or job paused)


rows = []
for f in glob.glob("runs/sched/*.log.mem"):
    j = Path(f).name[:-8]
    m = re.match(r"^(p4|t2|xw)_(\w\w?)_(.+)$", j)
    key, unit = None, None
    if j in W231:
        key, unit = W231[j], "gse1_231w"
    elif m and m[2] in MODEL:
        key = m[2]
        rest = m[3] if m[1] != "xw" else f"gse1_{m[3]}w"
        rest = re.sub(r"^g3_", "gse303070_", rest)
        rest = re.sub(r"^a_(animal\d)", lambda x: next(k for k in ncell if k.startswith(f"allen_{x[1]}_")), rest)
        rest = re.sub(r"^a_yao638850", "allen_yao638850_C57BL6J-638850.36", rest)
        cand = [k for k in ncell if rest == k or rest.startswith(k + "_")]
        unit = max(cand, key=len) if cand else None
    elif (m2 := re.match(r"^r2_(\w\w?)\d*_(gse1_231w|.+)$", j)) and m2[1] in MODEL:
        key, unit = m2[1], m2[2]
    if key is None or unit not in ncell or ncell[unit] is None:
        continue
    if not (Path(f"runs/sched/{j}.done").exists()):
        continue  # finished runs only
    p = parse(f)
    if p and p[0] > 0.2:
        rows.append(dict(job=j, model=MODEL[key], n=ncell[unit], minutes=p[0], peak_gb=p[1]))
out = [
    "# Feasibility review: memory and time against cells (finished runs only)",
    "",
    "| Model | Runs | Cells (range) | Peak GB (range) | GB per 100k cells (fit) | Minutes (range) | Min per 100k (fit) | Predicted at 741k: GB / min | Six slides (3.04M cells): h |",
    "|---|---|---|---|---|---|---|---|---|",
]
pred = {}
for mdl in MODEL.values():
    R = [r for r in rows if r["model"] == mdl]
    if len(R) < 3:
        continue
    n = np.array([r["n"] for r in R], float)
    g = np.array([r["peak_gb"] for r in R])
    t = np.array([r["minutes"] for r in R])
    fit = lambda y: np.polyfit(n / 1e5, y, 1) if np.ptp(n) > 2e4 else (0.0, y.mean())
    bg, ag = fit(g)
    bt, at = fit(t)
    g741, t741 = ag + bg * 7.41, at + bt * 7.41
    tall = sum(at + bt * x / 1e5 for x in (741155, 694553, 617498, 474471, 370074, 140542)) / 60
    pred[mdl] = dict(
        n_runs=len(R),
        cells=[int(n.min()), int(n.max())],
        gb_per_100k=float(bg),
        min_per_100k=float(bt),
        gb_741k=float(g741),
        min_741k=float(t741),
        hours_six=float(tall),
    )
    out.append(
        f"| {mdl} | {len(R)} | {int(n.min()):,}–{int(n.max()):,} | {g.min():.1f}–{g.max():.1f} | {bg:+.2f} | {t.min():.0f}–{t.max():.0f} | {bt:+.1f} | {g741:.1f} / {t741:.0f} | {tall:.1f} |"
    )
Path("runs/xw").mkdir(exist_ok=True)
Path("runs/xw/XW_SCALING.md").write_text("\n".join(out) + "\n")
Path("runs/xw/XW_SCALING.json").write_text(json.dumps(dict(rows=rows, fits=pred), indent=1))
print("\n".join(out))
big = sorted([r for r in rows if r["n"] > 100000], key=lambda r: (r["model"], r["n"]))
print("\nruns above 100k cells:")
[
    print(f"  {r['job']:28s} {r['model']:18s} {r['n']:>8,} cells {r['peak_gb']:5.2f} GB {r['minutes']:6.1f} min")
    for r in big
]
