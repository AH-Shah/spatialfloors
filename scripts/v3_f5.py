"""F5 floor v3 on one unit (PLAN §19; development units for the choice, validation units once with the chosen floor).
Inputs = the v2 floor's own InSituType inputs (runs/landscape/ist/<u>/all_const05_none: x.mtx, all-Pelka profiles, neg
.05), so v2 is re-run bit-for-bit (checked against its out.csv) with log-likelihoods saved (scripts/insitutype_ml.R).
Candidates (v2 = R insitutypeML; spill = the same likelihood in Python (nb_logliks, checked against R on v2); the spatial prior acts on the log-likelihoods):
  v2            locked floor
  v2+spatial    + tempered mean-field neighbour prior (15-NN, beta 1, 5 sweeps; log-likelihoods / median top-1 - top-2 margin)
  spill         background matrix = R's own per-cell background (estimateBackground: lm(neg ~ library - 1)) + gamma
                (.15) x library x the 6 nearest cells' normalised profile (neighbour transcripts explained as spill)
  spill+spatial both
  spill_x       spill without the neighbours of the candidate type's class (pass-1 labels = v2); + spatial (beta 1, .5)
Scored as scripts/probe_f5.py: held-out AUROC, marker-positive recall, class shares, gated macro-F1 (GSE1).
Usage (env floors, OMP_NUM_THREADS=3): [SF_OUT=runs/v3] python scripts/v3_f5.py <unit> -> <out>/f5_<unit>.json"""
import json, os, shutil, subprocess, sys, time
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp
from scipy.spatial import cKDTree
from spatialfloors.data.genes import filter_reason
from spatialfloors.floors import f5
from spatialfloors.metrics.f5 import heldout_auroc, if_gate2, macro_f1_gated, CLASSES
from spatialfloors.provenance import stamp

u = sys.argv[1]; t0 = time.time(); log = lambda *s: print(f"[{time.time() - t0:6.0f}s]", *s, flush=True)
OUT = Path(os.environ.get("SF_OUT", "runs/v3")); W = OUT / "ist" / u; SRC = Path(f"runs/landscape/ist/{u}/all_const05_none"); GAMMA, BETA = .15, 1.
RS = str(Path(".pixi/envs/r/bin/Rscript").resolve())
a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); MK = f5.for_unit(u, json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"])
keep = filter_reason(a.var_names, "human") == ""; Xf, vf = sp.csr_matrix(a.X[:, keep]), list(map(str, a.var_names[keep])); libf = np.maximum(np.asarray(Xf.sum(1)).ravel(), 1)
hs = {c: np.log1p(np.asarray(Xf[:, [vf.index(g) for g in MK[c] if g in vf]].sum(1)).ravel() / libf * 1e4) for c in CLASSES}
gate = if_gate2(a.obs.IF_PanCK.values, a.obs.IF_CD68_CK8_18.values)[0] if u.startswith("gse1_") else None
xy = a.obs[["x_um", "y_um"]].values.astype(float); del a
top = {c: hs[c] >= np.quantile(hs[c], .9) for c in CLASSES}; mpos = {c: top[c] & ~np.any([top[o] for o in CLASSES if o != c], 0) & (hs[c] > 0) for c in CLASSES}


def summary(c):
    h, per = heldout_auroc(c, Xf, vf, MK); r = dict(heldout_auroc=h, per_class={k: v["auroc"] for k, v in per.items()}, share={k: float((c == k).mean()) for k in (*CLASSES, "other")})
    r["marker_pos_recall"] = {k: float((c[mpos[k]] == k).mean()) for k in CLASSES}
    if gate is not None: r["gated_macro_f1"] = macro_f1_gated(gate, c, ("epithelial", "myeloid"))[0]
    return r


def run_r(tag):
    d = W / tag; d.mkdir(parents=True, exist_ok=True)
    for f in ("x.mtx", "profiles.csv", "neg.txt"):
        if not (d / f).exists(): shutil.copy(SRC / f, d / f)
    subprocess.run([RS, "scripts/insitutype_ml.R", str(d)], check=True)
    return pd.read_csv(d / "out.csv").clust.astype(str).values, pd.read_csv(d / "logliks.csv")


def spatial(L, nb, beta=BETA):
    ll = L.values.astype(float); s2 = np.sort(ll, 1); z0 = ll / max(float(np.median(s2[:, -1] - s2[:, -2])), 1e-6)
    R = np.exp(z0 - z0.max(1, keepdims=True)); R /= R.sum(1, keepdims=True)
    for _ in range(5):
        z = z0 + beta * R[nb].mean(1); z -= z.max(1, keepdims=True); R = np.exp(z); R /= R.sum(1, keepdims=True)
    return np.asarray(L.columns)[R.argmax(1)].astype(str)


def nb_logliks(B, names, chunk=4000, bg_fn=None):
    """insitutypeML in Python (InSituType 2.0: lldist + update_logliks_with_cohort_freqs, one cohort): NB size 10,
    mu = s * profile + bg, s = rowSums(max(x - bg, 0)) / sum(profile); bg a per-cell vector (lldist's vector branch, as
    R) or a cells x genes matrix (its matrix branch, which insitutypeML cannot call with several profiles), or per
    candidate type bg_fn(rows, k); logliks rounded to 2 digits, + log(cohort frequency + 1e-4) (nbaselinecells 100);
    cells without counts -> 'undefined'."""
    from scipy.special import gammaln
    n, G = x.shape; Kp = Pm.shape[1]; ll = np.zeros((n, Kp)); ps = Pm.sum(0); rr = 10.
    for s0 in range(0, n, chunk):
        r = np.arange(s0, min(n, s0 + chunk)); xd = x[r].toarray(); const = (gammaln(xd + rr) - gammaln(rr) - gammaln(xd + 1)).sum(1)
        if bg_fn is None: b = np.maximum(B[r][:, None] if B.ndim == 1 else B[r], 1e-5); sub = np.maximum(xd - b, 0).sum(1)  # R floors bg at 1e-5 (estimateBackground)
        for k in range(Kp):
            if bg_fn is not None: b = np.maximum(bg_fn(r, k), 1e-5); sub = np.maximum(xd - b, 0).sum(1)
            s = sub / ps[k]; s = np.where(s <= 0, xd.sum(1) / ps[k], s); mu = s[:, None] * Pm[None, :, k] + b
            ll[r, k] = const + (rr * np.log(rr / (rr + mu)) + xd * np.log(mu / (rr + mu))).sum(1)
    ll = np.round(ll, 2); cl = ll.argmax(1); n_k = np.bincount(cl, minlength=Kp); base = n_k / n
    f = (base * 100 + n_k) / (100 + n); ll = ll + np.log(f + 1e-4)[None]
    L = pd.DataFrame(ll, columns=names); lab = np.asarray(names)[ll.argmax(1)].astype(object); lab[np.asarray(x.sum(1)).ravel() == 0] = "undefined"
    return lab.astype(str), L


x = sio.mmread(str(SRC / "x.mtx")).tocsr(); assert x.shape[0] == len(xy), "x.mtx rows != cache cells"
Pdf = pd.read_csv(SRC / "profiles.csv"); Pm = Pdf.values.astype(float); pnames = list(Pdf.columns)
neg = np.loadtxt(SRC / "neg.txt"); lib = np.asarray(x.sum(1)).ravel(); bg0 = lib * float((neg * lib).sum() / (lib ** 2).sum())  # InSituType estimateBackground
nb6 = cKDTree(xy).query(xy, 7)[1][:, 1:]; nb15 = cKDTree(xy).query(xy, 16)[1][:, 1:]
lab = {}; lab["v2"], L2 = run_r("v2"); ref = pd.read_csv(SRC / "out.csv").clust.astype(str).values
rep = dict(unit=u, gamma=GAMMA, beta=BETA, v2_reproduced=float((lab["v2"] == ref).mean())); log("v2 rerun agreement with the locked run", rep["v2_reproduced"])
lab["v2+spatial"] = spatial(L2, nb15)
py_v2, _ = nb_logliks(bg0, pnames); rep["python_replica_agreement"] = float((py_v2 == lab["v2"]).mean()); log("Python insitutypeML replica vs R", rep["python_replica_agreement"])
S = np.zeros(x.shape, np.float32)
for s0 in range(0, x.shape[0], 10000):
    r = np.arange(s0, min(x.shape[0], s0 + 10000)); nbm = x[nb6[r].ravel()].toarray().reshape(len(r), 6, -1).sum(1); nbm /= np.maximum(nbm.sum(1, keepdims=True), 1)
    S[r] = bg0[r, None] + GAMMA * lib[r, None] * nbm
lab["spill"], LS = nb_logliks(S, pnames); del S; lab["spill+spatial"] = spatial(LS, nb15)
# spill_x: neighbours of the candidate type's class (pass-1 = v2 labels; F5 class, or the fine type outside the 3 classes)
# are not spill: a myeloid cell among macrophages keeps its myeloid evidence (spill cut myeloid recall .04-.05 on 278N / 118N)
grp = lambda l: np.array([c if c != "other" else str(t) for t, c in zip(l, f5.coarse(l))]); gk = grp(np.asarray(pnames)); nbg = grp(lab["v2"])[nb6]
NBX = {}


def bg_x(r, k):
    if r[0] not in NBX: NBX.clear(); NBX[r[0]] = x[nb6[r].ravel()].toarray().reshape(len(r), 6, -1)
    keep = (nbg[r] != gk[k]).astype(np.float32); nbm = (NBX[r[0]] * keep[..., None]).sum(1); nbm /= np.maximum(nbm.sum(1, keepdims=True), 1)
    return bg0[r, None] + GAMMA * lib[r, None] * nbm


lab["spill_x"], LX = nb_logliks(None, pnames, bg_fn=bg_x); lab["spill_x+spatial"] = spatial(LX, nb15); lab["spill_x+spatial_b05"] = spatial(LX, nb15, .5)
rep["summaries"] = {k: dict(**summary(f5.coarse(v)), agree_v2_fine=float((v == lab["v2"]).mean())) for k, v in lab.items()}
for k, s in rep["summaries"].items():
    log(k, round(s["heldout_auroc"], 4), {c: round(v, 3) for c, v in s["share"].items()}, "recall", {c: round(v, 3) for c, v in s["marker_pos_recall"].items()}, round(s.get("gated_macro_f1", np.nan), 3))
np.savez_compressed(OUT / f"f5_{u}_labels.npz", **{k: v for k, v in lab.items()})
(OUT / f"f5_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float)); print("V3F5_DONE", u, flush=True)
