"""P0: validate the GSE1 IF gate on the development slide 232 (SPEC §F5). Thresholds and pass criteria fixed here before
any number: confident-tail quantiles of log1p mean intensity over all cells of the slide.
  PanCK  +: >= q.60   -: <= q.40
  CD45   +: >= q.80   -: <= q.50
  CD68_CK8_18 +: >= q.85   -: <= q.50
3-class: epithelial = PanCK+ & CD45- ; myeloid = CD45+ & CD68ck+ & PanCK- ; immune_nonmyeloid = CD45+ & CD68ck- & PanCK-.
PASS (3-class) iff (a) gated myeloid mean CP10k of each held-out marker CD68, CD163, C1QC >= 2x the max over the other two
classes, (b) AUROC(myeloid vs immune_nonmyeloid; log1p CP10k sum of the 3 markers) >= .75, (c) AUROC(immune_nonmyeloid vs
epithelial; log1p CP10k sum of PTPRC, CD3E, IGKC, MS4A1) >= .75, (d) each class >= 200 cells.
Else fallback = 2-class gate (epithelial = PanCK+ ; myeloid = CD68ck+ & PanCK-), checked with (a) vs epithelial and
AUROC(myeloid vs epithelial; 3 markers) >= .75, reported separately (PLAN §5)."""
import json, numpy as np, anndata as ad, scipy.sparse as sp
from sklearn.metrics import roc_auc_score
from spatialfloors.provenance import stamp

a = ad.read_h5ad(__import__("spatialfloors.paths", fromlist=["RAW"]).RAW / "GSE1/cosmx_wt_232.h5ad")
L = {c: np.log1p(a.obs[f"Mean.{c}"].astype(float).values) for c in ["PanCK", "CD45", "CD68_CK8_18"]}
Q = dict(PanCK=(.40, .60), CD45=(.50, .80), CD68_CK8_18=(.50, .85))
thr = {c: (float(np.quantile(L[c], lo)), float(np.quantile(L[c], hi))) for c, (lo, hi) in Q.items()}
neg = {c: L[c] <= thr[c][0] for c in L}; pos = {c: L[c] >= thr[c][1] for c in L}
lib = np.asarray(a.X.sum(1)).ravel(); g = list(a.var_names); Xc = sp.csc_matrix(a.X)
cp = lambda n: np.asarray(Xc[:, g.index(n)].todense()).ravel() / np.maximum(lib, 1) * 1e4
MY, IM = ["CD68", "CD163", "C1QC"], ["PTPRC", "CD3E", "IGKC", "MS4A1"]
mys = np.log1p(np.sum([cp(n) for n in MY], 0)); ims = np.log1p(np.sum([cp(n) for n in IM], 0)); mk = {n: cp(n) for n in MY}
ist = a.obs["ist"].astype(str).values


def summarise(lab, classes):
    out = {}
    for k in classes:
        m = lab == k; out[k] = dict(n=int(m.sum()), coverage=round(float(m.mean()), 4), **{f"mean_cp10k_{n}": round(float(mk[n][m].mean()), 3) for n in MY},
                                    ist_top={s: int(c) for s, c in zip(*np.unique(ist[m], return_counts=True)) if c >= .05 * m.sum()})
    return out


auc = lambda y1, y0, s: float(roc_auc_score(np.r_[np.ones(y1.sum()), np.zeros(y0.sum())], np.r_[s[y1], s[y0]]))
lab3 = np.full(a.n_obs, "", object)
lab3[pos["PanCK"] & neg["CD45"]] = "epithelial"; lab3[pos["CD45"] & pos["CD68_CK8_18"] & neg["PanCK"]] = "myeloid"
lab3[pos["CD45"] & neg["CD68_CK8_18"] & neg["PanCK"]] = "immune_nonmyeloid"
C3 = ["epithelial", "immune_nonmyeloid", "myeloid"]; S3 = summarise(lab3, C3)
ok_n = all(S3[k]["n"] >= 200 for k in C3)
fold = {n: S3["myeloid"][f"mean_cp10k_{n}"] / max(max(S3[k][f"mean_cp10k_{n}"] for k in C3 if k != "myeloid"), 1e-9) for n in MY} if ok_n else {}
auc_my = auc(lab3 == "myeloid", lab3 == "immune_nonmyeloid", mys) if ok_n else None
auc_im = auc(lab3 == "immune_nonmyeloid", lab3 == "epithelial", ims) if ok_n else None
pass3 = bool(ok_n and min(fold.values()) >= 2 and auc_my >= .75 and auc_im >= .75)
lab2 = np.full(a.n_obs, "", object); lab2[pos["PanCK"]] = "epithelial"; lab2[pos["CD68_CK8_18"] & neg["PanCK"]] = "myeloid"
C2 = ["epithelial", "myeloid"]; S2 = summarise(lab2, C2)
fold2 = {n: S2["myeloid"][f"mean_cp10k_{n}"] / max(S2["epithelial"][f"mean_cp10k_{n}"], 1e-9) for n in MY}
auc2 = auc(lab2 == "myeloid", lab2 == "epithelial", mys); pass2 = bool(min(fold2.values()) >= 2 and auc2 >= .75)
# channel diagnosis: does each IF channel separate the RNA-labelled (InSituType) classes at all?
diag = {c: {k: round(float(np.median(L[c][np.isin(ist, v)])), 3) for k, v in dict(epi=["epi1", "epi2", "epi3", "epi4"], T=["TC"], plasma=["PC_IgA", "PC_IgG", "PC_IgM"],
                                                                                     myeloid=["mye1", "mye2"], fibro=["fib1", "fib2"]).items()} for c in L}
rep = stamp(dict(slide="GSE1 232", quantiles=Q, thresholds_log1p=thr, three_class=dict(classes=S3, fold_vs_other_max=fold, auroc_myeloid_vs_immune=auc_my,
           auroc_immune_vs_epithelial=auc_im, passed=pass3), two_class=dict(classes=S2, fold_vs_epithelial=fold2, auroc_myeloid_vs_epithelial=auc2, passed=pass2),
           channel_medians_by_insitutype=diag, decision="3-class" if pass3 else ("2-class fallback" if pass2 else "gate failed: F5 on GSE1 not testable")))
open("manifests/if_gate_gse1_232.json", "w").write(json.dumps(rep, indent=1, default=float)); print(json.dumps(rep, indent=1, default=float))
