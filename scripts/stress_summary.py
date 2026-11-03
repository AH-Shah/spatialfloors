"""Compact tables of the PLAN §13 development stress panels (runs/p3/f{1,2,3,5}_stress_<unit>.json).
Usage: python scripts/stress_summary.py [family ...] -> stdout (one row per method x unit)"""
import glob, json, sys
from pathlib import Path

FAM = sys.argv[1:] or ["f1", "f2", "f3", "f5"]
COLS = {"f1": ("embeddings", ["excess_k5", "excess_k15", "excess_k50", "top10_gene_share", "genes_improved", "depth_leak_r2",
                              "knn_macro_recall", "knn_rare_recall", "kmeans_ami_types", "IF_PanCK_nbr_spearman", "IF_CD45_nbr_spearman"]),
        "f2": ("partitions", ["ari", "ami", "matched_f1", "min_region_recall", "ari_boundary", "ari_interior", "ari_rotated_null",
                              "excess_over_blob", "size_entropy", "largest_share", "seed_stability", "ami_fine_regions", "ami_celltype", "nesting_homogeneity"]),
        "f3": ("variants", ["top10", "pc1", "top3", "gene_median", "r2", "worst_top10", "above_null"]),
        "f5": ("predictions", ["heldout_auroc", "worst_class", "depth_strat_mean", "threshold_mean", "gated_macro_f1"])}
short = lambda u: u.replace("gse303070_148191A_", "").replace("allen_animal4_C57BL6J-4.003", "a4")
fmt = lambda v: f"{v:.3f}" if isinstance(v, (int, float)) and v == v else "-"
for fam in FAM:
    key, cols = COLS[fam]; print(f"## {fam}  " + " ".join(cols))
    for f in sorted(glob.glob(f"runs/p3/{fam}_stress_*.json")):
        r = json.loads(Path(f).read_text()); u = short(r["unit"])
        if fam == "f2": print(f"# {u} k={r['k']} blob_ari={r['blob_ari']:.3f}")
        for name, v in r[key].items():
            v = v.get("summary", v); print(u, name, *(fmt(v.get(c)) for c in cols))
