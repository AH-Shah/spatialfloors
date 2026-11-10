"""P3 positive control for Novae (PLAN §7; decision 2026-10-28): the authors' tutorial slide Xenium_V1_Human_Colon_Cancer_P5
(HF dataset prism-oncology/novae, human/colon; sha256 in manifests/models.json) with their published annotations
(annotations/*.parquet: novae_domains_res0.8, 7 domains). Tutorial path: spatial_neighbors(radius=80 um),
pinned novae-human-0, compute_representations(zero_shot=True), assign_domains(level=7). Rule fixed before running:
pass iff ARI(ours, novae_domains_res0.8) >= .50 and the collapse gate holds (effective rank of the representation >= 5,
largest domain <= 50% of cells). Usage (env novae): python scripts/p3_posctl_novae.py [cpu|mps] -> runs/p3/posctl_novae.json"""
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score
from spatialfloors.provenance import stamp

D = __import__("spatialfloors.paths", fromlist=["RAW"]).RAW / "novae_p5"; acc = sys.argv[1] if len(sys.argv) > 1 else "cpu"
import novae
a = ad.read_h5ad(D / "slide.h5ad"); ann = pd.read_csv(D / "annotations.csv.gz", index_col=0).reindex(a.obs_names); t = time.time()  # parquet converted in env floors (no pyarrow here)
novae.utils.spatial_neighbors(a, radius=80)
m = novae.Novae.from_pretrained("prism-oncology/novae-human-0", revision="b8c0a5d7612bac6bc719ab57ed3cd16ad814728c")
m.compute_representations(a, zero_shot=True, accelerator=acc); key = m.assign_domains(a, level=7)
Z = np.asarray(a.obsm["novae_latent"]); s = np.linalg.svd(Z - Z.mean(0), compute_uv=False); p = s ** 2 / (s ** 2).sum(); erank = float(np.exp(-(p * np.log(p + 1e-12)).sum()))
ok = ann["novae_domains_res0.8"].notna().values; lab = a.obs[key].astype(str).values; share = float(pd.Series(lab).value_counts(normalize=True).iloc[0])
rep = dict(slide="Xenium_V1_Human_Colon_Cancer_P5_CRC_Add_on_FFPE_outs", n_cells=int(a.n_obs), n_annotated=int(ok.sum()), domain_key=key,
           ari_vs_published=float(adjusted_rand_score(ann["novae_domains_res0.8"].values[ok].astype(str), lab[ok])), effective_rank=erank,
           largest_domain_share=share, n_domains=int(len(set(lab))), s=round(time.time() - t, 1), accelerator=acc)
rep["pass"] = rep["ari_vs_published"] >= .5 and erank >= 5 and share <= .5
print(json.dumps(rep), flush=True); Path("runs/p3/posctl_novae.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
