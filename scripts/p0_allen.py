"""P0: Allen design check + section choice + crops (coordinates + region labels only). Rule fixed before any model:
a section passes if Thalamus and Isocortex hold >= 50 cells and Fiber_tracts >= 20 cells of each of the 5 Cellina types;
sections ranked by the minimum of those 15 counts; F4 uses the top 3 (animal 2: the camera-ready 036/039/041); F1-F3 use the
top-1 section, cropped to ~50k by spatialfloors.data.crops. Animal 4 (sagittal, lateral sections, no thalamus) fails ->
Allen (Yao) MERFISH-C57BL6J-638850 (500-gene panel; regions from CCF parcellation: TH, Isocortex, category 'fiber tracts')."""
import json, numpy as np, pandas as pd, anndata as ad
from spatialfloors.data.crops import choose_crop
from spatialfloors.provenance import stamp
T = ["GABAergic neuron", "astrocyte", "endothelial cell", "glutamatergic neuron", "oligodendrocyte"]; R = ["Thalamus", "Isocortex", "Fiber_tracts"]
from spatialfloors.paths import RAW, DERIVED
Z = f"{RAW}/allen_zhuang"; Y = f"{RAW}/allen_yao_638850"
SRC = {"animal1": f"{Z}/WB_MERFISH_animal1_coronal.h5ad", "animal2": f"{DERIVED}/WB_MERFISH_animal2_coronal.h5ad",
       "animal3": f"{Z}/WB_MERFISH_animal3_sagittal.h5ad", "animal4": f"{Z}/WB_MERFISH_animal4_sagittal.h5ad"}


def table(obs):
    rows = []
    for s, g in obs.groupby("section", observed=True):
        r = dict(section=s, n=len(g))
        for reg in R:
            gg = g[g.region == reg]; r[reg] = len(gg); r[reg + "_min5"] = int(min((gg.t5 == t).sum() for t in T))
        rows.append(r)
    D = pd.DataFrame(rows); D["pass"] = (D.Thalamus_min5 >= 50) & (D.Isocortex_min5 >= 50) & (D.Fiber_tracts_min5 >= 20)
    D["score"] = D[[r + "_min5" for r in R]].min(1); return D.sort_values(["pass", "score"], ascending=False)


out = {}; crops = json.load(open("manifests/crops/index.json"))
for an, p in SRC.items():
    A = ad.read_h5ad(p, backed="r"); o = A.obs
    obs = pd.DataFrame(dict(section=o.brain_section_label.astype(str).values, region=o.major_brain_region.astype(str).values, t5=o.cell_type.astype(str).values), index=o.index)
    D = table(obs); passing = D[D["pass"]].section.tolist()
    f4 = ["C57BL6J-2.036", "C57BL6J-2.039", "C57BL6J-2.041"] if an == "animal2" else passing[:3]
    rec = dict(source=p, title=str(A.uns.get("title")), n_cells=int(A.n_obs), n_genes=int(A.n_vars), orientation="sagittal" if "sagittal" in p else "coronal",
               n_sections=int(D.shape[0]), n_passing=len(passing), usable=len(passing) >= 3, f4_sections=f4 if len(passing) >= 3 else [],
               f1_f3_section=passing[0] if passing else None, table=D.head(15).to_dict("records"))
    if passing:
        s = passing[0]; m = (obs.section == s).values; xy = np.asarray(A.obsm["X_spatial_coords"][np.flatnonzero(m)], float)
        idx, cr = choose_crop(xy, obs.region.values[m]); name = f"allen_{an}_{s}"
        pd.Series(A.obs_names[m][idx], name="cell_id").to_csv(f"manifests/crops/{name}.csv.gz", index=False)
        crops[name] = dict(cr, dataset=f"Allen Zhuang {an}", role="F1-F3", region_key="major_brain_region"); rec["f1_f3_crop"] = name
    out[an] = rec; A.file.close(); print(an, rec["orientation"], "passing", len(passing), "F4", rec["f4_sections"], "F1-F3", rec["f1_f3_section"], flush=True)

m = pd.read_csv(f"{Y}/cell_metadata_with_parcellation_annotation.csv", usecols=["cell_label", "brain_section_label", "class", "subclass", "parcellation_division", "parcellation_category", "x_section", "y_section"])
c, s_ = m["class"].astype(str), m["subclass"].astype(str)
t5 = np.select([c.str.contains("Glut"), c.str.contains("GABA"), s_.str.contains(r"\bAstro"), s_.str.contains(r"\bEndo\b"), s_.str.contains(r"\bOligo NN")], T[3:4] + T[0:1] + T[1:3] + T[4:], default="other")
reg = np.select([m.parcellation_division == "TH", m.parcellation_division == "Isocortex", m.parcellation_category == "fiber tracts"], R, default=m.parcellation_division.astype(str))
obs = pd.DataFrame(dict(section=m.brain_section_label.values, region=reg, t5=t5), index=m.cell_label.values)
assert (obs.t5 == "glutamatergic neuron").sum() > 1e6  # Glut classes map to glutamatergic
D = table(obs); passing = D[D["pass"]].section.tolist(); s = passing[0]; mm = (obs.section == s).values
idx, cr = choose_crop(m[["x_section", "y_section"]].values[mm].astype(float), obs.region.values[mm]); name = f"allen_yao638850_{s}"
pd.Series(obs.index[mm][idx], name="cell_id").to_csv(f"manifests/crops/{name}.csv.gz", index=False); crops[name] = dict(cr, dataset="Allen Yao 638850", role="F1-F3", region_key="CCF parcellation")
out["yao638850"] = dict(source=f"{Y} (ABC atlas 20230830 raw h5ad, metadata 20231215 parcellation view)", n_cells=len(m), panel="500 genes (Yao), vs 1,120 (Zhuang)",
                        orientation="coronal", type_map="class contains Glut/GABA; subclass Astro*/Endo/Oligo NN", region_map="TH; Isocortex; parcellation_category 'fiber tracts'",
                        n_sections=int(D.shape[0]), n_passing=len(passing), f4_sections=passing[:3], f1_f3_section=s, f1_f3_crop=name, table=D.head(15).to_dict("records"))
out["decision"] = dict(prospective_allen=["animal1", "animal3", "yao638850"], replaced=dict(animal4="yao638850"),
                       reason="animal 4: 3 lateral sagittal sections, no thalamus (control domain) in any section; animal 3 sagittal passes (12 medial sections)")
json.dump(stamp(out), open("manifests/allen.json", "w"), indent=1, default=str); json.dump(crops, open("manifests/crops/index.json", "w"), indent=1, default=str)
print("yao", "passing", len(passing), "F4", passing[:3], "F1-F3", s, cr["n_crop"], cr.get("crop_shares"))
