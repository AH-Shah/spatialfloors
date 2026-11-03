"""P0: crop manifests (coordinates + region labels only) and the GSE303070 donor map. Output: manifests/crops/*.csv.gz
(cell ids) + manifests/crops/index.json, manifests/gse303070_donors.json."""
import json, numpy as np, pandas as pd, anndata as ad
from pathlib import Path
from spatialfloors.data.crops import choose_crop
from spatialfloors.provenance import stamp

RAW = Path("/Users/Hemato/project_3/data/raw"); OUT = Path("manifests/crops"); OUT.mkdir(parents=True, exist_ok=True); index = {}


def save(name, ids, rec, extra=None):
    pd.Series(np.asarray(ids), name="cell_id").to_csv(OUT / f"{name}.csv.gz", index=False); rec = dict(rec, **(extra or {})); index[name] = rec
    print(name, {k: rec.get(k) for k in ("n_slide", "n_crop", "whole_slide", "min_share", "meets_min_share", "crop_shares")}, flush=True)


# GSE1 (231 excluded: same patient as 232)
for s in [110, 120, 210, 221, 222, 232, 242]:
    a = ad.read_h5ad(RAW / f"GSE1/cosmx_wt_{s}.h5ad", backed="r"); o = a.obs
    xy = o[["CenterX_global_px", "CenterY_global_px"]].values.astype(float); reg = o["typ"].astype(str).str.extract(r"(REF|TVA|CRC)", expand=False).values
    idx, rec = choose_crop(xy, reg); save(f"gse1_{s}", o.index[idx], rec, dict(dataset="GSE1", role="development" if s == 232 else "confirmation", region_key="typ")); a.file.close()

# GSE0 S0 (slide B; no annotated regions -> densest contiguous square)
m = pd.read_csv(RAW / "GSE0/S0_metadata_file", usecols=["cell_id", "fov", "CenterX_global_px", "CenterY_global_px"])
idx, rec = choose_crop(m[["CenterX_global_px", "CenterY_global_px"]].values.astype(float)); save("gse0_S0", m.cell_id.values[idx], rec, dict(dataset="GSE0", role="development", region_key=None))

# GSE303070: unit = donor; dev slide 148191A, its donors excluded downstream; NC (no donor) FOVs dropped
M = pd.read_csv(RAW / "GSE303070/GSE303070_All_cells_metadata.csv", index_col=False,
                usecols=["Barcode", "fov", "CenterX_global_px", "CenterY_global_px", "Slide", "Sample", "Donor", "Sample_type", "niches", "Mean.PanCK", "Mean.CD45", "Mean.CD68", "Mean.CD3"])
M[["Slide", "Sample"]] = M[["Slide", "Sample"]].astype(str); M.loc[M.Donor.astype(str).isin(["DonorC", "nan"]), "Donor"] = np.nan  # mixed int/str slide ids; NC FOVs have no donor
fm = pd.read_csv("/Users/Hemato/project_3/antig/N2/trials/fap_spp1/gse303070_fov_meta.csv"); fm[["slide", "pref"]] = fm[["slide", "pref"]].astype(str)
mc = M.groupby(["Slide", "Sample"]).agg(cells=("Barcode", "size"), fovs=("fov", "nunique"), donor=("Donor", lambda x: ",".join(sorted(set(map(str, x.dropna())))))).reset_index()
fc = fm.groupby(["slide", "pref"]).agg(cells=("total cells", "sum"), fovs=("fov", "nunique"), donor=("donor", lambda x: ",".join(sorted(set(map(str, x.dropna())))))).reset_index()
J = mc.merge(fc, left_on=["Slide", "Sample"], right_on=["slide", "pref"], how="outer", suffixes=("_cells", "_fovcsv"))
J["match"] = (J.cells_cells == J.cells_fovcsv) & (J.fovs_cells == J.fovs_fovcsv) & (J.donor_cells == J.donor_fovcsv)
DEV, CONF = ["148191A"], ["10501A", "120131", "3461G", "Tumour_A", "Tumour_B"]
dev_donors = sorted(set(M.loc[M.Slide.isin(DEV), "Donor"].dropna()))
units = []
for (sl, sa), g in M[M.Slide.isin(DEV + CONF)].groupby(["Slide", "Sample"]):
    d = g.Donor.dropna().unique()
    if len(d) == 0: continue  # NC FOVs
    d = d[0]; role = "development" if sl in DEV else ("excluded_dev_donor" if d in dev_donors else "confirmation")
    idx, rec = choose_crop(g[["CenterX_global_px", "CenterY_global_px"]].values.astype(float))
    units.append(dict(slide=sl, sample=sa, donor=d, role=role, cells=len(g), fovs=int(g.fov.nunique()), crop=int(rec["n_crop"]),
                      if_available={c: bool(g[f"Mean.{c}"].notna().mean() > .9) for c in ["PanCK", "CD45", "CD68", "CD3"]}, niches_labelled=round(float(g.niches.notna().mean()), 3)))
    if role != "excluded_dev_donor": save(f"gse303070_{sl}_{sa}", g.Barcode.values[idx], rec, dict(dataset="GSE303070", role=role, donor=d, region_key=None))
conf_donors = sorted({u["donor"] for u in units if u["role"] == "confirmation"})
multi = {d: sorted({u["slide"] for u in units if u["donor"] == d and u["role"] == "confirmation"}) for d in conf_donors}
json.dump(stamp(dict(check_vs_fov_csv=J.to_dict("records"), all_match=bool(J.match.all()), dev_slides=DEV, dev_donors=dev_donors, conf_slides=CONF,
                     units=units, confirmation_donors=conf_donors, n_confirmation_donors=len(conf_donors),
                     confirmation_donors_on_several_conf_slides={d: v for d, v in multi.items() if len(v) > 1})), open("manifests/gse303070_donors.json", "w"), indent=1, default=str)
print("GSE303070 csv match:", bool(J.match.all()), "| dev donors", dev_donors, "| confirmation donors", conf_donors)
print(J[~J.match].to_string() if (~J.match).any() else "")

# SPATCH Xenium 5K (F1 only)
a = ad.read_h5ad(RAW / "spatch_xenium_5k/spatch_xenium_5k.h5ad", backed="r"); idx, rec = choose_crop(np.asarray(a.obsm["spatial"], float))
save("spatch_xenium5k", a.obs_names[idx], rec, dict(dataset="SPATCH", role="confirmation_F1", region_key=None)); a.file.close()

# Allen animal 2 (on disk): the three F4 sections, regions = major_brain_region
A = ad.read_h5ad("/Users/Hemato/project_3/antig/N2/trials/perturb/WB_MERFISH_animal2_coronal.h5ad", backed="r")
for sid in ["C57BL6J-2.036", "C57BL6J-2.039", "C57BL6J-2.041"]:
    m = (A.obs.brain_section_label == sid).values; xy = np.asarray(A.obsm["X_spatial_coords"][np.flatnonzero(m)], float)
    idx, rec = choose_crop(xy, A.obs.major_brain_region.values[m].astype(str)); save(f"allen2_{sid}", A.obs_names[m][idx], rec, dict(dataset="Allen Zhuang animal 2", role="retrospective", region_key="major_brain_region"))
json.dump(stamp(index), open(OUT / "index.json", "w"), indent=1, default=str)
