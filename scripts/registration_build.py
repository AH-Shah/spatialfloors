"""Registration package (registration/). Two steps:
  copy          locked configs -> configs/, crop / panel manifests -> units/, units/units.tsv (role, families, unit,
                stratum, cells, SHA-256 of each crop manifest), the sealed confirmation job list -> p4_jobs.txt
  seal <commit> code archive spatialfloors-<commit>.tar.gz (git archive; .gitattributes export-ignore applies) and
                SHA256SUMS over every package file and the archive
Usage (env floors): python scripts/registration_build.py copy | seal <commit>"""
import hashlib, json, shutil, subprocess, sys
from pathlib import Path
from spatialfloors.stats import stratum

R = Path("registration"); sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
CONFIGS = ["manifests/locked_v3.json", "manifests/locked_v3_ncem.json", "manifests/tier2.json", "manifests/power_v3.json", "manifests/locked_v2.json",
           "manifests/locked.json", "manifests/models.json", "manifests/allen.json", "manifests/gse303070_donors.json", "manifests/if_gate_gse1_232.json",
           "manifests/if_qc_gse303070.json", "manifests/panels/index.json"] + sorted(map(str, Path("manifests/refs").glob("*.json")))
ALLEN = ["allen_animal1_C57BL6J-1.086", "allen_animal2_C57BL6J-2.039", "allen_animal3_C57BL6J-3.003", "allen_yao638850_C57BL6J-638850.36"]

if sys.argv[1] == "copy":
    for d in ("configs", "units/crops", "units/panels"): (R / d).mkdir(parents=True, exist_ok=True)
    for f in CONFIGS: shutil.copy2(f, R / "configs" / f.replace("manifests/", "").replace("/", "_"))
    idx = json.loads(Path("manifests/crops/index.json").read_text()); shutil.copy2("manifests/crops/index.json", R / "units/crops/index.json")
    donor = {f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"] for d in json.loads(Path("manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]}
    conf = {u for u, v in idx.items() if isinstance(v, dict) and v.get("role") == "confirmation"}
    gse1, g3, bru = {f"gse1_{s}" for s in (110, 120, 210, 221, 222, 242)}, {u for u in conf if u.startswith("gse303070_")}, {"bruker_mobreast", "bruker_skin", "bruker_kidney"}
    FAM = (("F1", gse1 | g3 | bru | {"spatch_xenium5k"} | set(ALLEN)), ("F2", gse1 | set(ALLEN)), ("F3", gse1 | g3 | bru | set(ALLEN)), ("F4", {"gse1_110"}), ("F5", gse1 | g3 | bru))
    rows = ["crop\tunit\trole\tfamilies\tstratum\tn_cells\tcrop_manifest_sha256"]
    for u, v in sorted(idx.items(), key=lambda kv: (str(kv[1].get("role")) if isinstance(kv[1], dict) else "", kv[0])):
        if not isinstance(v, dict): continue
        f = Path(f"manifests/crops/{u}.csv.gz"); shutil.copy2(f, R / "units/crops" / f.name)
        p = Path(f"manifests/panels/{u}.json")
        if p.exists(): shutil.copy2(p, R / "units/panels" / p.name)
        role = {"confirmation_F1": "confirmation", "F1-F3": "confirmation", "retrospective": "confirmation (F4 retrospective)"}.get(v["role"], v["role"])
        fams = ",".join(k for k, s in FAM if u in s) if role == "confirmation" else ("F4" if "retrospective" in role else "")
        unit = donor.get(u, u.split("_C57BL6J")[0] if u.startswith("allen") else u)
        rows.append(f"{u}\t{unit}\t{role}\t{fams}\t{stratum(u.replace('allen2', 'allen'))}\t{v.get('n_crop', '')}\t{sha(f)}")
    (R / "units/units.tsv").write_text("\n".join(rows) + "\n"); shutil.copy2("runs/sched/p4_jobs.txt", R / "p4_jobs.txt")
    print(len(rows) - 1, "crops ->", R / "units/units.tsv", "| units.tsv", sha(R / "units/units.tsv"), "| p4_jobs.txt", sha(R / "p4_jobs.txt"))
elif sys.argv[1] == "seal":
    c = subprocess.run(["git", "rev-parse", "--short=12", sys.argv[2]], capture_output=True, text=True, check=True).stdout.strip(); arc = R / f"spatialfloors-{c}.tar.gz"
    subprocess.run(["git", "archive", "--format=tar.gz", f"--prefix=spatialfloors-{c}/", "-o", str(arc), c], check=True)
    files = sorted(p for p in R.rglob("*") if p.is_file() and p.name not in ("SHA256SUMS", ".DS_Store"))
    (R / "SHA256SUMS").write_text("".join(f"{sha(p)}  {p.relative_to(R)}\n" for p in files)); print(len(files), "files ->", R / "SHA256SUMS", "| archive", arc, sha(arc))
