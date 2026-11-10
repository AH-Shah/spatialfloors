# Deviations and additions after registration

Registration: DOI [10.5281/zenodo.23931264](https://doi.org/10.5281/zenodo.23931264). The full run log of the confirmatory
analysis is archived with `results-v1` (DOI [10.5281/zenodo.24044400](https://doi.org/10.5281/zenodo.24044400)).

| Date | Event | Effect on the results |
|---|---|---|
| 2026-11-04 | F3: on three GSE303070 crops (`120131_170T`, `120131_280T`, `3461G_7T`) every receiver identity falls into one lineage, so no model has an other-lineage context. They are not scored (PROTOCOL §6, model failure); their tier-2 jobs were not run. | F3 has 22 of 25 units: donors 7 and 170 are missing, donor 280 is scored on 280N only. |
| 2026-11-05 – 11-06 | Confirmation jobs that crashed were rerun after fixes to output paths and memory use. GSE1 120's L6 floor uses a memory-lean loader that reproduces the registered loader exactly (`scripts/l6_lean_check.py`). | None: no fix changes what a job computes. |
| 2026-11-08 | Before release, the tier-1 analysis was rerun once: the first run read GSE1 222's F4 scores with the wrong held-out type and so left that unit out of F4. The released report (`P5_corrected`) includes every registered unit; the first run is archived with `results-v1`. | Verdicts identical. |
| 2026-11-08 – 11-10 | Exploratory analyses, no decision weight. Listed in PROTOCOL §6: interior-only re-score, unit size as a moderator. Added after `results-v1`: F5 at a matched call rate, F1 with Euclidean k-NN, F4 by cancer and provenance, replication levels, positive controls, development effort. Not run: stress panels on the confirmation units, the protein shift-null, the F1 depth-matched reference. | Reported in `RESULTS.md` as exploratory. |
| 2026-11-09 – 11-10 | Whole-slide arm on the six GSE1 confirmation slides, planned in `results/WHOLESLIDE_PLAN.md` before any whole-slide output existed. Before the run the stop rule was changed so that a failed job would be set aside rather than stop the arm; no job failed. A diagnostic of the F5 slide-scale result followed (`X_f5_scale`). | Exploratory; reported in `RESULTS.md`. |
| 2026-11-10 | SPATCH Xenium 5K tissue, unverified at registration, identified as colon adenocarcinoma (DATA.md). The registered leave-one-cancer-out keeps it as its own group. | With it counted as CRC, F1 without CRC stays non-inferior (upper bounds: scVI −.0061, Nicheformer −.0269). |
