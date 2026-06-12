# Automating the AMReX → CSV → CP-HIFI pipeline (work log)

This is a record of the work done with Claude to automate preprocessing of
the LBNL AMReX plotfile datasets and wire them into CP-HIFI. It captures the
goal, what was built, the key findings (including two bugs the manual
workflow was hiding), how everything was validated, and how to use it.

For the clean, standalone recipe see
[`../README_Preprocessing.md`](../README_Preprocessing.md); for dataset
specifics see [`README.md`](README.md). This file is the narrative.

---

## The goal

The original workflow turned AMReX plotfiles into CP-HIFI input by hand in
the ParaView GUI: open the AMReX reader, check the field arrays, apply a
**Slice**, apply **Cell Centers**, then **Save Data** as CSV — one timestep
at a time. That does not scale to the oscillating-cylinder dataset:

- `osc_cyl_all.zip` is **~9.6 GB** (Zip64), **234 timesteps**, ~476 MB each
  uncompressed (~111 GB unpacked) — too big to extract in full.
- Multiple fields and many timesteps make manual clicking impractical.

The ask: automate the export, stream the big archive within the available
disk, and mirror the vortex folder's MATLAB code so the cylinder data can be
loaded and run through CP-HIFI the same way.

---

## What was built

| File | Location | Role |
|---|---|---|
| `amrex_to_csv.py` | repo root | Headless ParaView exporter (reader → Slice → Cell Centers → CSV). Two input modes: stream from a zip, or read unpacked `plt*` dirs in place. Resumable. |
| `README_Preprocessing.md` | repo root | Standalone recipe: prerequisites, commands, CSV format, gotchas. |
| `load_osc_cyl.m` | `lbnl_data/` | MATLAB loader: header-aware column mapping, auto-detects the constant slice axis, de-duplicates coincident points, returns `(c1,c2,t,v)`. |
| `xytv_to_tensor_unaligned.m` | `lbnl_data/` | Generic packer (copied verbatim from the vortex folder) so `lbnl_data` is self-contained. |
| `run_cp_hifi.m` | `lbnl_data/` | End-to-end: paths → load → mask solid cells → tensor → CP-HIFI → viz → reconstruction → storage. Mirrors the vortex script, adapted to `(y,z,t)`. |

The 24 stride-10 slices are exported under `lbnl_data/csv/`
(`osc_cyl__0.csv … osc_cyl__2300.csv`, ~70 MB each).

---

## The pipeline

For each selected timestep `amrex_to_csv.py` runs the exact GUI recipe
headlessly under ParaView's bundled Python (`pvpython`):

```
AMReX reader (all fields, Level≤99)  →  Slice (plane)  →  Cell Centers  →  CSV
```

For the cylinder it streams **one timestep at a time straight out of the
zip** using Python's `zipfile` (selective, Zip64-safe extraction): extract
`plt<NNNNN>/`, export, delete, move on. **Peak disk ≈ one timestep (~0.5 GB)**
instead of ~111 GB. The run is resumable — a timestep whose CSV already
exists is skipped — so you can do a coarse subset now and fill in later.

The command that produced the current `csv/` (run from `lbnl_data/`):

```powershell
& "C:\Program Files\ParaView 6.1.0\bin\pvpython.exe" `
  ..\amrex_to_csv.py `
  --zip osc_cyl_all.zip --out-dir csv --prefix osc_cyl `
  --fields density,y_velocity,z_velocity,vorticity,vfrac `
  --normal 1,0,0 --origin 0,0,0 --level 99 --stride 10
```

24 timesteps, 879,360 points each, ~50–75 s/step, 0 failures, ~23 min total.

---

## Key findings

Three things surfaced while building and validating this. The first two were
latent problems in the *manual* workflow; the third would have corrupted the
CP-HIFI fit on the new dataset.

### 1. The manual ParaView export was under-resolved

ParaView's AMReX reader defaults to `Level=1`, loading only the two coarsest
AMR levels and silently dropping finer cells. On the vortex this lost ~13%
of the slice points (19,968 vs 22,872 at full resolution) and under-resolved
the vortex core. The exporter caps the level high (`--level 99`) so the
**finest mesh is always loaded**. (Validated: at `Level=1` the script
reproduces the old reference CSVs byte-for-byte; the extra points at full
resolution are exactly the finest-level cells.)

### 2. The system `unzip` silently drops the data

This machine's `unzip` wildcard does **not** span `/`, so
`unzip "plt00010/*"` extracts only the `Header` and misses every `Level_*`
data file — you would feed ParaView a header with no data. The script
sidesteps this entirely with Python's `zipfile` and exact-prefix matching.
(`unzip -l` for listing is fine; it only reads the table of contents.)

### 3. The x=0 slice duplicates points, and the tensor builder *sums* them

The cylinder's x=0 slice lands exactly on a cell face (8 cells, symmetric
about x=0), so ParaView emits **two coincident cell-centers per (y,z)
location** with identical values: each CSV has 879,360 rows but only 565,888
distinct points. That is harmless *only if de-duplicated first*, because
`tensor_unaligned` (and Tensor Toolbox `sptensor`) **sum duplicate
subscripts**. Verified directly: at t=0, where `z_velocity` is uniformly 50,
the raw tensor came out with range **[50, 100]** — 313,472 locations doubled.

**Fix:** `load_osc_cyl.m` drops exact `(c1,c2,t)` duplicates (lossless, since
the values are identical) and warns if any coincident points ever disagree.
After the fix the t=0 tensor is a clean **[50, 50]**.

---

## Validation performed

- **Exporter vs reference (vortex):** reproducing the manual recipe
  (`Level=1`, z=0 slice, `y_velocity`) gave identical point count (19,968)
  and, at the reference's stored precision, **zero coordinate/value
  discrepancy** across all points.
- **Streaming mode (cylinder):** one timestep end-to-end through the script
  — extract → export → auto-delete — confirmed before the full run.
- **MATLAB load chain (R2024b):** load → de-dup → vfrac mask → tensor build
  on real CSVs; sizes, ranges, field mapping, and the dedup fix all checked
  (t=0 tensor values `[50,50]`, `nDupesDropped` as expected).
- The full CP-HIFI **solve was not run** here (heavy: ~21M nonzeros over 24
  steps) — that is meant to be run interactively via `run_cp_hifi.m`.

---

## Quick start

Export (or extend) the CSVs — see the command above, or
`pvpython amrex_to_csv.py --help`. Smaller `--stride` adds more timesteps;
existing files are skipped.

Then in MATLAB:

```matlab
cd lbnl_data
run('run_cp_hifi.m')          % load → mask → tensor → CP-HIFI → plots
```

Or just the data:

```matlab
XYTV = load_osc_cyl('z_velocity');   % N x 4: [y, z, t, z_velocity]
TU   = xytv_to_tensor_unaligned(XYTV);
```

---

## Open items / optional follow-ups

- **CSV size:** the face-slice duplication makes the CSVs ~1.5× larger than
  necessary. Re-exporting on a cell-center plane (`--origin 0.01953125,0,0`,
  physically identical in this thin-x domain) yields clean, half-size files
  and removes the need for loader de-dup. ~20 min; not required.
- **CP-HIFI tuning:** `run_cp_hifi.m` carries the vortex's kernel/
  regularization (`kernfunc_gaussian(1)`, `lambda=1e-1`, `rho=1e-6`, `R=30`)
  as a starting point. These will likely need retuning for this dataset's
  `(y ∈ [-5,5], z ∈ [0,20])` scales and much larger size.
- **Field choice:** `run_cp_hifi.m` defaults to `z_velocity`; `vorticity`
  may be a more interesting target for the wake dynamics. Change the `field`
  variable at the top.
- **Re-exporting the vortex** at full resolution (`--level 99`) if you want
  its CP-HIFI input to include the finest level it was missing.
