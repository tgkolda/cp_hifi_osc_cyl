# Oscillating Cylinder Dataset from LBNL

Data from Mahesh Natarajan, provided via

Natarajan, M. (2026). Oscillating cylinder (AMReX plotfile) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.20091893

## What the simulation is

A 3D compressible flow over an oscillating cylinder, run in an AMReX-based
solver with embedded-boundary (cut-cell) geometry and block-structured
adaptive mesh refinement. The domain is a tall, thin channel with the
cylinder near the inflow; flow is along **z** at roughly Mach-scale speed
(`z_velocity` is ~50 in the uniform inflow). This is the larger, "more
interesting" follow-on to the isentropic-vortex sanity check in
[`../lbnl_data_isentropic_vortex/`](../lbnl_data_isentropic_vortex/).

Specifics, read straight from the plotfile headers:

- **Format**: AMReX/BoxLib plotfile (`HyperCLaw-V1.1`), 3D.
- **Domain**: `[-0.15625, 0.15625] x [-5, 5] x [0, 20]`. The **x** extent is
  only 8 cells wide at the coarsest level -- the physics is effectively 2D
  in the **y-z** plane, with x as the thin direction.
- **AMR**: 4 levels (0-3), refinement ratio 2 in each direction. The finest
  mesh tracks the cylinder and its wake.
- **Fields** (5): `density`, `y_velocity`, `z_velocity`, `vorticity`,
  `vfrac`. Note there is **no `x_velocity`**; `vfrac` is the embedded-
  boundary volume fraction (0 inside the solid cylinder, 1 in open flow).
- **Snapshots**: 234 plotfiles `plt00000, plt00010, ..., plt02330`
  (stride 10 in timestep count, dt-uniform).

Everything ships as a single ~9.6 GB archive `osc_cyl_all.zip`; each
timestep is ~476 MB uncompressed.

## Directory layout

```
lbnl_data/
├── README.md             <- this file
├── osc_cyl_all.zip       <- raw data from Zenodo (Zip64, ~9.6 GB; do not modify)
├── csv/                  <- exported x=0 slices, one CSV per timestep
│   ├── osc_cyl__0.csv
│   ├── osc_cyl__100.csv
│   └── ...  (osc_cyl__<plt index>.csv)
├── load_osc_cyl.m        <- MATLAB loader (csv/ -> coordinate format)
├── xytv_to_tensor_unaligned.m  <- packs [c1 c2 t v] rows into a tensor_unaligned
├── run_cp_hifi.m         <- end-to-end: load -> mask -> tensor -> CP-HIFI
```

The headless ParaView exporter lives one level up:

```
../amrex_to_csv.py        <- AMReX plotfile -> CSV (used by both datasets)
```

## From raw data to CP-HIFI input

The manual ParaView GUI recipe used for the vortex (open reader, check Cell
Arrays, Slice, Cell Centers, Save Data) does not scale to 234 timesteps of
a 9.6 GB archive that will not fit on disk if fully unpacked. The pipeline
below is fully scripted instead.

### 1. Export slices with `amrex_to_csv.py`

[`../amrex_to_csv.py`](../amrex_to_csv.py) drives ParaView headlessly:
for each selected timestep it runs **AMReX reader -> Slice -> Cell Centers
-> CSV writer**, exactly the GUI recipe, and writes one CSV per timestep
with all requested fields as columns:

```
Points:0, Points:1, Points:2, density, y_velocity, z_velocity, vorticity, vfrac
```

For this dataset it reads timesteps **straight out of the zip, one at a
time**, using Python's `zipfile` for selective Zip64 extraction (the
system `unzip` is unreliable here: its `*` wildcard does not span `/`, so
`unzip "plt00010/*"` silently grabs only the `Header` and misses the
`Level_*` data). Each timestep is extracted to `_tmp_extract/`, exported,
then deleted before the next, so **peak disk stays around one timestep
(~0.5 GB)** instead of ~111 GB for all 234.

It must run under ParaView's bundled Python (`pvpython`). On this machine,
to export **all 234 timesteps** (clear any stale CSVs and the cached MATLAB
load first, so old pre-de-dup files are not left behind):

```powershell
cd "C:\Users\TammyKolda\Documents\GIT Repos\CP-HIFI\lbnl_data"
Remove-Item csv\osc_cyl__*.csv -ErrorAction SilentlyContinue   # drop stale CSVs
Remove-Item D.mat -ErrorAction SilentlyContinue                # drop cached load

& "C:\Program Files\ParaView 6.1.0\bin\pvpython.exe" `
  ..\amrex_to_csv.py `
  --zip osc_cyl_all.zip `
  --out-dir csv `
  --prefix osc_cyl `
  --fields density,y_velocity,z_velocity,vorticity,vfrac `
  --normal 1,0,0 --origin 0,0,0 `
  --level 99 `
  --stride 1 `
  --gzip
```

`--gzip` writes `osc_cyl__<n>.csv.gz` (~5-12x smaller; these low-precision,
repetitive slices compress extremely well). `load_osc_cyl.m` reads `.csv.gz`
transparently, decompressing one timestep at a time. Omit it for plain CSVs;
to compress an existing plain export after the fact, `gzip csv\osc_cyl__*.csv`
(or PowerShell `Get-ChildItem csv\osc_cyl__*.csv | ForEach-Object { & gzip $_ }`).
Do not leave both `osc_cyl__5.csv` and `osc_cyl__5.csv.gz` for the same index
-- the loader would read both and double that timestep.

- `--normal 1,0,0 --origin 0,0,0` takes the **x=0 slice** (the y-z plane),
  dropping the thin x direction. The raw slice has 879,360 cell-centers; the
  exporter de-dups the slice-on-cell-face copies down to **565,888** per CSV
  (see the duplicate-points note below). Pass `--no-dedup` to keep the raw
  copies.
- `--level 99` caps the AMR level high so the **finest mesh is loaded**.
  This matters: ParaView's GUI default is `Level=1`, which silently drops
  the finest cells. (On the vortex dataset that default lost ~13% of the
  points and under-resolved the core; see the validation note below.)
- `--stride 1` processes **every** plotfile (all 234:
  `plt00000, plt00010, ..., plt02330`). Use `--stride 10` for a coarse
  24-step subset, or `--start/--stop` (in plt-number units) and
  `--steps plt00500,plt01000` to pick an explicit subset.

The run is **resumable**: a timestep whose CSV already exists is skipped,
so you can export a coarse subset now and fill in more later without
redoing work. (When refreshing existing CSVs after a code change, either
delete them first as above or add `--overwrite`, since a skip keeps the old
file.) `pvpython amrex_to_csv.py --help` lists every option.

The same script also reads already-extracted plotfiles in place via
`--input-dir DIR` (used to re-export the vortex dataset; see that folder).

### 2. Load into MATLAB

```matlab
cd lbnl_data
XYTV = load_osc_cyl('z_velocity');     % N x 4: [y, z, t, z_velocity]
size(XYTV)

% Sanity plot: first timestep, y-z plane.
m = XYTV(:,3) == min(XYTV(:,3));
scatter(XYTV(m,1), XYTV(m,2), 5, XYTV(m,4), 'filled');
axis equal; colorbar; xlabel('y'); ylabel('z');
title('z\_velocity, first snapshot');
```

[`load_osc_cyl.m`](load_osc_cyl.m) reads the per-timestep CSVs, maps
columns by **name** from the header (so column order is irrelevant), and
auto-detects the constant slice axis (x here) -- returning the two in-plane
coordinates plus time and the chosen field, in CP-HIFI's `(c1, c2, t, v)`
coordinate format. Pass any exported field name; default is `z_velocity`.

To mask cells inside the cylinder, also load `vfrac` and drop the solid
cells:

```matlab
V  = load_osc_cyl('z_velocity');
VF = load_osc_cyl('vfrac');
keep = VF(:,4) > 0.5;     % open-flow cells only
V = V(keep,:);
```

Then pack into a tensor with the (generic) packer in this folder:

```matlab
TU = xytv_to_tensor_unaligned(V);
```

Or just run [`run_cp_hifi.m`](run_cp_hifi.m), which does the load, masks the
solid cells, builds the tensor, and fits CP-HIFI end to end (mirrored from
the vortex folder's script, adapted to this dataset's (y, z, t) coordinates).

## Validation

The exporter was checked against the vortex dataset's manually-produced
reference CSVs. Reproducing the manual recipe (`Level=1`, z=0 slice,
`y_velocity`) gave **byte-for-byte agreement**: identical point count
(19,968) and values matching to the reference's stored precision. Run at
full resolution (`--level 99`) the same slice yields 22,872 points -- the
extra ~2,900 are the finest level-2 cells the manual GUI export had been
dropping. The automation is faithful to the recipe and strictly more
accurate at the AMR boundaries.

## Notes

- **Time coordinate**: the loader uses the integer plt index (`0, 100, 200,
  ...`) as `t`. For CP-HIFI only relative spacing matters and the stride is
  uniform, so this is a fine time axis.
- **Duplicate slice points**: the x=0 slice lies on a cell face, so the raw
  slice has two coincident cell-centers per (y, z) with identical values
  (879,360 rows, 565,888 distinct). `amrex_to_csv.py` now drops these at
  export time (keeping the first of each (y, z)), so the CSVs are already
  duplicate-free (565,888 rows) and `load_osc_cyl.m` does no de-dup pass --
  it only sanity-checks the first file. Pass `--no-dedup` to keep the raw
  copies; if you do, collapse on (y, z, t) yourself before building the
  tensor, which *sums* duplicate subscripts and would otherwise double the
  values. CSVs exported before this change are stale: the loader will error
  and ask you to re-export.
- **Provenance warning**: `csv/` and `_tmp_extract/` are derived from
  `osc_cyl_all.zip`. If the archive is replaced, delete `csv/` and re-run
  the export.
- **Why `pvpython`**: the AMReX/BoxLib reader and `paraview.simple` live
  inside ParaView's own Python. The system/conda Python cannot import them,
  so the exporter is always launched via `pvpython.exe`.
```
