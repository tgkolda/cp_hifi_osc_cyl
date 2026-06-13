# CP-HIFI on the Oscillating Cylinder

Fits a CP-HIFI model to a 3D tensor of `(y, z, t)` flow data from an
oscillating-cylinder simulation, and visualizes the reconstruction.

## Quick start

```matlab
cd lbnl_data
run_cp_hifi
```

That script does everything end to end: load CSVs (cached to `D.mat`), mask
to a spatial/time window, randomly subsample to a target nnz, build the
tensor, fit CP-HIFI, and draw the figures. First run builds `D.mat` (slow);
later runs reload it.

Tensor modes: **1 = y**, **2 = z**, **3 = t**.

## Knobs (top of `run_cp_hifi.m`)

| Variable | What it does |
|---|---|
| `field` | which field to fit (`density`, `y_velocity`, `z_velocity`, `vorticity`) |
| `ylim_keep`, `zlim_keep` | spatial window kept for the tensor |
| `val_scale` | multiplier on field values before fitting |
| `mask_vfrac` | drop solid-cylinder cells (`vfrac<=0.5`) when `true` |
| `target_nnz` | random-subsample cap on nonzeros (`Inf` to disable) |
| `gauss_width{1,2,3}` | RKHS Gaussian kernel sigma per mode (also drawn on the fiber plots) |
| `R` | CP-HIFI rank |

## Figures

| # | Shows |
|---|---|
| 1–3 | mode-1/2/3 fibers, with the per-mode Gaussian kernel overlaid |
| 4 | static comparison grid: full data / sampled data / reconstruction, a few timesteps |
| 5 | `viz` of the first few CP components |
| 6 | side-by-side **movie** (full \| sampled \| reconstruction) animated over all timesteps |

Movie speed: set `fps` in the figure-6 section (default `10`).

## Regenerating the data

The CSVs and the loader are derived from the raw Zenodo archive
(`osc_cyl_all.zip`). That export pipeline (`amrex_to_csv.py` →
`load_osc_cyl.m`) is documented separately — see
[`README_from_CLAUDE.md`](README_from_CLAUDE.md). You only need it if the CSVs
are missing or the archive changes; for normal runs, just run the script.

Source: Natarajan, M. (2026). *Oscillating cylinder (AMReX plotfile)* [Data
set]. Zenodo. https://doi.org/10.5281/zenodo.20091893
