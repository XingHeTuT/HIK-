# Changelog

## 2026-05-13 — Cold Start Experiment Complete

### Results
- Cold-start sweep (9α × 50 iterations) completed in ~12 hours
- U-shaped Err_high curve confirmed: optimal α≈1.20 (theory: π/2=1.57)
- Small α (≤0.5): Err_high=33-52 (4-6× worse than optimal)
- RMSE shows first evidence of α-dependence (1.571→1.561)

### Files Created
- `sweep_cold_start.m` — Cold-start sweep script
- `cold_start_final.mat` — Final results
- `cold_start_intermediate.mat` — Intermediate checkpoint

---

## 2026-05-12 — Hot Start Experiment Complete

### Results
- Pipeline (txt→xlsx→MDPR) completed in 8.2 hours
- Hot-start RMSE constant at ~1.58 rad across all 9 α values
- Err_high varies (1.2→25.5) but dominated by zero-coefficient artifacts
- Phase profile figures generated (9 FIG + 9 PNG)

### Key Finding
Hot start masks sampling strategy differences. MDPR converges to the same solution regardless of α because the ideal-coefficient initial guess is too close to the true solution.

### Files Created
- `pipeline_final_results.mat` — Hot-start results
- `pipeline_intermediate.mat` — Intermediate checkpoint
- `figures/phase_profile_alpha_*.fig` — Phase profile figures

---

## 2026-05-11 — Zemax PSF Batch Export

### Results
- MATLAB ZOS-API successfully connected via `ConnectAsExtension`
- 590 PSF files exported across 9 alpha configurations
- Surface 6 (user-added defocus surface, thickness=0=no defocus) verified working
- Peak ratio dz=1mm/dz=0 = 0.38, confirming defocus effect

### Technical Details
- Defocus surface: S6 (added by user in Zemax)
- PSF config: 512×512, 0.425 μm pixel, Huygens PSF
- Connection: `ZOSAPI_Connection()` → `ConnectAsExtension(0)`
- API pitfalls: `ConnectToApplication()` fails for manually-opened Zemax;
  `Tools.OpenHuygensPsf()` → correct: `Analyses.New_HuygensPsf()`

### Files Created
- `Master_Zemax_Batch_Export.m` — MATLAB ZOS-API batch export
- `alpha_configs.mat` — Defocus position lists for 9α
- `zemax_export/` — 590 PSF txt files (4.7 GB, gitignored)
- `export_psf_defocus.ZPL` — ZPL backup (NOT USED)
- `run_zemax_batch.py` — Python backup (NOT USED)

---

## 2026-05-10 — Framework Setup

### Environment
- Working directory: `D:\LiLinhan_2026\CC-多距离强度相位复原\HIK相位测量`
- MATLAB R2023a, Python 3.9.2, Zemax OpticStudio 2023 R1
- pythonnet 3.0.5 installed via manual wheel download (pip SSL blocked by proxy)
- Dependencies: clr_loader 0.2.9, cffi 2.0.0, pycparser 2.22

### Files Created
- `compute_theoretical_params.m` — Theoretical parameter calculator
- `run_mdpr_validation.m` — Parameterized MDPR validation function
- `generate_alpha_configs.m` — Defocus position generator
- `pipeline_zemax_to_validation.m` — Full hot-start pipeline
- `sweep_alpha_quick.m`, `sweep_alpha_v2.m`, `sweep_alpha_validation.m` — Earlier versions (OBSOLETE)
- `.gitignore` — Excludes large binary/data files

### Git
- Commit `d72bb44`: 26 files, 7608 lines
- Pushed to `github.com:XingHeTuT/HIK-.git` (master)

---

## Experiment Evolution

### v1 (Quick Test)
- 4 α values, 15 iterations, fixed Z_span = [-1.0, +0.5] mm
- **Problem**: K varied dramatically (30→3), confounding results
- **Status**: ABANDONED

### v2 (Fixed K)
- 6 α values, K=9 fixed, 30 iterations
- **Problem**: α=0.10 with K=9 covers only ±0.16 mm — too small a range
- **Status**: ABANDONED

### v3 (Fixed Z_span, Full K — FINAL)
- 9 α values, Z_span = ±5 mm, K varies naturally (277→10)
- Fresh Zemax PSF export for each α configuration
- Two experiments: hot-start (30 iter) and cold-start (50 iter)
- **Status**: COMPLETE
