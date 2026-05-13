# Multi-Distance Phase Retrieval — Theory Validation Framework

## Project Overview

Validating the defocus sampling theory for large-aperture refractive-meta hybrid systems (Paper Section 2.2-2.3). The theory predicts an optimal dimensionless defocus step **α_opt = π/2 ≈ 1.57** for high-frequency phase recovery, based on Fisher Information Matrix analysis.

**System under test**: D50.4f200+ML1 (D=44.6mm, f≈200mm, λ=10.6μm, z_c=0.363mm)

## Quick Reference

| Item | Path |
|------|------|
| Paper | `../多距离强度相位复原论文/初稿R2.3V2.docx` |
| Main MDPR code | `../phi_re_Gemini_final_v5forHIK_v5.m` |
| Zemax system | `../D50.4f200+ML1_test_v1/D50.4f200+ML1.zos` |
| PSF export data | `zemax_export/` (590 files, 4.7 GB, gitignored) |
| MATLAB ZOS-API script | `zemax_automation/Master_Zemax_Batch_Export.m` |
| Git remote | `git@github.com:XingHeTuT/HIK-.git` |

## Directory Structure

```
validate_theory/
├── compute_theoretical_params.m    # Calculate z_c, α, SBP for all systems
├── run_mdpr_validation.m           # Core: parameterized MDPR function
├── sweep_alpha_quick.m             # v1: quick test (4 alpha, 15 iter) — OBSOLETE
├── sweep_alpha_v2.m                # v2: fixed K=9 design — OBSOLETE
├── sweep_alpha_validation.m        # v1 sweep script — OBSOLETE
├── sweep_cold_start.m              # v3: COLD START (zero init), 9α×50iter ← CURRENT
├── cold_start_final.mat            # Cold start results
├── cold_start_intermediate.mat     # Cold start crash-recovery checkpoint
├── figures/                        # Phase profile plots (9α, FIG+PNG)
└── zemax_automation/
    ├── Master_Zemax_Batch_Export.m # MATLAB ZOS-API: batch PSF export
    ├── generate_alpha_configs.m    # Generate defocus position lists
    ├── alpha_configs.mat           # Generated: 9α configs with defocus positions
    ├── pipeline_zemax_to_validation.m  # Full pipeline: txt→xlsx→MDPR (HOT START)
    ├── pipeline_final_results.mat  # Hot-start results (9α×30iter)
    ├── pipeline_intermediate.mat   # Hot-start checkpoint
    ├── export_psf_defocus.ZPL      # ZPL macro backup (NOT USED)
    └── run_zemax_batch.py          # Python ZOS-API backup (NOT USED)
```

## Key Theoretical Parameters (System A)

```
NA = 0.0965
z_c = λ/(π·NA²) = 0.363 mm
α_opt = π/2 = 1.571  →  Δz_opt = 0.570 mm
SBP = D·NA/λ ≈ 406
```

## Experiment Design (v3, Final)

- **Fixed scan range**: Z_span = ±5 mm (β ≈ 28)
- **9 alpha values**: [0.10, 0.20, 0.50, 0.80, 1.20, π/2, 2.00, 2.50, 3.00]
- **Variable K (not capped)**: K = Z_span / Δz, ranging from 277 to 10
- **Zemax export**: Surface 6 controls defocus (thickness=0 = no defocus)
- **PSF settings**: 512×512, 0.425 μm pixel, Huygens PSF

## Results Summary

### Hot Start (phi_init = ideal Zemax coefficients, 30 iterations)

| α | K | RMSE | Err_high |
|---|----|------|----------|
| 0.10→3.00 | 277→10 | ~1.58 (constant) | 1.2→25.5 |

**Finding**: RMSE constant across all α — hot start masks sampling differences. MDPR converges to same solution regardless of measurement quality.

### Cold Start (phi_init = 0, 50 iterations)

| α | K | RMSE | Err_high |
|---|----|------|----------|
| 0.10 | 277 | 1.571 | 52.5 |
| 0.20 | 138 | 1.572 | 44.9 |
| 0.50 | 57 | 1.569 | 33.6 |
| 0.80 | 36 | 1.569 | 12.1 |
| **1.20** | 24 | 1.570 | **8.7** ← BEST |
| 1.571 | 19 | 1.569 | 22.6 |
| 2.00 | 15 | 1.569 | 21.7 |
| 2.50 | 13 | 1.566 | 60.5 |
| 3.00 | 11 | 1.561 | 31.9 |

**Finding**: U-shaped Err_high curve confirmed. Optimal α ≈ 1.20 (theory: π/2 = 1.57, deviation 23%). Small α (≤0.5) gives 4-6× worse Err_high than optimal. RMSE begins to show slight variation.

## How to Reproduce

### 1. Zemax PSF Export (already done)
```matlab
% In MATLAB:
cd('validate_theory/zemax_automation');
Master_Zemax_Batch_Export
% Requires: Zemax open with D50.4f200+ML1.zos, Interactive Extension enabled
```

### 2. Hot-Start Validation (already done)
```matlab
pipeline_zemax_to_validation
% txt→xlsx→MDPR with ideal-coefficient initial guess, 30 iterations
```

### 3. Cold-Start Validation (already done)
```matlab
sweep_cold_start
% Zero-phase initial guess, 50 iterations
```

### 4. View Results
```matlab
load('cold_start_final.mat');
% Variables: alpha_list, rmse_cold, err_high_cold, results_cold
```

## Known Issues

1. **RMSE nearly constant**: MDPR converges to same fixed point regardless of α. The constraint-based algorithm limits sensitivity to measurement information.
2. **Err_high metric**: Only meaningful for terms with |A_ideal| > 1e-12. Zero-coefficient terms (r^10+) inflate raw Err_high.
3. **Phase unit**: Zemax Binary 2 coefficients are in WAVES. MATLAB code in `phi_re_Gemini_final_v5forHIK_v5.m` may be missing 2π factor.
4. **zemax_export/ is 4.7 GB**: Gitignored. PSF data stored as xlsx files, one per defocus plane.

## Next Steps

- [ ] Low-frequency validation: cross-system β_min vs SBP (Task #3)
- [ ] Fix 2π factor in main MDPR code
- [ ] Generate paper-quality figures from cold start results
- [ ] Add Gaussian noise to PSF data for robustness testing
