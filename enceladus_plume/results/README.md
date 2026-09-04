# Committed fit caches

Small, version-controlled products of the diurnal-profile fit (`../run_fit.py`),
so that the manuscript fit figures can be redrawn and the fit numbers
re-reported without re-running the (slow) fit or rebuilding the lookup tables.

| File | Content | Produced by |
|---|---|---|
| `diurnal_fit.json` | Adopted best fit (mode B): Δδ, L, δ_eff*, α, φ₂, σ_φ, φ₀, A, χ² | `run_fit.py --refine` (local Nelder–Mead polish; deterministic) |
| `diurnal_fit_modeA.json` | Secondary optimum (mode A, φ₂ ≈ 84°) | `analysis/two_modes.py` |
| `gas_lut.npz` | Gas-column lookup table (r-table) used by every flux evaluation | `run_fit.py` (rebuilt only if missing) |
| `weff_grid.npz` | On-attractor closure width δ_eff*(Δδ, L) on a coarse grid | `run_fit.build_weff_interp` (rebuilt only if missing) |

Keys in the JSON files use the code's internal names: `dw` = Δδ [m], `L` [m],
`w_eff` = δ_eff* [m], `harm_scale` = α (2f/1f amplitude ratio), `harm_phase` = φ₂ [deg],
`sigma` = σ_φ [deg], `phi0` [deg], `A` (slab density per unit model flux), `chi2`, `dof`, `chi2_red`.

Redraw the manuscript figures (written into the paper repo, or set `ENC_PAPER_FIGS`):

```bash
cd enceladus_plume
PYTHONPATH=. python run_fit.py --plot-only    # Fig: diurnal_fit.pdf
PYTHONPATH=. python run_fit.py --ensemble     # Fig: diurnal_fit_ensemble.pdf
PYTHONPATH=. python run_fit.py --refine       # re-polish the fit from the stored result, then redraw
```
