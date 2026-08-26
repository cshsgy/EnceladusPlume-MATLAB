# Diurnal-fit statistical analyses

Driver scripts for the statistical analyses behind the diurnal-profile fit in the
paper (the maximum-likelihood fit, the two fit modes, the exact-model MCMC
posterior, and the digitization sensitivity test). They build on the fit
machinery in [`../run_fit.py`](../run_fit.py) (`_neg_loglike`, `_flux_curve`,
`build_weff_interp`, `_ensemble_smooth`, `plot_corner`, `plot_ensemble`,
`plot_overlay`, `MLE_BOUNDS`, ...).

## Running

Run from the package directory (`enceladus_plume/`, the one containing
`run_fit.py`) so that `import run_fit` resolves, single-threaded per worker:

```bash
cd enceladus_plume
OMP_NUM_THREADS=1 PYTHONPATH=. python analysis/run_mcmc_exact.py
```

They cache the gas lookup table and write intermediate `.npz`/PDF products under
the system temp dir (`/tmp/...`); the lookup table is rebuilt on demand
(`GasLookupTable(..., clean=True)`), so no manual setup is required. The observed
target curve is the digitized Ingersoll+ (2020) slab-density profile referenced
by `run_fit._DATA`.

## Scripts

| Script | Purpose | Key outputs |
|---|---|---|
| `run_mcmc_exact.py` | Rigorous affine-invariant MCMC sampling the **exact** forward model (no emulator), 64 walkers, parallelized. Produces the posterior behind Figure 11. | `/tmp/diurnal_posterior_exact.npz`, `/tmp/mcmc_exact_chain.npy` |
| `two_modes.py` | Local optimization around each of the two fit optima (mode A: φ₂≈84°; mode B: φ₂≈123°), with fit-overlay figures for visual comparison. | `/tmp/modeA.npz`, `/tmp/modeB.npz`, `/tmp/mode{A,B}_fit.pdf` |
| `modeB_numbers.py` | Reports the adopted (mode B) fit numbers: ensemble vs single-crack reduced χ², and the exact-posterior depth interval. | stdout |
| `gen_figs.py` | Regenerates the ensemble figure (Fig 10) and the corner figure (Fig 11) for the adopted fit, with the truth-cross at mode B. | `/tmp/diurnal_fit_ensemble.pdf`, `/tmp/diurnal_fit_posterior.pdf` |
| `sensitivity.py` | Digitization sensitivity test: Monte-Carlo the digitized points within their per-bin scatter (N=30) and refit from both modes. Reports depth/α/φ₂/σ spread and how often each mode wins. | `/tmp/sensitivity.npz` |
| `sens_fig.py` | Renders the sensitivity figure (SI Figure S1) from `sensitivity.npz`. | `sensitivity.pdf` |

The adopted best fit is mode B (Δw≈7 mm, L≈8 km, φ₂≈123°, α≈0.81, σ_φ≈32°,
reduced χ²≈2.87); mode A (φ₂≈84°, σ_φ≈22°, χ²≈3.27) is reported as a secondary
optimum. See the paper's fit section and Supporting Information Text S4.
