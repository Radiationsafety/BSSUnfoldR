# BSSUnfoldR

[![R CMD check](https://github.com/Radiationsafety/BSSUnfoldR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Radiationsafety/BSSUnfoldR/actions/workflows/R-CMD-check.yaml)
[![Tests](https://github.com/Radiationsafety/BSSUnfoldR/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/Radiationsafety/BSSUnfoldR/actions/workflows/test-coverage.yaml)
[![pkgdown](https://github.com/Radiationsafety/BSSUnfoldR/actions/workflows/pkgdown.yaml/badge.svg)](https://radiationsafety.github.io/BSSUnfoldR/)
[![DOI](https://zenodo.org/badge/1372098194.svg)](https://doi.org/10.5281/zenodo.22790648)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.2.1-blue.svg)](https://github.com/Radiationsafety/BSSUnfoldR/releases)
[![Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

Neutron spectrum unfolding for Bonner Sphere Spectrometers — an R port of the
Python package [bssunfold](https://github.com/Radiationsafety/bssunfold).

## What's included

- **Detector R6 class** — bundles together an energy grid, per-sphere
  sensitivities, ICRP-116 conversion coefficients and a result-history list.
  One-line construction from any of the seven built-in response functions
  (`GSF`, `PTB`, `LANL`, `JINR`, `FERMILAB`, `EURADOS`, `IHEP`).
- **Sixty-plus classic unfolding algorithms** (split into eight batches):
  - **Batch 1 (v0.1.0)**: MLEM, GRAVEL, MAXED, SAND-II, BUNKI, FERDOR,
    STAY'SL, Landweber, CGLS, TSVD.
  - **Batch 2 (v0.1.1)**: OSEM, MAP-EM, BSREM, SART, Kaczmarz,
    Randomized Kaczmarz, Lanczos, Tikhonov-Legendre, ReBUNKI, Doroshenko.
  - **Batch 3 (v0.1.2)**: Bayes, Directed divergence, Express, Iterative
    refinement, BUNKI-UT, MLEM-STOP, StatReg, Tikhonov-TV, GKS, Crystal Ball.
  - **Batch 4 (v0.1.3)**: IMAXED, AMAXED, FISTA, Bayes-spline, NSDUAZ,
    MLEM-BS, NSpline, MCMC, Reconst, Ensemble / Cascade / Composite.
  - **Batch 5 (v0.1.4)**: Scipy direct, EKI, RFSP-JUL, FRUIT-like,
    AMAXED-Reg, CS, Bayesian parametric, Hybrid parametric, Hybrid GMRES,
    Binned.
  - **Batch 6 (v0.1.5)**: Parametric (FRUIT 3-component with P constraint),
    Parametric2 (BON95 4-component with grid search), EPIC (Equal Posterior
    Information Condition Tikhonov), NN-KSVD (non-negative dictionary
    learning + sparse inversion), Genetic (simulated annealing on log
    spectrum with Landweber warm-start), Mystic (differential evolution),
    QUBO (quantum-inspired binary-encoded simulated annealing), LMfit
    (Levenberg-Marquardt via optim L-BFGS-B), QPsolvers (Tikhonov-NNLS via
    augmented matrix), CVXPY (convex optimization with L1/L2, IRLS for L1).
  - **Batch 7 (v0.2.0)**: PDHG (Primal-Dual Hybrid Gradient with TV),
    Douglas-Rachford splitting, AMG (algebraic multigrid preconditioned
    Krylov), SSR (Sign-Simplicity Regression with minimum-statistic
    scoring), P-spline REML (mixed-model smoothing), SMT (exact KKT/L2
    solver), Combined (named-method pipeline), Mystic-hybrid
    (differential-evolution global + L-BFGS-B local), MAEO
    (Multi-Algorithm Evolutionary Optimization), MLEM-ODL (operator
    interface), Zfit (Poisson-likelihood Minuit analogue), plus
    `max_neutron_energy` cutoff on every method, ggplot2 visualization
    utilities, comparison/benchmarking helpers and interpretation report.
  - **Batch 8 (v0.2.1)**: GEE (Generalized Estimating Equations with
    independence / exchangeable / AR-1 working correlations and gaussian /
    poisson / gamma quasi-likelihood families; robust and naive sandwich
    covariance estimators of Liang & Zeger), Uno (constrained quadratic
    unfolding solved with Lagrange-Newton iterations — `filter_sqp` preset
    with exact-Hessian SQP and the Vanaret-Leyffer filter, `ipopt_like`
    preset with primal-dual interior point, exact or BFGS Hessian),
    Douglas-Rachford Detector-facing wrapper, and a complete rewrite of
    `solve_nnksvd` matching the Python reference (proper dictionary
    orientation, equivalent-dictionary normalization, scale-fix step).
- **Monte-Carlo uncertainty estimation** — add Gaussian noise to readings,
  run the unfolding N times, return per-bin mean / std / median / 5-95
  percentiles / all-samples matrix.
- **Dose-rate calculation** — uses ICRP-74, ICRP-116 and NRB99-2009
  conversion-coefficient datasets shipped with the package (numeric data
  is taken verbatim from the upstream `bssunfold.constants` module).
- **Unified `run_unfolding` pipeline** — every method shares the same
  validation / system-build / standardize-output / Monte-Carlo / save-result
  workflow, mirroring `bssunfold/core/_base_unfolder.py`.

Each algorithm has two entry points:
- `solve_<method>(A, b, x0, ...)` — the raw solver that operates directly
  on a response matrix `A` and measurement vector `b`;
- `unfold_<method>(detector_names, n_energy_bins, ...)` — the wrapper that
  plugs into `run_unfolding` and is called from the `Detector` class.

## Installation

```r
# Install the source tarball
install.packages("BSSUnfoldR_0.2.1.tar.gz", repos = NULL, type = "source")

# Or, if the package directory is on disk:
install.packages("path/to/bssunfoldr", repos = NULL, type = "source")
```

### System dependencies

- R >= 4.0
- From CRAN: `Matrix` (>= 1.6), `lsei` (>= 1.3), `R6` (>= 2.5)

## Quick start

```r
library(BSSUnfoldR)

# 1) Construct a Detector from the built-in PTB response function
det <- Detector$new(
    response_function = "PTB",
    detector_names    = c("0in", "3in", "5in", "8in", "12in")
)

# 2) Synthesize some readings from a known spectrum
E  <- det$E_MeV
phi_true <- 0.5 / E +
            exp(-((log10(E) + 7)^2) / 2) +      # thermal peak
            exp(-((log10(E))^2) / 0.4)          # fission bump
phi_true <- phi_true / sum(phi_true) * 1000 /
            sum(det$sensitivities[["0in"]] * phi_true)
readings <- sapply(det$detector_names, function(n)
    sum(det$sensitivities[[n]] * phi_true)
)

# 3) Unfold with each algorithm (each call also saves to det$history)
r_mlem   <- det$unfold_mlem(readings, max_iterations = 200)
r_gravel <- det$unfold_gravel(readings, max_iterations = 200)
r_maxed  <- det$unfold_maxed(readings, max_iterations = 200)
r_ferdor <- det$unfold_ferdor(readings, max_iterations = 100)
r_tsvd   <- det$unfold_tsvd(readings, method = "gcv")

# 4) Each result is a standardised list
str(r_gravel, max.level = 1)
# $ energy, $ spectrum, $ effective_readings, $ residual, $ residual_norm,
# $ method, $ doserates, $ iterations, $ converged, $ tolerance, $ regularization

# 5) Monte-Carlo uncertainty
r_mc <- det$unfold_mlem(
    readings, max_iterations = 200,
    calculate_errors = TRUE, n_montecarlo = 50, random_state = 42
)
head(r_mc$spectrum_uncert_std)
```

## Status

This port covers the core infrastructure and **sixty-plus** classic
unfolding algorithms (eight batches) from the upstream Python package,
plus Generalized Estimating Equations (GEE) and the Uno NLP solver.

R CMD check status: **OK** (0 errors, 0 warnings, 1 note — offline
future-timestamp verification, an environment artefact).
Tests: 595+ testthat tests covering every solver, the Detector R6 class,
validators, Monte-Carlo, dose rates, B-spline basis (partition of unity),
GEE / Uno / Douglas-Rachford end-to-end, and PTB integration tests.

## Citation

If you use BSSUnfoldR in your research, please cite:

```bibtex
@software{Chizhov2026BSSUnfoldR,
  author       = {Chizhov, Konstantin},
  title        = {{BSSUnfoldR}: Neutron Spectrum Unfolding for Bonner Sphere
                  Spectrometers. {R} port of the Python package 'bssunfold'
                  ({Radiationsafety}/bssunfold)},
  year         = {2026},
  version      = {0.2.1},
  url          = {https://github.com/Radiationsafety/BSSUnfoldR},
  orcid        = {0000-0003-1591-4289},
  affiliation = {Joint Institute for Nuclear Research}
}
```

See `inst/CITATION` for the R-native citation format and `CITATION.cff`
for the CFF format that GitHub renders in the "Cite this repository"
button.
