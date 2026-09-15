# BSSUnfoldR

Neutron spectrum unfolding for Bonner Sphere Spectrometers — an R port of the
Python package [bssunfold](https://github.com/Radiationsafety/bssunfold).

## What's included

- **Detector R6 class** — bundles together an energy grid, per-sphere
  sensitivities, ICRP-116 conversion coefficients and a result-history list.
  One-line construction from any of the seven built-in response functions
  (`GSF`, `PTB`, `LANL`, `JINR`, `FERMILAB`, `EURADOS`, `IHEP`).
- **Ten classic unfolding algorithms**:
  - MLEM (Maximum Likelihood Expectation Maximization)
  - GRAVEL
  - MAXED (Maximum Entropy Deconvolution)
  - SAND-II
  - BUNKI (SPUNIT)
  - FERDOR (constrained least-squares with second-difference smoothing)
  - STAY'SL (linear Bayesian update)
  - Landweber
  - CGLS (Conjugate Gradient Least Squares)
  - TSVD (Truncated SVD, with automatic k-selection via discrepancy
    principle, L-curve, GCV, etc.)
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
install.packages("BSSUnfoldR_0.1.0.tar.gz", repos = NULL, type = "source")

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

This port covers the core infrastructure and ten classic algorithms from
the upstream Python package. The remaining ~70 specialised methods
(MLEM-BS, NSDUAZ, NSpline, MCMC, FISTA, Kaczmarz variants, parametric /
hybrid methods, etc.) can be ported on top of the same `run_unfolding`
pipeline in the same style as the ten included here.

R CMD check status: **OK** (no ERRORs, WARNINGs, or NOTEs).
