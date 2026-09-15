# Contributing to BSSUnfoldR

First of all: thank you for taking the time to contribute! 🎉

This document outlines how to propose changes, report bugs, and add new
unfolding algorithms to BSSUnfoldR.

## Table of contents

- [Reporting a bug](#reporting-a-bug)
- [Suggesting a new feature](#suggesting-a-new-feature)
- [Porting a new unfolding method from Python bssunfold](#porting-a-new-unfolding-method)
- [Pull-request workflow](#pull-request-workflow)
- [Code style](#code-style)
- [Running checks locally](#running-checks-locally)

## Reporting a bug

Open a [Bug report issue](https://github.com/Radiationsafety/BSSUnfoldR/issues/new?template=bug_report.yml)
and include a reproducible R example plus the output of `packageVersion("BSSUnfoldR")`
and `R.version.string`.

## Suggesting a new feature

Open a [Feature request issue](https://github.com/Radiationsafety/BSSUnfoldR/issues/new?template=feature_request.yml).
When suggesting a port from the upstream Python package, please mention the
specific module (e.g. `core/unfold_mlem_bs.py`) so the work can be tracked
against the original implementation.

## Porting a new unfolding method

BSSUnfoldR is designed so that every unfolding algorithm follows exactly the
same three-file pattern. To add a new method `<name>` ported from Python's
`bssunfold/core/unfold_<name>.py`:

1. **Create `R/unfold_<name>.R`** with:
   - `solve_<name>(A, b, x0, ...)` — the raw solver returning either a
     numeric spectrum or `list(spectrum, iterations, converged)`.
   - `unfold_<name>(detector_names, n_energy_bins, E_MeV, sensitivities,
     cc_icrp116, save_result_callback, readings, ...)` — the wrapper that
     calls `run_unfolding(...)`.
   - A file header comment referencing the Python source module, e.g.:
     `# Port of bssunfold/core/unfold_<name>.py (Radiationsafety/bssunfold).`
2. **Add a one-line method to the `Detector` R6 class** in `R/Detector.R`:
   ```r
   unfold_<name> = function(readings, ...) {
       unfold_<name>(self$detector_names, self$n_energy_bins,
                      self$E_MeV, self$sensitivities, self$cc_icrp116,
                      function(out) self$save_result(out),
                      readings, ...)
   },
   ```
3. **Document with Roxygen2** — every exported function needs `@export`,
   `@param`, `@return`, and a runnable `@examples` block.
4. **Regenerate NAMESPACE** by running `roxygen2::roxygenize(".")` from
   the package root.

Numeric constants (ICRP coefficients, response functions) must be added to
the upstream `bssunfold/src/bssunfold/constants.py` first; the R data file
`R/constants_data.R` is **generated** by
`scripts/extract_constants_to_r.py` and must not be hand-edited.

## Pull-request workflow

1. Fork the repository and create a feature branch from `main`.
2. Make your changes following the [Code style](#code-style) below.
3. Run `R CMD check --as-cran` locally — your PR will be blocked by CI if
   this fails.
4. Open a pull request against `main` and fill in the PR template.
5. Make sure all CI checks pass before requesting review.

## Code style

- 4-space indentation, no tabs.
- `snake_case` for function names and `PascalCase` for R6 classes.
- Every public function has Roxygen2 docs with a runnable example.
- Forbid `library()` inside package code — always use the `pkg::fun()`
  form or `@importFrom pkg fun` in `R/zzz.R`.
- Use `tryCatch()` rather than `try(..., silent = TRUE)` for error
  handling in production code paths.
- Always return `list(spectrum, iterations, converged)` from iterative
  solvers so that `run_unfolding` can extract the metadata.

## Running checks locally

```r
# From the parent directory of bssunfoldr/

# Regenerate NAMESPACE + man/*.Rd
roxygen2::roxygenize("bssunfoldr")

# Full check
devtools::check("bssunfoldr")
# Or equivalently
library(tools)
check_packages_in_dir("bssunfoldr")

# Just run examples
devtools::run_examples("bssunfoldr")
```

A full end-to-end smoke test covering all 10 unfolding methods + Monte-Carlo
+ Detector R6 lives in `scripts/smoke_test.R` at the project root (outside
the package tarball).

## Contact

For anything not covered here, contact **Konstantin Chizhov**
<kchizhov@jinr.ru>, or open a [Discussion](https://github.com/Radiationsafety/BSSUnfoldR/discussions).
