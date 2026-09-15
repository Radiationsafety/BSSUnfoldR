#' Monte-Carlo uncertainty estimation for unfolding methods
#'
#' Adds Gaussian noise to the detector readings, runs the supplied unfolding
#' function for each noisy sample, and computes summary statistics of the
#' resulting spectra. This is the R port of
#' \code{core/_montecarlo.py:monte_carlo_uncertainty}.
#'
#' @param func Function with signature \code{func(noisy_readings, ...)} returning
#'   a numeric spectrum vector.
#' @param readings Named numeric vector of detector readings.
#' @param noise_level Numeric; relative Gaussian noise std as a fraction of value.
#' @param n_samples Integer; number of Monte-Carlo samples.
#' @param n_energy_bins Integer; length of the unfolded spectrum.
#' @param random_state Optional integer seed for reproducibility.
#' @param ... Additional keyword arguments forwarded to \code{func}.
#' @return A list with the following components:
#' \describe{
#'   \item{spectrum_uncert_mean}{mean spectrum (length n_energy_bins).}
#'   \item{spectrum_uncert_std}{per-bin standard deviation.}
#'   \item{spectrum_uncert_min}{per-bin minimum.}
#'   \item{spectrum_uncert_max}{per-bin maximum.}
#'   \item{spectrum_uncert_median}{per-bin median.}
#'   \item{spectrum_uncert_percentile_5}{5th percentile per bin.}
#'   \item{spectrum_uncert_percentile_95}{95th percentile per bin.}
#'   \item{spectrum_uncert_all}{n_samples x n_energy_bins matrix of all samples.}
#' }
#' @export
#' @examples
#' # Tiny synthetic test: identity-like response matrix
#' A <- matrix(c(1, 0.5, 0.2, 0.1), nrow = 2)
#' b <- c(1, 0.6)
#' readings <- c(d1 = 1, d2 = 0.6)
#' # Trivial "solver": return b projected onto n=2 bins
#' toy_solver <- function(rds, ...) as.numeric(rds)
#' mc <- monte_carlo_uncertainty(toy_solver, readings, 0.05, 10, 2, 42)
#' names(mc)
monte_carlo_uncertainty <- function(func, readings, noise_level, n_samples,
                                    n_energy_bins, random_state = NULL, ...) {
    if (!is.numeric(readings) || length(readings) == 0L) {
        stop("'readings' must be a non-empty numeric vector")
    }
    if (!is.function(func)) stop("'func' must be a function")
    if (!is.numeric(noise_level) || length(noise_level) != 1L ||
        noise_level <= 0 || noise_level > 1) {
        stop("'noise_level' must be a number in (0, 1]")
    }
    if (!is.numeric(n_samples) || length(n_samples) != 1L ||
        n_samples < 0 || n_samples != as.integer(n_samples)) {
        stop("'n_samples' must be a non-negative integer")
    }
    n_samples <- as.integer(n_samples)
    if (!is.numeric(n_energy_bins) || length(n_energy_bins) != 1L ||
        n_energy_bins <= 0 || n_energy_bins != as.integer(n_energy_bins)) {
        stop("'n_energy_bins' must be a positive integer")
    }
    n_energy_bins <- as.integer(n_energy_bins)

    if (!is.null(random_state)) {
        if (!is.numeric(random_state) || length(random_state) != 1L ||
            random_state < 0) {
            stop("'random_state' must be a non-negative integer or NULL")
        }
        set.seed(as.integer(random_state))
    }

    keys <- names(readings)
    values <- as.numeric(readings)
    n_readings <- length(keys)

    # (n_samples x n_readings) Gaussian noise factors
    noise_factors <- matrix(1.0 + rnorm(n_samples * n_readings,
                                        mean = 0, sd = noise_level),
                            nrow = n_samples, ncol = n_readings)

    spectra <- matrix(0.0, nrow = n_samples, ncol = n_energy_bins)
    for (i in seq_len(n_samples)) {
        noisy_values <- values * noise_factors[i, ]
        names(noisy_values) <- keys
        spectrum <- tryCatch(func(noisy_values, ...),
                             error = function(e) rep(NA_real_, n_energy_bins))
        spectra[i, ] <- as.numeric(spectrum)
    }
    # Replace any all-NA rows (failed solves) with the row mean so summary stats
    # don't blow up -- this matches the spirit of the Python implementation.
    ok <- rowSums(!is.na(spectra)) > 0
    if (!all(ok)) {
        if (!any(ok)) {
            warning("All Monte-Carlo samples failed; returning zero spectra.")
            spectra <- matrix(0.0, nrow = n_samples, ncol = n_energy_bins)
        } else {
            spectra[!ok, ] <- rep(colMeans(spectra[ok, , drop = FALSE],
                                           na.rm = TRUE),
                                  each = sum(!ok))
        }
    }

    list(
        spectrum_uncert_mean           = colMeans(spectra),
        spectrum_uncert_std            = matrixStats_col_sd(spectra),
        spectrum_uncert_min            = apply(spectra, 2L, min),
        spectrum_uncert_max           = apply(spectra, 2L, max),
        spectrum_uncert_median         = apply(spectra, 2L, median),
        spectrum_uncert_percentile_5  = apply(spectra, 2L, quantile, probs = 0.05,
                                               names = FALSE),
        spectrum_uncert_percentile_95 = apply(spectra, 2L, quantile, probs = 0.95,
                                               names = FALSE),
        spectrum_uncert_all           = spectra
    )
}

# Column-wise sd that does not require matrixStats; identical to
# apply(X, 2, sd) but keeps NA handling explicit.
matrixStats_col_sd <- function(x) {
    n <- nrow(x)
    if (n <= 1L) return(rep(NA_real_, ncol(x)))
    sx <- colSums(x)
    sx2 <- colSums(x * x)
    sqrt((sx2 - sx * sx / n) / (n - 1L))
}
