#' BUNKI-UT (BON31G) unfolding
#'
#' R port of \code{bssunfold/core/unfold_bunkiut.py}. BUNKI-UT is the
#' University of Texas modernisation of the BUNKI code (K. A. Miller et al.)
#' and implements the BON31G variant of the SPUNIT iterative unfolding
#' algorithm.
#'
#' BON31G operates on the lethargy-weighted response matrix \code{aleth}
#' (built by the Detector class) transformed by the initial spectrum
#' \code{alethnew = aleth * x0}, with the starting spectrum \code{spl = 1}:
#' \deqn{bk_{jm} = \sum_i alethnew_{ij} * alethnew_{im}}
#' \deqn{vect_j = \sum_i alethnew_{ij} * b_i}
#' \deqn{ax_j = \sum_m spl_m * bk_{jm}}
#' \deqn{spll_j = spl_j * vect_j / ax_j}
#' Then \code{spl <- 3-point-smoothed spll} (bins 0, 1 kept verbatim), and
#' the final spectrum is the inverse transform \code{x = spl * x0}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param smoothing Numeric; 3-point smoothing factor. Default 0.05.
#' @param max_iterations Positive integer; default 1000.
#' @param tolerance Positive numeric; default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_bunkiut(A, b, rep(1, 3), max_iterations = 50)
solve_bunkiut <- function(A, b, x0, smoothing = 0.05,
                            max_iterations = 1000L, tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    n <- ncol(A)
    if (any(b < 0)) stop("BUNKI-UT requires strictly positive measurements")
    if (any(b == 0)) {
        keep <- b > 0
        A <- A[keep, , drop = FALSE]
        b <- b[keep]
        if (length(b) == 0L) {
            stop("BUNKI-UT requires strictly positive measurements")
        }
    }
    x0_safe <- pmax(x0, 0.0)
    aleth <- A * matrix(x0_safe, nrow = nrow(A), ncol = n, byrow = TRUE)
    spl <- rep(1.0, n)
    bk <- crossprod(aleth)            # n x n
    vect <- as.numeric(t(aleth) %*% b)  # n
    denom_s <- 1.0 + 2.0 * smoothing
    converged <- FALSE; iterations <- 0L

    for (k in seq_len(max_iterations)) {
        iterations <- k
        spl_old <- spl
        # ax[j] = sum_m spl[m] * bk[j, m]
        ax <- pmax(as.numeric(spl %*% t(bk)), 1e-37)
        spll <- spl * vect / ax
        spll <- ifelse(spll < 1e-37, 0.0, spll)
        new_spl <- spll
        if (n > 2L) {
            for (j in 2:(n - 1L)) {
                hi <- if (j + 1L <= n) spll[j + 1L] else 0.0
                new_spl[j] <- (spll[j - 1L] * smoothing + spll[j] +
                                hi * smoothing) / denom_s
            }
        }
        new_spl[1L] <- spll[1L]
        if (n > 1L) new_spl[2L] <- spll[2L]
        spl <- new_spl
        rel <- sqrt(sum((spl - spl_old)^2)) / (sqrt(sum(spl_old^2)) + 1e-12)
        if (rel < tolerance) { converged <- TRUE; break }
    }
    spectrum <- spl * x0_safe
    list(spectrum = as.numeric(spectrum), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_bunkiut}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_bunkiut
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_bunkiut <- function(detector_names, n_energy_bins, E_MeV,
                             sensitivities, cc_icrp116, save_result_callback,
                             readings, initial_spectrum = NULL,
                             smoothing = 0.05, max_iterations = 1000L,
                             tolerance = 1e-6,
                             calculate_errors = FALSE,
                             noise_level = 0.01, n_montecarlo = 100L,
                             save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    x0_default[1L] <- 0.0
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_bunkiut,
                                         smoothing = smoothing,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "BUNKI-UT",
        extra_output = list(smoothing = smoothing),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
