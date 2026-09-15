#' ReBUNKI (SPUNIT) unfolding
#'
#' R port of \code{bssunfold/core/unfold_rebunki.py}. ReBUNKI (Lacerda et al.,
#' 2018) is a modern, open reimplementation of the BUNKI Bonner-sphere
#' unfolding code; the Python version (and this R port) supports the SPUNIT
#' iterative algorithm with the default settings recommended by the ReBUNKI
#' documentation: iterations run until the relative change of the solution
#' falls below a ~1% tolerance (bounded by \code{max_iterations}), using the
#' three-point SPUNIT smoothing.
#'
#' This is a thin wrapper around \code{\link{solve_bunki}} implementing the
#' same SPUNIT scheme; the Detector-facing entry point is
#' \code{Detector$unfold_rebunki()}.
#'
#' @inheritParams solve_bunki
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_rebunki(A, b, rep(1, 3), max_iterations = 50)
solve_rebunki <- function(A, b, x0, smoothing = 0.1, max_iterations = 1000L,
                           tolerance = 0.01) {
    solve_bunki(A = A, b = b, x0 = x0, smoothing = smoothing,
                max_iterations = max_iterations, tolerance = tolerance)
}

#' Wrapper around \code{\link{solve_rebunki}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_rebunki
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_rebunki <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            smoothing = 0.1, max_iterations = 1000L,
                            tolerance = 0.01,
                            calculate_errors = FALSE,
                            noise_level = 0.01, n_montecarlo = 100L,
                            save_result = FALSE, random_state = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_rebunki,
                                         smoothing = smoothing,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "ReBUNKI",
        extra_output = list(smoothing = smoothing),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
