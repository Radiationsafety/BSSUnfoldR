#' Kaczmarz (ART) unfolding
#'
#' R port of \code{bssunfold/core/unfold_kaczmarz.py}. Cyclic row-action
#' algorithm that projects the current solution onto each hyperplane defined
#' by a single detector row.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial guess (length n).
#' @param max_iterations Positive integer; default 1000.
#' @param omega Numeric relaxation in (0, 2]; default 1.0.
#' @param tolerance Positive numeric; default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_kaczmarz(A, b, rep(0, 3), max_iterations = 30)
solve_kaczmarz <- function(A, b, x0, max_iterations = 1000L, omega = 1.0,
                           tolerance = 1e-6) {
    v <- validate_system(A, b, x0 = x0,
                        max_iterations = max_iterations,
                        tolerance = tolerance)
    A <- v$A; b <- v$b; x0 <- v$x0
    m <- nrow(A)
    x <- x0
    if (omega <= 0 || omega > 2) {
        warning("omega=", omega, " outside recommended range (0, 2]")
    }
    row_norms_sq <- as.numeric(rowSums(A * A))
    converged <- FALSE
    iterations <- 0L
    x_old <- x
    for (k in seq_len(max_iterations)) {
        i <- ((k - 1L) %% m) + 1L
        if (row_norms_sq[i] > 0) {
            Ai <- A[i, ]
            update <- (b[i] - sum(Ai * x)) / row_norms_sq[i]
            x <- pmax(x + omega * update * Ai, 0.0)
        }
        if (k %% m == 0L) {
            if (sqrt(sum((x - x_old)^2)) < tolerance) {
                converged <- TRUE
                iterations <- k
                break
            }
            x_old <- x
        }
    }
    if (!converged) iterations <- max_iterations
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_kaczmarz}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_kaczmarz
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_kaczmarz <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              max_iterations = 1000L, omega = 1.0,
                              tolerance = 1e-6,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_kaczmarz,
                                         max_iterations = max_iterations,
                                         omega = omega,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "Kaczmarz",
        extra_output = list(tolerance = tolerance, omega = omega),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
