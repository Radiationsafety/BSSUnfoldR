#' MLEM with J-factor early stopping (MLEM-STOP)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_mlem_stop.py}. Based
#' on:
#' Montgomery et al., "A novel MLEM stopping criterion for unfolding neutron
#' fluence spectra in radiation therapy", Nucl. Instrum. Meth. A 957 (2020)
#' 163400. \url{https://doi.org/10.1016/j.nima.2020.163400}
#'
#' The J-factor is \eqn{J = \sum (b_i - (A x)_i)^2 / \sum (A x)_i}. Iterations
#' stop as soon as \eqn{J \le J_{threshold}}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial guess (length n).
#' @param max_iterations Positive integer; default 15000.
#' @param cps_crossover Numeric; crossover CPS value for automatic J threshold.
#'   Default 30000.
#' @param j_threshold Optional numeric; if \code{NULL} computed as
#'   \code{mean(b) / cps_crossover}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_mlem_stop(A, b, rep(1, 3), max_iterations = 100)
solve_mlem_stop <- function(A, b, x0, max_iterations = 15000L,
                              cps_crossover = 30000.0, j_threshold = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    if (is.null(j_threshold)) j_threshold <- mean(b) / cps_crossover
    x <- pmax(x0, 1e-10)
    AT <- t(A)
    converged <- FALSE
    iterations <- 0L
    for (i in seq_len(max_iterations)) {
        iterations <- i
        Ax <- as.numeric(A %*% x)
        j_factor <- sum((b - Ax)^2) / max(sum(Ax), .Machine$double.eps)
        if (j_factor <= j_threshold) {
            return(list(spectrum = as.numeric(x), iterations = iterations,
                        converged = TRUE))
        }
        Ax <- pmax(Ax, 1e-10)
        ratio <- b / Ax
        correction <- as.numeric(AT %*% ratio)
        x <- pmax(x * correction, 0)
    }
    j_final <- sum((b - as.numeric(A %*% x))^2) /
               max(sum(as.numeric(A %*% x)), .Machine$double.eps)
    converged <- j_final <= j_threshold
    list(spectrum = as.numeric(x), iterations = as.integer(max_iterations),
         converged = converged)
}

#' Wrapper around \code{\link{solve_mlem_stop}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_mlem_stop
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_mlem_stop <- function(detector_names, n_energy_bins, E_MeV,
                               sensitivities, cc_icrp116, save_result_callback,
                               readings, initial_spectrum = NULL,
                               max_iterations = 15000L,
                               cps_crossover = 30000.0,
                               j_threshold = NULL,
                               calculate_errors = FALSE,
                               noise_level = 0.01, n_montecarlo = 100L,
                               save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.5, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_mlem_stop,
                                         max_iterations = max_iterations,
                                         cps_crossover = cps_crossover,
                                         j_threshold = j_threshold),
        solve_kwargs = list(),
        method_name = "MLEM-STOP",
        extra_output = list(),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
