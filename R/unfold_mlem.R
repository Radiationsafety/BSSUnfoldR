#' MLEM (Maximum Likelihood Expectation Maximization) unfolding
#'
#' Iteratively solves \eqn{A x = b} with the multiplicative update
#' \eqn{x_{k+1} = x_k \cdot (A^T (b / (A x_k)))}{x_{k+1} = x_k * (A' * (b / (A x_k)))}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum guess (length n).
#' @param max_iterations Positive integer; default 1000.
#' @param tolerance Positive numeric; convergence tolerance on the relative
#'   spectrum change. Default 1e-6.
#' @return A list \code{list(spectrum = ..., iterations = ..., converged = ...)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_mlem(A, b, rep(1, 3), max_iterations = 50)
#' r$converged
solve_mlem <- function(A, b, x0, max_iterations = 1000L, tolerance = 1e-6) {
    v <- validate_system(A, b, x0 = x0,
                        max_iterations = max_iterations,
                        tolerance = tolerance)
    A <- v$A; b <- v$b; x0 <- v$x0
    AT <- t(A)
    x <- pmax(x0, 1e-10)
    converged <- FALSE
    iterations <- 0L
    for (i in seq_len(max_iterations)) {
        Ax <- A %*% x
        Ax <- pmax(as.numeric(Ax), 1e-10)
        ratio <- b / Ax
        correction <- as.numeric(AT %*% ratio)
        x_new <- x * correction
        diff <- sqrt(sum((x_new - x)^2)) / (sqrt(sum(x^2)) + 1e-10)
        x <- pmax(x_new, 0)
        iterations <- i
        if (diff < tolerance) { converged <- TRUE; break }
    }
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_mlem}} that plugs into the unified
#' \code{\link{run_unfolding}} workflow.
#'
#' @inheritParams run_unfolding
#' @param max_iterations Positive integer; default 1000.
#' @param tolerance Positive numeric; default 1e-6.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' det_names <- c("d1", "d2", "d3")
#' sens <- lapply(setNames(det_names, det_names), function(n)
#'     A[which(det_names == n), ])
#' E <- c(1e-9, 1e-6, 1e-3)
#' res <- unfold_mlem(
#'     detector_names = det_names, n_energy_bins = 3L, E_MeV = E,
#'     sensitivities = sens, cc_icrp116 = NULL, save_result_callback = NULL,
#'     readings = c(d1 = 1, d2 = 0.6, d3 = 0.4), max_iterations = 50
#' )
#' res$method
unfold_mlem <- function(detector_names, n_energy_bins, E_MeV,
                        sensitivities, cc_icrp116, save_result_callback,
                        readings, initial_spectrum = NULL,
                        max_iterations = 1000L, tolerance = 1e-6,
                        calculate_errors = FALSE,
                        noise_level = 0.01, n_montecarlo = 100L,
                        save_result = FALSE, random_state = NULL) {
    x0_default <- rep(0.5, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_mlem,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "MLEM",
        extra_output = list(),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
