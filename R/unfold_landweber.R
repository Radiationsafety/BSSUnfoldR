#' Landweber iteration unfolding
#'
#' Iterative Landweber method for solving \eqn{A x = b}:
#' \eqn{x_{k+1} = P_+(x_k - tau * A^T (A x_k - b))}{x_{k+1} = P+(x_k - tau * A'(A x_k - b))}
#' with step size \eqn{tau = 1 / sigma_max^2}{tau = 1 / sigma_max^2}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum guess (length n).
#' @param max_iterations Positive integer; default 1000.
#' @param tolerance Positive numeric; convergence tolerance on the residual
#'   norm. Default 1e-6.
#' @return A list \code{list(spectrum = ..., iterations = ..., converged = ...)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_landweber(A, b, rep(0, 3), max_iterations = 50)
#' r$converged
solve_landweber <- function(A, b, x0, max_iterations = 1000L, tolerance = 1e-6) {
    v <- validate_system(A, b, x0 = x0,
                        max_iterations = max_iterations,
                        tolerance = tolerance)
    A <- v$A; b <- v$b; x0 <- v$x0
    x <- x0
    sigma_max <- norm(A, type = "2")
    if (sigma_max == 0) {
        warning("Response matrix has zero norm; returning initial guess.")
        return(list(spectrum = x, iterations = 0L, converged = FALSE))
    }
    step_size <- 1.0 / (sigma_max^2)
    AT <- t(A)
    ATb <- as.numeric(AT %*% b)
    converged <- FALSE
    iterations <- 0L
    for (i in seq_len(max_iterations)) {
        Ax <- as.numeric(A %*% x)
        grad <- as.numeric(AT %*% Ax) - ATb
        residual_norm <- sqrt(sum((Ax - b)^2))
        if (residual_norm < tolerance) {
            converged <- TRUE; iterations <- i - 1L; break
        }
        x <- pmax(x - step_size * grad, 0)
        iterations <- i
    }
    if (!converged) iterations <- max_iterations
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_landweber}} that plugs into the unified
#' \code{\link{run_unfolding}} workflow.
#'
#' @inheritParams run_unfolding
#' @param max_iterations Positive integer; default 1000.
#' @param tolerance Positive numeric; default 1e-6.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_landweber <- function(detector_names, n_energy_bins, E_MeV,
                             sensitivities, cc_icrp116, save_result_callback,
                             readings, initial_spectrum = NULL,
                             max_iterations = 1000L, tolerance = 1e-6,
                             calculate_errors = FALSE,
                             noise_level = 0.01, n_montecarlo = 100L,
                             save_result = FALSE, random_state = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_landweber,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "Landweber",
        extra_output = list(),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
