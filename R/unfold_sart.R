#' SART (simultaneous algebraic reconstruction technique) unfolding
#'
#' R port of \code{bssunfold/core/unfold_sart.py}. A relaxed weighted
#' least-squares algebraic reconstruction:
#' \deqn{x^{n+1} = x^n + alpha(n)/(A^T 1 + eps) * A^T ((b - A x^n) / (A 1 + eps))}
#' The residual is normalised by the forward-projected unit image (\code{A 1})
#' and the update by the back-projected unit image (\code{A^T 1}).
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param max_iterations Positive integer; default 50.
#' @param tolerance Positive numeric; default 1e-6.
#' @param relaxation Numeric constant (default 0.8) or a function
#'   \code{function(n)} returning the relaxation at iteration n.
#'   Default \code{NULL} = constant 0.8.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_sart(A, b, rep(0, 3), max_iterations = 50)
solve_sart <- function(A, b, x0, max_iterations = 50L, tolerance = 1e-6,
                       relaxation = NULL) {
    v <- validate_system(A, b, x0 = x0,
                        max_iterations = max_iterations,
                        tolerance = tolerance)
    A <- v$A; b <- v$b; x0 <- v$x0
    if (is.null(relaxation)) {
        relax_seq <- function(n) 0.8
    } else if (is.function(relaxation)) {
        relax_seq <- relaxation
    } else {
        relax_val <- as.numeric(relaxation)
        relax_seq <- function(n) relax_val
    }
    eps <- 1e-11
    x <- pmax(x0, 0)
    x0_first <- x[1L]
    norm_back <- colSums(A)            # A^T 1
    norm_forward <- rowSums(A)         # A 1
    converged <- FALSE
    iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        x_old <- x
        alpha <- relax_seq(it)
        Ax <- as.numeric(A %*% x)
        residual <- (b - Ax) / (norm_forward + eps)
        update <- as.numeric(t(A) %*% residual)
        x <- x + alpha * update / (norm_back + eps)
        x <- pmax(x, 0.0)
        x[1L] <- x0_first
        rel <- sqrt(sum((x - x_old)^2)) / (sqrt(sum(x_old^2)) + eps)
        if (rel < tolerance) { converged <- TRUE; break }
    }
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_sart}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_sart
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_sart <- function(detector_names, n_energy_bins, E_MeV,
                         sensitivities, cc_icrp116, save_result_callback,
                         readings, initial_spectrum = NULL,
                         max_iterations = 50L, tolerance = 1e-6,
                         relaxation = NULL,
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
        solve_func = make_solve_wrapper(solve_sart,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         relaxation = relaxation),
        solve_kwargs = list(),
        method_name = "SART",
        extra_output = list(
            relaxation = if (is.numeric(relaxation)) relaxation else NULL
        ),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
