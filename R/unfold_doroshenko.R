#' Doroshenko coordinate-update unfolding
#'
#' R port of \code{bssunfold/core/unfold_doroshenko.py}. Uses incremental
#' residual update for O(n) per-coordinate complexity instead of O(n^2) from
#' full matrix-vector products.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial guess (length n).
#' @param max_iterations Positive integer; default 1000.
#' @param tolerance Positive numeric; default 1e-6.
#' @param regularization Non-negative numeric; regularization strength to
#'   prevent division by zero. Default 0.0.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_doroshenko(A, b, rep(1, 3), max_iterations = 30)
solve_doroshenko <- function(A, b, x0, max_iterations = 1000L,
                              tolerance = 1e-6, regularization = 0.0) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    n <- ncol(A)
    x <- x0
    denom_cache <- as.numeric(colSums(A * A)) + regularization

    residual <- b - as.numeric(A %*% x)
    converged <- FALSE
    iterations <- 0L

    for (i in seq_len(max_iterations)) {
        x_old <- x
        for (j in seq_len(n)) {
            if (denom_cache[j] <= 0) next
            Aj <- A[, j]
            old_xj <- x[j]
            numerator <- sum(Aj * residual) + denom_cache[j] * old_xj
            new_xj <- max(0.0, numerator / denom_cache[j])
            delta <- new_xj - old_xj
            if (delta != 0) {
                residual <- residual - delta * Aj
                x[j] <- new_xj
            }
        }
        if (sqrt(sum((x - x_old)^2)) < tolerance) {
            converged <- TRUE
            iterations <- i
            break
        }
    }
    if (!converged) iterations <- max_iterations
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_doroshenko}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_doroshenko
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_doroshenko <- function(detector_names, n_energy_bins, E_MeV,
                                sensitivities, cc_icrp116, save_result_callback,
                                readings, initial_spectrum = NULL,
                                max_iterations = 1000L, tolerance = 1e-6,
                                regularization = 0.0,
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
        solve_func = make_solve_wrapper(solve_doroshenko,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         regularization = regularization),
        solve_kwargs = list(),
        method_name = "Doroshenko",
        extra_output = list(tolerance = tolerance,
                            regularization = regularization),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
