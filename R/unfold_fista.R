#' FISTA unfolding (Fast Iterative Shrinkage-Thresholding Algorithm)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_fista.py}.
#' Accelerated proximal gradient method achieving O(1/k^2) convergence.
#' Supports L1, TV, Tikhonov regularization, and box constraints.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial spectrum (length n). Default \code{NULL} = scaled ones.
#' @param max_iterations Positive integer; default 500.
#' @param tolerance Positive numeric; default 1e-8.
#' @param regularization Numeric Tikhonov regularization parameter. Default 0.
#' @param l1_penalty Numeric L1 penalty for sparsity. Default 0.
#' @param tv_penalty Numeric total-variation penalty. Default 0.
#' @param nonnegativity Logical; enforce non-negativity. Default \code{TRUE}.
#' @param x_min,x_max Numeric box bounds. Default \code{0} and \code{Inf}.
#' @param noise_level Optional numeric; relative noise level for discrepancy
#'   principle stopping. Default \code{NULL}.
#' @param eta Numeric; safety factor for discrepancy principle. Default 1.01.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_fista(A, b, rep(1, 3), max_iterations = 50,
#'                  regularization = 0.01)
solve_fista <- function(A, b, x0 = NULL, max_iterations = 500L,
                         tolerance = 1e-8, regularization = 0.0,
                         l1_penalty = 0.0, tv_penalty = 0.0,
                         nonnegativity = TRUE,
                         x_min = 0.0, x_max = Inf,
                         noise_level = NULL, eta = 1.01) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (is.null(x0)) {
        x <- rep(mean(b) / max(mean(A), 1e-10), n)
    } else {
        x <- as.numeric(x0)
    }
    if (nonnegativity) x <- pmax(x, 0)
    y <- x
    t <- 1.0
    L <- tryCatch(norm(A, type = "2")^2, error = function(e) {
        # Power iteration fallback
        v <- rnorm(n); v <- v / sqrt(sum(v^2))
        for (i in 1:20) {
            u <- as.numeric(A %*% v)
            v_new <- as.numeric(t(A) %*% u)
            v <- v_new / sqrt(sum(v_new^2))
        }
        sum((as.numeric(A %*% v))^2) / sum(v^2)
    })
    L <- max(L, 1e-10)
    step_size <- 1.0 / L
    D <- NULL
    if (tv_penalty > 0 && n > 1L) {
        D <- as.matrix(create_derivative_matrix(n, 1L))
    }
    discrepancy_threshold <- if (!is.null(noise_level)) {
        eta * noise_level * sqrt(sum(b^2))
    } else NULL

    converged <- FALSE; iterations <- 0L
    for (k in seq_len(max_iterations)) {
        iterations <- k
        x_old <- x
        residual <- as.numeric(A %*% y) - b
        gradient <- as.numeric(t(A) %*% residual)
        if (regularization > 0) gradient <- gradient + regularization * y
        if (tv_penalty > 0 && !is.null(D)) {
            gradient <- gradient + tv_penalty * as.numeric(t(D) %*% (D %*% y))
        }
        x_temp <- y - step_size * gradient
        if (l1_penalty > 0) {
            # Soft thresholding
            x_temp <- sign(x_temp) * pmax(abs(x_temp) - step_size * l1_penalty, 0)
        }
        if (nonnegativity) x_temp <- pmax(x_temp, x_min)
        if (is.finite(x_max)) x_temp <- pmin(x_temp, x_max)
        x <- x_temp
        t_new <- (1.0 + sqrt(1.0 + 4.0 * t * t)) / 2.0
        y <- x + ((t - 1.0) / t_new) * (x - x_old)
        t <- t_new
        rel_change <- sqrt(sum((x - x_old)^2)) / max(sqrt(sum(x_old^2)), 1e-10)
        if (rel_change < tolerance) { converged <- TRUE; break }
        if (!is.null(discrepancy_threshold)) {
            current_res <- sqrt(sum((as.numeric(A %*% x) - b)^2))
            if (current_res <= discrepancy_threshold) {
                converged <- TRUE; break
            }
        }
    }
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_fista}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_fista
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_fista <- function(detector_names, n_energy_bins, E_MeV,
                           sensitivities, cc_icrp116, save_result_callback,
                           readings, initial_spectrum = NULL,
                           max_iterations = 500L, tolerance = 1e-8,
                           regularization = 0.0, l1_penalty = 0.0,
                           tv_penalty = 0.0, nonnegativity = TRUE,
                           x_min = 0.0, x_max = Inf,
                           noise_level = NULL, eta = 1.01,
                           calculate_errors = FALSE,
                           n_montecarlo = 100L,
                           save_result = FALSE, random_state = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_fista,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         regularization = regularization,
                                         l1_penalty = l1_penalty,
                                         tv_penalty = tv_penalty,
                                         nonnegativity = nonnegativity,
                                         x_min = x_min, x_max = x_max,
                                         noise_level = noise_level, eta = eta),
        solve_kwargs = list(),
        method_name = "FISTA",
        extra_output = list(regularization = regularization,
                            l1_penalty = l1_penalty,
                            tv_penalty = tv_penalty),
        calculate_errors = calculate_errors,
        noise_level = if (is.null(noise_level)) 0.01 else noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
