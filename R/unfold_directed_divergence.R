#' Directed-divergence unfolding
#'
#' R port of \code{bssunfold/core/unfold_directed_divergence.py}. Minimizes
#' the Poisson/I-divergence data term with multiplicative updates. Optional
#' first- or second-order Tikhonov smoothing is applied as a non-negative
#' proximal step after each update.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum guess (length n).
#' @param max_iterations Positive integer; default 200.
#' @param tol_chi2 Numeric; chi-squared stopping threshold. Default 1.0.
#' @param tol_rel Numeric; relative-change stopping threshold. Default 1e-6.
#' @param relative_uncertainty Numeric; relative measurement uncertainty used
#'   to derive sigma when \code{sigma} is not given. Default 0.05.
#' @param sigma Optional explicit per-detector uncertainties.
#' @param smoothness_order Integer; 0, 1, or 2. Default 0.
#' @param smoothness_weight Numeric non-negative penalty weight. Default 0.0.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_directed_divergence(A, b, rep(1, 3), max_iterations = 50)
solve_directed_divergence <- function(A, b, x0, max_iterations = 200L,
                                        tol_chi2 = 1.0, tol_rel = 1e-6,
                                        relative_uncertainty = 0.05,
                                        sigma = NULL,
                                        smoothness_order = 0L,
                                        smoothness_weight = 0.0) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    x <- pmax(x0, 1e-30)
    if (any(b < 0)) stop("Directed divergence requires non-negative measurements")
    if (!(smoothness_order %in% c(0L, 1L, 2L))) {
        stop("smoothness_order must be 0, 1, or 2")
    }
    if (smoothness_weight < 0) stop("smoothness_weight must be non-negative")

    if (is.null(sigma)) {
        sigma <- pmax(relative_uncertainty * pmax(b, 1e-30), 1e-30)
    } else {
        sigma <- as.numeric(sigma)
        if (length(sigma) != length(b) || any(sigma <= 0)) {
            stop("sigma must be positive and match b")
        }
    }
    weights <- 1.0 / sigma^2
    denominator <- pmax(colSums(A), 1e-30)
    penalty <- NULL
    if (smoothness_order != 0L && smoothness_weight > 0) {
        L <- as.matrix(create_derivative_matrix(length(x), smoothness_order))
        penalty <- smoothness_weight * crossprod(L)
    }
    converged <- FALSE; iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        predicted <- pmax(as.numeric(A %*% x), 1e-30)
        chi2 <- mean(((predicted - b)^2) * weights)
        if (chi2 <= tol_chi2) { converged <- TRUE; break }
        update <- as.numeric(t(A) %*% (b / predicted))
        new_x <- x * update / denominator
        if (!is.null(penalty)) {
            new_x <- as.numeric(solve(diag(length(x)) + penalty, new_x))
        }
        new_x <- pmax(new_x, 1e-30)
        rel <- max(abs(new_x - x) / pmax(x, 1e-30))
        x <- new_x
        if (rel <= tol_rel) { converged <- TRUE; break }
    }
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_directed_divergence}} for the unified
#' workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_directed_divergence
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_directed_divergence <- function(detector_names, n_energy_bins, E_MeV,
                                         sensitivities, cc_icrp116,
                                         save_result_callback, readings,
                                         initial_spectrum = NULL,
                                         max_iterations = 200L,
                                         tol_chi2 = 1.0, tol_rel = 1e-6,
                                         relative_uncertainty = 0.05,
                                         smoothness_order = 0L,
                                         smoothness_weight = 0.0,
                                         calculate_errors = FALSE,
                                         noise_level = 0.01,
                                         n_montecarlo = 100L,
                                         save_result = FALSE,
                                         random_state = NULL) {
    steps <- compute_log_steps(E_MeV)
    default <- rep(1.0, n_energy_bins)
    A <- do.call(rbind, lapply(detector_names, function(n) as.numeric(sensitivities[[n]])))
    b <- as.numeric(readings[detector_names])
    scale <- mean(b) / max(mean(rowSums(A)), 1e-30)
    default <- default * scale

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = default,
        solve_func = make_solve_wrapper(solve_directed_divergence,
                                         max_iterations = max_iterations,
                                         tol_chi2 = tol_chi2,
                                         tol_rel = tol_rel,
                                         relative_uncertainty = relative_uncertainty,
                                         smoothness_order = smoothness_order,
                                         smoothness_weight = smoothness_weight),
        solve_kwargs = list(),
        method_name = "Directed divergence",
        extra_output = list(
            tol_chi2 = tol_chi2, tol_rel = tol_rel,
            smoothness_order = as.integer(smoothness_order),
            smoothness_weight = smoothness_weight,
            log_steps = steps
        ),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
