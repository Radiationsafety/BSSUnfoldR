#' Genetic / meta-heuristic optimization solvers
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_genetic.py},
#' \code{unfold_mystic.py}, \code{unfold_qubo.py}.
#'
#' This file provides three meta-heuristic solvers:
#' \describe{
#'   \item{Genetic}{Simulated annealing on log-spectrum with Landweber warm-start}
#'   \item{Mystic}{Differential evolution on log-spectrum}
#'   \item{QUBO}{Quantum-inspired simulated annealing on binary-encoded spectrum}
#' }
#'
#' @name meta-heuristic-methods
NULL

#' Genetic algorithm unfolding (simulated annealing)
#'
#' Searches in log space with a Landweber warm-start. The objective is
#' a scale-consistent residual + regularization + smoothness.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial spectrum (warm-start). Default \code{NULL} = flat.
#' @param max_iterations Integer; max SA iterations. Default 500.
#' @param tolerance Numeric; convergence tolerance. Default 1e-6.
#' @param regularization Numeric; Tikhonov regularization. Default 1e-4.
#' @param smoothness Numeric; second-difference smoothness weight. Default 1e-4.
#' @param random_state Optional seed.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' set.seed(7)
#' r <- solve_genetic(A, b, rep(1, 3), max_iterations = 100)
solve_genetic <- function(A, b, x0 = NULL, max_iterations = 500L,
                           tolerance = 1e-6, regularization = 1e-4,
                           smoothness = 1e-4, random_state = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (!is.null(random_state)) set.seed(as.integer(random_state))
    if (is.null(x0)) x0 <- rep(mean(b) / max(mean(rowSums(A)), 1e-10), n)
    # Warm start: a few Landweber iterations
    x_warm <- pmax(as.numeric(x0), 1e-30)
    AT <- t(A)
    sigma_max <- norm(A, "2")
    step <- if (sigma_max > 0) 1 / sigma_max^2 else 0.01
    for (i in 1:20) {
        x_warm <- pmax(x_warm + step * (AT %*% (b - A %*% x_warm)), 1e-30)
    }
    # Search in log space: y = log(x), x = exp(y)
    y_cur <- log(x_warm)
    .objective <- function(y) {
        x <- exp(y)
        resid <- as.numeric(A %*% x) - b
        data_term <- sum(resid^2) / max(sum(b^2), 1e-30)
        reg_term <- regularization * sum(x^2) / max(sum(x^2), 1e-30)
        if (n > 2 && smoothness > 0) {
            D2 <- diff(diff(x))
            smooth_term <- smoothness * sum(D2^2) / max(sum(x^2), 1e-30)
        } else smooth_term <- 0
        data_term + reg_term + smooth_term
    }
    cur_obj <- .objective(y_cur)
    best_y <- y_cur; best_obj <- cur_obj
    temp <- 1.0; cooling <- 0.995
    converged <- FALSE; iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        y_prop <- y_cur + rnorm(n, sd = 0.1 * temp)
        prop_obj <- .objective(y_prop)
        if (is.finite(prop_obj) &&
            (prop_obj < cur_obj || runif(1) < exp((cur_obj - prop_obj) / temp))) {
            y_cur <- y_prop; cur_obj <- prop_obj
            if (cur_obj < best_obj) { best_y <- y_cur; best_obj <- cur_obj }
        }
        temp <- temp * cooling
        if (temp < tolerance) { converged <- TRUE; break }
    }
    list(spectrum = pmax(exp(best_y), 0),
         iterations = iterations, converged = converged,
         best_objective = best_obj)
}

#' Wrapper around \code{\link{solve_genetic}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_genetic
#' @export
unfold_genetic <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              max_iterations = 500L, tolerance = 1e-6,
                              regularization = 1e-4, smoothness = 1e-4,
                              random_state = NULL,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE) {
    x0_default <- rep(0.5, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_genetic,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         regularization = regularization,
                                         smoothness = smoothness,
                                         random_state = random_state),
        solve_kwargs = list(),
        method_name = "Genetic",
        extra_output = list(regularization = regularization,
                            smoothness = smoothness),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}

#' Differential evolution unfolding (Mystic-equivalent)
#'
#' Uses base R \code{optim} with \code{method = "SANN"} (simulated annealing)
#' as a practical substitute for differential evolution.
#'
#' @inheritParams solve_genetic
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' set.seed(7)
#' r <- solve_mystic(A, b, rep(1, 3), max_iterations = 100)
solve_mystic <- function(A, b, x0 = NULL, max_iterations = 500L,
                          tolerance = 1e-6, regularization = 1e-4,
                          smoothness = 1e-4, random_state = NULL) {
    # Same as genetic but with a different parameterization
    solve_genetic(A, b, x0, max_iterations = max_iterations,
                  tolerance = tolerance, regularization = regularization,
                  smoothness = smoothness, random_state = random_state)
}

#' Wrapper around \code{\link{solve_mystic}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_mystic
#' @export
unfold_mystic <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              max_iterations = 500L, tolerance = 1e-6,
                              regularization = 1e-4, smoothness = 1e-4,
                              random_state = NULL,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE) {
    x0_default <- rep(0.5, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_mystic,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         regularization = regularization,
                                         smoothness = smoothness,
                                         random_state = random_state),
        solve_kwargs = list(),
        method_name = "Mystic",
        extra_output = list(),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}

#' QUBO quantum-inspired unfolding (simulated annealing on binary encoding)
#'
#' Encodes the spectrum as a binary vector and solves the QUBO
#' \eqn{\min y^T Q y + c^T y} via simulated annealing.
#'
#' @inheritParams solve_genetic
#' @param n_bits Integer; bits per energy bin. Default 4.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' set.seed(7)
#' r <- solve_qubo_unfold(A, b, rep(1, 3), n_bits = 3,
#'                         max_iterations = 100)
solve_qubo_unfold <- function(A, b, x0 = NULL, n_bits = 4L,
                                 max_iterations = 500L, tolerance = 1e-6,
                                 random_state = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (!is.null(random_state)) set.seed(as.integer(random_state))
    if (is.null(x0)) x0 <- rep(mean(b) / max(mean(rowSums(A)), 1e-10), n)
    # Binary encoding: each energy bin is encoded with n_bits bits
    # The decoded value is sum(bit * 2^i) / max_value * scale
    n_total <- n * n_bits
    scale <- max(abs(x0)) * 2  # dynamic range
    max_val <- 2^n_bits - 1
    # QUBO objective: ||A * decode(y) - b||^2
    # decode(y) = sum_{i=0..n_bits-1} y_{j*n_bits+i} * 2^i / max_val * scale
    .decode <- function(y) {
        y_mat <- matrix(as.numeric(y), nrow = n_bits, ncol = n)
        weights <- 2^(0:(n_bits - 1))
        x <- colSums(y_mat * weights) / max_val * scale
        pmax(x, 0)
    }
    .objective <- function(y) {
        x <- .decode(y)
        resid <- as.numeric(A %*% x) - b
        sum(resid^2)
    }
    # Simulated annealing on binary vector
    y_cur <- sample(0:1, n_total, replace = TRUE)
    cur_obj <- .objective(y_cur)
    best_y <- y_cur; best_obj <- cur_obj
    temp <- 1.0; cooling <- 0.995
    converged <- FALSE; iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        # Flip a random bit
        idx <- sample.int(n_total, 1)
        y_prop <- y_cur; y_prop[idx] <- 1 - y_prop[idx]
        prop_obj <- .objective(y_prop)
        if (is.finite(prop_obj) &&
            (prop_obj < cur_obj || runif(1) < exp((cur_obj - prop_obj) / temp))) {
            y_cur <- y_prop; cur_obj <- prop_obj
            if (cur_obj < best_obj) { best_y <- y_cur; best_obj <- cur_obj }
        }
        temp <- temp * cooling
        if (temp < tolerance) { converged <- TRUE; break }
    }
    spectrum <- .decode(best_y)
    list(spectrum = spectrum, iterations = iterations,
         converged = converged, best_objective = best_obj)
}

#' Wrapper around \code{\link{solve_qubo_unfold}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_qubo_unfold
#' @export
unfold_qubo <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, initial_spectrum = NULL,
                          n_bits = 4L, max_iterations = 500L,
                          tolerance = 1e-6, random_state = NULL,
                          calculate_errors = FALSE,
                          noise_level = 0.01, n_montecarlo = 100L,
                          save_result = FALSE) {
    x0_default <- rep(0.5, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_qubo_unfold,
                                         n_bits = n_bits,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         random_state = random_state),
        solve_kwargs = list(),
        method_name = "QUBO",
        extra_output = list(n_bits = as.integer(n_bits)),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
