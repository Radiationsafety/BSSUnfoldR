#' LMfit / QP / CVXPY-equivalent unfolding methods
#'
#' R ports of \code{bssunfold/src/bssunfold/core/unfold_lmfit.py},
#' \code{unfold_qpsolvers.py}, \code{unfold_cvxpy.py}.
#'
#' @name optimization-methods
NULL

#' Levenberg-Marquardt unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_lmfit.py}.
#' Minimises \eqn{\|A x - b\|^2 + \alpha \|L x\|^2} using
#' \code{\link[stats]{optim}} with \code{L-BFGS-B} method and box constraints.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial guess. Default \code{NULL} = flat.
#' @param alpha Numeric; regularization parameter. Default 0.01.
#' @param smoothness_order Integer; 0 (identity), 1, or 2. Default 0.
#' @param max_iterations Integer; default 1000.
#' @param tolerance Numeric; default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_lmfit(A, b, rep(1, 3), alpha = 0.01, max_iterations = 50)
solve_lmfit <- function(A, b, x0 = NULL, alpha = 0.01,
                          smoothness_order = 0L, max_iterations = 1000L,
                          tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    L <- make_regularization_operator(n, smoothness_order,
                                       identity_for_zero = FALSE)
    .objective <- function(x) {
        x <- pmax(x, 0)
        resid <- as.numeric(A %*% x) - b
        data <- sum(resid^2)
        reg <- if (!is.null(L)) alpha * sum(as.numeric(L %*% x)^2) else 0
        data + reg
    }
    .gradient <- function(x) {
        x <- pmax(x, 0)
        resid <- as.numeric(A %*% x) - b
        grad <- as.numeric(t(A) %*% resid)
        if (!is.null(L)) grad <- grad + alpha * as.numeric(t(L) %*% (L %*% x))
        grad
    }
    if (is.null(x0)) x0 <- rep(mean(b) / max(mean(rowSums(A)), 1e-10), n)
    result <- stats::optim(x0, .objective, .gradient, method = "L-BFGS-B",
                            lower = rep(0, n), upper = rep(Inf, n),
                            control = list(maxit = max_iterations))
    list(spectrum = pmax(as.numeric(result$par), 0),
         iterations = as.integer(result$counts[1]),
         converged = (result$convergence == 0))
}

#' Wrapper around \code{\link{solve_lmfit}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_lmfit
#' @export
unfold_lmfit <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            alpha = 0.01, smoothness_order = 0L,
                            max_iterations = 1000L, tolerance = 1e-6,
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
        solve_func = make_solve_wrapper(solve_lmfit,
                                         alpha = alpha,
                                         smoothness_order = smoothness_order,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "LMfit",
        extra_output = list(alpha = alpha,
                            smoothness_order = as.integer(smoothness_order)),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}

#' QP solver unfolding (Tikhonov-regularized NNLS)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_qpsolvers.py}.
#' Solves \eqn{\min \|A x - b\|^2 + \alpha \|L x\|^2} with \eqn{x \geq 0}
#' via the augmented-matrix NNLS approach.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused.
#' @param alpha Numeric; regularization. Default 0.01.
#' @param smoothness_order Integer; 0, 1, or 2. Default 0.
#' @return A list \code{list(spectrum, iterations = 1, converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_qpsolvers(A, b, NULL, alpha = 0.01)
solve_qpsolvers <- function(A, b, x0 = NULL, alpha = 0.01,
                              smoothness_order = 0L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    L <- make_regularization_operator(n, smoothness_order,
                                       identity_for_zero = FALSE)
    if (!is.null(L) && alpha > 0) {
        Aw <- rbind(A, sqrt(alpha) * L)
        bw <- c(b, rep(0, nrow(L)))
    } else {
        Aw <- A; bw <- b
    }
    x <- tryCatch(as.numeric(lsei::nnls(Aw, bw)$x),
                  error = function(e) as.numeric(qr.solve(A, b)))
    list(spectrum = pmax(x, 0), iterations = 1L, converged = TRUE)
}

#' Wrapper around \code{\link{solve_qpsolvers}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_qpsolvers
#' @export
unfold_qpsolvers <- function(detector_names, n_energy_bins, E_MeV,
                                 sensitivities, cc_icrp116, save_result_callback,
                                 readings, initial_spectrum = NULL,
                                 alpha = 0.01, smoothness_order = 0L,
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
        solve_func = make_solve_wrapper(solve_qpsolvers,
                                         alpha = alpha,
                                         smoothness_order = smoothness_order),
        solve_kwargs = list(),
        method_name = "QPsolvers",
        extra_output = list(alpha = alpha,
                            smoothness_order = as.integer(smoothness_order)),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}

#' CVXPY-equivalent convex optimization unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_cvxpy.py}.
#' Solves \eqn{\min \|A x - b\|_2 + \alpha \|x\|_p} with \eqn{x \geq 0}
#' using NNLS (for L2 norm) or iterative reweighting (for L1 norm).
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused.
#' @param alpha Numeric; regularization. Default 0.01.
#' @param norm Integer; 1 (L1) or 2 (L2). Default 2.
#' @return A list \code{list(spectrum, iterations = 1, converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_cvxpy(A, b, NULL, alpha = 0.01, norm = 2)
solve_cvxpy <- function(A, b, x0 = NULL, alpha = 0.01, norm = 2L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (norm == 2L) {
        # L2: augmented NNLS
        Aw <- rbind(A, sqrt(alpha) * diag(n))
        bw <- c(b, rep(0, n))
        x <- tryCatch(as.numeric(lsei::nnls(Aw, bw)$x),
                      error = function(e) as.numeric(qr.solve(A, b)))
    } else {
        # L1: iterative reweighted least squares (IRLS)
        x <- tryCatch(as.numeric(lsei::nnls(A, b)$x),
                      error = function(e) as.numeric(qr.solve(A, b)))
        for (iter in 1:10) {
            w <- 1 / pmax(abs(x), 1e-8)
            Aw <- rbind(A, sqrt(alpha * w) * diag(n))
            bw <- c(b, rep(0, n))
            x <- tryCatch(as.numeric(lsei::nnls(Aw, bw)$x),
                          error = function(e) x)
        }
    }
    list(spectrum = pmax(x, 0), iterations = 1L, converged = TRUE)
}

#' Wrapper around \code{\link{solve_cvxpy}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_cvxpy
#' @export
unfold_cvxpy <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            alpha = 0.01, norm = 2L,
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
        solve_func = make_solve_wrapper(solve_cvxpy,
                                         alpha = alpha, norm = norm),
        solve_kwargs = list(),
        method_name = "CVXPY",
        extra_output = list(alpha = alpha, norm = as.integer(norm)),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
