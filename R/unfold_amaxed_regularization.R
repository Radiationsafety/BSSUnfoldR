#' AMAXED with Tikhonov regularization
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_amaxed_regularization.py}.
#' Implements AMAXED with Tikhonov-style regularization (Wong 2024). Instead
#' of fixing chi-squared and minimizing cross-entropy (as in AMAXED), this
#' method simultaneously minimizes both chi-squared and the regularizing
#' function, providing more stable convergence.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric reference (prior) spectrum (length n).
#' @param sigma_factor Numeric; relative measurement uncertainty. Default 0.1.
#' @param tau Numeric; regularization parameter. Larger values favor solutions
#'   closer to the prior. Default 1.0.
#' @param max_iterations Positive integer; default 5000.
#' @param tolerance Positive numeric; gradient convergence tolerance. Default 1e-8.
#' @param line_search_tol Numeric; Armijo c1. Default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_amaxed_regularization(A, b, rep(1, 3), max_iterations = 50)
solve_amaxed_regularization <- function(A, b, x0, sigma_factor = 0.1,
                                          tau = 1.0, max_iterations = 5000L,
                                          tolerance = 1e-8,
                                          line_search_tol = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    m <- nrow(A); n <- ncol(A)
    phi_floor <- 1e-12
    phi_0 <- pmax(x0, 1e-300)
    phi_0_sum <- sum(phi_0)
    phi_0_norm <- phi_0 / phi_0_sum
    b_work <- b / phi_0_sum
    b_safe <- pmax(b_work, 1e-300)
    sigma <- sigma_factor * b_safe
    S_b_diag <- 1.0 / sigma^2
    A_weighted <- A * matrix(S_b_diag, nrow = m, ncol = n, byrow = FALSE)
    At_Sb_A <- t(A) %*% A_weighted
    ATb <- as.numeric(t(A) %*% (b_work * S_b_diag))

    # Objective: 0.5 * chi^2 + tau * sum(phi_0_norm / phi - 1 + log(phi/phi_0))
    # Gradient: AT_Sb_A phi - ATb + tau * (-phi_0_norm / phi^2 + 1/phi)
    # Hessian: AT_Sb_A + tau * diag(2 * phi_0_norm / phi^3 - 1/phi^2)
    objective <- function(phi) {
        p <- pmax(phi, phi_floor)
        residual <- as.numeric(A %*% p) - b_work
        chi2 <- 0.5 * sum(residual^2 * S_b_diag)
        kl <- tau * sum(phi_0_norm / p - 1 + log(p / phi_0_norm))
        as.numeric(chi2 + kl)
    }
    gradient <- function(phi) {
        p <- pmax(phi, phi_floor)
        residual <- as.numeric(A %*% p) - b_work
        as.numeric(t(A) %*% (residual * S_b_diag) +
                    tau * (-phi_0_norm / p^2 + 1 / p))
    }
    hessian <- function(phi) {
        p <- pmax(phi, phi_floor)
        At_Sb_A + tau * diag(2 * phi_0_norm / p^3 - 1 / p^2, nrow = n)
    }
    phi <- phi_0_norm
    c1 <- min(max(line_search_tol, 1e-12), 0.5)
    grad_norm <- Inf
    iteration <- 0L
    for (iteration in seq_len(max_iterations)) {
        grad <- gradient(phi)
        grad_norm <- sqrt(sum(grad^2))
        if (grad_norm < tolerance) break
        Hess <- hessian(phi)
        delta <- tryCatch(as.numeric(qr.solve(Hess, -grad)),
                          error = function(e) {
            reg <- 1e-6 * max(abs(diag(Hess)))
            if (reg <= 0) reg <- 1e-12
            as.numeric(qr.solve(Hess + reg * diag(n), -grad))
        })
        slope <- sum(grad * delta)
        if (slope >= 0) { delta <- -grad; slope <- sum(grad * delta) }
        base <- objective(phi)
        beta <- 1.0
        accepted <- FALSE
        for (j in 1:30) {
            trial <- pmax(phi + beta * delta, phi_floor)
            if (objective(trial) <= base + c1 * beta * slope) {
                accepted <- TRUE; break
            }
            beta <- beta * 0.5
        }
        if (!accepted) beta <- 0.01
        phi <- pmax(phi + beta * delta, phi_floor)
    }
    if (grad_norm >= tolerance) {
        grad <- gradient(phi)
        grad_norm <- sqrt(sum(grad^2))
    }
    phi <- phi * phi_0_sum
    list(spectrum = as.numeric(phi), iterations = as.integer(iteration),
         converged = (grad_norm < tolerance))
}

#' Wrapper around \code{\link{solve_amaxed_regularization}} for the unified
#' workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_amaxed_regularization
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_amaxed_regularization <- function(detector_names, n_energy_bins, E_MeV,
                                             sensitivities, cc_icrp116,
                                             save_result_callback, readings,
                                             initial_spectrum = NULL,
                                             sigma_factor = 0.1, tau = 1.0,
                                             max_iterations = 5000L,
                                             tolerance = 1e-8,
                                             line_search_tol = 1e-6,
                                             calculate_errors = FALSE,
                                             noise_level = 0.01,
                                             n_montecarlo = 100L,
                                             save_result = FALSE,
                                             random_state = NULL) {
    x0_ref <- if (!is.null(initial_spectrum)) as.numeric(initial_spectrum)
               else rep(1.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = x0_ref,
        default_initial = rep(1.0, n_energy_bins),
        solve_func = make_solve_wrapper(solve_amaxed_regularization,
                                         sigma_factor = sigma_factor,
                                         tau = tau,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         line_search_tol = line_search_tol),
        solve_kwargs = list(),
        method_name = "AMAXED-Reg",
        extra_output = list(sigma_factor = sigma_factor, tau = tau),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
