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
#' @param line_search_tol Numeric; kept for signature compatibility with the
#'   Python implementation, which does not use it here (the Armijo constant is
#'   hard-coded to \code{1e-4}). Default 1e-6.
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

    # Objective L = tau * D_KL(phi || phi_0_norm) + chi^2 with
    #   D_KL(phi || phi_0) = sum_i phi_i log(phi_i / phi_0_i) - phi_i + phi_0_i
    #   chi^2             = residual' S_b residual,  residual = A phi - b_work
    # so the gradient is tau * log(phi / phi_0_norm) + 2 A' S_b residual and
    # the Hessian is tau * diag(1 / phi) + 2 A' S_b A.
    objective <- function(phi) {
        p <- pmax(phi, phi_floor)
        residual <- as.numeric(A %*% p) - b_work
        chi2 <- as.numeric(residual %*% (S_b_diag * residual))
        kl <- sum(p * log(p / phi_0_norm + 1e-300) - p + phi_0_norm)
        as.numeric(tau * kl + chi2)
    }
    gradient <- function(phi) {
        p <- pmax(phi, phi_floor)
        residual <- as.numeric(A %*% p) - b_work
        kl_grad <- log(p / phi_0_norm + 1e-300)
        chi2_grad <- as.numeric(t(A) %*% (S_b_diag * residual))
        tau * kl_grad + 2 * chi2_grad
    }
    hessian <- function(phi) {
        p <- pmax(phi, phi_floor)
        kl_hess <- diag(1.0 / (p + 1e-300), nrow = n)
        tau * kl_hess + 2 * At_Sb_A
    }
    phi <- phi_0_norm
    grad_norm <- Inf
    iteration <- 0L
    for (iteration in seq_len(max_iterations)) {
        grad <- gradient(phi)
        grad_norm <- sqrt(sum(grad^2))
        if (grad_norm < tolerance) break
        Hess <- hessian(phi)
        delta <- tryCatch(as.numeric(solve(Hess, -grad)),
                          error = function(e) {
            reg <- 1e-6 * max(abs(diag(Hess)))
            if (!is.finite(reg) || reg <= 0) reg <- 1e-12
            as.numeric(solve(Hess + reg * diag(n), -grad))
        })
        base_obj <- objective(phi)
        slope <- sum(grad * delta)
        beta <- 1.0
        accepted <- FALSE
        for (j in 1:30) {
            new_phi <- pmax(phi + beta * delta, phi_floor)
            if (objective(new_phi) <= base_obj + 1e-4 * beta * slope) {
                accepted <- TRUE; break
            }
            beta <- beta * 0.5
        }
        if (!accepted) beta <- 0.01
        phi <- pmax(phi + beta * delta, phi_floor)
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
                                             random_state = NULL,
                              max_neutron_energy = NULL) {
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
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
