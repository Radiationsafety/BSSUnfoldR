#' AMAXED (Alternative MAXED) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_amaxed.py}.
#' Implements Alternative Maximum Entropy Deconvolution (Wong 2024) with
#' reversed cross-entropy definition compared to MAXED, using Newton's
#' method with backtracking line search on the KKT residual.
#'
#' Minimises KL divergence while adhering to the chi-squared constraint
#' using Lagrangian multipliers. The work happens on the normalized simplex
#' (phi sums to 1) and the result is rescaled at the end.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric reference (prior) spectrum (length n).
#' @param sigma_factor Numeric; relative measurement uncertainty. Default 0.1.
#' @param target_chi2 Optional numeric; target chi-squared. Default \code{NULL}
#'   = \code{m} (number of measurements).
#' @param max_iterations Positive integer; default 5000.
#' @param tolerance Positive numeric; gradient convergence tolerance. Default 1e-8.
#' @param line_search_tol Numeric; line search tolerance. Default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_amaxed(A, b, rep(1, 3), max_iterations = 50)
solve_amaxed <- function(A, b, x0, sigma_factor = 0.1, target_chi2 = NULL,
                          max_iterations = 5000L, tolerance = 1e-8,
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
    if (is.null(target_chi2)) target_chi2 <- as.numeric(m)
    Omega <- as.numeric(target_chi2)
    phi_sol <- phi_0_norm
    mu <- 1.0
    A_weighted <- A * matrix(S_b_diag, nrow = m, ncol = n, byrow = FALSE)
    a_vec <- as.numeric(b_work %*% A_weighted)  # length n
    At_Sb_A <- t(A) %*% A_weighted  # n x n

    lagrangian_grads <- function(phi, mu) {
        phi_safe <- pmax(phi, phi_floor)
        phi_sum <- sum(phi_safe)
        R_phi <- as.numeric(A %*% phi_safe)
        residual <- R_phi - b_work
        b_term <- as.numeric(phi_safe %*% At_Sb_A)  # length n
        grad_phi <- (rep(1.0, n) / phi_sum - phi_0_norm / phi_safe +
                      2 * mu * (b_term - a_vec))
        grad_mu <- sum(residual^2 * S_b_diag) - Omega
        list(grad_phi = grad_phi, grad_mu = grad_mu)
    }
    hessian <- function(phi, mu) {
        phi_safe <- pmax(phi, phi_floor)
        phi_sum <- sum(phi_safe)
        b_term <- as.numeric(phi_safe %*% At_Sb_A)
        ones_outer <- matrix(1.0, n, n)
        diag_term <- diag(phi_0_norm / phi_safe^2, n, n)
        H_phi_phi <- -ones_outer / phi_sum^2 + diag_term + 2 * mu * At_Sb_A
        H_phi_mu <- 2 * (b_term - a_vec)
        H <- matrix(0.0, n + 1L, n + 1L)
        H[1:n, 1:n] <- H_phi_phi
        H[1:n, n + 1L] <- H_phi_mu
        H[n + 1L, 1:n] <- H_phi_mu
        H
    }

    grad_norm <- Inf
    iteration <- 0L
    for (iteration in seq_len(max_iterations)) {
        g <- lagrangian_grads(phi_sol, mu)
        state_grad <- c(g$grad_phi, g$grad_mu)
        grad_norm <- sqrt(sum(state_grad^2))
        if (grad_norm < tolerance) break
        Hess <- hessian(phi_sol, mu)
        delta_state <- tryCatch(as.numeric(qr.solve(Hess, -state_grad)),
                                error = function(e) {
            reg <- 1e-6 * max(abs(diag(Hess)))
            if (reg <= 0) reg <- 1e-12
            as.numeric(qr.solve(Hess + reg * diag(n + 1L), -state_grad))
        })
        delta_phi <- delta_state[1:n]
        delta_mu <- delta_state[n + 1L]
        kkt_residual <- function(beta, phi_cur, d_phi, mu_cur, d_mu) {
            new_phi <- pmax(phi_cur + beta * d_phi, phi_floor)
            new_mu <- mu_cur + beta * d_mu
            g_new <- lagrangian_grads(new_phi, new_mu)
            sqrt(sum(g_new$grad_phi^2) + g_new$grad_mu^2)
        }
        beta <- 1.0
        accepted <- FALSE
        for (j in 1:30) {
            if (kkt_residual(beta, phi_sol, delta_phi, mu, delta_mu) <=
                (1.0 - 1e-4 * beta) * grad_norm) {
                accepted <- TRUE; break
            }
            beta <- beta * 0.5
        }
        if (!accepted) beta <- 0.01
        phi_sol <- pmax(phi_sol + beta * delta_phi, phi_floor)
        mu <- mu + beta * delta_mu
    }
    phi_sol <- phi_sol * phi_0_sum
    list(spectrum = as.numeric(phi_sol), iterations = as.integer(iteration),
         converged = (grad_norm < tolerance))
}

#' Wrapper around \code{\link{solve_amaxed}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_amaxed
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_amaxed <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            sigma_factor = 0.1, target_chi2 = NULL,
                            max_iterations = 5000L, tolerance = 1e-8,
                            line_search_tol = 1e-6,
                            calculate_errors = FALSE,
                            noise_level = 0.01, n_montecarlo = 100L,
                            save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_ref <- if (!is.null(initial_spectrum)) as.numeric(initial_spectrum)
              else rep(1.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = x0_ref,
        default_initial = rep(1.0, n_energy_bins),
        solve_func = make_solve_wrapper(solve_amaxed,
                                         sigma_factor = sigma_factor,
                                         target_chi2 = target_chi2,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         line_search_tol = line_search_tol),
        solve_kwargs = list(),
        method_name = "AMAXED",
        extra_output = list(sigma_factor = sigma_factor,
                            target_chi2 = target_chi2),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
