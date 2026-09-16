#' IMAXED (Improved MAXED) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_imaxed.py}.
#' Implements Improved Maximum Entropy Deconvolution (Wong 2024) using
#' Newton's method with Armijo backtracking line search instead of L-BFGS-B.
#'
#' Minimises
#' \deqn{f(\phi) = 0.5 (A \phi - b)^T S_b (A \phi - b) + \sum_i \phi_i \ln(\phi_i / \phi_{0,i}) - \phi_i + \phi_{0,i}}
#' where \eqn{S_b = diag(1/\sigma^2)} and \eqn{\sigma = \sigma_factor * b}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric reference (prior) spectrum (length n).
#' @param sigma_factor Numeric; relative measurement uncertainty. Default 0.1.
#' @param max_iterations Positive integer; default 5000.
#' @param tolerance Positive numeric; gradient convergence tolerance. Default 1e-8.
#' @param line_search_tol Numeric; Armijo c1 constant. Default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_imaxed(A, b, rep(1, 3), max_iterations = 50)
solve_imaxed <- function(A, b, x0, sigma_factor = 0.1, max_iterations = 5000L,
                          tolerance = 1e-8, line_search_tol = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    n <- ncol(A)
    phi_floor <- 1e-12
    b_safe <- pmax(b, 1e-300)
    sigma <- sigma_factor * b_safe
    S_b_diag <- 1.0 / sigma^2
    phi_0 <- pmax(x0, 1e-300)
    if (length(phi_0) != n) stop("x0 length ", length(phi_0), " != n ", n)
    log_phi_0 <- log(phi_0)
    A_weighted <- A * matrix(S_b_diag, nrow = nrow(A), ncol = ncol(A), byrow = FALSE)
    At_Sb_A <- t(A) %*% A_weighted

    objective <- function(phi) {
        p <- pmax(phi, phi_floor)
        residual <- as.numeric(A %*% p) - b
        chi2 <- 0.5 * sum(residual^2 * S_b_diag)
        kl <- sum(p * (log(p + 1e-300) - log_phi_0) - p + phi_0)
        as.numeric(chi2 + kl)
    }
    gradient <- function(phi) {
        p <- pmax(phi, phi_floor)
        residual <- as.numeric(A %*% p) - b
        as.numeric(t(A) %*% (residual * S_b_diag)) + log(p + 1e-300) - log_phi_0
    }
    hessian <- function(phi) {
        p <- pmax(phi, phi_floor)
        At_Sb_A + diag(1.0 / (p + 1e-300))
    }

    phi <- phi_0
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
    list(spectrum = as.numeric(phi), iterations = as.integer(iteration),
         converged = (grad_norm < tolerance))
}

#' Wrapper around \code{\link{solve_imaxed}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_imaxed
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_imaxed <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            sigma_factor = 0.1, max_iterations = 5000L,
                            tolerance = 1e-8, line_search_tol = 1e-6,
                            calculate_errors = FALSE,
                            noise_level = 0.01, n_montecarlo = 100L,
                            save_result = FALSE, random_state = NULL) {
    x0_ref <- if (!is.null(initial_spectrum)) as.numeric(initial_spectrum)
              else rep(1.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = x0_ref,
        default_initial = rep(1.0, n_energy_bins),
        solve_func = make_solve_wrapper(solve_imaxed,
                                         sigma_factor = sigma_factor,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         line_search_tol = line_search_tol),
        solve_kwargs = list(),
        method_name = "IMAXED",
        extra_output = list(sigma_factor = sigma_factor),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
