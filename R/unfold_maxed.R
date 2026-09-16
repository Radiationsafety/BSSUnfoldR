#' MAXED unfolding (Maximum Entropy Deconvolution)
#'
#' Implements the MAXED algorithm of Reginatto & Goldhagen (1999) by
#' minimising, in log-space \eqn{y_i = ln(x_i)}, the primal
#' \deqn{f(x) = -S(x) + 0.5 * sum_j (b_j - (A x)_j)^2 / sigma_j^2}
#' where \eqn{S(x) = -sum_i x_i ln(x_i / x0_i) + sum_i x_i - sum_i x0_i} is
#' the Shannon entropy relative to the reference spectrum \code{x0}. The
#' L-BFGS-B optimiser of \code{optim} is used.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric reference (prior) spectrum (length n).
#' @param sigma_factor Numeric; relative measurement uncertainty. Default 0.1.
#' @param max_iterations Positive integer; default 5000.
#' @param tolerance Positive numeric; gradient convergence tolerance. Default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_maxed(A, b, rep(1, 3))
solve_maxed <- function(A, b, x0, sigma_factor = 0.1, max_iterations = 5000L,
                        tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    b_safe <- pmax(b, 1e-300)
    sigma <- sigma_factor * b_safe
    sigma2_inv <- 1.0 / (sigma^2)
    phi_0 <- pmax(x0, 1e-300)
    log_phi_0 <- log(phi_0)
    AT <- t(A)

    fg <- function(y) {
        x <- exp(y)
        folded <- as.numeric(A %*% x)
        residual <- b - folded
        f_ent <- sum(x * (y - log_phi_0) - x + phi_0)
        f_chi <- 0.5 * sum(residual^2 * sigma2_inv)
        AT_resid_over_sigma2 <- as.numeric(AT %*% (residual * sigma2_inv))
        grad <- x * (y - log_phi_0 - AT_resid_over_sigma2)
        list(f = f_ent + f_chi, g = grad)
    }

    f_eval <- function(y) fg(y)$f
    g_eval <- function(y) fg(y)$g

    y0 <- log(phi_0)
    result <- stats::optim(
        par = y0, fn = f_eval, gr = g_eval,
        method = "L-BFGS-B",
        control = list(maxit = max_iterations,
                       ndeps = rep(1e-8, length(y0)),
                       pgtol = tolerance)
    )
    x_opt <- exp(result$par)
    iterations <- if (!is.null(result$counts[1L])) result$counts[1L] else 0L
    converged <- (result$convergence == 0L)
    list(spectrum = as.numeric(x_opt), iterations = as.integer(iterations),
         converged = converged)
}

#' Wrapper around \code{\link{solve_maxed}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_maxed
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_maxed <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, initial_spectrum = NULL,
                          sigma_factor = 0.1, max_iterations = 5000L,
                          tolerance = 1e-6,
                          calculate_errors = FALSE,
                          noise_level = 0.01, n_montecarlo = 100L,
                          save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_ref <- if (!is.null(initial_spectrum)) {
        as.numeric(initial_spectrum)
    } else {
        rep(1.0, n_energy_bins)
    }
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = x0_ref,
        default_initial = rep(1.0, n_energy_bins),
        solve_func = make_solve_wrapper(solve_maxed,
                                         sigma_factor = sigma_factor,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "MAXED",
        extra_output = list(sigma_factor = sigma_factor),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
