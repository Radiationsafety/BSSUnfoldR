#' STAY'SL unfolding (linear Bayesian update)
#'
#' A single-step (one-shot) linear Bayesian update:
#' \deqn{x = x0 + Cx A^T (Cb + A Cx A^T)^{-1} (b - A x0)}
#' Default covariances are diagonal: \eqn{Cb = diag((rel * b)^2)},
#' \eqn{Cx = diag((prior * x0)^2)}. A small Tikhonov term is added to the
#' bracket for numerical stability.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric prior spectrum (length n).
#' @param relative_uncertainty Numeric; default 0.1.
#' @param prior_uncertainty Numeric; default 1.0.
#' @param Cb Optional explicit measurement covariance (m x m). Default \code{NULL}.
#' @param Cx Optional explicit prior covariance (n x n). Default \code{NULL}.
#' @param regularization Numeric Tikhonov term. Default 1e-12.
#' @return A list \code{list(spectrum = ..., iterations = 1L, converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_staysl(A, b, rep(1, 3))
solve_staysl <- function(A, b, x0, relative_uncertainty = 0.1,
                         prior_uncertainty = 1.0, Cb = NULL, Cx = NULL,
                         regularization = 1e-12) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    if (nrow(A) == 0L || length(b) == 0L) {
        stop("Response matrix and measurements must be non-empty")
    }
    m <- nrow(A); n <- ncol(A)
    AT <- t(A)
    if (is.null(Cb)) {
        b_safe <- pmax(abs(b), 1e-12)
        Cb <- diag((relative_uncertainty * b_safe)^2, nrow = m, ncol = m)
    } else {
        Cb <- as.matrix(Cb)
    }
    if (is.null(Cx)) {
        x_safe <- pmax(abs(x0), 1e-12)
        Cx <- diag((prior_uncertainty * x_safe)^2, nrow = n, ncol = n)
    } else {
        Cx <- as.matrix(Cx)
    }
    bracket <- Cb + A %*% Cx %*% AT + regularization * diag(m)
    gain <- Cx %*% AT %*% solve(bracket)
    spectrum <- x0 + as.numeric(gain %*% (b - as.numeric(A %*% x0)))
    list(spectrum = pmax(spectrum, 0.0), iterations = 1L, converged = TRUE)
}

#' Wrapper around \code{\link{solve_staysl}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_staysl
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_staysl <- function(detector_names, n_energy_bins, E_MeV,
                           sensitivities, cc_icrp116, save_result_callback,
                           readings, initial_spectrum = NULL,
                           relative_uncertainty = 0.1,
                           prior_uncertainty = 1.0,
                           calculate_errors = FALSE,
                           noise_level = 0.01, n_montecarlo = 100L,
                           save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_staysl,
                                         relative_uncertainty = relative_uncertainty,
                                         prior_uncertainty = prior_uncertainty),
        solve_kwargs = list(),
        method_name = "STAY'SL",
        extra_output = list(relative_uncertainty = relative_uncertainty,
                            prior_uncertainty = prior_uncertainty),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
