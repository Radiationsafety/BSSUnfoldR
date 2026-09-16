#' Ensemble Kalman Inversion (EKI) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_eki.py}.
#' Implements the EKI method of Iglesias et al. (2013) for Bayesian posterior
#' approximation without MCMC. The ensemble of particles is propagated through
#' the forward model and updated via the Kalman gain equation with optional
#' regularization for stability.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial guess (length n), used as the centre of the
#'   initial ensemble.
#' @param n_ensemble Integer; number of ensemble members. Default 50.
#' @param n_iterations Integer; number of EKI iterations. Default 50.
#' @param regularization Numeric; Tikhonov-style regularization added to the
#'   covariance diagonal for numerical stability. Default 1e-4.
#' @param inflation Numeric; covariance inflation factor to prevent ensemble
#'   collapse. Default 1.02.
#' @param noise_std Optional numeric; standard deviation of measurement noise.
#'   If \code{NULL}, estimated as 5% of \code{||b|| / sqrt(m)}.
#' @param random_state Optional integer seed.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' set.seed(7)
#' r <- solve_eki(A, b, rep(1, 3), n_ensemble = 20, n_iterations = 30)
solve_eki <- function(A, b, x0, n_ensemble = 50L, n_iterations = 50L,
                        regularization = 1e-4, inflation = 1.02,
                        noise_std = NULL, random_state = NULL) {
    v <- validate_system(A, b, x0 = x0)
    A <- v$A; b <- v$b; x0 <- v$x0
    m <- nrow(A); n <- ncol(A)
    if (!is.null(random_state)) set.seed(as.integer(random_state))
    if (is.null(noise_std)) {
        noise_std <- if (m > 0) 0.05 * sqrt(sum(b^2)) / sqrt(m) else 1e-6
    }
    noise_var <- noise_std^2
    sigma_prior <- abs(x0) + 1e-6
    ensemble <- matrix(rnorm(n * n_ensemble, mean = x0, sd = sigma_prior),
                        nrow = n, ncol = n_ensemble)
    for (it in seq_len(n_iterations)) {
        predictions <- A %*% ensemble
        pred_mean <- rowMeans(predictions)
        state_mean <- rowMeans(ensemble)
        pred_pert <- predictions - pred_mean
        state_pert <- ensemble - state_mean
        C_dd <- (pred_pert %*% t(pred_pert)) / max(n_ensemble - 1L, 1L)
        C_dd <- C_dd + (noise_var + regularization) * diag(m)
        C_md <- (state_pert %*% t(pred_pert)) / max(n_ensemble - 1L, 1L)
        C_d_inv <- tryCatch(solve(C_dd, diag(m)), error = function(e) {
            # Pseudo-inverse fallback
            sv <- svd(C_dd, nu = m, nv = 0)
            s_inv <- ifelse(sv$d > 1e-10 * max(sv$d), 1 / sv$d, 0)
            sv$u %*% (s_inv * (t(sv$u)))
        })
        # Innovation: b + noise*perturbation - predictions
        innovation <- b + noise_std * matrix(rnorm(m * n_ensemble),
                                              nrow = m, ncol = n_ensemble) -
                      predictions
        ensemble <- ensemble + C_md %*% C_d_inv %*% innovation
        ensemble <- ensemble * inflation
        ensemble <- pmax(ensemble, 0)
    }
    mean_spectrum <- rowMeans(ensemble)
    list(spectrum = pmax(as.numeric(mean_spectrum), 0),
         iterations = as.integer(n_iterations),
         converged = TRUE,
         ensemble = ensemble)
}

#' Wrapper around \code{\link{solve_eki}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_eki
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_eki <- function(detector_names, n_energy_bins, E_MeV,
                        sensitivities, cc_icrp116, save_result_callback,
                        readings, initial_spectrum = NULL,
                        n_ensemble = 50L, n_iterations = 50L,
                        regularization = 1e-4, inflation = 1.02,
                        noise_std = NULL,
                        calculate_errors = FALSE,
                        noise_level = 0.01, n_montecarlo = 100L,
                        save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(1.0 / max(n_energy_bins, 1L), n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_eki,
                                         n_ensemble = n_ensemble,
                                         n_iterations = n_iterations,
                                         regularization = regularization,
                                         inflation = inflation,
                                         noise_std = noise_std,
                                         random_state = random_state),
        solve_kwargs = list(),
        method_name = "EKI",
        extra_output = list(
            n_ensemble = as.integer(n_ensemble),
            n_iterations = as.integer(n_iterations),
            regularization = regularization,
            inflation = inflation
        ),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
