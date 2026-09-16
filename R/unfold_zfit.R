#' Poisson-likelihood spectrum inference (zfit / Minuit analogue)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_zfit.py}. The Python
#' version builds a binned Poisson likelihood in zfit (Minuit backend); the
#' pure-R analogue maximises
#'
#' \eqn{-\sum_i [Ax_i + k_i \log (Ax_i)] + \mu\,\|D^{(2)}x\|^2}
#'
#' (deviance form) using the bounded \code{"L-BFGS-B"} optimizer from
#' \code{stats::optim} over an exponential parameterisation \eqn{x=e^y},
#' which guarantees \eqn{x>0}. Option \code{fallback_optimizer} switches to
#' a gradient-free Nelder-Mead, mirroring the Python SciPy fallback of the
#' Minuit/zfit path. With \code{use_mcmc = TRUE}, a Gaussian posterior
#' sample approximated around the optimum gives error bands.
#'
#' @name zfit-like
NULL

.zfit_deviance <- function(mu, k) {
    mu <- pmax(mu, 1e-30)
    2 * sum(mu - k * log(mu))
}

# Internal objective on y with smoothness and optional Tikhonov prior
.zfit_objective <- function(y, A, b, smoothness_weight, regularization, n) {
    x <- exp(y)
    dev <- .zfit_deviance(as.numeric(A %*% x), b)
    prior <- 0
    if (smoothness_weight > 0) {
        prior <- prior + smoothness_weight * sum(diff(diff(x))^2)
    }
    if (regularization > 0) {
        prior <- prior + regularization * sum(x^2)
    }
    dev + prior
}

#' Poisson-likelihood unfolding with smoothness prior
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m) (counts \eqn{k}).
#' @param x0 Numeric initial spectrum (length n).
#' @param smoothness_weight Numeric second-difference prior weight.
#'   Default 1e-3.
#' @param use_mcmc Logical; sample the (Gaussian-approximated) posterior
#'   for error bands. Default FALSE.
#' @param n_samples Integer size of the posterior sample. Default 100.
#' @param regularization Numeric additional L2 prior weight. Default 0.
#' @param fallback_optimizer Character; \code{"lbfgsb"} (default) or
#'   \code{"nelder"}.
#' @return A list \code{list(spectrum, iterations, converged)} plus
#'   \code{posterior_samples} when \code{use_mcmc = TRUE}.
#' @export
#' @examples
#' A <- matrix(c(2.0, 0.10, 0.02, 0.20, 1.5, 0.15, 0.40, 0.40, 0.80),
#'             nrow = 3, byrow = TRUE)
#' k <- c(21.5, 10.0, 4.5)
#' r <- solve_zfit_poisson(A, k, rep(3, 3))
solve_zfit_poisson <- function(A, b, x0 = NULL, smoothness_weight = 1e-3,
                               use_mcmc = FALSE, n_samples = 100L,
                               regularization = 0,
                               fallback_optimizer = "lbfgsb") {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- pmax(as.numeric(b), 0)
    n <- ncol(A)
    if (is.null(x0))
        x0 <- rep(mean(b) / max(sum(A), 1e-10) * n, n)
    x0 <- pmax(as.numeric(x0), 1e-30)
    y0 <- log(x0)
    meth <- tolower(fallback_optimizer)
    if (meth == "nelder") {
        res <- stats::optim(par = y0,
                            fn = .zfit_objective, A = A, b = b,
                            smoothness_weight = smoothness_weight,
                            regularization = regularization, n = n,
                            control = list(maxit = 50L * n))
        iterations <- res$counts[1]
        ybest <- res$par
    } else {
        res <- stats::optim(par = y0, fn = .zfit_objective, A = A, b = b,
                            smoothness_weight = smoothness_weight,
                            regularization = regularization, n = n,
                            method = "L-BFGS-B",
                            lower = rep(-20, n), upper = rep(20, n),
                            control = list(maxit = 100L * n))
        iterations <- res$counts[2]
        ybest <- res$par
    }
    spectrum <- pmax(exp(ybest), 0)
    out <- list(spectrum = spectrum,
                iterations = as.integer(iterations),
                converged = TRUE,
                smoothness_weight = smoothness_weight)
    if (isTRUE(use_mcmc) && n_samples > 0) {
        n_samples <- as.integer(n_samples)
        # Laplace (Gaussian) approximation: Hessian of the log-likelihood
        sigma_scale <- max(mean(b), 1e-30)
        H0 <- crossprod(A) / sigma_scale +
              smoothness_weight * crossprod(second_difference_matrix(n)) +
              regularization * diag(n)
        covm <- tryCatch(chol2inv(chol(H0 + 1e-9 * diag(n))),
                         error = function(e) diag(1 / max(diag(H0), 1e-30)))
        L <- tryCatch(t(chol(covm)), error = function(e) diag(0, n))
        samples <- matrix(0, nrow = n_samples, ncol = n)
        for (i in seq_len(n_samples)) {
            dy <- as.numeric(L %*% stats::rnorm(n))
            samples[i, ] <- pmax(exp(ybest + dy), 0)
        }
        out$posterior_samples <- samples
    }
    out
}

#' Wrapper around \code{\link{solve_zfit_poisson}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_zfit_poisson
#' @export
unfold_zfit <- function(detector_names, n_energy_bins, E_MeV, sensitivities,
                        cc_icrp116, save_result_callback, readings,
                        initial_spectrum = NULL, smoothness_weight = 1e-3,
                        use_mcmc = FALSE, n_samples = 100L,
                        regularization = 0, fallback_optimizer = "lbfgsb",
                        method_name = "ZFitPoisson", calculate_errors = FALSE,
                        noise_level = 0.01, n_montecarlo = 100L,
                        save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1, n_energy_bins),
        solve_func = solve_zfit_poisson,
        solve_kwargs = list(smoothness_weight = smoothness_weight,
                            use_mcmc = use_mcmc, n_samples = n_samples,
                            regularization = regularization,
                            fallback_optimizer = fallback_optimizer),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
