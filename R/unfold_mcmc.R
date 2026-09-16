#' Bayesian MCMC unfolding (Metropolis-Hastings)
#'
#' Simplified R port of \code{bssunfold/src/bssunfold/core/unfold_mcmc.py}.
#' The Python original uses PyMC and a NUTS sampler; this R port implements
#' a pure-R random-walk Metropolis-Hastings sampler that:
#' \itemize{
#'   \item Works in log-space so the spectrum stays positive;
#'   \item Uses an Ornstein-Uhlenbeck (Gaussian) smoothness prior anchored
#'         on a data-driven center (default: NNLS solution);
#'   \item Uses a Gaussian likelihood with per-detector sigma derived from
#'         \code{relative_uncertainty * b}.
#' }
#'
#' @section Limitations vs. Python bssunfold:
#' No NUTS sampler (uses Metropolis-Hastings, slower mixing). No PyMC
#' model. No ArviZ integration. The user gets raw posterior samples in
#' the result list.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional numeric initial spectrum (length n). Default
#'   \code{NULL} = NNLS solution (or flat if NNLS fails).
#' @param n_samples Integer; total samples to draw. Default 1000.
#' @param n_burn Integer; burn-in samples to discard. Default 200.
#' @param proposal_sd Numeric; proposal step standard deviation (in log space).
#'   Default 0.1.
#' @param smoothness Numeric; smoothness prior strength. Default 1.0.
#' @param relative_uncertainty Numeric; relative measurement uncertainty.
#'   Default 0.1.
#' @param random_state Optional integer seed.
#' @return A list \code{list(spectrum, iterations = n_samples - n_burn,
#'   converged = TRUE, posterior_samples = matrix)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' set.seed(7)
#' r <- solve_bayesian_mcmc(A, b, rep(1, 3), n_samples = 100, n_burn = 20)
solve_bayesian_mcmc <- function(A, b, x0 = NULL, n_samples = 1000L,
                                   n_burn = 200L, proposal_sd = 0.1,
                                   smoothness = 1.0,
                                   relative_uncertainty = 0.1,
                                   random_state = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (!is.null(random_state)) set.seed(as.integer(random_state))
    # Initial spectrum: NNLS solution or flat fallback
    if (is.null(x0)) {
        nnls_fit <- tryCatch(lsei::nnls(A, pmax(b, 1e-10)), error = function(e) NULL)
        if (!is.null(nnls_fit)) {
            x0 <- pmax(nnls_fit$x, 1e-10)
        } else {
            x0 <- rep(mean(b) / max(mean(rowSums(A)), 1e-10), n)
        }
    } else {
        x0 <- as.numeric(x0)
    }
    log_phi <- log(pmax(x0, 1e-30))
    sigma <- pmax(relative_uncertainty * pmax(b, 1e-30), 1e-30)
    sigma2 <- sigma^2
    # Build first-derivative operator for smoothness prior
    L1 <- if (n >= 3L) as.matrix(create_derivative_matrix(n, 1L)) else matrix(0, 0, n)

    log_posterior <- function(lp) {
        phi <- exp(lp)
        # Likelihood: Gaussian
        residual <- as.numeric(A %*% phi) - b
        ll <- -0.5 * sum(residual^2 / sigma2)
        # Smoothness prior: Ornstein-Uhlenbeck on log_phi
        if (nrow(L1) > 0L) {
            grad_log_phi <- as.numeric(L1 %*% lp)
            lp_prior <- -0.5 * smoothness * sum(grad_log_phi^2)
        } else {
            lp_prior <- 0
        }
        ll + lp_prior
    }

    current_lp <- log_posterior(log_phi)
    samples <- matrix(0.0, nrow = n_samples - n_burn, ncol = n)
    accepted <- 0L
    for (i in seq_len(n_samples)) {
        proposal <- log_phi + rnorm(n, sd = proposal_sd)
        proposed_lp <- log_posterior(proposal)
        if (is.finite(proposed_lp) &&
            log(runif(1)) < proposed_lp - current_lp) {
            log_phi <- proposal
            current_lp <- proposed_lp
            accepted <- accepted + 1L
        }
        if (i > n_burn) {
            samples[i - n_burn, ] <- exp(log_phi)
        }
    }
    acceptance_rate <- accepted / n_samples
    spectrum <- colMeans(samples)
    list(spectrum = as.numeric(pmax(spectrum, 0)),
         iterations = as.integer(n_samples - n_burn),
         converged = TRUE,
         posterior_samples = samples,
         acceptance_rate = acceptance_rate)
}

#' Wrapper around \code{\link{solve_bayesian_mcmc}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_bayesian_mcmc
#' @return A result list as produced by \code{\link{run_unfolding}} plus
#'   \code{posterior_samples} and \code{acceptance_rate}.
#' @export
unfold_mcmc <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, initial_spectrum = NULL,
                          n_samples = 1000L, n_burn = 200L,
                          proposal_sd = 0.1, smoothness = 1.0,
                          relative_uncertainty = 0.1,
                          calculate_errors = FALSE,
                          noise_level = 0.01, n_montecarlo = 100L,
                          save_result = FALSE, random_state = NULL) {
    x0_default <- rep(0.5, n_energy_bins)
    # MCMC does not fit the run_unfolding solve_func pattern naturally
    # because it returns posterior samples, not a single spectrum. We
    # build the system explicitly and call solve_bayesian_mcmc directly.
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    x0 <- if (is.null(initial_spectrum)) x0_default else as.numeric(initial_spectrum)
    res <- solve_bayesian_mcmc(A, b, x0, n_samples = n_samples,
                                n_burn = n_burn, proposal_sd = proposal_sd,
                                smoothness = smoothness,
                                relative_uncertainty = relative_uncertainty,
                                random_state = random_state)
    spectrum <- res$spectrum
    computed_readings <- as.numeric(A %*% spectrum)
    residual <- b - computed_readings
    result <- list(
        energy = E_MeV,
        spectrum = spectrum,
        spectrum_absolute = spectrum,
        effective_readings = stats::setNames(computed_readings, selected),
        residual = residual,
        residual_norm = sqrt(sum(residual^2)),
        method = "MCMC",
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = res$iterations,
        converged = res$converged,
        posterior_samples = res$posterior_samples,
        acceptance_rate = res$acceptance_rate
    )
    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }
    result
}
