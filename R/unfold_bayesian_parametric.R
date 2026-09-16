#' Bayesian parametric unfolding (FRUIT-like with Metropolis-Hastings)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_bayesian_parametric.py}.
#' Bayesian inference on the 5 parameters of the FRUIT-like parametric model
#' (Maxwellian + 1/E + evaporation) using random-walk Metropolis-Hastings.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid (length n).
#' @param log_steps Numeric log-energy bin widths (length n).
#' @param n_samples Integer; total MCMC samples. Default 1000.
#' @param n_burn Integer; burn-in samples. Default 200.
#' @param proposal_sd Numeric; proposal step std (in log space). Default 0.3.
#' @param sigma Numeric; measurement noise std. Default 0.05 * mean(b).
#' @param random_state Optional integer seed.
#' @return A list \code{list(spectrum, n_iterations, converged, posterior_samples)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' set.seed(7)
#' r <- solve_bayesian_parametric(A, b, E, compute_log_steps(E) * log(10),
#'                                 n_samples = 50, n_burn = 10)
solve_bayesian_parametric <- function(A, b, E, log_steps, n_samples = 1000L,
                                        n_burn = 200L, proposal_sd = 0.3,
                                        sigma = NULL, random_state = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E); log_steps <- as.numeric(log_steps)
    if (!is.null(random_state)) set.seed(as.integer(random_state))
    if (is.null(sigma)) sigma <- 0.05 * mean(b)
    # Parameter bounds (matches Python original)
    lower <- c(0, 1e-9, 0, 0, 0.1)
    upper <- c(1e-3, 1e-3, 1e-3, 1e-3, 20.0)
    p_init <- c(1e-6, 0.025e-6, 1e-6, 1e-6, 2.0)
    log_posterior <- function(p) {
        if (any(p < lower) || any(p > upper)) return(-Inf)
        spectrum <- .parametric_model(E, p[1], p[2], p[3], p[4], p[5])
        computed <- as.numeric(A %*% (spectrum * log_steps))
        residual <- b - computed
        ll <- -0.5 * sum((residual / sigma)^2)
        ll
    }
    log_p <- log_posterior(p_init)
    samples <- matrix(0.0, nrow = n_samples - n_burn, ncol = 5L)
    accepted <- 0L
    p_cur <- p_init
    for (i in seq_len(n_samples)) {
        proposal <- p_cur + rnorm(5, sd = proposal_sd)
        proposal <- pmax(pmin(proposal, upper), lower)
        log_p_new <- log_posterior(proposal)
        if (is.finite(log_p_new) && log(runif(1)) < log_p_new - log_p) {
            p_cur <- proposal
            log_p <- log_p_new
            accepted <- accepted + 1L
        }
        if (i > n_burn) samples[i - n_burn, ] <- p_cur
    }
    # Mean spectrum from posterior samples
    spectra <- matrix(0.0, nrow = n_samples - n_burn, ncol = length(E))
    for (i in seq_len(n_samples - n_burn)) {
        p <- samples[i, ]
        sp <- .parametric_model(E, p[1], p[2], p[3], p[4], p[5]) * log_steps
        spectra[i, ] <- sp
    }
    spectrum <- colMeans(spectra)
    list(spectrum = pmax(as.numeric(spectrum), 0),
         n_iterations = as.integer(n_samples - n_burn),
         converged = TRUE,
         posterior_samples = samples,
         acceptance_rate = accepted / n_samples)
}

#' Wrapper around \code{\link{solve_bayesian_parametric}} for the unified
#' workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_bayesian_parametric
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_bayesian_parametric <- function(detector_names, n_energy_bins, E_MeV,
                                          sensitivities, cc_icrp116,
                                          save_result_callback, readings,
                                          initial_spectrum = NULL,
                                          n_samples = 1000L, n_burn = 200L,
                                          proposal_sd = 0.3, sigma = NULL,
                                          random_state = NULL,
                                          calculate_errors = FALSE,
                                          noise_level = 0.01,
                                          n_montecarlo = 100L,
                                          save_result = FALSE) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    log_steps <- compute_log_steps(E_MeV) * log(10.0)
    res <- solve_bayesian_parametric(A, b, E_MeV, log_steps,
                                       n_samples = n_samples, n_burn = n_burn,
                                       proposal_sd = proposal_sd, sigma = sigma,
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
        method = "BayesianParametric",
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = res$n_iterations,
        converged = res$converged,
        posterior_samples = res$posterior_samples,
        acceptance_rate = res$acceptance_rate
    )
    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }
    result
}
