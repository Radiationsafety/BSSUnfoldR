#' Bayesian parametric unfolding (FRUIT-like model with Metropolis-Hastings)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_bayesian_parametric.py}.
#' Random-walk Metropolis-Hastings over the five parameters of the FRUIT-like
#' spectral model (Maxwellian thermal + \code{1/E} epithermal + evaporation
#' fast).  Note that this method, unlike \code{\link{unfold_fruit_like}}, is fed
#' the base-10 lethargy widths \code{compute_log_steps(E_MeV)}, not their
#' natural-log equivalents.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid (length n).
#' @param log_steps Numeric base-10 log-energy bin widths (length n).
#' @param sigma Numeric; measurement uncertainty on the readings. Default 0.02.
#' @param initial_params Optional named numeric vector of starting values for
#'   \code{A_th, T_th, A_epi, A_f, T_ev}.
#' @param n_samples Integer; number of post-burn-in samples. Default 1000.
#' @param burn_in Integer; samples discarded before the chain is recorded.
#'   Default 200.  The chain runs \code{n_samples + burn_in} steps in total.
#' @param proposal_scale Numeric; the proposal is a relative random walk,
#'   \code{v + rnorm(0, proposal_scale * |v|)}. Default 0.1.
#' @param random_state Optional integer seed.
#' @return A list \code{list(spectrum, success, message, n_samples, means,
#'   samples, acceptance_rate)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5) * 0.2
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_bayesian_parametric(A, b, E, compute_log_steps(E),
#'                                 n_samples = 60, burn_in = 20)
solve_bayesian_parametric <- function(A, b, E, log_steps, sigma = 0.02,
                                      initial_params = NULL,
                                      n_samples = 1000L, burn_in = 200L,
                                      proposal_scale = 0.1,
                                      random_state = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E); log_steps <- as.numeric(log_steps)
    n_samples <- as.integer(n_samples); burn_in <- as.integer(burn_in)
    if (is.null(initial_params)) {
        initial_params <- c(A_th = 1e-6, T_th = 0.025e-6, A_epi = 1e-6,
                            A_f = 1e-6, T_ev = 2.0)
    }
    init <- initial_params[c("A_th", "T_th", "A_epi", "A_f", "T_ev")]
    chain <- .bayes_param_mh(A, b, E, log_steps, sigma, init,
                             n_samples, burn_in, proposal_scale, random_state)
    m <- chain$means
    spectrum <- .parametric_model(E, m[["A_th"]], m[["T_th"]], m[["A_epi"]],
                                  m[["A_f"]], m[["T_ev"]]) * log_steps
    list(spectrum = as.numeric(spectrum), success = TRUE,
         message = "Bayesian estimation completed",
         n_samples = as.integer(n_samples),
         means = m, samples = chain$samples,
         acceptance_rate = chain$acceptance_rate)
}

#' Log prior of the parametric model: independent uniform densities over the
#' admissible range of each parameter, exactly as in the Python original.
#' @keywords internal
#' @noRd
.bayes_param_log_prior <- function(p) {
    if (p[["A_th"]] < 0 || p[["A_th"]] > 1e-3) return(-Inf)
    lp <- -log(1e-3)
    if (p[["T_th"]] < 1e-9 || p[["T_th"]] > 1e-3) return(-Inf)
    lp <- lp - log(1e-3 - 1e-9)
    if (p[["A_epi"]] < 0 || p[["A_epi"]] > 1e-3) return(-Inf)
    lp <- lp - log(1e-3)
    if (p[["A_f"]] < 0 || p[["A_f"]] > 1e-3) return(-Inf)
    lp <- lp - log(1e-3)
    if (p[["T_ev"]] < 0.1 || p[["T_ev"]] > 20) return(-Inf)
    lp - log(19.9)
}

#' @keywords internal
#' @noRd
.bayes_param_log_posterior <- function(p, A, b, E, log_steps, sigma) {
    lp <- .bayes_param_log_prior(p)
    if (!is.finite(lp)) return(-Inf)
    spectrum <- .parametric_model(E, p[["A_th"]], p[["T_th"]], p[["A_epi"]],
                                  p[["A_f"]], p[["T_ev"]])
    residual <- b - as.numeric(A %*% (spectrum * log_steps))
    ll <- -0.5 * sum((residual / sigma)^2)
    if (!is.finite(ll)) return(-Inf)
    lp + ll
}

#' Random-walk Metropolis-Hastings matching the Python chain step for step:
#' relative Gaussian proposal, a floor (never a ceiling) applied per parameter,
#' \code{n_samples + burn_in} iterations.
#' @keywords internal
#' @noRd
.bayes_param_mh <- function(A, b, E, log_steps, sigma, init,
                            n_samples, burn_in, proposal_scale, random_state) {
    if (!is.null(random_state)) set.seed(as.integer(random_state))
    current <- init
    current_lp <- .bayes_param_log_posterior(current, A, b, E, log_steps, sigma)
    samples <- matrix(0.0, nrow = n_samples, ncol = 5L,
                      dimnames = list(NULL, names(init)))
    accepted <- 0L
    temp_idx <- c(2L, 5L)          # T_th, T_ev
    amp_idx <- c(1L, 3L, 4L)       # A_th, A_epi, A_f
    for (i in seq_len(n_samples + burn_in)) {
        proposed <- current + rnorm(5L, sd = proposal_scale * abs(current) +
                                        1e-15)
        proposed[temp_idx] <- pmax(proposed[temp_idx], 1e-9)
        proposed[amp_idx] <- pmax(proposed[amp_idx], 0)
        proposed_lp <- .bayes_param_log_posterior(proposed, A, b, E,
                                                 log_steps, sigma)
        if (log(runif(1)) < proposed_lp - current_lp) {
            current <- proposed
            current_lp <- proposed_lp
            accepted <- accepted + 1L
        }
        if (i > burn_in) samples[i - burn_in, ] <- current
    }
    list(means = stats::setNames(colMeans(samples), names(init)),
         samples = samples, acceptance_rate = accepted / (n_samples + burn_in))
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
                                       sigma = 0.02, n_samples = 1000L,
                                       burn_in = 200L, proposal_scale = 0.1,
                                       calculate_errors = FALSE,
                                       noise_level = 0.01, n_montecarlo = 100L,
                                       save_result = FALSE,
                                       random_state = NULL,
                                       max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b
    # base-10 lethargy widths, deliberately *not* converted to natural log
    log_steps <- compute_log_steps(E_MeV)
    solver <- function(A, b, x0 = NULL, ...) {
        res <- solve_bayesian_parametric(A, b, E_MeV, log_steps,
                                         sigma = sigma, n_samples = n_samples,
                                         burn_in = burn_in,
                                         proposal_scale = proposal_scale,
                                         random_state = random_state)
        list(spectrum = res$spectrum, iterations = res$n_samples,
             converged = res$success)
    }
    x0_default <- rep(mean(b) / max(mean(rowSums(A)), 1e-30), n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver, solve_kwargs = list(),
        method_name = "bayesian_parametric",
        extra_output = list(sigma = sigma, n_samples = n_samples,
                            burn_in = burn_in,
                            proposal_scale = proposal_scale),
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result, max_neutron_energy = max_neutron_energy)
}
