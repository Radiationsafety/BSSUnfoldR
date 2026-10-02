#' Stochastic parametric unfolding with the Fission model (BonnerFinder)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_fission_ga.py},
#' itself a port of the two-stage \code{BonnerFinder()} algorithm of
#' I. N. Ogorodnikov, "Inverse problems of spectroscopy and spectrometry in
#' applied research", No. 2 (10), pp. 42-83 (2024), sections 4-5.  Stage 1
#' is a stochastic global search of the seven-parameter Fission-model
#' hypercube (article: SciLab \code{optim_ga}; Python and this port use
#' differential evolution) minimizing the L1 discrepancy of the folded
#' readings (article eq. 4.32).  Stage 2 refines the best point with
#' nonlinear least squares (article: SciLab \code{leastsq}).  An optional
#' free overall scale \code{phi_scale} (log10-parameterized) is fitted so
#' the model matches absolutely calibrated readings; \code{fit_scale =
#' FALSE} recovers the exact 7-parameter normalized formulation.
#'
#' @name fission-ga-methods
NULL

# Fixed constants of the Fission model (article eq. 4.29); identical to the
# FRUIT constants (.fruit_T0 / .fruit_Ed) already in unfold_parametric.R.
.fga_T0 <- 2.53e-8   # Thermal peak energy (MeV)
.fga_Ed <- 7.07e-8   # Epithermal cutoff parameter (MeV)

.fga_param_names <- c("a1", "a2", "a3", "b", "beta", "alpha", "TF")
.fga_param_bounds <- list(a1 = c(0, 1), a2 = c(0, 1), a3 = c(0, 1),
                          b = c(-0.5, 0.5), beta = c(1e-4, 1),
                          alpha = c(0, 1), TF = c(1, 2))

.fga_theta_bounds <- function(fit_scale) {
    lo <- vapply(.fga_param_names, function(p) .fga_param_bounds[[p]][1],
                 numeric(1))
    hi <- vapply(.fga_param_names, function(p) .fga_param_bounds[[p]][2],
                 numeric(1))
    if (isTRUE(fit_scale)) {
        lo <- c(lo, -12)
        hi <- c(hi, 12)
    }
    list(lo = lo, hi = hi)
}

.fga_default_theta <- function(fit_scale) {
    theta <- c(0.3, 0.3, 0.4, 0, 0.1, 0.5, 1.5)
    if (isTRUE(fit_scale)) theta <- c(theta, 0)   # phi_scale = 1
    theta
}

#' Three-fraction Fission model spectrum (article eq. 4.29), fluence per
#' unit lethargy.  Port of fission_model().
.fga_model <- function(E, a1, a2, a3, b, beta, alpha, TF) {
    E <- as.numeric(E)
    thermal <- (E / (.fga_T0^2)) * exp(-E / .fga_T0)
    epithermal <- (1 - exp(-((E / .fga_Ed)^2))) * E^(b - 1) * exp(-E / beta)
    fast <- E^alpha * exp(-E / TF)
    a1 * thermal + a2 * epithermal + a3 * fast
}

#' Model shape for the flat parameter vector, including phi_scale when the
#' 8th element is present.  Port of _model_shape().
.fga_shape <- function(theta, E) {
    phi_scale <- if (length(theta) > 7L) 10^theta[8L] else 1
    phi_scale * .fga_model(E, theta[1], theta[2], theta[3], theta[4],
                           theta[5], theta[6], theta[7])
}

#' Fold the model through the response matrix.  Port of _folded_readings().
.fga_fold <- function(theta, A, b_vec, E, ln_steps) {
    as.numeric(A %*% (.fga_shape(theta, E) * ln_steps))
}

.fga_residual <- function(theta, A, b_vec, E, ln_steps) {
    .fga_fold(theta, A, b_vec, E, ln_steps) - b_vec
}

#' Article eq. 4.32: L1 discrepancy of the folded readings.
.fga_target <- function(theta, A, b_vec, E, ln_steps) {
    sum(abs(.fga_residual(theta, A, b_vec, E, ln_steps)))
}

#' Convert a user parameter list into the flat theta vector, or NULL.
#' Port of _theta_from_initial().
.fga_theta_from_initial <- function(initial_params, fit_scale) {
    if (!is.list(initial_params) || length(initial_params) == 0L) {
        return(NULL)
    }
    theta <- .fga_default_theta(fit_scale)
    for (i in seq_along(.fga_param_names)) {
        nm <- .fga_param_names[i]
        if (!is.null(initial_params[[nm]])) theta[i] <- as.numeric(
            initial_params[[nm]])
    }
    if (isTRUE(fit_scale) && !is.null(initial_params$phi_scale)) {
        theta[8L] <- log10(max(as.numeric(initial_params$phi_scale), 1e-30))
    }
    bb <- .fga_theta_bounds(fit_scale)
    pmax(bb$lo, pmin(bb$hi, theta))
}

#' Differential evolution, best1bin strategy, reproducing the scipy
#' defaults used by Python (popsize multiplier, mutation U(0.5,1),
#' recombination 0.7, latinhypercube init, no polish, workers = 1).  The
#' RNG stream cannot match scipy bit-for-bit; given the same seed the R
#' run is self-reproducible.
.fga_differential_evolution <- function(fn, lo, hi, popsize, maxiter, tol,
                                         ...) {
    d <- length(lo)
    np <- max(as.integer(popsize) * d, 4L)
    # latinhypercube initial population: one point per stratum per dim
    pop <- matrix(0, nrow = np, ncol = d)
    for (j in seq_len(d)) {
        perm <- sample.int(np)
        pop[, j] <- (perm - 1 + stats::runif(np)) / np
    }
    pop <- sweep(pop, 2, hi - lo, `*`)
    pop <- sweep(pop, 2, lo, `+`)

    fpop <- vapply(seq_len(np), function(i) fn(pop[i, ], ...), numeric(1))
    best_i <- which.min(fpop)
    prev_best <- Inf
    hist <- numeric(0L)
    n_evals <- np

    for (it in seq_len(max(1L, as.integer(maxiter)))) {
        for (i in seq_len(np)) {
            idx <- sample(np, 3L)
            while (idx[1] == i) idx[1] <- sample(np, 1L)
            while (idx[2] == i || idx[2] == idx[1]) idx[2] <- sample(np, 1L)
            while (idx[3] == i || idx[3] %in% idx[1:2]) {
                idx[3] <- sample(np, 1L)
            }
            f <- stats::runif(1, 0.5, 1.0)
            donor <- pop[best_i, ] +
                f * (pop[idx[1], ] - pop[idx[2], ])
            cross <- stats::runif(d) < 0.7
            jrand <- sample.int(d, 1L)
            cross[jrand] <- TRUE
            trial <- ifelse(cross, donor, pop[i, ])
            trial <- pmax(lo, pmin(hi, trial))
            ft <- fn(trial, ...)
            n_evals <- n_evals + 1L
            if (ft <= fpop[i]) {
                pop[i, ] <- trial
                fpop[i] <- ft
                if (ft < fpop[best_i]) best_i <- i
            }
        }
        cur_best <- fpop[best_i]
        hist <- c(hist, cur_best)
        if (length(hist) > 10L) hist <- hist[seq_len(10L)]
        if (is.finite(cur_best) && is.finite(prev_best)) {
            if (abs(prev_best - cur_best) < tol &&
                (max(hist) - min(hist)) < tol) {
                break
            }
        }
        prev_best <- cur_best
    }
    list(x = pop[best_i, ], fun = fpop[best_i], nfev = n_evals,
         nit = as.integer(length(hist)))
}

#' Levenberg-Marquardt style nonlinear least-squares refinement with a
#' finite-difference Jacobian.  scipy.optimize.least_squares has no base-R
#' equivalent; the bounded ("trf") mode clips every trial step to the
#' parameter box, the unbounded ("lm") mode relies on the caller to clip
#' afterwards, mirroring _run_lm().
.fga_least_squares <- function(res_fn, theta0, lo, hi, bounded,
                               max_nfev) {
    theta <- as.numeric(theta0)
    r <- res_fn(theta)
    nfev <- 1L
    cost <- sum(r^2)
    lambda <- 1e-3
    p <- length(theta)
    converged <- FALSE
    for (it in seq_len(200L)) {
        if (nfev >= max_nfev) break
        J <- matrix(0, nrow = length(r), ncol = p)
        for (j in seq_len(p)) {
            h <- 1e-7 * max(1, abs(theta[j]))
            tp <- theta
            tp[j] <- tp[j] + h
            rp <- res_fn(tp)
            nfev <- nfev + 1L
            J[, j] <- (rp - r) / h
            if (nfev >= max_nfev) break
        }
        g <- as.numeric(crossprod(J, r))
        if (max(abs(g)) < 1e-14) {
            converged <- TRUE
            break
        }
        improved <- FALSE
        for (try_step in seq_len(30L)) {
            if (nfev >= max_nfev) break
            JJ <- crossprod(J) + lambda * diag(p)
            delta <- tryCatch(as.numeric(solve(JJ, -g)),
                              error = function(e) NULL)
            if (is.null(delta)) {
                lambda <- lambda * 10
                next
            }
            trial <- theta + delta
            if (bounded) trial <- pmax(lo, pmin(hi, trial))
            rt <- res_fn(trial)
            nfev <- nfev + 1L
            ct <- sum(rt^2)
            if (ct < cost) {
                theta <- trial
                r <- rt
                cost <- ct
                lambda <- max(lambda * 0.5, 1e-12)
                improved <- TRUE
                break
            }
            lambda <- lambda * 2
        }
        if (!improved) break
    }
    list(theta = theta, cost = sqrt(cost), nfev = nfev,
         success = converged || is.finite(cost),
         message = if (converged) {
             "Optimization terminated successfully: gradient converged."
         } else if (nfev >= max_nfev) {
             "The maximum number of function evaluations is exceeded."
         } else {
             "Optimization terminated successfully: step size underflow."
         })
}

#' Solve by Fission-model GA (differential evolution + least squares)
#'
#' Low-level entry point, port of \code{solve_fission_ga}
#' (\code{unfold_fission_ga.py:258}).
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid in MeV (length n).
#' @param ln_steps Numeric natural-logarithmic bin widths (d ln E,
#'   length n).
#' @param initial_params Optional named list of starting parameter values
#'   (keys \code{a1, a2, a3, b, beta, alpha, TF, phi_scale}); refined as an
#'   extra stage-2 start.
#' @param fit_scale Logical; fit a free overall scale factor
#'   \code{phi_scale} (default TRUE).
#' @param ga_popsize Population multiplier of the differential-evolution
#'   stage (default 15, a scipy-style multiplier: population size is
#'   \code{ga_popsize * n_params}).
#' @param ga_maxiter Maximum generations of the stochastic stage.
#' @param ga_tol Convergence tolerance of the stochastic stage.
#' @param lm_method Stage-2 least-squares method, \code{"trf"} (bounded,
#'   default) or \code{"lm"} (unbounded with clipping, as SciLab
#'   \code{leastsq}).
#' @param lm_max_nfev Maximum function evaluations of one stage-2 run.
#' @param random_state Optional integer seed for the stochastic stage.
#' @return A list \code{list(spectrum, iterations, converged, message,
#'   params)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' E <- c(0.025, 0.5, 2.0)
#' r <- solve_fission_ga(A, b, E, log(E) * 0 + 1,
#'                        ga_maxiter = 10L, random_state = 1L)
solve_fission_ga <- function(A, b, E, ln_steps, initial_params = NULL,
                             fit_scale = TRUE, ga_popsize = 15L,
                             ga_maxiter = 100L, ga_tol = 1e-10,
                             lm_method = "trf", lm_max_nfev = 2000L,
                             random_state = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    E <- as.numeric(E)
    ln_steps <- as.numeric(ln_steps)
    bb <- .fga_theta_bounds(fit_scale)
    lo <- bb$lo; hi <- bb$hi
    total_nfev <- 0L

    if (!is.null(random_state)) set.seed(as.integer(random_state))

    # ---- Stage 1: stochastic (differential-evolution) global search ----
    ga <- .fga_differential_evolution(
        .fga_target, lo, hi, popsize = ga_popsize,
        maxiter = as.integer(ga_maxiter), tol = ga_tol,
        A = A, b_vec = b, E = E, ln_steps = ln_steps)
    total_nfev <- total_nfev + as.integer(ga$nfev)

    candidates <- list(pmax(lo, pmin(hi, ga$x)))

    theta_user <- .fga_theta_from_initial(initial_params, fit_scale)
    if (!is.null(theta_user)) candidates[[length(candidates) + 1L]] <-
        theta_user

    # ---- Stage 2: nonlinear least-squares refinement ----
    res_fn <- function(theta) .fga_residual(theta, A, b, E, ln_steps)
    best <- NULL
    for (theta0 in candidates) {
        if (identical(lm_method, "lm")) {
            # Levenberg-Marquardt (SciLab leastsq) does not support
            # bounds; project the start inside and clip after the fit.
            start <- pmax(lo + 1e-9 * (hi - lo),
                          pmin(hi - 1e-9 * (hi - lo), theta0))
            fit <- .fga_least_squares(res_fn, start, lo, hi,
                                      bounded = FALSE,
                                      max_nfev = lm_max_nfev)
            fit$theta <- pmax(lo, pmin(hi, fit$theta))
            fit$cost <- sqrt(sum(res_fn(fit$theta)^2))
        } else {
            start <- pmax(lo, pmin(hi, theta0))
            fit <- .fga_least_squares(res_fn, start, lo, hi,
                                      bounded = TRUE,
                                      max_nfev = lm_max_nfev)
        }
        total_nfev <- total_nfev + as.integer(fit$nfev)
        if (is.null(best) || fit$cost < best$cost) best <- fit
    }

    spectrum <- .fga_shape(best$theta, E) * ln_steps

    params <- as.list(stats::setNames(as.numeric(best$theta[1:7]),
                                      .fga_param_names))
    if (isTRUE(fit_scale)) {
        params$phi_scale <- as.numeric(10^best$theta[8L])
    }
    weight_sum <- sum(best$theta[1:3])
    params$weight_fractions <- if (weight_sum > 0) {
        as.list(stats::setNames(as.numeric(best$theta[1:3] / weight_sum),
                                c("a1", "a2", "a3")))
    } else {
        list(a1 = 0, a2 = 0, a3 = 0)
    }
    params$cost <- best$cost

    list(spectrum = spectrum, iterations = total_nfev,
         converged = isTRUE(best$success), message = best$message,
         params = params)
}

#' Article validation criteria ("Validatsiya rascheta"): per-sphere
#' relative uncertainties, sign alternation, fit-quality FOM and the
#' normalized-spectrum norm check.  Port of _validate_fit().
.fga_validate_fit <- function(computed, measured, spectrum_bins,
                              eps_threshold = 0.05, norm_range = NULL) {
    measured <- as.numeric(measured)
    computed <- as.numeric(computed)
    eps <- ifelse(measured != 0, (computed - measured) / measured, NA_real_)
    finite_eps <- eps[is.finite(eps)]
    max_eps <- if (length(finite_eps)) max(abs(finite_eps)) else Inf
    fom <- if (length(finite_eps)) {
        100 * sqrt(mean(finite_eps^2))
    } else Inf
    signs <- sign(finite_eps)
    sign_changes <- if (length(signs) > 1L) {
        sum(signs[-1] * signs[-length(signs)] < 0)
    } else 0L
    signs_mixed <- any(signs > 0) && any(signs < 0)
    spectrum_norm <- sum(as.numeric(spectrum_bins))
    residuals_ok <- is.finite(max_eps) && max_eps <= eps_threshold
    norm_ok <- if (!is.null(norm_range)) {
        norm_range[1] <= spectrum_norm && spectrum_norm <= norm_range[2]
    } else NULL
    passed <- if (is.null(norm_ok)) residuals_ok else residuals_ok && norm_ok
    list(fom_percent = fom,
         max_relative_uncertainty = max_eps,
         relative_uncertainties = eps,
         residual_sign_changes = as.integer(sign_changes),
         signs_mixed = signs_mixed,
         spectrum_norm = spectrum_norm,
         residuals_ok = residuals_ok,
         norm_ok = norm_ok,
         eps_threshold = as.numeric(eps_threshold),
         passed = passed)
}

#' Unfold with the Fission-model GA + least-squares algorithm
#'
#' Detector-level workflow, R port of \code{unfold_fission_ga} from
#' \code{bssunfold/src/bssunfold/core/unfold_fission_ga.py:466}.  The
#' article's validation criteria are attached to the result as the
#' \code{validation} entry, with the fitted parameters in
#' \code{model_params}.
#'
#' @inheritParams run_unfolding
#' @param ln_steps Optional natural-logarithmic bin widths; default
#'   \code{compute_log_steps(E_MeV) * log(10)}.
#' @param initial_params Optional named list of starting parameter values
#'   (see \code{\link{solve_fission_ga}}).
#' @param fit_scale,ga_popsize,ga_maxiter,ga_tol,lm_method,lm_max_nfev
#'   See \code{\link{solve_fission_ga}}.
#' @param eps_threshold Threshold on the per-sphere relative uncertainty
#'   used by the validation criteria.
#' @param reading_uncertainties,reading_covariance,noise_model,measurement_time
#'   Monte-Carlo options accepted for API parity with the Python original.
#' @return A result list as produced by \code{\link{run_unfolding}} with
#'   extra keys \code{model_params} and \code{validation}.
#' @export
unfold_fission_ga <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116,
                              save_result_callback, readings,
                              ln_steps = NULL,
                              initial_spectrum = NULL,
                              initial_params = NULL,
                              fit_scale = TRUE,
                              ga_popsize = 15L,
                              ga_maxiter = 100L,
                              ga_tol = 1e-10,
                              lm_method = "trf",
                              lm_max_nfev = 2000L,
                              eps_threshold = 0.05,
                              calculate_errors = FALSE,
                              noise_level = 0.01,
                              n_montecarlo = 100L,
                              save_result = FALSE,
                              random_state = NULL,
                              reading_uncertainties = NULL,
                              reading_covariance = NULL,
                              noise_model = "gaussian",
                              measurement_time = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b_clean <- sys$b
    if (is.null(ln_steps)) {
        ln_steps <- compute_log_steps(E_MeV) * log(10)
    }
    ln_steps <- as.numeric(ln_steps)

    holder <- new.env(parent = emptyenv())

    solve_wrapper <- function(A_mat, b_vec, x0 = NULL, ...) {
        out <- solve_fission_ga(A_mat, b_vec, E_MeV, ln_steps,
                                 initial_params = initial_params,
                                 fit_scale = fit_scale,
                                 ga_popsize = ga_popsize,
                                 ga_maxiter = ga_maxiter,
                                 ga_tol = ga_tol,
                                 lm_method = lm_method,
                                 lm_max_nfev = lm_max_nfev,
                                 random_state = random_state)
        # Record the clean-fit artifacts only: Monte-Carlo replicates
        # receive perturbed readings and would overwrite them otherwise.
        if (identical(as.numeric(b_vec), as.numeric(b_clean))) {
            computed <- as.numeric(A_mat %*% out$spectrum)
            norm_range <- if (!isTRUE(fit_scale)) c(0.6, 1.2) else NULL
            assign("model_params", out$params, envir = holder)
            assign("validation",
                   .fga_validate_fit(computed, b_vec, out$spectrum,
                                     eps_threshold = eps_threshold,
                                     norm_range = norm_range),
                   envir = holder)
        }
        out
    }

    x0_default <- rep(1, n_energy_bins) *
        (mean(b_clean) / mean(rowSums(A)))

    result <- run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solve_wrapper,
        solve_kwargs = list(),
        method_name = "fission_ga",
        extra_output = list(initial_params = initial_params,
                            fit_scale = fit_scale,
                            ga_popsize = ga_popsize,
                            ga_maxiter = ga_maxiter,
                            lm_method = lm_method,
                            T0 = .fga_T0, Ed = .fga_Ed),
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)

    if (!is.null(holder$model_params)) result$model_params <-
        holder$model_params
    if (!is.null(holder$validation)) result$validation <-
        holder$validation
    result
}
