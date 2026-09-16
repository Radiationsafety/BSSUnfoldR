#' Mystic hybrid two-stage (global + local) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_mystic_hybrid.py}.
#' Stage 1 performs global exploration of the penalized least-squares
#' objective in log-space with differential evolution (\code{diffev2}
#' analogue); stage 2 refines the best individual with a bounded
#' \code{"L-BFGS-B"} (\code{fmin_powell} analogue) descent.
#'
#' @name mystic-hybrid
NULL

# classical rand/1/bin differential evolution (internal helper)
.mystic_de <- function(fun, n, lower, upper, pop_size, maxiter) {
    lower <- rep(lower, length.out = n)
    upper <- rep(upper, length.out = n)
    pop <- matrix(runif(pop_size * n, lower, upper), ncol = n)
    fit <- apply(pop, 1, fun)
    best_i <- which.min(fit)
    best <- pop[best_i, ]
    best_f <- fit[best_i]
    for (it in seq_len(maxiter)) {
        for (i in seq_len(pop_size)) {
            others <- sample.int(pop_size, 3L)
            others <- others[others != i][1:min(3, length(others))]
            if (length(others) < 3L) next
            a <- others[1]; b <- others[2]; c <- others[3]
            trial <- pop[a, ] + 0.8 * (pop[b, ] - pop[c, ])
            mask <- runif(n) < 0.7
            j <- sample.int(n, 1)
            mask[j] <- TRUE
            candidate <- pop[i, ]
            candidate[mask] <- trial[mask]
            candidate_f <- tryCatch(fun(candidate), error = function(e) Inf)
            if (is.finite(candidate_f) && candidate_f < fit[i]) {
                pop[i, ] <- candidate
                fit[i] <- candidate_f
                if (candidate_f < best_f) {
                    best_f <- candidate_f
                    best <- candidate
                }
            }
        }
    }
    list(par = best, value = best_f)
}

#' Solve by mystic hybrid (differential-evolution global + L-BFGS-B local)
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (warm start).
#' @param global_maxiter Integer DE generations. Default 50.
#' @param global_maxfun Integer; upper bound on global stage evaluations
#'   (Python \code{global_maxfun}). Default 2000.
#' @param local_maxiter Integer L-BFGS-B iterations. Default 200.
#' @param npop Integer population size. Default \code{max(10F, 4n)}.
#' @param regularization Numeric Tikhonov weight. Default 1e-4.
#' @param smoothness Numeric second-difference weight. Default 1e-4.
#' @param random_state Optional seed.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' set.seed(3)
#' r <- solve_mystic_hybrid(A, b, rep(1, 3), global_maxiter = 10L)
solve_mystic_hybrid <- function(A, b, x0 = NULL,
                                global_maxiter = 50L,
                                global_maxfun = 2000L,
                                local_maxiter = 200L,
                                npop = NULL,
                                regularization = 1e-4, smoothness = 1e-4,
                                random_state = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (!is.null(random_state)) set.seed(as.integer(random_state))
    if (is.null(x0)) x0 <- rep(mean(b) / max(sum(A), 1e-10) * n, n)
    scale <- max(sum(x0), 1e-30)
    .objective <- function(y) {
        x <- pmax(exp(y) * scale, 0)
        resid <- as.numeric(A %*% x) - b
        r <- sum(resid^2) / max(sum(b^2), 1e-30)
        tik <- regularization * sum(x^2) / scale^2
        sm <- if (n > 2 && smoothness > 0)
                  smoothness * sum(diff(diff(x))^2) / scale^2
              else 0
        r + tik + sm
    }
    # DE global stage on log-spectrum in bounded box (warm start bounds the
    # box around zero log-scale deviation like the Python diffev2 setup).
    npop <- max(as.integer(npop), 4L * min(n, 50L), 10L)
    iters <- min(as.integer(global_maxiter),
                 max(as.integer(global_maxfun) %/% max(npop, 1L), 1L))
    de <- .mystic_de(.objective, n, rep(-6.0, n), rep(3.0, n), npop,
                     max(iters, 1L))
    # Local refinement
    res_local <- tryCatch(stats::optim(par = de$par, fn = .objective,
                                       method = "L-BFGS-B",
                                       control = list(maxit =
                                           as.integer(local_maxiter))),
                          error = function(e) list(par = de$par,
                                                   value = de$value))
    spectrum <- pmax(exp(res_local$par) * scale, 0)
    list(spectrum = spectrum,
         iterations = as.integer(iters + res_local$counts[[1]]),
         converged = TRUE)
}

#' Wrapper around \code{\link{solve_mystic_hybrid}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_mystic_hybrid
#' @export
unfold_mystic_hybrid <- function(detector_names, n_energy_bins, E_MeV,
                                 sensitivities, cc_icrp116,
                                 save_result_callback, readings,
                                 initial_spectrum = NULL,
                                 global_maxiter = 50L, global_maxfun = 2000L,
                                 local_maxiter = 200L, npop = NULL,
                                 regularization = 1e-4, smoothness = 1e-4,
                                 method_name = "MysticHybrid",
                                 calculate_errors = FALSE, noise_level = 0.01,
                                 n_montecarlo = 100L, save_result = FALSE,
                                 random_state = NULL,
                              max_neutron_energy = NULL) {
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1, n_energy_bins),
        solve_func = solve_mystic_hybrid,
        solve_kwargs = list(global_maxiter = global_maxiter,
                            global_maxfun = global_maxfun,
                            local_maxiter = local_maxiter, npop = npop,
                            regularization = regularization,
                            smoothness = smoothness),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
