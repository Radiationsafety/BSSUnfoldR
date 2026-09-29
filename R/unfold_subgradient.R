#' Projected subgradient descent unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_subgradient.py}.
#' Solves the nonsmooth unfolding problem
#' \deqn{\min_{x\ge0} \frac12\|Ax-b\|^2 + \lambda_1\|x\|_1 + \lambda_{tv}\|Dx\|_1}
#' by projected subgradient descent with a subgradient
#' \deqn{g = A^T(Ax-b) + \lambda_1\,\mathrm{sign}(x) + \lambda_{tv}\,D^T\mathrm{sign}(Dx).}
#' Step-size policies are \code{"polyak"}, \code{"diminishing"} (default) or
#' \code{"fixed"}.  Because the iterate sequence need not converge, the best
#' iterate by objective value is returned.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Initial spectrum (length n); projected onto the nonnegative orthant.
#' @param max_iterations Positive integer; default 3000L.
#' @param tolerance Positive numeric relative change tolerance; default 1e-8.
#' @param l1_penalty Numeric L1 (sparsity) penalty weight; default 0.
#' @param tv_penalty Numeric total-variation penalty weight; default 0.
#' @param step_policy Character; \code{"polyak"}, \code{"diminishing"} or
#'   \code{"fixed"}; default \code{"diminishing"}.
#' @param step_size Numeric base step for the fixed/diminishing policies; default 1.
#' @param decay Numeric decay rate of the diminishing step; default 1.
#' @param polyak_margin Numeric relative shrink of the running-best objective
#'   used as the \eqn{f^*} estimate in the Polyak rule; default 0.05.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
solve_subgradient <- function(A, b, x0, max_iterations = 3000L, tolerance = 1e-8,
                              l1_penalty = 0.0, tv_penalty = 0.0,
                              step_policy = "diminishing", step_size = 1.0,
                              decay = 1.0, polyak_margin = 0.05) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)

    steps <- c("polyak", "diminishing", "fixed")
    if (!step_policy %in% steps) {
        stop("step_policy must be one of polyak, diminishing, fixed")
    }

    n <- ncol(A)
    l1_penalty <- max(as.numeric(l1_penalty), 0.0)
    tv_penalty <- max(as.numeric(tv_penalty), 0.0)
    D <- if (tv_penalty > 0) .sg_difference_matrix(n) else NULL

    objective <- function(z) {
        r <- as.numeric(A %*% z) - b
        val <- 0.5 * sum(r * r)
        if (l1_penalty) val <- val + l1_penalty * sum(abs(z))
        if (!is.null(D)) val <- val + tv_penalty * sum(abs(as.numeric(D %*% z)))
        val
    }

    subgradient <- function(z) {
        g <- as.numeric(t(A) %*% (as.numeric(A %*% z) - b))
        if (l1_penalty) g <- g + l1_penalty * sign(z)
        if (!is.null(D)) g <- g + tv_penalty * as.numeric(t(D) %*% sign(as.numeric(D %*% z)))
        g
    }

    # Scale-aware base step: t0 * ||g(x0)|| ~ ||x||_ref makes one step move the
    # iterate by a reference solution magnitude (||b|| / ||A||_2).
    x_scale <- sqrt(sum(b^2)) / max(.bss_spectral_norm(A), 1e-30)
    g0 <- subgradient(pmax(as.numeric(x0), 0.0))
    g0_norm <- sqrt(sum(g0^2))
    if (is.null(step_size)) step_size <- 1.0
    t0 <- as.numeric(step_size) * max(x_scale, 1e-300) / max(g0_norm, 1e-300)

    x <- pmax(as.numeric(x0), 0.0)
    best_x <- x
    best_f <- objective(x)
    f_star_est <- best_f * max(1.0 - polyak_margin, 0.0)

    converged <- FALSE
    iterations <- 0L
    for (k in seq_len(max_iterations)) {
        kk <- k - 1L                            # Python's 0-based loop index
        g <- subgradient(x)
        g_norm_sq <- sum(g * g)

        if (step_policy == "polyak") {
            f_cur <- objective(x)
            if (f_cur < best_f) {
                best_f <- f_cur
                best_x <- x
                f_star_est <- best_f * max(1.0 - polyak_margin, 0.0)
            }
            t <- if (g_norm_sq > 0) max((f_cur - f_star_est) / g_norm_sq, 0.0) else 0.0
        } else if (step_policy == "diminishing") {
            t <- t0 / (1.0 + decay * kk)
        } else {
            t <- t0
        }

        x_new <- pmax(x - t * g, 0.0)
        rel_change <- sqrt(sum((x_new - x)^2)) / max(sqrt(sum(x^2)), 1e-30)
        x <- x_new
        iterations <- k

        f_new <- objective(x)
        if (f_new < best_f) {
            best_f <- f_new
            best_x <- x
        }

        if (rel_change < tolerance) {
            converged <- TRUE
            break
        }
    }

    list(spectrum = as.numeric(best_x), iterations = iterations,
         converged = converged)
}

# First-order difference operator D with (D x)_i = x_{i+1} - x_i.
.sg_difference_matrix <- function(n) {
    if (n < 2L) return(matrix(numeric(0), nrow = 0L, ncol = n))
    D <- matrix(0.0, nrow = n - 1L, ncol = n)
    idx <- seq_len(n - 1L)
    D[cbind(idx, idx)] <- -1.0
    D[cbind(idx, idx + 1L)] <- 1.0
    D
}

#' @rdname solve_subgradient
#' @inheritParams run_unfolding
#' @param initial_spectrum Optional initial spectrum guess.
#' @param step_policy Character; \code{"polyak"}, \code{"diminishing"} or \code{"fixed"}.
#' @param step_size Numeric base step size.
#' @param decay Numeric diminishing-step decay rate.
#' @param polyak_margin Numeric Polyak optimal-value margin.
#' @param calculate_errors Logical; run Monte-Carlo uncertainty.
#' @param noise_level Numeric; MC noise level.
#' @param n_montecarlo Integer; MC samples.
#' @param save_result Logical; call the save callback.
#' @param random_state Optional integer seed.
#' @param max_neutron_energy Optional energy cutoff in MeV.
#' @return A result list produced by \code{\link{run_unfolding}}.
#' @export
unfold_subgradient <- function(detector_names, n_energy_bins, E_MeV,
                               sensitivities, cc_icrp116,
                               save_result_callback, readings,
                               initial_spectrum = NULL,
                               max_iterations = 3000L, tolerance = 1e-8,
                               l1_penalty = 0.0, tv_penalty = 0.0,
                               step_policy = "diminishing", step_size = 1.0,
                               decay = 1.0, polyak_margin = 0.05,
                               calculate_errors = FALSE, noise_level = 0.01,
                               n_montecarlo = 100L, save_result = FALSE,
                               random_state = NULL,
                               max_neutron_energy = NULL) {
    steps <- c("polyak", "diminishing", "fixed")
    if (!step_policy %in% steps) {
        stop("step_policy must be one of polyak, diminishing, fixed")
    }

    # scale = mean(all readings) / mean(mean of each sensitivity vector)
    sens_mean <- mean(vapply(sensitivities, function(v) mean(as.numeric(v)),
                             numeric(1)))
    scale <- max(mean(as.numeric(readings)), 1e-30) / max(sens_mean, 1e-30)
    x0_default <- rep(scale / max(n_energy_bins, 1), n_energy_bins)

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_subgradient,
                                        max_iterations = max_iterations,
                                        tolerance = tolerance,
                                        l1_penalty = l1_penalty,
                                        tv_penalty = tv_penalty,
                                        step_policy = step_policy,
                                        step_size = step_size,
                                        decay = decay,
                                        polyak_margin = polyak_margin),
        solve_kwargs = list(),
        method_name = "Subgradient Descent",
        extra_output = list(l1_penalty = l1_penalty,
                            tv_penalty = tv_penalty,
                            step_policy = step_policy),
        calculate_errors = calculate_errors,
        noise_level = noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
