#' Frank--Wolfe (conditional gradient) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_frank_wolfe.py}.
#' Solves the fluence-constrained unfolding problem
#' \deqn{\min_x \frac12\|Ax-b\|^2 \quad \text{s.t. } x \ge 0,\ \text{sum}(x) = F}
#' by the Frank--Wolfe algorithm.  Every iteration linearizes the objective and
#' solves the linear minimization oracle over the simplex (picking the "vertex"
#' bin descent wants most), keeping the iterate a convex combination of simplex
#' vertices so the total fluence is preserved exactly.  Optional Wolfe
#' away-steps accelerate local convergence.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Initial spectrum (length n); projected onto the simplex.
#' @param total_fluence Positive numeric simplex level \eqn{F}.
#' @param max_iterations Positive integer; default 1000L.
#' @param tolerance Positive numeric Frank--Wolfe gap tolerance; default 1e-8.
#' @param away_steps Logical; enable Wolfe away-steps; default \code{TRUE}.
#' @param line_search Character; \code{"exact"} or \code{"backtracking"};
#'   default \code{"exact"}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
solve_frank_wolfe <- function(A, b, x0, total_fluence, max_iterations = 1000L,
                              tolerance = 1e-8, away_steps = TRUE,
                              line_search = "exact") {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)

    if (as.numeric(total_fluence) <= 0) stop("total_fluence must be positive")
    F <- as.numeric(total_fluence)

    x <- pmax(as.numeric(x0), 0.0)
    s <- sum(x)
    x <- if (s > 0) F * x / s else rep(F / n, n)

    ATb <- as.numeric(t(A) %*% b)
    G <- crossprod(A)                     # n x n Gram matrix, precomputed once

    converged <- FALSE
    iterations <- 0L
    for (k in seq_len(max_iterations)) {
        grad <- as.numeric(G %*% x) - ATb

        # Frank-Wolfe (duality) gap: <grad, x - s> >= f(x) - f*
        s_v <- .fw_lmo_simplex(grad, F)
        fw_gap <- sum(grad * (x - s_v))
        iterations <- k
        if (fw_gap <= tolerance * max(1.0, abs(sum(grad * x)))) {
            converged <- TRUE
            break
        }

        # Regular FW step: move from x towards the LMO vertex
        d <- s_v - x
        gamma_max <- 1.0

        # Away-step candidate
        if (isTRUE(away_steps)) {
            mask <- x > 0
            if (any(mask) && sum(mask) > 1) {
                masked <- ifelse(mask, grad, -Inf)
                i_away <- which.max(masked)
                x_rest <- x
                x_rest[i_away] <- 0.0
                x_rest <- F * x_rest / sum(x_rest)
                away_gap <- sum(grad * (x - x_rest))
                if (away_gap > fw_gap) {
                    d <- x_rest - x
                    gamma_max <- 1.0
                }
            }
        }

        if (line_search == "exact") {
            g_d <- sum(grad * d)
            gd <- sum(d * as.numeric(G %*% d))
            gamma <- if (gd <= 0) gamma_max else min(max(-g_d / gd, 0.0), gamma_max)
        } else {
            phi <- function(t) {
                r <- as.numeric(A %*% (x + t * d)) - b
                0.5 * sum(r * r)
            }
            res <- .bss_golden_section(phi, 0.0, gamma_max)
            gamma <- res$t_opt
        }

        if (gamma <= 1e-16) {
            converged <- TRUE
            break
        }
        x <- pmax(x + gamma * d, 0.0)
        x <- F * x / sum(x)
    }

    list(spectrum = as.numeric(x), iterations = iterations, converged = converged)
}

# Linear minimization oracle over the simplex: put all mass on argmin.
.fw_lmo_simplex <- function(gradient, total_fluence) {
    s <- numeric(length(gradient))
    s[which.min(gradient)] <- total_fluence
    s
}

# Total-fluence estimate from an unconstrained NNLS fit, mirroring
# core/_matrix_utils.py::estimate_total_fluence (with ln_steps = NULL, so the
# returned value is sum(x_nnls), the lethargy integral).
.bss_estimate_total_fluence <- function(A, b) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    total <- tryCatch({
        x_nnls <- .bss_nnls(A, b)
        sum(x_nnls)
    }, error = function(e) 0.0)
    if (is.finite(total) && total > 0.0) return(total)
    mean_response <- max(mean(A), 1e-30)
    mean(b) / mean_response * n
}

#' @rdname solve_frank_wolfe
#' @inheritParams run_unfolding
#' @param initial_spectrum Optional initial spectrum guess.
#' @param total_fluence Optional simplex level \eqn{F}; estimated from a
#'   uniform NNLS fit when \code{NULL}.
#' @param away_steps Logical; use Wolfe away-steps.
#' @param line_search Character; \code{"exact"} or \code{"backtracking"}.
#' @param calculate_errors Logical; run Monte-Carlo uncertainty.
#' @param noise_level Numeric; MC noise level.
#' @param n_montecarlo Integer; MC samples.
#' @param save_result Logical; call the save callback.
#' @param random_state Optional integer seed.
#' @param max_neutron_energy Optional energy cutoff in MeV.
#' @return A result list produced by \code{\link{run_unfolding}}.
#' @export
unfold_frank_wolfe <- function(detector_names, n_energy_bins, E_MeV,
                               sensitivities, cc_icrp116,
                               save_result_callback, readings,
                               initial_spectrum = NULL, total_fluence = NULL,
                               max_iterations = 1000L, tolerance = 1e-8,
                               away_steps = TRUE, line_search = "exact",
                               calculate_errors = FALSE, noise_level = 0.01,
                               n_montecarlo = 100L, save_result = FALSE,
                               random_state = NULL,
                               max_neutron_energy = NULL) {
    if (is.null(total_fluence)) {
        est <- .bss_select_system(detector_names, readings, sensitivities)
        total_fluence <- .bss_estimate_total_fluence(est$A, est$b)
    }
    x0_default <- rep(as.numeric(total_fluence) / n_energy_bins, n_energy_bins)

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_frank_wolfe,
                                        total_fluence = total_fluence,
                                        max_iterations = max_iterations,
                                        tolerance = tolerance,
                                        away_steps = away_steps,
                                        line_search = line_search),
        solve_kwargs = list(),
        method_name = "Frank-Wolfe",
        extra_output = list(total_fluence = total_fluence,
                            away_steps = away_steps,
                            line_search = line_search),
        calculate_errors = calculate_errors,
        noise_level = noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}

# Build the (A, b) system from the FULL (weighted) sensitivities and the
# selected readings, mirroring run_unfolding's system construction so the
# fluence estimate matches Python's estimate_total_fluence(A_est, b_est).
.bss_select_system <- function(detector_names, readings, sensitivities) {
    selected <- detector_names[detector_names %in% names(readings)]
    b <- as.numeric(readings[selected])
    A <- do.call(rbind, lapply(selected, function(nm) as.numeric(sensitivities[[nm]])))
    list(A = A, b = b)
}
