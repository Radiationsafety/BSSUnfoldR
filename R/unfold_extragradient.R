#' Extragradient (Korpelevich) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_extragradient.py}.
#' Solves the robust unfolding problem
#' \deqn{\min_{x\ge0} \frac12\|Ax-b\|^2 + \delta\|Ax-b\|_2}
#' (with \eqn{\delta = \text{noise\_level}\,\|b\|_2}) through its bilinear
#' saddle form using Korpelevich's two-step extragradient scheme.  The
#' prediction and correction steps each use the gradient
#' \eqn{F_x(x, y) = A^T(Ax - b) + \delta A^T y}{Fx(x,y) = A'(Ax-b) + delta A'y},
#' with \eqn{x} projected onto the nonnegative orthant and \eqn{y} onto the unit
#' L2 ball.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Initial spectrum (length n); projected onto the nonnegative orthant.
#' @param max_iterations Positive integer; default 2000L.
#' @param tolerance Positive numeric relative change tolerance; default 1e-8.
#' @param noise_level Numeric relative noise-ball radius, \eqn{\delta = \text{noise\_level}\|b\|_2}; default 0.02.
#' @param step_size Optional numeric extragradient step \eqn{\eta}; when
#'   \code{NULL} (default) set to \code{0.9 / L} with
#'   \code{L = ||A||_2^2 + delta*||A||_2 + 1}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
solve_extragradient <- function(A, b, x0, max_iterations = 2000L,
                                tolerance = 1e-8, noise_level = 0.02,
                                step_size = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A)

    noise_level <- max(as.numeric(noise_level), 0.0)
    delta <- noise_level * sqrt(sum(b^2))

    x <- pmax(as.numeric(x0), 0.0)
    y <- numeric(m)                      # dual noise direction

    norm_A <- .bss_spectral_norm(A)       # ||A||_2
    eta <- if (is.null(step_size)) {
        L <- norm_A^2 + delta * norm_A + 1.0
        0.9 / L
    } else {
        as.numeric(step_size)
    }

    F_x <- function(xv, yv) {
        as.numeric(t(A) %*% (as.numeric(A %*% xv) - b)) +
            delta * as.numeric(t(A) %*% yv)
    }

    converged <- FALSE
    iterations <- 0L
    for (k in seq_len(max_iterations)) {
        residual <- as.numeric(A %*% x) - b

        # ---- prediction step ------------------------------------------------
        x_tilde <- pmax(x - eta * F_x(x, y), 0.0)
        y_tilde <- .xg_project_ball(y + eta * residual, 1.0)

        # ---- correction step (extra gradient evaluation) --------------------
        residual_t <- as.numeric(A %*% x_tilde) - b
        x_new <- pmax(x - eta * F_x(x_tilde, y_tilde), 0.0)
        y_new <- .xg_project_ball(y + eta * residual_t, 1.0)

        rel_change <- sqrt(sum((x_new - x)^2)) / max(sqrt(sum(x^2)), 1e-30)
        x <- x_new; y <- y_new
        iterations <- k
        if (rel_change < tolerance) {
            converged <- TRUE
            break
        }
    }

    list(spectrum = as.numeric(x), iterations = iterations, converged = converged)
}

# Euclidean projection onto the L2 ball of the given radius.
.xg_project_ball <- function(v, radius) {
    nv <- sqrt(sum(v^2))
    if (nv > radius) v * (radius / max(nv, 1e-300)) else v
}

#' @rdname solve_extragradient
#' @inheritParams run_unfolding
#' @param initial_spectrum Optional initial spectrum guess.
#' @param noise_level Numeric relative noise-ball radius (default 0.02).
#' @param step_size Optional extragradient step.
#' @param mc_noise_level Numeric; Monte-Carlo noise level (default 0.01).
#' @param calculate_errors Logical; run Monte-Carlo uncertainty.
#' @param n_montecarlo Integer; MC samples.
#' @param save_result Logical; call the save callback.
#' @param random_state Optional integer seed.
#' @param max_neutron_energy Optional energy cutoff in MeV.
#' @return A result list produced by \code{\link{run_unfolding}}.
#' @export
unfold_extragradient <- function(detector_names, n_energy_bins, E_MeV,
                                 sensitivities, cc_icrp116,
                                 save_result_callback, readings,
                                 initial_spectrum = NULL,
                                 max_iterations = 2000L, tolerance = 1e-8,
                                 noise_level = 0.02, step_size = NULL,
                                 calculate_errors = FALSE,
                                 mc_noise_level = 0.01, n_montecarlo = 100L,
                                 save_result = FALSE, random_state = NULL,
                                 max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_extragradient,
                                        max_iterations = max_iterations,
                                        tolerance = tolerance,
                                        noise_level = noise_level,
                                        step_size = step_size),
        solve_kwargs = list(),
        method_name = "Extragradient",
        extra_output = list(noise_ball_radius = noise_level),
        calculate_errors = calculate_errors,
        noise_level = mc_noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
