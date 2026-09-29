#' Mirror descent unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_mirror_descent.py}.
#' Mirror descent replaces the Euclidean geometry of gradient descent with a
#' Bregman geometry induced by a mirror map \eqn{\psi}.  The \emph{entropy}
#' map produces multiplicative updates renormalized onto the simplex
#' \eqn{\{x \ge 0, \text{sum}(x) = F\}} (generalizing MLEM / GRAVEL / SAND-II);
#' the \emph{log} (log-barrier), \emph{l2} (projected gradient) and \emph{pnorm}
#' maps give other nonnegative geometries.  When \code{step_size} is \code{NULL}
#' (default) the step is selected each iteration by golden-section line search
#' along the mirror trajectory.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Initial spectrum (length n); must be strictly positive for the
#'   \code{"entropy"} and \code{"log"} maps.
#' @param max_iterations Positive integer; default 1000L.
#' @param tolerance Positive numeric relative change tolerance; default 1e-8.
#' @param mirror_map Character; \code{"entropy"}, \code{"log"}, \code{"l2"} or
#'   \code{"pnorm"}; default \code{"entropy"}.
#' @param step_size Optional numeric mirror step; \code{NULL} (default) uses a
#'   golden-section line search each iteration.
#' @param total_fluence Optional simplex level for the entropy map; defaults to
#'   \code{sum(x0)}.
#' @param regularization Numeric Tikhonov (L2) strength; default 0.
#' @param p Numeric order of the p-norm mirror map (\code{p > 1}); default 3.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @keywords internal
solve_mirror_descent <- function(A, b, x0, max_iterations = 1000L,
                                 tolerance = 1e-8, mirror_map = "entropy",
                                 step_size = NULL, total_fluence = NULL,
                                 regularization = 0.0, p = 3.0) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)

    mirror_maps <- c("entropy", "log", "l2", "pnorm")
    if (!mirror_map %in% mirror_maps) {
        stop("mirror_map must be one of entropy, log, l2, pnorm")
    }
    if (mirror_map == "pnorm" && p <= 1.0) {
        stop("p-norm mirror map requires p > 1")
    }

    x0 <- as.numeric(x0)
    if (mirror_map == "entropy") {
        F <- if (!is.null(total_fluence)) as.numeric(total_fluence) else sum(x0)
        if (F <= 0) stop("entropy mirror map requires positive total fluence")
        x <- pmax(x0, 1e-300)
        x <- F * x / sum(x)
    } else {
        F <- NULL
        max_A <- max(max(A), 1e-30)
        floor_val <- max(max(b) / max_A, 1e-12) / max(ncol(A), 1)
        x <- pmax(x0, floor_val)
    }

    L <- .bss_spectral_norm(A)^2 + max(as.numeric(regularization), 0.0)

    objective <- function(z) {
        r <- as.numeric(A %*% z) - b
        val <- 0.5 * sum(r * r)
        if (regularization) val <- val + 0.5 * regularization * sum(z * z)
        val
    }

    converged <- FALSE
    iterations <- 0L
    for (k in seq_len(max_iterations)) {
        gradient <- as.numeric(t(A) %*% (as.numeric(A %*% x) - b))
        if (regularization) gradient <- gradient + regularization * x

        if (is.null(step_size)) {
            hi <- .md_eta_max(gradient, L, mirror_map, x)
            if (hi <= 0) break
            phi <- function(tt) objective(.md_mirror_next(x, gradient, tt,
                                                          mirror_map, F, p))
            eta_k <- .bss_golden_section(phi, 0.0, hi, tolerance = hi * 1e-4)$t_opt
        } else {
            eta_k <- as.numeric(step_size)
        }

        x_new <- .md_mirror_next(x, gradient, eta_k, mirror_map, F, p)

        rel_change <- sqrt(sum((x_new - x)^2)) / max(sqrt(sum(x^2)), 1e-300)
        x <- x_new
        iterations <- k
        if (rel_change < tolerance) {
            converged <- TRUE
            break
        }
    }

    list(spectrum = as.numeric(x), iterations = iterations, converged = converged)
}

# Mirror-descent proximal map applied with step `eta`.
.md_mirror_next <- function(x, gradient, eta, mirror_map, total_fluence, p) {
    if (mirror_map == "entropy") {
        log_step <- -eta * gradient
        shift <- max(log_step)                 # overflow guard; cancels in norm
        x_new <- x * exp(log_step - shift)
        s <- sum(x_new)
        if (!is.finite(s) || s <= 0) return(x)
        return(total_fluence * x_new / s)
    }
    if (mirror_map == "log") {
        inv <- 1.0 / pmax(x, 1e-300) + eta * gradient
        if (!all(is.finite(inv)) || any(inv <= 0)) return(x)
        return(1.0 / inv)
    }
    if (mirror_map == "l2") {
        return(pmax(x - eta * gradient, 0.0))
    }
    if (mirror_map == "pnorm") {
        xp <- pmax(x, 1e-300)^(p - 1.0) - eta * gradient
        return(pmax(xp, 0.0)^(1.0 / (p - 1.0)))
    }
    stop("Unknown mirror map '", mirror_map, "'")
}

# Upper bracket for the per-iteration line search.
.md_eta_max <- function(gradient, L, mirror_map, x) {
    g_max <- max(abs(gradient))
    if (mirror_map %in% c("entropy", "log")) {
        hi <- 4.0 / max(g_max, 1e-300)         # scale-free multiplicative bracket
        if (mirror_map == "log") {
            neg <- gradient < 0
            if (any(neg)) {
                inv_x <- 1.0 / pmax(x, 1e-300)
                limit <- min(inv_x[neg] / abs(gradient[neg]))
                hi <- min(hi, 0.25 * limit)
            }
        }
        return(hi)
    }
    if (L > 0) 2.0 / L else 1.0
}

# Golden-section minimization over [lo, hi]; port of
# core/_line_search.py::golden_section_minimize.  Returns list(t_opt, f_opt).
.bss_golden_section <- function(func, lo, hi, tolerance = 1e-8,
                                max_iterations = 200L) {
    if (!is.finite(lo) || !is.finite(hi)) {
        stop("golden_section_minimize requires finite bounds")
    }
    if (hi <= lo) return(list(t_opt = lo, f_opt = as.numeric(func(lo))))

    inv_phi <- (sqrt(5.0) - 1.0) / 2.0
    a <- lo; bb <- hi
    cc <- bb - inv_phi * (bb - a)
    dd <- a + inv_phi * (bb - a)
    fc <- as.numeric(func(cc))
    fd <- as.numeric(func(dd))

    for (i in seq_len(max_iterations)) {
        if (bb - a <= tolerance) break
        if (fc < fd) {
            bb <- dd; dd <- cc; fd <- fc
            cc <- bb - inv_phi * (bb - a)
            fc <- as.numeric(func(cc))
        } else {
            a <- cc; cc <- dd; fc <- fd
            dd <- a + inv_phi * (bb - a)
            fd <- as.numeric(func(dd))
        }
    }
    t_opt <- 0.5 * (a + bb)
    list(t_opt = t_opt, f_opt = as.numeric(func(t_opt)))
}

#' @rdname solve_mirror_descent
#' @inheritParams run_unfolding
#' @param initial_spectrum Optional initial spectrum guess.
#' @param mirror_map Character; mirror map name.
#' @param step_size Optional fixed mirror step.
#' @param total_fluence Optional simplex level for the entropy map.
#' @param p Numeric p-norm order.
#' @param calculate_errors Logical; run Monte-Carlo uncertainty.
#' @param noise_level Numeric; MC noise level.
#' @param n_montecarlo Integer; MC samples.
#' @param save_result Logical; call the save callback.
#' @param random_state Optional integer seed.
#' @param max_neutron_energy Optional energy cutoff in MeV.
#' @return A result list produced by \code{\link{run_unfolding}}.
#' @export
unfold_mirror_descent <- function(detector_names, n_energy_bins, E_MeV,
                                  sensitivities, cc_icrp116,
                                  save_result_callback, readings,
                                  initial_spectrum = NULL,
                                  max_iterations = 2000L, tolerance = 1e-8,
                                  mirror_map = "entropy", step_size = NULL,
                                  total_fluence = NULL, regularization = 0.0,
                                  p = 3.0, calculate_errors = FALSE,
                                  noise_level = 0.01, n_montecarlo = 100L,
                                  save_result = FALSE, random_state = NULL,
                                  max_neutron_energy = NULL) {
    mirror_maps <- c("entropy", "log", "l2", "pnorm")
    if (!mirror_map %in% mirror_maps) {
        stop("mirror_map must be one of entropy, log, l2, pnorm")
    }

    est <- .bss_select_system(detector_names, readings, sensitivities)
    A_est <- est$A; b_est <- est$b

    if (mirror_map == "entropy") {
        if (is.null(total_fluence)) {
            total_fluence <- .bss_estimate_total_fluence(A_est, b_est)
        }
        x0_default <- rep(as.numeric(total_fluence) / n_energy_bins, n_energy_bins)
    } else {
        scale <- max(mean(b_est), 1e-30) / max(mean(A_est), 1e-30)
        x0_default <- rep(scale / max(n_energy_bins, 1), n_energy_bins)
    }

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_mirror_descent,
                                        max_iterations = max_iterations,
                                        tolerance = tolerance,
                                        mirror_map = mirror_map,
                                        step_size = step_size,
                                        total_fluence = total_fluence,
                                        regularization = regularization,
                                        p = p),
        solve_kwargs = list(),
        method_name = "Mirror Descent",
        extra_output = list(mirror_map = mirror_map,
                            total_fluence = if (mirror_map == "entropy")
                                total_fluence else NULL),
        calculate_errors = calculate_errors,
        noise_level = noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
