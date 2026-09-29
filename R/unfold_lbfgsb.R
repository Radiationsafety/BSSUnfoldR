#' L-BFGS-B (quasi-Newton) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_lbfgsb.py}. Minimises
#' the smooth Tikhonov-type objective
#' \deqn{\frac12 \|A x - b\|^2 + \frac{\lambda}{2}\|x\|^2
#'       + \frac{\mu}{2}\|D_2 x\|^2}{0.5*||Ax-b||^2 + reg/2*||x||^2 + smooth/2*||D2 x||^2}
#' over the box \eqn{x_{min} \le x \le x_{max}{x_min <= x <= x_max}} with the
#' bound-constrained limited-memory BFGS method (Byrd, Lu, Nocedal & Zhu, 1995)
#' of \code{\link[stats]{optim}}, which is the same Fortran implementation that
#' \code{scipy.optimize.minimize(method = "L-BFGS-B")} drives. Analytic
#' gradients are supplied, so \code{tolerance} maps to scipy's \code{gtol}
#' (\code{pgtol}) and \code{lbfgs_history} to \code{maxcor} (\code{lmm}).
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial guess (length n). Default \code{NULL} = zeros.
#' @param max_iterations Positive integer; maximum iterations
#'   (scipy \code{maxiter}). Default 500.
#' @param tolerance Positive numeric; projected-gradient stopping tolerance
#'   (scipy \code{gtol}). Default 1e-8.
#' @param regularization Numeric Tikhonov (L2) strength. Default 0.0.
#' @param smoothness Numeric second-difference (curvature) penalty. Default 0.0.
#' @param x_min Numeric lower box bound. Default 0.0.
#' @param x_max Numeric upper box bound. Default \code{Inf}.
#' @param lbfgs_history Integer; number of correction pairs (scipy
#'   \code{maxcor}). Default 10.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_lbfgsb(A, b, rep(0, 3), max_iterations = 100L)
solve_lbfgsb <- function(A, b, x0 = NULL, max_iterations = 500L,
                         tolerance = 1e-8, regularization = 0.0,
                         smoothness = 0.0, x_min = 0.0, x_max = Inf,
                         lbfgs_history = 10L) {
    v <- validate_system(A, b, x0 = x0, max_iterations = max_iterations,
                         tolerance = tolerance)
    A <- v$A; b <- v$b
    n <- ncol(A)
    if (is.null(x0)) x0 <- rep(0.0, n)
    regularization <- max(as.numeric(regularization), 0.0)
    smoothness <- max(as.numeric(smoothness), 0.0)
    D2 <- if (smoothness > 0) second_difference_matrix(n) else NULL
    AT <- t(A)

    ## scipy's ``jac = True`` callback returns (f, g); both are needed on every
    ## evaluation, so they are computed together and cached per iterate to keep
    ## R's separate fn/gr calls to one matrix-product pair per point.
    cache <- new.env(parent = emptyenv())
    .fg <- function(x) {
        if (!is.null(cache$x) && identical(cache$x, x)) return(cache$v)
        r <- as.numeric(A %*% x) - b
        f <- 0.5 * sum(r * r)
        g <- as.numeric(AT %*% r)
        if (regularization > 0) {
            f <- f + 0.5 * regularization * sum(x * x)
            g <- g + regularization * x
        }
        if (!is.null(D2)) {
            Dx <- as.numeric(D2 %*% x)
            f <- f + 0.5 * smoothness * sum(Dx * Dx)
            g <- g + smoothness * as.numeric(t(D2) %*% Dx)
        }
        cache$x <- x; cache$v <- list(f = f, g = g)
        cache$v
    }
    .objective <- function(x) .fg(x)$f
    .gradient <- function(x) .fg(x)$g

    ## scipy clips the start into the box before the first iteration.
    start <- pmax(as.numeric(x0), as.numeric(x_min))
    upper <- rep(if (is.finite(x_max)) as.numeric(x_max) else Inf, n)
    res <- stats::optim(par = start, fn = .objective, gr = .gradient,
                        method = "L-BFGS-B",
                        lower = rep(as.numeric(x_min), n), upper = upper,
                        control = list(maxit = as.integer(max_iterations),
                                       pgtol = as.numeric(tolerance),
                                       lmm = as.integer(lbfgs_history)))
    x <- as.numeric(res$par)
    ## scipy reports the number of L-BFGS-B iterations (result.nit); R only
    ## exposes the evaluation counts, whose second entry is the gradient (i.e.
    ## major-iteration) count.
    iterations <- as.integer(res$counts[2L])
    converged <- isTRUE(res$convergence == 0L) || iterations < max_iterations
    list(spectrum = x, iterations = iterations, converged = converged)
}

#' Wrapper around \code{\link{solve_lbfgsb}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_lbfgsb
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_lbfgsb <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, initial_spectrum = NULL,
                          max_iterations = 500L, tolerance = 1e-8,
                          regularization = 0.0, smoothness = 0.0,
                          x_min = 0.0, x_max = Inf, lbfgs_history = 10L,
                          calculate_errors = FALSE,
                          noise_level = 0.01, n_montecarlo = 100L,
                          save_result = FALSE, random_state = NULL,
                          max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_lbfgsb,
                                        max_iterations = max_iterations,
                                        tolerance = tolerance,
                                        regularization = regularization,
                                        smoothness = smoothness,
                                        x_min = x_min, x_max = x_max,
                                        lbfgs_history = lbfgs_history),
        solve_kwargs = list(),
        method_name = "L-BFGS-B",
        extra_output = list(regularization = regularization,
                            smoothness = smoothness,
                            x_min = x_min, x_max = x_max),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
