#' Express piecewise-exponential unfolding
#'
#' R port of \code{bssunfold/core/unfold_express.py}. Fits a
#' piecewise-exponential spectrum (linear in log-space between boundaries)
#' directly to sphere readings via nonlinear least-squares.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid (length n).
#' @param x0 Optional initial spectrum guess. Default \code{NULL} = flat.
#' @param n_groups Integer; number of piecewise-exponential groups. Default 6.
#' @param interval_boundaries Optional numeric vector of group boundary
#'   energies. Default \code{NULL} = evenly spaced.
#' @param max_iterations Integer; maximum number of function evaluations
#'   (will be multiplied by 100 internally). Default 3.
#' @param tol_iteration Numeric; convergence tolerance on relative residual.
#'   Default 0.05.
#' @param relative_uncertainty Numeric; relative measurement uncertainty.
#'   Default 0.05.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 10)
#' A <- matrix(runif(3 * 10), nrow = 3)
#' b <- as.numeric(A %*% rep(0.5, 10))
#' r <- solve_express(A, b, E, n_groups = 4)
solve_express <- function(A, b, E, x0 = NULL, n_groups = 6L,
                            interval_boundaries = NULL,
                            max_iterations = 3L, tol_iteration = 0.05,
                            relative_uncertainty = 0.05) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E)
    if (any(b < 0)) stop("Express requires non-negative readings")
    if (any(diff(E) <= 0)) stop("Express requires strictly increasing E")
    if (is.null(interval_boundaries)) {
        if (n_groups < 2L) stop("n_groups must be at least 2")
        boundaries <- seq(E[1L], E[length(E)], length.out = n_groups + 1L)
    } else {
        boundaries <- as.numeric(interval_boundaries)
        if (length(boundaries) < 2L || any(diff(boundaries) <= 0)) {
            stop("interval_boundaries must be strictly increasing")
        }
    }
    centers <- boundaries
    initial <- if (is.null(x0)) rep(1.0, length(E)) else pmax(as.numeric(x0), 1e-30)
    guess <- approx(x = E, y = log(initial), xout = centers, rule = 2)$y
    sigma <- pmax(relative_uncertainty * pmax(b, 1e-30), 1e-30)

    .express_model <- function(log_values) {
        spectrum <- exp(approx(x = boundaries, y = log_values,
                                xout = E, rule = 2)$y)
        as.numeric(A %*% spectrum)
    }
    .express_resid <- function(log_values) {
        (.express_model(log_values) - b) / sigma
    }
    # Levenberg-Marquardt via stats::nls.lm would need minpack.lm; use a
    # simple Levenberg-Marquardt loop with numerical Jacobian.
    p <- guess
    lambda <- 1e-3
    nfev <- 0L
    converged <- FALSE
    for (iter in seq_len(max(1L, max_iterations) * 100L)) {
        nfev <- iter
        r0 <- .express_resid(p)
        # Numerical Jacobian: finite differences
        J <- matrix(0.0, nrow = length(r0), ncol = length(p))
        eps <- 1e-8
        for (j in seq_along(p)) {
            dp <- rep(0.0, length(p)); dp[j] <- eps * max(1.0, abs(p[j]))
            r1 <- .express_resid(p + dp)
            J[, j] <- (r1 - r0) / dp[j]
        }
        # Levenberg-Marquardt update: (J'J + lambda*diag(J'J)) dp = -J'r0
        JtJ <- crossprod(J)
        Jtr <- as.numeric(t(J) %*% r0)
        diag_JtJ <- diag(JtJ)
        diag_JtJ[diag_JtJ == 0] <- 1e-12
        A_lm <- JtJ + lambda * diag(diag_JtJ)
        dp <- tryCatch(as.numeric(qr.solve(A_lm, -Jtr)),
                       error = function(e) rep(0.0, length(p)))
        p_new <- p + dp
        r_new <- .express_resid(p_new)
        if (sum(r_new^2) < sum(r0^2)) {
            p <- p_new
            lambda <- max(lambda * 0.5, 1e-12)
        } else {
            lambda <- min(lambda * 2.0, 1e8)
        }
        if (max(abs(dp)) < 1e-8) { converged <- TRUE; break }
    }
    spectrum <- exp(approx(x = boundaries, y = p, xout = E, rule = 2)$y)
    rel_change <- sqrt(sum((as.numeric(A %*% spectrum) - b)^2)) /
                  (sqrt(sum(b^2)) + 1e-30)
    converged <- converged || rel_change <= tol_iteration
    list(spectrum = pmax(spectrum, 0.0), iterations = nfev,
         converged = converged)
}

#' Wrapper around \code{\link{solve_express}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_express
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_express <- function(detector_names, n_energy_bins, E_MeV,
                             sensitivities, cc_icrp116, save_result_callback,
                             readings, initial_spectrum = NULL,
                             n_groups = 6L, interval_boundaries = NULL,
                             max_iterations = 3L, tol_iteration = 0.05,
                             relative_uncertainty = 0.05,
                             calculate_errors = FALSE,
                             noise_level = 0.01, n_montecarlo = 100L,
                             save_result = FALSE, random_state = NULL) {
    boundaries <- interval_boundaries
    if (is.null(boundaries)) {
        boundaries <- seq(E_MeV[1L], E_MeV[length(E_MeV)],
                           length.out = n_groups + 1L)
    }
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1.0, n_energy_bins),
        solve_func = make_solve_wrapper(
            function(A, b, x0 = NULL) solve_express(
                A, b, E_MeV, x0 = x0, n_groups = n_groups,
                interval_boundaries = boundaries,
                max_iterations = max_iterations,
                tol_iteration = tol_iteration,
                relative_uncertainty = relative_uncertainty
            )
        ),
        solve_kwargs = list(),
        method_name = "Express",
        extra_output = list(
            n_groups = as.integer(length(boundaries) - 1L),
            interval_boundaries = as.numeric(boundaries),
            tol_iteration = tol_iteration,
            log_steps = compute_log_steps(E_MeV)
        ),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
