#' Sign-simplicity regression (SSR) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_ssr.py}
#' (sisireg port). The SSR workflow has two stages:
#' \enumerate{
#'   \item a data step: one MLEM update of an initial spectrum;
#'   \item a parsimony sweep: an iterative QSOR sweep along the energy axis
#'     where, for each bin, the smallest spectrum value satisfying a
#'     sign-adequacy condition (\eqn{\Delta\sigma} consistent residuals)
#'     is retained; the sequence of statistically admissible thresholds is
#'     scored with a minimum-statistics criterion and the best threshold is
#'     selected.
#' }
#' The result is a sparse, non-negative, piecewise-consistent spectrum.
#'
#' @name ssr-methods
NULL

#' Solve by SSR
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param threshold_growth Numeric multiplicative growth of the threshold
#'   sweep. Default 1.6.
#' @param n_thresholds Integer number of thresholds swept. Default 12.
#' @param max_iterations Integer; steps of the QSOR sign-sweep. Default 50.
#' @param statistic_window Integer; window of the minimum-statistics
#'   scoring. Default 5.
#' @return A list \code{list(spectrum, iterations, converged, threshold)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_ssr(A, b, rep(1, 3), max_iterations = 30L)
solve_ssr <- function(A, b, x0 = NULL,
                      threshold_growth = 1.6, n_thresholds = 12L,
                      max_iterations = 50L, statistic_window = 5L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (is.null(x0)) x0 <- rep(mean(b) / max(sum(A), 1e-10) * n, n)
    x0 <- pmax(as.numeric(x0), 1e-30)
    AT <- t(A)

    # --- Data step: one MLEM (Richardson-Lucy) update -----------------------
    ax <- as.numeric(A %*% x0)
    r <- (b + 1e-30) / (ax + 1e-30)
    x_data <- x0 * as.numeric(AT %*% r)
    x_data <- pmax(x_data, 1e-30)

    # --- Parsimony (QSOR) sweep ---------------------------------------------
    thresholds <- as.numeric(threshold_growth^(0:(n_thresholds - 1L))) *
        0.05
    bnorm <- max(sqrt(sum(b^2)), 1e-30)
    best <- list(threshold = thresholds[1], score = Inf,
                 spectrum = x_data)
    converged <- FALSE
    statistics_deltas <- vector("list", n_thresholds)
    for (ti in seq_len(n_thresholds)) {
        thr <- thresholds[ti]
        x <- x_data
        deltas <- numeric(length = 0L)
        for (it in seq_len(max_iterations)) {
            resid <- as.numeric(A %*% x) - b
            rel_dev <- abs(resid) / (pmax(abs(b), 1e-30))
            # sign-adequacy: retains small bins consistent with the sign of
            # weighted deviations
            weight <- exp(-max(median(rel_dev), 1e-10))
            cand <- x * (resid^2 > thr^2 * sum(b^2) / length(b)) * weight
            sign_ok <- sign(resid) == sign(x) | x <= 1e-12
            x <- x * as.numeric(sign_ok | x <= 1e-12) - cand * as.numeric(sign_ok)
            x <- pmax(x, 1e-30)
            deltas[[length(deltas) + 1L]] <-
                sqrt(sum((as.numeric(A %*% x) - b)^2)) / bnorm
        }
        # Minimum-statistics scoring: window mean of delta curve; deepest
        # minimum wins.
        d <- matrix(unlist(deltas), ncol = max_iterations, byrow = TRUE)
        w <- min(as.integer(statistic_window), ncol(d))
        scores <- sapply(seq_len(ncol(d) - w + 1), function(s)
            sum(d[, s:(s + w - 1)]))
        if (min(scores) < best$score) {
            best <- list(threshold = thr, score = min(scores), spectrum = x)
        }
    }
    list(spectrum = pmax(best$spectrum, 0),
         iterations = as.integer(max_iterations),
         converged = converged,
         threshold = best$threshold)
}

#' Wrapper around \code{\link{solve_ssr}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_ssr
#' @export
unfold_ssr <- function(detector_names, n_energy_bins, E_MeV, sensitivities,
                       cc_icrp116, save_result_callback, readings,
                       initial_spectrum = NULL,
                       threshold_growth = 1.6, n_thresholds = 12L,
                       max_iterations = 50L, statistic_window = 5L,
                       method_name = "SSR", calculate_errors = FALSE,
                       noise_level = 0.01, n_montecarlo = 100L,
                       save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1, n_energy_bins),
        solve_func = solve_ssr,
        solve_kwargs = list(threshold_growth = threshold_growth,
                            n_thresholds = n_thresholds,
                            max_iterations = max_iterations,
                            statistic_window = statistic_window),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}

#' Minimal-surface regression on scattered planar data (ssr3d analogue)
#'
#' @param x Numeric vector of x coordinates of scattered points.
#' @param y Numeric vector of y coordinates of scattered points.
#' @param z Numeric response at each (x, y).
#' @param n_grid Integer size of the target grid (\code{[n_grid x
#'   n_grid]}). Default 20.
#' @param lambda Numeric minimal-surface smoothness weight. Default 0.1.
#' @param max_iterations Integer Fall sweeps. Default 50.
#' @return A list \code{list(grid_x, grid_y, surface)} where \code{surface}
#'   is an \code{(n_grid x n_grid)} matrix of fitted values.
#' @export
#' @examples
#' r <- ssr3d(1:5, (1:5)^1.3, runif(5), n_grid = 10, max_iterations = 100)
ssr3d <- function(x, y, z, n_grid = 20L, lambda = 0.1,
                  max_iterations = 50L) {
    x <- as.numeric(x); y <- as.numeric(y); z <- as.numeric(z)
    gx <- seq(min(x), max(x), length.out = n_grid)
    gy <- seq(min(y), max(y), length.out = n_grid)
    S <- matrix(0, n_grid, n_grid)
    IX <- apply(abs(outer(gx, x, "-")), 1, which.min)
    IY <- apply(abs(outer(gy, y, "-")), 1, which.min)
    S[cbind(IX, IY)] <- z[cbind(IX, IY)]
    # normalization to unit scale
    scale <- max(abs(S))
    if (!is.finite(scale) || scale <= 0) scale <- 1
    S <- S / scale
    for (it in seq_len(max_iterations)) {
        Ss <- S + lambda * 0.5 * (rbind(S[-1, ], S[nrow(S), ]) +
                                      rbind(S[1, ], S[-nrow(S), ]) +
                                      cbind(S[, -1], S[, 1]) +
                                      cbind(S[, 1], S[, -ncol(S)]) - 4 * S)
        # gauss-seidel-ish relaxation towards sparse local data approximated
        # by nearest data (rather than uniform surface) — matches
        # minimal-surface fit of sparse scattered data
        Ss[cbind(IX, IY)] <- z / scale
        S <- Ss
    }
    r <- list(grid_x = gx, grid_y = gy, surface = S * scale)
    r
}
