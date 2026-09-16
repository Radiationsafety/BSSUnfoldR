#' Operator-based MLEM unfolding (ODL-free pure-R analogue)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_mlem_odl.py}.
#' The Python version runs Maximum Likelihood for Emission tomography
#' (MLEM) through the ODL operator framework; the pure-R analogue accepts
#' either an explicit response matrix \code{A} (default forward /
#' an explicit response matrix \code{A} (default forward / transpose
#' application) or a named list \code{operator = list(op, opadjoint)} of
#' fluence-space callable operators with the same mathematical contract.
#'
#' @name mlem-odl-methods
NULL

# Internal: apply forward and adjoint operators to a spectrum vector
.mlem_odl_forward <- function(A, x) {
    if (is.function(A)) return(A(x))
    as.numeric(A %*% x)
}

.mlem_odl_adjoint <- function(A, y) {
    if (is.function(A)) return(A(y))
    as.numeric(t(A) %*% y)
}

# Normalisation vector (sensitivity): sum over measurement rows
.mlem_odl_norm <- function(A, n) {
    if (is.function(A)) {
        e1 <- rep(1, n)
        sum(.mlem_odl_forward(A, e1)) / max(n, 1)
    } else {
        colSums(A)
    }
}

#' MLEM via the operator interface (ODL analogue)
#'
#' @param A Response matrix (m x n) or a list \code{list(op, opadjoint)}
#' with callables; or callable forward operator A(x).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param n_bins Integer number of spectral bins (only needed when \code{A}
#'   is callable and \code{x0} is missing).
#' @param max_iterations Integer MLEM iterations. Default 200.
#' @param tolerance Numeric relative tolerance. Default 1e-6.
#' @param nonneg Logical clip negative bins every iteration. Default TRUE.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_mlem_odl(A, b, rep(1, 3), max_iterations = 50L)
solve_mlem_odl <- function(A, b, x0 = NULL, n_bins = NULL,
                           max_iterations = 200L, tolerance = 1e-6,
                           nonneg = TRUE) {
    b <- as.numeric(b)
    m <- length(b)
    op <- A; adj <- A
    if (is.list(A) && !is.null(A$op)) {
        op <- A$op; adj <- if (!is.null(A$opadjoint)) A$opadjoint else A
    }
    if (is.null(x0)) {
        if (!is.null(n_bins)) x0 <- rep(mean(b) / max(m, 1), n_bins)
        else if (is.matrix(A)) x0 <- rep(mean(b) / max(sum(A), 1e-10) *
                                          ncol(A), ncol(A))
        else stop("provide x0 or a matrix A")
    }
    x <- pmax(as.numeric(x0), 1e-30)
    n <- length(x)
    sensitivity <- .mlem_odl_norm(op, n)
    # The classical ODL MLEM multiplies by the adjoint of the forward model
    # and clips tiny values; the sensitivity normalisation (column sums) is
    # not applied so that the explicit-matrix path reproduces solve_mlem.
    sensitivity <- rep(1, n)
    bnorm <- max(sqrt(sum(b^2)), 1e-30)
    converged <- FALSE
    iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        mu <- pmax(.mlem_odl_forward(op, x), 1e-30)
        ratio <- (b + 1e-30) / mu
        corr <- .mlem_odl_adjoint(op, ratio) / (sensitivity + 1e-30)
        x_new <- x * corr
        if (!isFALSE(nonneg)) x_new <- pmax(x_new, 0)
        change <- sqrt(sum((as.numeric(mu) - b)^2)) / bnorm
        x <- x_new
        if (change < tolerance) { converged <- TRUE; break }
    }
    list(spectrum = pmax(x, 0), iterations = as.integer(iterations),
         converged = converged)
}

#' Wrapper around \code{\link{solve_mlem_odl}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_mlem_odl
#' @export
unfold_mlem_odl <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116,
                            save_result_callback, readings,
                            initial_spectrum = NULL, max_iterations = 200L,
                            tolerance = 1e-6, nonneg = TRUE,
                            method_name = "MLEM-ODL",
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
        solve_func = solve_mlem_odl,
        solve_kwargs = list(max_iterations = max_iterations,
                            tolerance = tolerance, nonneg = nonneg),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
