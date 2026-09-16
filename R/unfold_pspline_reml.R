#' P-spline REML unfolding (mixed-model smoothing)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_pspline_reml.py}.
#' Represents the unfolded spectrum in a B-spline basis and estimates the
#' smoothing parameter by restricted maximum likelihood (REML), as in
#' mixed-model smoothing (Eilers & Marx 1996; Currie et al. 2004).
#' Iterates between (a) fitting the spline coefficients for fixed
#' \code{lambda}, and (b) re-estimating \code{lambda} on a REML grid around
#' the previous value.
#'
#' @name pspline-reml-methods
NULL

# ---- internal helpers ------------------------------------------------------

# Compact-cardinal B-spline basis matrix: n rows (bin positions on [0,1])
# x n_basis columns; positive Epanechnikov-compact kernels on interior knots.
.pspline_basis <- function(n, n_basis) {
    t <- seq(0, 1, length.out = n)
    u <- seq(0, 1, length.out = n_basis)
    h <- max(1 / max(n_basis - 1, 1), 1e-6)
    B <- outer(t, u, function(ti, uj) {
        x <- (ti - uj) / h
        (1 - x^2)^2 * (abs(x) < 1)
    })
    B + matrix(1e-6, nrow = nrow(B), ncol = ncol(B))
}

# Penalized coefficient fit: minimises ||RB c - b||^2 / sig2 + lambda ||D2 c||^2
.pspline_fit <- function(A, b, n_basis, lambda, DTD, sig2) {
    B <- .pspline_basis(ncol(A), n_basis)
    RB <- A %*% B
    lhs <- crossprod(RB) / sig2 + (lambda / sig2) * DTD
    rhs <- crossprod(RB, b) / sig2
    tryCatch(as.numeric(solve(lhs + 1e-10 * diag(ncol(lhs)), rhs)),
             error = function(e) rep(0, n_basis))
}

# Simplified REML criterion: residual deviance + log|H|
.pspline_reml_criterion <- function(A, b, n_basis, lambda, DTD, sig2) {
    B <- .pspline_basis(ncol(A), n_basis)
    RB <- A %*% B
    lhs <- crossprod(RB) / sig2 + (lambda / sig2) * DTD
    rhs <- crossprod(RB, b) / sig2
    pd <- tryCatch(solve(lhs + 1e-10 * diag(ncol(lhs)), rhs),
                   error = function(e) NULL)
    if (is.null(pd)) return(Inf)
    coef <- as.numeric(pd)
    resid <- as.numeric(RB %*% coef) - b
    dev <- sum(resid^2) / sig2
    logdet <- determinant(lhs, logarithm = TRUE)$modulus[[1]]
    dev + logdet
}

#' Solve by P-spline REML
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n); reserved for API
#'   compatibility (the fit is linear in the basis).
#' @param n_basis Integer; number of B-spline basis functions. Default 15.
#' @param lambda_init Numeric; starting smoothness parameter. Default 1.0.
#' @param max_iterations Integer; outer REML iterations. Default 50.
#' @param tolerance Numeric; convergence on \code{log10(lambda)}. Default
#'   1e-3.
#' @param n_grid Integer; number of \code{lambda} values tried per iteration.
#'   Default 7.
#' @param nonneg Logical; clip final spectrum to be non-negative. Default
#'   TRUE.
#' @return A list \code{list(spectrum, iterations, converged, lambda,
#'   coefficients)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_pspline_reml(A, b, NULL, n_basis = 5, max_iterations = 5L)
solve_pspline_reml <- function(A, b, x0 = NULL, n_basis = 15L,
                               lambda_init = 1.0, max_iterations = 50L,
                               tolerance = 1e-3, n_grid = 7L,
                               nonneg = TRUE) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    n_basis <- max(as.integer(n_basis), 4L)
    sigma2 <- max(mean(b^2), 1e-30)
    D2 <- second_difference_matrix(n_basis)
    DTD <- crossprod(D2)
    lam <- log10(max(lambda_init, 1e-12))
    converged <- FALSE
    iterations <- 0L
    coef <- rep(0, n_basis)
    for (it in seq_len(max_iterations)) {
        iterations <- it
        lam_old <- lam
        lams <- lam + seq(-1.0, 1.0,
                          length.out = max(as.integer(n_grid), 3L))
        reml <- vapply(lams, function(l)
            .pspline_reml_criterion(A, b, n_basis, 10.0^l, DTD, sigma2),
            numeric(1))
        lam <- lams[which.min(reml)]
        if (abs(lam - lam_old) < tolerance) {
            converged <- TRUE
            break
        }
    }
    coef <- .pspline_fit(A, b, n_basis, 10.0^lam, DTD, sigma2)
    B <- .pspline_basis(n, n_basis)
    spectrum <- as.numeric(B %*% coef)
    if (isTRUE(nonneg)) spectrum <- pmax(spectrum, 0)
    list(spectrum = spectrum,
         iterations = as.integer(iterations),
         converged = converged,
         lambda = 10.0^lam,
         coefficients = coef)
}

#' Wrapper around \code{\link{solve_pspline_reml}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_pspline_reml
#' @export
unfold_pspline_reml <- function(detector_names, n_energy_bins, E_MeV,
                                sensitivities, cc_icrp116,
                                save_result_callback, readings,
                                initial_spectrum = NULL, n_basis = 15L,
                                lambda_init = 1.0, max_iterations = 50L,
                                tolerance = 1e-3, n_grid = 7L, nonneg = TRUE,
                                method_name = "PSplineREML",
                                calculate_errors = FALSE,
                                noise_level = 0.01, n_montecarlo = 100L,
                                save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    n_basis <- max(as.integer(n_basis), 4L)
    if (n_basis > n_energy_bins + 2L) n_basis <- n_energy_bins + 2L
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1, n_energy_bins),
        solve_func = solve_pspline_reml,
        solve_kwargs = list(n_basis = n_basis, lambda_init = lambda_init,
                            max_iterations = max_iterations,
                            tolerance = tolerance, n_grid = n_grid,
                            nonneg = nonneg),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
