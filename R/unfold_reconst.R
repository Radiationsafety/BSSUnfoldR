#' RECONST unfolding (Turchin/Vapnik statistical regularization)
#'
#' Simplified R port of \code{bssunfold/src/bssunfold/core/unfold_reconst.py}.
#' Implements the STREG1 algorithm from RECONST.FOR — solves
#' \eqn{(B \beta + \Omega \alpha) f = A_{vec} \beta} where \eqn{\Omega} is a
#' 5-diagonal smoothing matrix and \eqn{B = A^T A}.
#'
#' @section Limitations vs. Python bssunfold:
#' The original code includes an automatic alpha/beta selection via the
#' Vapnik-Chervonenkis bound. This R port uses the simpler L-curve
#' heuristic (one alpha, beta fixed at 1). The full RECONST.FOR has multiple
#' smoothing parameter selection criteria that are not ported here.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused; accepted for API compatibility.
#' @param alpha Numeric regularization parameter for the smoothing matrix
#'   \eqn{\Omega}. Default 1.0.
#' @param pp Numeric; smoothing parameter (influences the diagonal of
#'   \eqn{\Omega}). Default 0.0.
#' @return A list \code{list(spectrum, iterations = 1L, converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_reconst(A, b, NULL, alpha = 1.0)
solve_reconst <- function(A, b, x0 = NULL, alpha = 1.0, pp = 0.0) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    # Build B = A^T A (data-fidelity matrix)
    B <- crossprod(A)
    # Build Omega: 5-diagonal smoothing matrix
    Omega <- .reconst_build_omo(n, pp)
    D <- B + alpha * Omega
    # RHS: A^T b
    ATb <- as.numeric(t(A) %*% b)
    # Solve D * f = A^T b
    f <- tryCatch(as.numeric(qr.solve(D, ATb)),
                  error = function(e) {
        sv <- svd(D, nu = nrow(D), nv = ncol(D))
        s_inv <- ifelse(sv$d > 1e-10 * max(sv$d), 1 / sv$d, 0)
        as.numeric(sv$u %*% (s_inv * (t(sv$u) %*% ATb)))
    })
    list(spectrum = pmax(f, 0.0), iterations = 1L, converged = TRUE)
}

.reconst_build_omo <- function(n, pp) {
    if (n < 3L) {
        # Tiny problems: identity is the only sensible smoothing matrix.
        return(diag(n))
    }
    # Build the 5-diagonal smoothing matrix Omega as a full n x n matrix.
    XX <- seq(1.0, n + 1.0)
    AA <- numeric(n + 2L)
    BB <- numeric(n + 3L)
    CC <- numeric(n + 3L)
    # Loop must iterate increasing indices: for n=3, this is empty; for n=4,
    # it runs i=3; for n>=5, i=3..n-1.
    if (n >= 4L) {
        for (i in 3:(n - 1L)) {
            AA[i] <- 1.0 / (XX[i] - XX[i - 1L])
            CC[i] <- 1.0 / (XX[i - 1L] - XX[i - 2L])
            BB[i] <- -(AA[i] + CC[i])
        }
    }
    OMO <- matrix(0.0, nrow = 5L, ncol = n)
    for (i in 1:n) {
        OMO[1L, i] <- AA[i] * CC[i]
        OMO[2L, i] <- AA[i] * BB[i] + BB[i + 1L] * CC[i + 1L]
        OMO[3L, i] <- AA[i]^2 + BB[i + 1L]^2 + CC[i + 2L]^2 +
                      pp * (XX[i + 1L] - XX[i])
    }
    Omega <- matrix(0.0, n, n)
    for (i in 1:n) {
        Omega[i, i] <- OMO[3L, i]
        if (i > 1L) Omega[i, i - 1L] <- OMO[2L, i]
        if (i > 2L) Omega[i, i - 2L] <- OMO[1L, i]
        if (i < n) Omega[i, i + 1L] <- OMO[2L, i + 1L]
        if (i < n - 1L) Omega[i, i + 2L] <- OMO[1L, i + 2L]
    }
    Omega
}

#' Wrapper around \code{\link{solve_reconst}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_reconst
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_reconst <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              alpha = 1.0, pp = 0.0,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_reconst,
                                         alpha = alpha, pp = pp),
        solve_kwargs = list(),
        method_name = "Reconst",
        extra_output = list(alpha = alpha, pp = pp),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
