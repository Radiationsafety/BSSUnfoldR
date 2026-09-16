#' Crystal Ball direct unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_crystal_ball.py}.
#' CRYSTAL BALL unfolds the spectrum directly (without iteration) by
#' representing the unknown spectrum \eqn{\phi} as a linear combination of
#' the detector response functions (rows of \code{A}):
#' \deqn{\phi_j = \sum_i \alpha_i A_{ij}}
#' Substituting into the measurement equation \eqn{b_i = \sum_j A_{ij} \phi_j}
#' gives \eqn{b = (A A^T) \alpha \Rightarrow \alpha = (A A^T + \lambda I)^{-1} b}
#' and the recovered spectrum is \eqn{\phi = A^T \alpha}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused; accepted for API compatibility.
#' @param regularization Numeric Tikhonov regularization strength added to
#'   the diagonal of \eqn{A A^T}. Default 0.0.
#' @return A list \code{list(spectrum, iterations = 1L, converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_crystal_ball(A, b, NULL, regularization = 0.01)
solve_crystal_ball <- function(A, b, x0 = NULL, regularization = 0.0) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    if (nrow(A) == 0L || length(b) == 0L) {
        stop("Response matrix and measurements must be non-empty")
    }
    if (all(b <= 0)) stop("All measurements are zero or negative")
    m <- nrow(A)
    G <- tcrossprod(A)  # A A^T
    if (regularization > 0) {
        G <- G + regularization * diag(m)
    }
    alpha <- tryCatch(as.numeric(qr.solve(G, b)),
                      error = function(e) {
                          # Last-resort: pseudo-inverse via SVD.
                          sv <- svd(G, nu = nrow(G), nv = 0)
                          s_inv <- ifelse(sv$d > 1e-10 * max(sv$d), 1 / sv$d, 0)
                          as.numeric(sv$u %*% (s_inv * (t(sv$u) %*% b)))
                      })
    spectrum <- as.numeric(t(A) %*% alpha)
    list(spectrum = pmax(spectrum, 0.0), iterations = 1L, converged = TRUE)
}

#' Wrapper around \code{\link{solve_crystal_ball}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_crystal_ball
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_crystal_ball <- function(detector_names, n_energy_bins, E_MeV,
                                  sensitivities, cc_icrp116,
                                  save_result_callback, readings,
                                  initial_spectrum = NULL,
                                  regularization = 0.0,
                                  calculate_errors = FALSE,
                                  noise_level = 0.01, n_montecarlo = 100L,
                                  save_result = FALSE, random_state = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_crystal_ball,
                                         regularization = regularization),
        solve_kwargs = list(),
        method_name = "CRYSTAL_BALL",
        extra_output = list(regularization = regularization),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
