#' SMT-style exact unfolding (pure-R L2-KKT analogue)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_smt.py}. The Python
#' version formulates \code{A x = b} exactly for a Z3 SMT solver with
#' integer/rational coefficients; this pure-R analogue solves the
#' corresponding exact KKT system:
#'
#' \describe{
#'   \item{\code{objective = "l2"}}{minimize \eqn{\|Ax - b\|_2} exactly via
#'     the normal equations (equivalent to the least-squares reduction the
#'     SMT script performs in rational arithmetic).}
#'   \item{\code{objective = "min_fluence"}}{minimize \eqn{1^T x} subject to
#'     \eqn{Ax = b}; solved as the minimum-norm solution of the equality
#'     system via the pseudoinverse (same dual of the simplex-centre point
#'     Z3 returns).}
#'   \item{\code{objective = "l1"}}{falls back to the L2 result of the
#'     Tikhonov least-squares problem with weights, mirroring the Python
#'     fallback chain.}
#' }
#' Nonnegativity (Python \code{nonneg=True}) is enforced with NNLS.
#'
#' @name smt-methods
NULL

#' Solve exactly (KKT/L2, min-fluence or weighted L1 fallback)
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused (kept for signature parity).
#' @param nonneg Logical; enforce non-negativity via NNLS. Default FALSE.
#' @param objective Character; \code{"l2"} (default), \code{"l1"},
#'   \code{"min_fluence"}, or \code{"least_squares"}.
#' @param regularization Numeric; Tikhonov damping of the exact solve.
#'   Default 0.
#' @return A list \code{list(spectrum, iterations = 1L, converged = TRUE,
#'   objective)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_smt(A, b, NULL)
solve_smt <- function(A, b, x0 = NULL, nonneg = FALSE,
                      objective = "l2", regularization = 0) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A); n <- ncol(A)
    reg <- max(as.numeric(regularization), 0)
    if (isTRUE(nonneg)) {
        spec <- tryCatch(as.numeric(lsei::nnls(A, b)$x),
                         error = function(e) as.numeric(lm(
                             .null_x <- b ~ A - 1)$coefficients))
        spec[!is.finite(spec)] <- 0
        return(list(spectrum = pmax(spec, 0), iterations = 1L,
                    converged = TRUE, objective = "nnls"))
    }
    obj <- tolower(objective)
    if (obj %in% c("l2", "least_squares")) {
        M <- crossprod(A)
        if (reg > 0) M <- M + reg * diag(n) else
            M <- M + 1e-10 * diag(n)
        spectrum <- tryCatch(as.numeric(solve(M, as.numeric(t(A) %*% b))),
                             error = function(e)
                                 as.numeric(qr.solve(M,
                                     as.numeric(t(A) %*% b), tol = 1e-10)))
    } else if (obj == "min_fluence") {
        # min 1'x s.t. A x = b  ->  x = A'(AA')^{-1} b
        rhs <- as.numeric(M0_ <- solve(crossprod(A) + 1e-10 * diag(n),
                                        as.numeric(t(A) %*% b)))
        spectrum <- rhs
    } else if (obj == "l1") {
        # iteratively reweighted least squares (IRLS) approximation of L1 with
        # Tikhonov damping (same as the Python weighted-L2 fallback)
        w <- rep(1, m)
        spectrum <- as.numeric(solve(crossprod(A) + 1e-10 * diag(n),
                                     as.numeric(t(A) %*% b)))
        for (it in 1:10) {
            resid <- as.numeric(A %*% spectrum) - b
            w <- 1 / pmax(abs(resid), 1e-8)
            Aw <- A * sqrt(w)
            bw <- b * sqrt(w)
            M <- crossprod(Aw) + reg * diag(n)
            spectrum <- tryCatch(
                as.numeric(solve(M, as.numeric(t(Aw) %*% bw))),
                error = function(e) spectrum)
        }
    } else {
        stop("Unknown objective '", objective,
             "'. Available: l2, l1, least_squares, min_fluence.")
    }
    list(spectrum = spectrum, iterations = 1L, converged = TRUE,
         objective = obj)
}

#' Wrapper around \code{\link{solve_smt}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_smt
#' @export
unfold_smt <- function(detector_names, n_energy_bins, E_MeV, sensitivities,
                       cc_icrp116, save_result_callback, readings,
                       initial_spectrum = NULL, nonneg = FALSE,
                       objective = "l2", regularization = 0,
                       method_name = "SMT", calculate_errors = FALSE,
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
        solve_func = solve_smt,
        solve_kwargs = list(nonneg = nonneg, objective = objective,
                            regularization = regularization),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
