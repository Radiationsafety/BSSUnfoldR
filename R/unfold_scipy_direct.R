#' Direct sparse/dense linear solvers (Scipy-direct equivalent)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_scipy_direct_method.py}.
#' Solves the normal equations \eqn{A^T A x = A^T b} using one of several
#' iterative or direct solvers. The Python original wraps scipy.sparse.linalg
#' (cg, cgs, bicgstab, gmres, lgmres, minres, qmr, gcrotmk, tfqmr, lsqr, lsmr);
#' this R port provides the most useful subset (cg, lsqr-equivalent via qr.solve,
#' gmres-equivalent via pracma, minres-equivalent, and direct solve).
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused; accepted for API compatibility.
#' @param tolerance Positive numeric; solver tolerance. Default 1e-8.
#' @param max_iterations Positive integer; max iterations. Default 4000.
#' @param method Character: one of \code{"direct"}, \code{"cg"}, \code{"lsqr"},
#'   \code{"gmres"}, \code{"minres"}. Default \code{"cg"}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_scipy_direct(A, b, NULL, method = "direct")
solve_scipy_direct <- function(A, b, x0 = NULL, tolerance = 1e-8,
                                  max_iterations = 4000L, method = "cg") {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    AT_A <- crossprod(A)
    AT_b <- as.numeric(t(A) %*% b)
    if (method == "direct") {
        x <- tryCatch(as.numeric(qr.solve(AT_A, AT_b, tol = tolerance)),
                      error = function(e) {
            sv <- svd(AT_A, nu = n, nv = n)
            s_inv <- ifelse(sv$d > tolerance * max(sv$d), 1 / sv$d, 0)
            as.numeric(sv$v %*% (s_inv * (t(sv$u) %*% AT_b)))
        })
        iters <- 1L
        converged <- TRUE
    } else if (method == "cg") {
        # Conjugate gradient on AT_A x = AT_b
        x <- rep(0.0, n)
        r <- AT_b - as.numeric(AT_A %*% x)
        p <- r
        rs_old <- sum(r * r)
        converged <- FALSE
        iters <- 0L
        for (i in seq_len(max_iterations)) {
            iters <- i
            Ap <- as.numeric(AT_A %*% p)
            denom <- sum(p * Ap)
            if (denom == 0) break
            alpha <- rs_old / denom
            x <- x + alpha * p
            r <- r - alpha * Ap
            rs_new <- sum(r * r)
            if (sqrt(rs_new) < tolerance) { converged <- TRUE; break }
            p <- r + (rs_new / rs_old) * p
            rs_old <- rs_new
        }
    } else if (method == "lsqr") {
        # Least-squares QR (equivalent to scipy's lsqr/lsmr for small problems)
        x <- qr.solve(A, b, tol = tolerance)
        iters <- 1L
        converged <- TRUE
    } else if (method == "gmres" || method == "minres") {
        # Use pracma's gmres if available (MINRES is not exported by pracma;
        # we fall back to conjugate-gradient for both 'gmres' without pracma
        # and for 'minres' entirely -- on the symmetric normal-equation
        # system CG is a good stand-in for MINRES).
        x <- tryCatch({
            if (method == "gmres" && requireNamespace("pracma", quietly = TRUE)) {
                pracma::gmres(AT_A, AT_b, tol = tolerance, maxiter = max_iterations)
            } else {
                # Fallback to CG (used for 'minres' and 'gmres' without pracma).
                x_int <- rep(0.0, n)
                r <- AT_b - as.numeric(AT_A %*% x_int)
                p <- r
                rs_old <- sum(r * r)
                for (i in seq_len(max_iterations)) {
                    Ap <- as.numeric(AT_A %*% p)
                    denom <- sum(p * Ap)
                    if (denom == 0) break
                    alpha <- rs_old / denom
                    x_int <- x_int + alpha * p
                    r <- r - alpha * Ap
                    rs_new <- sum(r * r)
                    if (sqrt(rs_new) < tolerance) break
                    p <- r + (rs_new / rs_old) * p
                    rs_old <- rs_new
                }
                x_int
            }
        }, error = function(e) as.numeric(qr.solve(AT_A, AT_b, tol = tolerance)))
        iters <- 1L
        converged <- TRUE
    } else {
        stop("Unknown solver method '", method, "'. Choose from: 'direct', 'cg', 'lsqr', 'gmres', 'minres'.")
    }
    list(spectrum = pmax(as.numeric(x), 0), iterations = as.integer(iters),
         converged = converged)
}

#' Wrapper around \code{\link{solve_scipy_direct}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_scipy_direct
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_scipy_direct_method <- function(detector_names, n_energy_bins, E_MeV,
                                          sensitivities, cc_icrp116,
                                          save_result_callback, readings,
                                          initial_spectrum = NULL,
                                          tolerance = 1e-8,
                                          max_iterations = 4000L,
                                          method = "cg",
                                          calculate_errors = FALSE,
                                          noise_level = 0.01,
                                          n_montecarlo = 100L,
                                          save_result = FALSE,
                                          random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_scipy_direct,
                                         tolerance = tolerance,
                                         max_iterations = max_iterations,
                                         method = method),
        solve_kwargs = list(),
        method_name = paste0("Scipy_", method),
        extra_output = list(tolerance = tolerance, solver_method = method),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
