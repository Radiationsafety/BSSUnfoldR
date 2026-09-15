#' Tikhonov regularization with Legendre polynomial basis unfolding
#'
#' R port of \code{bssunfold/core/unfold_tikhonov_legendre.py}. Projects the
#' response matrix onto a Legendre polynomial basis, applies second-derivative
#' regularization on the polynomial coefficients, and solves the augmented
#' least-squares system.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Ignored; accepted for API compatibility.
#' @param delta Numeric regularization parameter. Default 0.05.
#' @param n_polynomials Integer; number of Legendre polynomials. Default 15.
#' @return A list \code{list(spectrum, iterations = 0L, converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_tikhonov_legendre(A, b, NULL, n_polynomials = 3)
solve_tikhonov_legendre <- function(A, b, x0 = NULL, delta = 0.05,
                                     n_polynomials = 15L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n_energy <- ncol(A)
    n_polynomials <- as.integer(n_polynomials)
    if (n_polynomials < 3L) {
        stop("n_polynomials must be >= 3 for second-derivative regularization")
    }
    basis <- .build_legendre_basis(n_energy, n_polynomials)
    A_proj <- A %*% basis
    L <- as.matrix(create_derivative_matrix(n_polynomials, 2L))
    M <- rbind(A_proj, delta * L)
    rhs <- c(b, rep(0, nrow(L)))
    c_coefs <- tryCatch(qr.solve(M, rhs),
                        error = function(e) qr.solve(M, rhs, tol = 1e-8))
    spectrum <- as.numeric(basis %*% c_coefs[seq_len(n_polynomials)])
    list(spectrum = pmax(spectrum, 0.0), iterations = 0L, converged = TRUE)
}

.build_legendre_basis <- function(n_energy, n_polynomials) {
    x <- seq(-1, 1, length.out = n_energy)
    basis <- matrix(0.0, nrow = n_energy, ncol = n_polynomials)
    for (i in seq_len(n_polynomials)) {
        # Legendre polynomial of degree (i-1) using pracma::legendre
        # or via recursion to avoid hard dependency
        # P_0 = 1, P_1 = x, P_{n+1} = ((2n+1) x P_n - n P_{n-1}) / (n+1)
        deg <- i - 1L
        if (deg == 0L) {
            P <- rep(1, n_energy)
        } else if (deg == 1L) {
            P <- x
        } else {
            P_prev <- rep(1, n_energy)
            P_curr <- x
            for (n in seq_len(deg - 1L)) {
                P_next <- ((2 * n + 1) * x * P_curr - n * P_prev) / (n + 1)
                P_prev <- P_curr
                P_curr <- P_next
            }
            P <- P_curr
        }
        basis[, i] <- P
    }
    basis
}

#' Wrapper around \code{\link{solve_tikhonov_legendre}} for the unified
#' workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_tikhonov_legendre
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_tikhonov_legendre <- function(detector_names, n_energy_bins, E_MeV,
                                        sensitivities, cc_icrp116,
                                        save_result_callback, readings,
                                        initial_spectrum = NULL,
                                        delta = 0.05,
                                        n_polynomials = 15L,
                                        calculate_errors = FALSE,
                                        noise_level = 0.01,
                                        n_montecarlo = 100L,
                                        save_result = FALSE,
                                        random_state = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_tikhonov_legendre,
                                         delta = delta,
                                         n_polynomials = n_polynomials),
        solve_kwargs = list(),
        method_name = "Tikhonov_Legendre",
        extra_output = list(delta = delta, n_polynomials = n_polynomials),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
