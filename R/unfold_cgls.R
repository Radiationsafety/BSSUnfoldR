#' CGLS (Conjugate Gradient Least Squares) unfolding
#'
#' Applies the conjugate-gradient algorithm implicitly to the normal equations
#' \eqn{A^T A x = A^T b}{A'A x = A'b}. A regularized solution is obtained by
#' early termination of the iterations (semi-convergence). Optionally a
#' Tikhonov term \eqn{lambda^2 ||L x||^2} may be added.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial guess (length n). Default \code{NULL} = zero vector.
#' @param max_iterations Positive integer; default 100.
#' @param tolerance Positive numeric; default 1e-12.
#' @param noise_level Optional relative noise level used for
#'   discrepancy-principle stopping. Default \code{NULL}.
#' @param regularization Numeric Tikhonov parameter. Default 0.0.
#' @param smoothness_order Integer; 0 (identity), 1 or 2. Default 0.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_cgls(A, b, rep(0, 3), max_iterations = 50)
solve_cgls <- function(A, b, x0 = NULL, max_iterations = 100L,
                       tolerance = 1e-12, noise_level = NULL,
                       regularization = 0.0, smoothness_order = 0L) {
    v <- validate_system(A, b, x0 = x0, max_iterations = max_iterations,
                        tolerance = tolerance)
    A <- v$A; b <- v$b; x0 <- v$x0
    m <- nrow(A); n <- ncol(A)
    if (is.null(x0)) x <- rep(0.0, n) else x <- as.numeric(x0)
    nrmb <- sqrt(sum(b^2))
    if (nrmb == 0.0) {
        return(list(spectrum = rep(0.0, n), iterations = 0L,
                    converged = TRUE))
    }
    L <- NULL
    if (regularization > 0.0) {
        L <- make_regularization_operator(n, smoothness_order,
                                          identity_for_zero = FALSE)
    }
    AT <- t(A)
    r <- b - as.numeric(A %*% x)
    s <- as.numeric(AT %*% r)
    if (!is.null(L)) s <- s - regularization * as.numeric(t(L) %*% (L %*% x))
    d <- s
    nrmAtb <- sqrt(sum(as.numeric(AT %*% b)^2))
    rho <- sum(s * s)
    rtol <- if (!is.null(noise_level) && noise_level >= 0) {
        1.01 * noise_level * nrmb
    } else NULL
    iterations <- 0L; converged <- FALSE
    for (k in seq_len(max_iterations)) {
        Ad <- as.numeric(A %*% d)
        if (!is.null(L)) {
            Ld <- as.numeric(L %*% d)
            normAd2 <- sum(Ad * Ad) + regularization^2 * sum(Ld * Ld)
        } else {
            normAd2 <- sum(Ad * Ad)
        }
        if (normAd2 <= 0) break
        alpha_k <- rho / normAd2
        x <- x + alpha_k * d
        r <- r - alpha_k * Ad
        if (!is.null(L)) {
            s <- as.numeric(AT %*% r)
            Lx <- as.numeric(L %*% x)
            s <- s - regularization * as.numeric(t(L) %*% Lx)
        } else {
            s <- as.numeric(AT %*% r)
        }
        rho_new <- sum(s * s)
        beta <- if (rho > 0) rho_new / rho else 0.0
        rho <- rho_new
        d <- s + beta * d
        iterations <- k
        ne_res <- sqrt(sum(s * s))
        if (!is.null(rtol)) {
            res <- sqrt(sum(r * r))
            if (res <= rtol) { converged <- TRUE; break }
        }
        if (nrmAtb > 0 && ne_res <= tolerance * nrmAtb) {
            converged <- TRUE; break
        }
        if (ne_res <= tolerance) { converged <- TRUE; break }
    }
    list(spectrum = pmax(x, 0), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_cgls}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_cgls
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_cgls <- function(detector_names, n_energy_bins, E_MeV,
                         sensitivities, cc_icrp116, save_result_callback,
                         readings, initial_spectrum = NULL,
                         max_iterations = 100L, tolerance = 1e-12,
                         noise_level = NULL, regularization = 0.0,
                         smoothness_order = 0L,
                         calculate_errors = FALSE,
                         n_montecarlo = 100L,
                         save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_cgls,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         noise_level = noise_level,
                                         regularization = regularization,
                                         smoothness_order = smoothness_order),
        solve_kwargs = list(),
        method_name = "CGLS",
        extra_output = list(max_iterations = max_iterations,
                            regularization = regularization,
                            smoothness_order = as.integer(smoothness_order)),
        calculate_errors = calculate_errors,
        noise_level = if (is.null(noise_level)) 0.01 else noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
