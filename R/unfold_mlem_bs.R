#' MLEM-BS (B-spline MLEM with sieve regularization)
#'
#' R port (simplified) of \code{bssunfold/src/bssunfold/core/unfold_mlem_bs.py}.
#' Implements the regularized MLEM iteration on B-spline coefficients of
#' Mazankova et al., CNDGS'2026, \url{https://doi.org/10.47459/cndcgs.2026.61}.
#'
#' The spectrum is represented as \eqn{x(E) = \sum_s b_s B_s(E)} so the
#' effective system matrix is \eqn{RB = R B}. The B-spline coefficients
#' are updated by the multiplicative MLEM iteration with a second-derivative
#' penalty \eqn{P(b) = \|D^{(2)} b\|^2_2} and a sieve restriction to
#' non-negative coefficients.
#'
#' @section Limitations vs. Python bssunfold:
#' This is a simplified port that does not include the Poisson-bootstrap
#' confidence intervals nor the auto-selection of \eqn{N_s} / \eqn{\beta}
#' by minimisation of the \eqn{K_S} statistic. Users must specify
#' \code{n_basis} and \code{beta} (or \code{beta_relative}) manually.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param E_MeV Numeric energy grid (length n).
#' @param n_basis Integer; number of B-spline basis functions. Default
#'   \code{NULL} = \code{max(5, n/4)}.
#' @param beta Numeric; absolute penalty strength. Exactly one of
#'   \code{beta} / \code{beta_relative} may be non-zero. Default 0.
#' @param beta_relative Numeric; penalty strength relative to the mean
#'   column sum of \eqn{RB}. Default 0.
#' @param max_iterations Positive integer; default 200.
#' @param tolerance Positive numeric; default 1e-6.
#' @param knot_spacing Character; \code{"auto"}, \code{"uniform"}, or
#'   \code{"log"}. Default \code{"auto"} = log when span > 2 decades.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_mlem_bs(A, b, rep(1, 60), E, n_basis = 10L,
#'                    beta_relative = 1e-3, max_iterations = 50L)
solve_mlem_bs <- function(A, b, x0, E_MeV, n_basis = NULL,
                            beta = 0.0, beta_relative = 0.0,
                            max_iterations = 200L, tolerance = 1e-6,
                            knot_spacing = "auto") {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0); E_MeV <- as.numeric(E_MeV)
    n <- length(E_MeV)
    if (is.null(n_basis)) n_basis <- max(5L, n %/% 4L)
    if (n_basis < 4L) {
        stop("n_basis must be at least 4 for cubic B-spline basis")
    }
    if (beta != 0 && beta_relative != 0) {
        stop("Specify either beta or beta_relative, not both")
    }
    spacing <- tolower(knot_spacing)
    if (spacing == "auto") {
        e_range <- log10(max(E_MeV)) - log10(min(E_MeV))
        spacing <- if (e_range > 2.0) "log" else "uniform"
    }
    B <- .mlem_bs_basis(E_MeV, n_basis, spacing)  # n x n_basis
    RB <- A %*% B  # m x n_basis
    if (beta_relative > 0) {
        s_bar <- mean(colSums(RB))
        beta_eff <- beta_relative * s_bar
    } else {
        beta_eff <- beta
    }
    D2 <- as.matrix(create_derivative_matrix(n_basis, 2L))  # (n_basis-2, n_basis)
    D2T_D2 <- crossprod(D2)  # n_basis x n_basis
    # Initial coefficients: positive least-squares fit to x0
    b_init <- tryCatch(
        as.numeric(lsei::nnls(B, pmax(x0, 1e-10))$x),
        error = function(e) rep(1.0, n_basis)
    )
    b_coef <- pmax(b_init, 1e-10)
    RB_colsums <- colSums(RB)
    converged <- FALSE; iterations <- 0L
    eps <- 1e-300
    for (k in seq_len(max_iterations)) {
        iterations <- k
        b_old <- b_coef
        Ax <- as.numeric(RB %*% b_coef)
        Ax_safe <- pmax(Ax, eps)
        ratio <- b / Ax_safe
        correction <- as.numeric(t(RB) %*% ratio)
        if (beta_eff > 0) {
            grad_P <- 2.0 * as.numeric(D2T_D2 %*% b_coef)
        } else {
            grad_P <- 0.0
        }
        denom <- RB_colsums + beta_eff * grad_P + eps
        b_coef <- pmax(b_coef * correction / denom, 0.0)
        rel <- sqrt(sum((b_coef - b_old)^2)) / (sqrt(sum(b_old^2)) + eps)
        if (rel < tolerance) { converged <- TRUE; break }
    }
    spectrum <- as.numeric(B %*% b_coef)
    list(spectrum = pmax(spectrum, 0.0), iterations = iterations,
         converged = converged)
}

.mlem_bs_basis <- function(E, n_basis, spacing) {
    # Clamped B-spline basis of degree 3 with n_basis basis functions.
    # Knot vector length = n_basis + degree + 1; clamped = (degree+1)
    # copies of the endpoints + (n_basis - degree - 1) interior knots.
    degree <- 3L
    n <- length(E)
    e_min <- min(E); e_max <- max(E)
    n_interior <- n_basis - degree - 1L
    if (n_interior > 0L) {
        if (spacing == "log") {
            t_seq <- 10^seq(log10(e_min), log10(e_max), length.out = n_interior + 2L)
        } else {
            t_seq <- seq(e_min, e_max, length.out = n_interior + 2L)
        }
        t_int <- t_seq[-c(1L, n_interior + 2L)]  # exclude endpoints
    } else {
        t_int <- numeric(0L)
    }
    knots <- c(rep(e_min, degree + 1L), t_int, rep(e_max, degree + 1L))
    B <- matrix(0.0, nrow = n, ncol = n_basis)
    for (i in seq_len(n)) {
        B[i, ] <- .mlem_bs_basis_eval(E[i], knots, degree, n_basis)
    }
    B
}

.mlem_bs_basis_eval <- function(x, knots, degree, n_basis) {
    # Clean Cox-de Boor recursion. Returns the n_basis-long vector of basis
    # function values N_i(x) at point x. Partition of unity holds for
    # clamped B-splines at every interior point.
    n_knots <- length(knots)
    # Find span index i (1-based) such that knots[i] <= x < knots[i+1],
    # and i in [degree+1, n_basis] (interior span).
    if (x >= knots[n_knots]) {
        i <- n_basis
    } else if (x <= knots[1L]) {
        i <- degree + 1L
    } else {
        # Binary search could be used; linear is fine for small knot vectors
        i <- degree + 1L
        for (k in (degree + 1L):(n_basis)) {
            if (knots[k] <= x && x < knots[k + 1L]) { i <- k; break }
        }
    }
    # Degree-0 basis functions (only the one containing x is non-zero)
    N0 <- numeric(degree + 1L)
    for (j in 0:degree) {
        kk <- i - degree + j  # 1-based knot index
        if (x >= knots[kk] && x < knots[kk + 1L]) N0[j + 1L] <- 1.0
    }
    # Special case: x at the right endpoint
    if (x == knots[n_knots]) N0[degree + 1L] <- 1.0
    # Cox-de Boor recursion. At degree 0, only N0[degree+1] (= N_{i, 0}) is
    # 1. At degree d, the non-zero offsets are j = (degree-d)..degree (the
    # LAST d+1 entries of N). We use a separate N_prev snapshot to avoid
    # overwriting d-1 values we still need.
    N <- N0
    for (d in 1:degree) {
        N_prev <- N  # snapshot of d-1 values
        for (j in (degree - d):degree) {
            kk <- i - degree + j  # 1-based knot index
            denom_a <- knots[kk + d] - knots[kk]
            denom_b <- knots[kk + d + 1L] - knots[kk + 1L]
            term_a <- if (denom_a > 0) (x - knots[kk]) / denom_a * N_prev[j + 1L] else 0
            right_N <- if (j + 2L <= degree + 1L) N_prev[j + 2L] else 0
            term_b <- if (denom_b > 0) (knots[kk + d + 1L] - x) / denom_b * right_N else 0
            N[j + 1L] <- term_a + term_b
        }
        # Zero out the FIRST (degree-d) entries — they fall outside the
        # active support at degree d.
        if (degree - d >= 1L) {
            for (j in 0:(degree - d - 1L)) {
                N[j + 1L] <- 0
            }
        }
    }
    # Map N[1..degree+1] to basis vector at positions (i-degree)..(i)
    basis <- numeric(n_basis)
    for (j in 0:degree) {
        idx <- i - degree + j  # 1-based position in basis vector
        # idx is in [i-degree, i] which should be in [1, n_basis]
        if (idx >= 1L && idx <= n_basis) basis[idx] <- N[j + 1L]
    }
    basis
}

#' Build B-spline basis matrix (clamped cubic)
#'
#' @param E_MeV Numeric energy grid.
#' @param n_basis Integer; number of basis functions.
#' @param knot_spacing Character; \code{"uniform"} or \code{"log"}.
#' @return Numeric matrix of shape \code{(length(E_MeV), n_basis)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' B <- build_bspline_basis(E, 10L, "log")
#' dim(B)
build_bspline_basis <- function(E_MeV, n_basis, knot_spacing = "auto") {
    E_MeV <- as.numeric(E_MeV); n_basis <- as.integer(n_basis)
    if (n_basis < 4L) stop("n_basis must be at least 4")
    spacing <- tolower(knot_spacing)
    if (spacing == "auto") {
        e_range <- log10(max(E_MeV)) - log10(min(E_MeV))
        spacing <- if (e_range > 2.0) "log" else "uniform"
    }
    .mlem_bs_basis(E_MeV, n_basis, spacing)
}

#' Second-difference matrix (for MLEM-BS penalty)
#'
#' Returns the second-difference operator of shape \eqn{(n-2, n)} with
#' stencil \code{[1, -2, 1]}.
#'
#' @param n Integer; length of the coefficient vector.
#' @return Numeric matrix of shape \eqn{(n-2, n)}.
#' @export
#' @examples
#' D2 <- second_difference_matrix(6L)
#' dim(D2)
second_difference_matrix <- function(n) {
    as.matrix(create_derivative_matrix(n, 2L))
}

#' K_S goodness-of-fit statistic
#'
#' \eqn{K_S = | sum_i (n_i - model_i)^2 / sum_i model_i - 1 |}
#'
#' @param measurements Numeric vector.
#' @param estimate Numeric vector (model prediction).
#' @return Numeric scalar.
#' @export
#' @examples
#' ks_statistic(c(1, 2, 3), c(1.1, 1.9, 3.1))
ks_statistic <- function(measurements, estimate) {
    num <- sum((as.numeric(measurements) - as.numeric(estimate))^2)
    denom <- sum(as.numeric(estimate))
    if (denom <= 0) return(Inf)
    abs(num / denom - 1.0)
}

#' Wrapper around \code{\link{solve_mlem_bs}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_mlem_bs
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_mlem_bs <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            n_basis = NULL, beta = 0.0,
                            beta_relative = 0.0,
                            max_iterations = 200L, tolerance = 1e-6,
                            knot_spacing = "auto",
                            calculate_errors = FALSE,
                            noise_level = 0.01, n_montecarlo = 100L,
                            save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.5, n_energy_bins)
    # Use closure to forward E_MeV
    solver_with_E <- function(A, b, x0 = NULL, ...) {
        solve_mlem_bs(A = A, b = b, x0 = if (is.null(x0)) x0_default else x0,
                      E_MeV = E_MeV, n_basis = n_basis,
                      beta = beta, beta_relative = beta_relative,
                      max_iterations = max_iterations, tolerance = tolerance,
                      knot_spacing = knot_spacing)
    }
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver_with_E,
        solve_kwargs = list(),
        method_name = "MLEM-BS",
        extra_output = list(n_basis = if (is.null(n_basis))
                                           max(5L, n_energy_bins %/% 4L)
                                       else as.integer(n_basis),
                            beta = beta, beta_relative = beta_relative,
                            knot_spacing = knot_spacing),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
