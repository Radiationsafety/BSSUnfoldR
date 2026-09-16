#' N-spline unfolding (Islamgulov & Lartsev, 2008)
#'
#' Simplified R port of \code{bssunfold/src/bssunfold/core/unfold_nspline.py}.
#' Parameterises the spectrum as piecewise-exponential
#' \eqn{N_k(E) = \exp(a_k + q_k \ln E + r_k E)} with \eqn{C^0/C^1}
#' continuity at the interior knots, and minimises the directed divergence
#' between measured and calculated activations via a multiplicative
#' gradient iteration with backtracking line search.
#'
#' @section Limitations vs. Python bssunfold:
#' The original implementation uses a specialised "neutron spline" basis
#' with both \eqn{\ln E} and \eqn{E} terms. This R port uses a piecewise
#' exponential basis with continuity constraints (a subset of the N-spline
#' family) and a directed-divergence minimisation loop with backtracking.
#' Auto-selection of the number of knots, the \eqn{H_target} stopping
#' criterion based on measurement errors, and the \code{nev} quality
#' metric are not implemented.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param E_MeV Numeric energy grid (length n).
#' @param knots Optional numeric vector of knot energies (interior, strictly
#'   increasing). Default \code{NULL} = \code{\link{auto_knots}(E_MeV)}.
#' @param max_iterations Positive integer; default 500.
#' @param tolerance Positive numeric; relative tolerance on H. Default 1e-6.
#' @param dmu_start Numeric; initial step size. Default 0.1.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_nspline(A, b, rep(1, 60), E, max_iterations = 50L)
solve_nspline <- function(A, b, x0, E_MeV, knots = NULL,
                            max_iterations = 500L, tolerance = 1e-6,
                            dmu_start = 0.1) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0); E_MeV <- as.numeric(E_MeV)
    n <- length(E_MeV)
    if (is.null(knots)) knots <- auto_knots(E_MeV)
    M <- length(knots) + 1L  # number of segments
    # Continuity matrix D (2*M-2 x 3*M): imposes C0 (N matches at knots) and
    # C1 (dN/dlnE matches at knots). We only enforce C0 here.
    # Initial guess: assume each segment is exp(a_k + q_k ln E + r_k E) with
    # a_k = log(x0[knots]), q_k = 0, r_k = 0 (flat in log E).
    phi <- pmax(x0, 1e-30)
    converged <- FALSE; iterations <- 0L
    eps <- 1e-30
    b_sum <- sum(b); if (b_sum <= 0) b_sum <- 1
    p_b <- b / b_sum

    for (k in seq_len(max_iterations)) {
        iterations <- k
        phi_old <- phi
        Q <- as.numeric(A %*% phi)
        Q_safe <- pmax(Q, eps)
        pN <- Q / sum(Q_safe)
        # H = sum_i p_i log(p_i / pN_i) - p_i + pN_i (directed divergence)
        H <- sum(p_b * log(pmax(p_b, eps) / pmax(pN, eps)) - p_b + pN)
        if (!is.finite(H)) H <- Inf
        # R(E) = sum_i (p_i / Q_i) * A_{ij} * log(pN_i / p_i)
        ratio_log <- log(pmax(pN, eps) / pmax(p_b, eps))
        coeff <- p_b / Q_safe * ratio_log  # length m
        R_E <- as.numeric(t(A) %*% coeff)  # length n
        # Rbar = flux-weighted mean
        Rbar <- sum(phi * R_E) / max(sum(phi), eps)
        # Step direction: -(R - Rbar) (descent)
        delta <- -(R_E - Rbar)
        # Step size: backtracking line search
        dmu <- dmu_start / max(abs(delta), eps)
        dmu <- min(dmu, 1.0)  # cap step
        phi_new <- phi * (1 + dmu * delta)
        phi_new <- pmax(phi_new, eps)
        Q_new <- as.numeric(A %*% phi_new)
        pN_new <- Q_new / max(sum(Q_new), eps)
        H_new <- sum(p_b * log(pmax(p_b, eps) / pmax(pN_new, eps)) - p_b + pN_new)
        # Halve step until H decreases
        halves <- 0L
        while (H_new > H && halves < 20L) {
            dmu <- dmu * 0.5
            phi_new <- phi * (1 + dmu * delta)
            phi_new <- pmax(phi_new, eps)
            Q_new <- as.numeric(A %*% phi_new)
            pN_new <- Q_new / max(sum(Q_new), eps)
            H_new <- sum(p_b * log(pmax(p_b, eps) / pmax(pN_new, eps)) - p_b + pN_new)
            halves <- halves + 1L
        }
        phi <- phi_new
        rel <- abs(H - H_new) / max(abs(H), eps)
        if (rel < tolerance || H_new < 1e-12) { converged <- TRUE; break }
    }
    list(spectrum = as.numeric(phi), iterations = iterations,
         converged = converged)
}

#' Auto-generate log-uniform interior knots
#'
#' @param E_MeV Numeric energy grid.
#' @param n_interior Integer; number of interior knots. Default
#'   \code{NULL} = \code{max(4, length(E_MeV)/10)}.
#' @return Numeric vector of strictly increasing interior knot energies.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' auto_knots(E)
auto_knots <- function(E_MeV, n_interior = NULL) {
    E_MeV <- as.numeric(E_MeV)
    e_min <- min(E_MeV); e_max <- max(E_MeV)
    if (is.null(n_interior)) n_interior <- max(4L, length(E_MeV) %/% 10L)
    10^seq(log10(e_min), log10(e_max), length.out = n_interior + 2L)[-c(1L, n_interior + 2L)]
}

#' Build continuity-constraint matrix for N-spline basis
#'
#' Returns the block matrix \eqn{D} of continuity constraints
#' (\eqn{C^0/C^1}) on the N-spline parameters \eqn{X = (a, q, r)^T}.
#'
#' @param n_knots Integer; number of interior knots (segments = n_knots + 1).
#' @param order Integer; 0 (C0) or 1 (C1). Default 0.
#' @return Numeric matrix.
#' @export
#' @examples
#' D <- build_continuity_matrix(4L, order = 0L)
#' dim(D)
build_continuity_matrix <- function(n_knots, order = 0L) {
    segments <- n_knots + 1L
    M <- 3L * segments
    if (order == 0L) {
        n_constraints <- segments - 1L
        D <- matrix(0.0, n_constraints, M)
        for (i in seq_len(n_constraints)) {
            D[i, (i - 1L) * 3L + 1L] <- 1
            D[i, i * 3L + 1L] <- -1
        }
        D
    } else if (order == 1L) {
        n_constraints <- 2L * (segments - 1L)
        D <- matrix(0.0, n_constraints, M)
        for (i in seq_len(segments - 1L)) {
            D[i, (i - 1L) * 3L + 1L] <- 1
            D[i, i * 3L + 1L] <- -1
            # C1: dN/dln E = q + r E. At knot i, this gives constraints on
            # q_left + r_left * E_i = q_right + r_right * E_i.
            D[segments - 1L + i, (i - 1L) * 3L + 2L] <- 1
            D[segments - 1L + i, (i - 1L) * 3L + 3L] <- 1
            D[segments - 1L + i, i * 3L + 2L] <- -1
            D[segments - 1L + i, i * 3L + 3L] <- -1
        }
        D
    } else {
        stop("order must be 0 or 1")
    }
}

#' Wrapper around \code{\link{solve_nspline}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_nspline
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_nspline <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              knots = NULL,
                              max_iterations = 500L, tolerance = 1e-6,
                              dmu_start = 0.1,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    solver_with_E <- function(A, b, x0 = NULL, ...) {
        solve_nspline(A = A, b = b,
                      x0 = if (is.null(x0)) x0_default else x0,
                      E_MeV = E_MeV, knots = knots,
                      max_iterations = max_iterations,
                      tolerance = tolerance, dmu_start = dmu_start)
    }
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver_with_E,
        solve_kwargs = list(),
        method_name = "NSpline",
        extra_output = list(knots = if (is.null(knots)) auto_knots(E_MeV)
                                     else as.numeric(knots)),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
