#' ODL-style PDHG and Douglas-Rachford unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_odl_advanced.py},
#' exposed through \code{Detector.unfold_odl_pdhg} and
#' \code{Detector.unfold_odl_douglas_rachford}.
#'
#' The Python module follows the Operator Discretization Library (ODL)
#' formulation of the two primal-dual / operator-splitting schemes, but is
#' implemented in pure NumPy (ODL 1.0's own solvers break on translated data
#' terms).  These are NOT the same algorithms as \code{\link{solve_pdhg}} and
#' \code{\link{solve_douglas_rachford}}: the ODL variants augment the operator
#' with the TV block (\eqn{K = [A;\; w D]}, \eqn{d = [b;\; 0]}) and use
#' Chambolle-Pock's simplified explicit update, a fixed relaxation
#' \eqn{\gamma = 1}, Chambolle's dual gradient-ascent TV prox and the
#' "residual did not get worse" stopping rule.
#'
#' @name odl-advanced-methods
NULL

# 1-D forward-difference operator D of shape (n - 1) x n, exactly as
# _forward_diff_matrix() in unfold_odl_advanced.py.
.odl_forward_diff <- function(n) {
    n <- as.integer(n)
    D <- matrix(0.0, nrow = max(n - 1L, 0L), ncol = n)
    if (n > 1L) {
        idx <- seq_len(n - 1L)
        D[cbind(idx, idx)] <- -1.0
        D[cbind(idx, idx + 1L)] <- 1.0
    }
    D
}

# D^T p for the forward-difference operator (numpy builds the same vector by
# grad[1:] = p; grad[:-1] -= p).
.odl_dt_apply <- function(p, n) {
    grad <- numeric(n)
    if (n > 1L) {
        grad[-1L] <- p
        grad[-n] <- grad[-n] - p
    }
    grad
}

# Proximal operator of lam * ||D x||_1 (1-D total-variation denoising),
# solved with Chambolle's dual gradient-ascent algorithm:
#   argmin_x 0.5 ||x - f||^2 + lam ||D x||_1
# Ported literally from _tv_prox(), including the fixed spectral bound L = 4,
# the step rho = 1 / (L lam^2), the dual clip to [-1, 1] and the 200 default
# iterations.
.odl_tv_prox <- function(f, lam, n_iter = 200L) {
    f <- as.numeric(f)
    n <- length(f)
    if (n <= 1L || !(lam > 0)) return(f)
    Df <- as.numeric(diff(f))
    L <- 4.0
    rho <- 1.0 / (L * lam^2)
    p <- numeric(n - 1L)
    for (k in seq_len(as.integer(n_iter))) {
        grad <- .odl_dt_apply(p, n)
        p <- p + rho * (lam * Df - lam^2 * (grad[-1L] - grad[-n]))
        p <- pmin(pmax(p, -1.0), 1.0)
    }
    f - lam * .odl_dt_apply(p, n)
}

# Largest eigenvalue of K' K, mirroring np.linalg.eigvalsh(K.T @ K).max().
.odl_max_eig_gram <- function(K) {
    vals <- eigen(crossprod(K), symmetric = TRUE, only.values = TRUE)$values
    max(as.numeric(vals))
}

#' Solve by the ODL formulation of PDHG (Chambolle-Pock)
#'
#' Solves \eqn{\min_x 0.5 \|A x - b\|^2 + tv\_weight \|D x\|_1}{min 0.5||Ax-b||^2 + tv_weight ||Dx||_1}
#' with the primal-dual hybrid gradient scheme used by
#' \code{bssunfold.core.unfold_odl_advanced.solve_odl_pdhg}.  The TV term is
#' folded into an augmented operator \eqn{K = [A;\; w D]}{K = [A; w D]} with
#' augmented data \eqn{d = [b;\; 0]}{d = [b; 0]}, so the dual update is the
#' closed-form prox of the quadratic data term; only the non-negativity
#' indicator is handled by a primal projection.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional numeric initial spectrum (length n).  Default \code{NULL}
#'   = \code{rep(0.5, n)}, matching Python's \code{np.ones(n) * 0.5}.
#' @param max_iterations Positive integer; the iteration count is not adaptive
#'   (the loop always runs to completion). Default \code{100L}.
#' @param tau Numeric primal step size, or \code{NULL} for the automatic
#'   \code{0.99 / ||K||}. Default \code{NULL}.
#' @param sigma Numeric dual step size, or \code{NULL} for the automatic
#'   \code{0.99 / ||K||}. Default \code{NULL}.
#' @param use_tv Logical; include the TV block in the augmented operator.
#'   Default \code{TRUE}.
#' @param tv_weight Numeric weight of the TV term. Default 0.1.
#' @param nonnegativity Logical; project onto \eqn{x \ge 0}{x >= 0} inside the
#'   loop and on the returned spectrum. Default \code{TRUE}.
#' @param tolerance Numeric; Python declares convergence when the relative
#'   residual did not worsen by more than \code{1 + tolerance * 100}.
#'   Default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_odl_pdhg(A, b, rep(0.5, 3), max_iterations = 50L)
solve_odl_pdhg <- function(A, b, x0 = NULL, max_iterations = 100L,
                           tau = NULL, sigma = NULL,
                           use_tv = TRUE, tv_weight = 0.1,
                           nonnegativity = TRUE, tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (is.null(x0)) x0 <- rep(0.5, n)
    x <- as.numeric(x0)

    # Augmented operator K = [A; w D] and data d = [b; 0] so that
    #     0.5 ||A x - b||^2 + tv_weight ||D x||_1  =  0.5 ||K x - d||^2
    # in the primal-dual sense.
    if (isTRUE(use_tv)) {
        D <- .odl_forward_diff(n)
        K <- rbind(A, as.numeric(tv_weight) * D)
        d <- c(b, numeric(max(n - 1L, 0L)))
    } else {
        K <- A
        d <- b
    }

    eig_max <- .odl_max_eig_gram(K)
    op_norm <- sqrt(max(eig_max, 1e-12))
    if (is.null(tau)) tau <- 0.99 / op_norm
    if (is.null(sigma)) sigma <- 0.99 / op_norm
    tau <- as.numeric(tau); sigma <- as.numeric(sigma)

    b_norm <- max(sqrt(sum(b^2)), 1e-300)
    res_before <- sqrt(sum((as.numeric(A %*% x) - b)^2)) / b_norm

    y <- numeric(nrow(K))
    x_bar <- x
    for (k in seq_len(as.integer(max_iterations))) {
        z <- y + sigma * as.numeric(K %*% x_bar)
        y <- (z - sigma * d) / (1.0 + sigma)
        x_new <- x - tau * as.numeric(t(K) %*% y)
        if (isTRUE(nonnegativity)) x_new <- pmax(x_new, 0.0)
        x_bar <- x_new + (x_new - x)
        x <- x_new
    }

    x_opt <- if (isTRUE(nonnegativity)) pmax(x, 0.0) else x
    finite <- all(is.finite(x_opt))
    res_after <- sqrt(sum((as.numeric(A %*% x_opt) - b)^2)) / b_norm
    converged <- finite && res_after <= res_before * (1.0 + tolerance * 100)

    list(spectrum = as.numeric(x_opt),
         iterations = as.integer(max_iterations),
         converged = converged)
}

#' Solve by the ODL formulation of Douglas-Rachford splitting
#'
#' Douglas-Rachford splitting for \eqn{\min_x \psi_1(x) + \psi_2(x)} with
#' \eqn{\psi_1(x) = 0.5 \|A x - b\|^2}{psi1 = 0.5||Ax-b||^2} (prox = precomputed
#' linear solve) and \eqn{\psi_2(x) = tv\_weight \|D x\|_1}{psi2 = tv_weight ||Dx||_1}
#' (prox = 1-D TV denoising, Chambolle's dual algorithm).  Ported literally
#' from \code{bssunfold.core.unfold_odl_advanced.solve_odl_douglas_rachford}:
#' the relaxation is fixed at \eqn{\gamma = 1}{gamma = 1}, non-negativity is
#' applied only as a final clamp of the \eqn{y}-iterate, and the loop always
#' runs for the full iteration count.
#'
#' @inheritParams solve_odl_pdhg
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_odl_douglas_rachford(A, b, rep(0.5, 3), max_iterations = 50L)
solve_odl_douglas_rachford <- function(A, b, x0 = NULL, max_iterations = 100L,
                                       use_tv = TRUE, tv_weight = 0.1,
                                       nonnegativity = TRUE,
                                       tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (is.null(x0)) x0 <- rep(0.5, n)
    x <- as.numeric(x0)

    gamma <- 1.0
    Atb <- as.numeric(t(A) %*% b)
    Mt <- diag(n) + gamma * crossprod(A)
    # prox_{psi1}(v) = (I + gamma A'A)^{-1} (v + gamma A'b); Python factors the
    # matrix once with a Cholesky decomposition and reuses it every iteration.
    R_chol <- tryCatch(chol(as.matrix(Mt)), error = function(e) NULL)
    prox_psi1 <- function(v) {
        rhs <- as.numeric(v) + gamma * Atb
        if (is.null(R_chol)) {
            return(as.numeric(solve(Mt, rhs, tol = 1e-12)))
        }
        as.numeric(backsolve(R_chol,
                             forwardsolve(t(R_chol), rhs,
                                          upper.tri = FALSE,
                                          transpose = FALSE),
                             transpose = FALSE))
    }
    prox_psi2 <- function(v) {
        if (isTRUE(use_tv)) .odl_tv_prox(v, gamma * as.numeric(tv_weight))
        else as.numeric(v)
    }

    b_norm <- max(sqrt(sum(b^2)), 1e-300)
    res_before <- sqrt(sum((as.numeric(A %*% x) - b)^2)) / b_norm

    y <- x
    z <- x
    for (k in seq_len(as.integer(max_iterations))) {
        u <- prox_psi1(2.0 * z - y)
        y <- prox_psi2(u)
        z <- z + y - u
    }

    x_opt <- if (isTRUE(nonnegativity)) pmax(y, 0.0) else y
    finite <- all(is.finite(x_opt))
    res_after <- sqrt(sum((as.numeric(A %*% x_opt) - b)^2)) / b_norm
    converged <- finite && res_after <= res_before * (1.0 + tolerance * 100)

    list(spectrum = as.numeric(x_opt),
         iterations = as.integer(max_iterations),
         converged = converged)
}

#' Wrapper around \code{\link{solve_odl_pdhg}} for the unified workflow
#'
#' Thin wrapper matching \code{bssunfold}'s
#' \code{Detector.unfold_odl_pdhg} / \code{unfold_odl_pdhg}: it builds the
#' response system through \code{\link{run_unfolding}} and hands it to the ODL
#' PDHG solver with a default initial spectrum of \code{0.5} in every bin.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_odl_pdhg
#' @param max_iterations Positive integer; default \code{100L}.
#' @param initial_spectrum Optional numeric initial spectrum (length
#'   \code{n_energy_bins}); default \code{NULL} = \code{rep(0.5, n_energy_bins)}.
#' @param calculate_errors Logical; run Monte-Carlo uncertainty estimation.
#'   Default \code{FALSE}.
#' @param noise_level Numeric; relative Gaussian noise for Monte-Carlo.
#'   Default 0.01.
#' @param n_montecarlo Integer; number of Monte-Carlo samples. Default 100.
#' @param save_result Logical; append the result to the detector history.
#'   Default \code{FALSE}.
#' @param random_state Optional integer seed for Monte-Carlo. Default \code{NULL}.
#' @param max_neutron_energy Optional energy cutoff in MeV. Default \code{NULL}.
#' @return A result list as produced by \code{\link{run_unfolding}} with
#'   \code{method = "ODL-PDHG"}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' det_names <- c("d1", "d2", "d3")
#' sens <- stats::setNames(lapply(seq_len(3), function(i) A[i, ]), det_names)
#' res <- unfold_odl_pdhg(det_names, 3L, c(1e-9, 1e-6, 1e-3), sens, NULL,
#'                        NULL, c(d1 = 1, d2 = 0.6, d3 = 0.4),
#'                        max_iterations = 50L)
unfold_odl_pdhg <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            max_iterations = 100L, tau = NULL, sigma = NULL,
                            use_tv = TRUE, tv_weight = 0.1,
                            nonnegativity = TRUE,
                            calculate_errors = FALSE, noise_level = 0.01,
                            n_montecarlo = 100L, save_result = FALSE,
                            random_state = NULL,
                            max_neutron_energy = NULL) {
    x0_default <- rep(1.0, n_energy_bins) * 0.5
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_odl_pdhg,
                                        max_iterations = max_iterations,
                                        tau = tau, sigma = sigma,
                                        use_tv = use_tv,
                                        tv_weight = tv_weight,
                                        nonnegativity = nonnegativity),
        solve_kwargs = list(),
        method_name = "ODL-PDHG",
        extra_output = list(max_iterations = as.integer(max_iterations),
                            use_tv = use_tv, tv_weight = tv_weight),
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}

#' Wrapper around \code{\link{solve_odl_douglas_rachford}} for the unified
#' workflow
#'
#' Thin wrapper matching \code{bssunfold}'s
#' \code{Detector.unfold_odl_douglas_rachford} /
#' \code{unfold_odl_douglas_rachford}.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_odl_douglas_rachford
#' @inheritParams unfold_odl_pdhg
#' @return A result list as produced by \code{\link{run_unfolding}} with
#'   \code{method = "ODL-DouglasRachford"}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' det_names <- c("d1", "d2", "d3")
#' sens <- stats::setNames(lapply(seq_len(3), function(i) A[i, ]), det_names)
#' res <- unfold_odl_douglas_rachford(det_names, 3L,
#'                                    c(1e-9, 1e-6, 1e-3), sens, NULL,
#'                                    NULL, c(d1 = 1, d2 = 0.6, d3 = 0.4),
#'                                    max_iterations = 50L)
unfold_odl_douglas_rachford <- function(detector_names, n_energy_bins, E_MeV,
                                        sensitivities, cc_icrp116,
                                        save_result_callback, readings,
                                        initial_spectrum = NULL,
                                        max_iterations = 100L,
                                        use_tv = TRUE, tv_weight = 0.1,
                                        nonnegativity = TRUE,
                                        calculate_errors = FALSE,
                                        noise_level = 0.01,
                                        n_montecarlo = 100L,
                                        save_result = FALSE,
                                        random_state = NULL,
                                        max_neutron_energy = NULL) {
    x0_default <- rep(1.0, n_energy_bins) * 0.5
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_odl_douglas_rachford,
                                        max_iterations = max_iterations,
                                        use_tv = use_tv,
                                        tv_weight = tv_weight,
                                        nonnegativity = nonnegativity),
        solve_kwargs = list(),
        method_name = "ODL-DouglasRachford",
        extra_output = list(max_iterations = as.integer(max_iterations),
                            use_tv = use_tv, tv_weight = tv_weight),
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
