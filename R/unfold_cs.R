#' Compressed Sensing (CS) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_cs.py}.
#' Implements a neutron spectrum unfolding method based on compressive
#' sensing. The spectrum \eqn{x} is represented sparsely in a dictionary
#' \eqn{D} as \eqn{x = D \alpha}. The measurement equation \eqn{b = A x}
#' becomes \eqn{b = (A D) \alpha}, which is solved for the sparse \eqn{\alpha}
#' using SL0. Finally the spectrum is reconstructed as \eqn{x = D \alpha}.
#'
#' This package provides the OMP and SL0 sub-solvers and the main CS wrapper.
#' K-SVD dictionary learning is also included but not used by the main wrapper
#' (which uses the trivial identity dictionary by default).
#'
#' @name cs-methods
#' @rdname cs-methods
NULL

#' Orthogonal Matching Pursuit (OMP)
#'
#' Finds a sparse coefficient vector \eqn{\alpha} such that \eqn{y \approx D \alpha}
#' with at most \code{sparsity} non-zero entries.
#'
#' @param D Numeric dictionary matrix (n x k).
#' @param y Numeric signal vector (length n).
#' @param sparsity Integer; maximum number of non-zero coefficients.
#' @param tolerance Numeric; residual tolerance for early stopping. Default 1e-6.
#' @return Numeric vector (length k).
#' @export
#' @examples
#' D <- matrix(runif(10 * 20), nrow = 10, ncol = 20)
#' y <- D[, 1] * 0.5 + D[, 5] * 0.3  # 2-sparse signal
#' alpha <- solve_omp(D, y, sparsity = 3)
solve_omp <- function(D, y, sparsity, tolerance = 1e-6) {
    D <- as.matrix(D); storage.mode(D) <- "double"
    y <- as.numeric(y)
    k <- ncol(D)
    alpha <- numeric(k)
    residual <- y
    support <- integer(0L)
    norms <- sqrt(colSums(D^2))
    norms <- ifelse(norms == 0, 1.0, norms)
    D_norm <- sweep(D, 2L, norms, "/")
    for (i in seq_len(min(sparsity, k))) {
        correlations <- abs(as.numeric(t(D_norm) %*% residual))
        if (length(support) > 0L) correlations[support] <- -1.0
        idx <- which.max(correlations)
        if (correlations[idx] <= 0) break
        support <- c(support, idx)
        D_s <- D[, support, drop = FALSE]
        coefs <- qr.solve(D_s, y, tol = 1e-12)
        residual <- y - as.numeric(D_s %*% coefs)
        if (sqrt(sum(residual^2)) < tolerance) break
    }
    if (length(support) > 0L) {
        D_s <- D[, support, drop = FALSE]
        coefs <- qr.solve(D_s, y, tol = 1e-12)
        alpha[support] <- coefs
    }
    alpha
}

#' K-SVD dictionary learning
#'
#' @param signals Numeric matrix (n x m), one column per training sample.
#' @param n_atoms Integer; number of dictionary atoms.
#' @param n_iterations Integer; number of K-SVD iterations. Default 20.
#' @param sparsity Integer; target sparsity for sparse coding. Default 5.
#' @param random_state Optional integer seed.
#' @return Numeric matrix (n x n_atoms).
#' @export
#' @examples
#' \dontrun{
#' signals <- matrix(runif(10 * 30), nrow = 10, ncol = 30)
#' D <- solve_ksvd(signals, n_atoms = 15, n_iterations = 5)
#' }
solve_ksvd <- function(signals, n_atoms, n_iterations = 20L, sparsity = 5L,
                         random_state = NULL) {
    if (!is.null(random_state)) set.seed(as.integer(random_state))
    signals <- as.matrix(signals); storage.mode(signals <- "double")
    n <- nrow(signals); m <- ncol(signals)
    n_atoms <- min(n_atoms, m)
    idx <- sample.int(m, size = n_atoms, replace = FALSE)
    D <- signals[, idx, drop = FALSE]
    norms <- sqrt(colSums(D^2))
    norms <- ifelse(norms == 0, 1.0, norms)
    D <- sweep(D, 2L, norms, "/")
    for (it in seq_len(n_iterations)) {
        coefficients <- matrix(0.0, nrow = n_atoms, ncol = m)
        for (j in seq_len(m)) {
            coefficients[, j] <- solve_omp(D, signals[, j], sparsity)
        }
        for (atom in seq_len(n_atoms)) {
            used <- which(coefficients[atom, ] != 0)
            if (length(used) == 0L) next
            D_restricted <- D
            D_restricted[, atom] <- 0.0
            E <- signals[, used, drop = FALSE] -
                 D_restricted %*% coefficients[, used, drop = FALSE]
            svdE <- svd(E, nu = min(dim(E)), nv = min(dim(E)))
            D[, atom] <- svdE$u[, 1L]
            coefficients[atom, used] <- svdE$d[1L] * svdE$v[, 1L]
        }
        norms <- sqrt(colSums(D^2))
        norms <- ifelse(norms == 0, 1.0, norms)
        D <- sweep(D, 2L, norms, "/")
    }
    D
}

#' Smoothed L0 (SL0) reconstruction
#'
#' Approximates the L0 norm by a smooth surrogate (Gaussian) and performs
#' a steepest-descent / projection iteration to find the sparsest solution
#' of the underdetermined linear system \eqn{b = A x}.
#'
#' @param A Numeric sensing matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param sigma_min Numeric; minimum sigma. Default 0.01.
#' @param sigma_decrease_factor Numeric; sigma decrease per outer iter. Default 0.5.
#' @param mu_0 Numeric; step-size factor. Default 1.0.
#' @param L Integer; inner steepest-descent iterations per sigma. Default 3.
#' @param max_iterations Integer; max outer iterations. Default 1000.
#' @param tolerance Numeric; convergence tolerance. Default 1e-6.
#' @return Numeric vector (length n).
#' @export
#' @examples
#' A <- matrix(runif(3 * 10), nrow = 3)
#' b <- as.numeric(A %*% c(1, 0, 0, 0, 2, 0, 0, 0, 0, 0))
#' x <- solve_sl0(A, b)
solve_sl0 <- function(A, b, sigma_min = 0.01, sigma_decrease_factor = 0.5,
                        mu_0 = 1.0, L = 3L, max_iterations = 1000L,
                        tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    AAT <- A %*% t(A)
    pinv_AAT <- tryCatch(solve(AAT), error = function(e) {
        sv <- svd(AAT, nu = nrow(AAT), nv = 0)
        s_inv <- ifelse(sv$d > 1e-10 * max(sv$d), 1 / sv$d, 0)
        sv$u %*% (s_inv * (t(sv$u)))
    })
    pinv_AT <- t(A) %*% pinv_AAT
    x <- as.numeric(pinv_AT %*% b)  # minimum-norm solution
    sigma <- 2.0 * max(abs(x))
    if (sigma == 0) sigma <- 1.0
    sigma <- max(sigma, sigma_min)
    for (iter in seq_len(max_iterations)) {
        x_prev <- x
        for (inner in seq_len(L)) {
            exp_term <- exp(-(x^2) / (2.0 * sigma^2))
            x <- x - mu_0 * x * exp_term
            x <- x - as.numeric(pinv_AT %*% (as.numeric(A %*% x) - b))
        }
        sigma <- sigma * sigma_decrease_factor
        if (sigma < sigma_min) break
        if (sqrt(sum((x - x_prev)^2)) < tolerance * max(1.0, sqrt(sum(x^2)))) break
    }
    x
}

#' Main CS solver
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused.
#' @param n_atoms Integer; number of dictionary atoms (for K-SVD). Default
#'   \code{NULL} = use identity dictionary.
#' @param sparsity Integer; sparsity for OMP/SL0. Default 5.
#' @param dictionary Optional pre-computed dictionary (n x n_atoms). Default
#'   \code{NULL} = identity (so \eqn{x = \alpha}).
#' @param sigma_min,sigma_decrease_factor,mu_0,L,max_iterations,tolerance
#'   SL0 parameters.
#' @param random_state Optional seed.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_cs(A, b, NULL, sparsity = 10L)
solve_cs <- function(A, b, x0 = NULL, n_atoms = NULL, sparsity = 5L,
                       dictionary = NULL, sigma_min = 0.01,
                       sigma_decrease_factor = 0.5, mu_0 = 1.0, L = 3L,
                       max_iterations = 1000L, tolerance = 1e-6,
                       random_state = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (is.null(dictionary)) {
        D <- diag(n)  # identity dictionary
    } else {
        D <- as.matrix(dictionary)
    }
    AD <- A %*% D
    alpha <- solve_sl0(AD, b, sigma_min = sigma_min,
                         sigma_decrease_factor = sigma_decrease_factor,
                         mu_0 = mu_0, L = L, max_iterations = max_iterations,
                         tolerance = tolerance)
    spectrum <- as.numeric(D %*% alpha)
    list(spectrum = pmax(spectrum, 0),
         iterations = as.integer(max_iterations),
         converged = TRUE,
         alpha = alpha)
}

#' Wrapper around \code{\link{solve_cs}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_cs
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_cs <- function(detector_names, n_energy_bins, E_MeV,
                        sensitivities, cc_icrp116, save_result_callback,
                        readings, initial_spectrum = NULL,
                        n_atoms = NULL, sparsity = 5L,
                        dictionary = NULL, sigma_min = 0.01,
                        sigma_decrease_factor = 0.5, mu_0 = 1.0, L = 3L,
                        max_iterations = 1000L, tolerance = 1e-6,
                        random_state = NULL,
                        calculate_errors = FALSE,
                        noise_level = 0.01, n_montecarlo = 100L,
                        save_result = FALSE) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_cs,
                                         n_atoms = n_atoms,
                                         sparsity = sparsity,
                                         dictionary = dictionary,
                                         sigma_min = sigma_min,
                                         sigma_decrease_factor = sigma_decrease_factor,
                                         mu_0 = mu_0, L = L,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         random_state = random_state),
        solve_kwargs = list(),
        method_name = "CS",
        extra_output = list(sparsity = sparsity, n_atoms = n_atoms),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
