#' Hybrid GMRES unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_hybrid_gmres.py}.
#' Combines the GMRES iterative solver with Tikhonov regularization applied
#' to the projected problem at each iteration. The regularization parameter
#' is selected automatically using GCV or discrepancy principle.
#'
#' @section Limitations:
#' This is a simplified port. The Python original uses GMRES with
#' reorthogonalization. The R port uses a similar approach via Lanczos /
#' Arnoldi-style iteration with optional Tikhonov regularization on the
#' projected problem (similar to Lanczos-hybrid but with a different
#' projection basis).
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused; accepted for API compatibility.
#' @param max_iterations Positive integer; max Krylov dimension. Default 100.
#' @param regularization_method Character; \code{"gcv"} (default), \code{"dp"},
#'   \code{"manual"}.
#' @param regularization Numeric; manual regularization. Default 0.0.
#' @param noise_level Optional numeric; relative noise for DP. Default \code{NULL}.
#' @param eta Numeric; safety factor for DP. Default 1.01.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_hybrid_gmres(A, b, NULL, max_iterations = 3)
solve_hybrid_gmres <- function(A, b, x0 = NULL, max_iterations = 100L,
                                  regularization_method = "gcv",
                                  regularization = 0.0,
                                  noise_level = NULL, eta = 1.01) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A); n <- ncol(A)
    max_k <- min(max_iterations, min(m, n))
    beta <- sqrt(sum(b^2))
    if (beta == 0) {
        return(list(spectrum = rep(0.0, n), iterations = 0L, converged = TRUE))
    }
    # Arnoldi iteration: build orthonormal Krylov basis V and Hessenberg H
    V <- matrix(0.0, nrow = m, ncol = max_k + 1L)
    H <- matrix(0.0, nrow = max_k + 1L, ncol = max_k)
    V[, 1L] <- b / beta
    best_x <- rep(0.0, n)
    converged <- FALSE; iterations <- 0L
    for (k in seq_len(max_k)) {
        # Arnoldi step
        w <- as.numeric(A %*% V[, k])
        for (j in seq_len(k)) {
            H[j, k] <- sum(w * V[, j])
            w <- w - H[j, k] * V[, j]
        }
        H[k + 1L, k] <- sqrt(sum(w^2))
        if (H[k + 1L, k] <= 1e-14) {
            converged <- TRUE
        } else {
            V[, k + 1L] <- w / H[k + 1L, k]
        }
        # Projected problem: min || H_k y - beta e_1 ||^2 + lambda ||y||^2
        H_k <- H[1:(k + 1L), 1:k, drop = FALSE]
        e1 <- c(beta, rep(0.0, k))
        # Select lambda
        if (regularization_method == "manual") {
            lam <- regularization
        } else if (regularization_method == "dp") {
            nl <- if (is.null(noise_level)) 0.01 else noise_level
            target <- eta * nl * sqrt(m)
            lambdas <- 10^seq(-12, 2, length.out = 50)
            res <- vapply(lambdas, function(lam) {
                lhs <- rbind(H_k, sqrt(lam) * diag(k))
                rhs <- c(e1, rep(0.0, k))
                y <- qr.solve(lhs, rhs)
                r <- as.numeric(A %*% (V[, 1:k] %*% y)) - b
                sqrt(sum(r^2))
            }, numeric(1))
            lam <- lambdas[which.min(abs(res - target))]
        } else {  # gcv
            sv <- svd(H_k, nu = min(dim(H_k)), nv = 0)
            s <- sv$d; s2 <- s^2
            c <- as.numeric(t(sv$u) %*% e1)
            lambdas <- 10^seq(-12, 2, length.out = 100)
            gcv <- vapply(lambdas, function(lam) {
                residual_sq <- sum((c * lam / (s2 + lam))^2)
                trace <- sum(s2 / (s2 + lam))
                denom <- (m - trace)^2
                if (denom > 0) residual_sq / denom else Inf
            }, numeric(1))
            lam <- lambdas[which.min(gcv)]
        }
        # Solve regularized least-squares
        lhs <- rbind(H_k, sqrt(lam) * diag(k))
        rhs <- c(e1, rep(0.0, k))
        y <- qr.solve(lhs, rhs)
        best_x <- as.numeric(V[, 1:k, drop = FALSE] %*% y)
        iterations <- k
        if (converged) break
    }
    list(spectrum = pmax(best_x, 0),
         iterations = as.integer(iterations),
         converged = converged)
}

#' Wrapper around \code{\link{solve_hybrid_gmres}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_hybrid_gmres
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_hybrid_gmres <- function(detector_names, n_energy_bins, E_MeV,
                                   sensitivities, cc_icrp116,
                                   save_result_callback, readings,
                                   initial_spectrum = NULL,
                                   max_iterations = 100L,
                                   regularization_method = "gcv",
                                   regularization = 0.0,
                                   noise_level = NULL, eta = 1.01,
                                   calculate_errors = FALSE,
                                   n_montecarlo = 100L,
                                   save_result = FALSE, random_state = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_hybrid_gmres,
                                         max_iterations = max_iterations,
                                         regularization_method = regularization_method,
                                         regularization = regularization,
                                         noise_level = noise_level,
                                         eta = eta),
        solve_kwargs = list(),
        method_name = "HybridGMRES",
        extra_output = list(regularization_method = regularization_method,
                            regularization = regularization),
        calculate_errors = calculate_errors,
        noise_level = if (is.null(noise_level)) 0.01 else noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
