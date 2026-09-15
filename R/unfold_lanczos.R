#' Lanczos-hybrid (Krylov-GCV) unfolding
#'
#' R port of \code{bssunfold/core/unfold_lanczos.py}. Performs Golub-Kahan
#' bidiagonalization of \code{A}, building a sequence of Krylov subspaces. On
#' the projected problem \code{min ||B_k y - bhat||^2} a Tikhonov term
#' \code{lambda * ||y||^2} is added, where \code{lambda} is selected
#' automatically by Generalized Cross Validation (GCV) at each iteration.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Ignored; accepted for API compatibility.
#' @param max_iterations Positive integer; Krylov dimension. Default
#'   \code{NULL} = \code{min(nrow(A), ncol(A))}.
#' @param regularization Numeric; fallback regularization used if GCV returns
#'   a degenerate value. Default 1e-8.
#' @param noise_level Optional numeric; relative noise level used for
#'   discrepancy-principle early stopping. Default \code{NULL}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_lanczos(A, b, NULL, max_iterations = 3)
solve_lanczos <- function(A, b, x0 = NULL, max_iterations = NULL,
                           regularization = 1e-8, noise_level = NULL) {
    v <- validate_system(A, b)
    A <- v$A; b <- v$b
    m <- nrow(A); n <- ncol(A)
    if (is.null(max_iterations)) max_iterations <- min(m, n)
    max_iterations <- max(1L, as.integer(max_iterations))

    beta <- sqrt(sum(b^2))
    if (beta == 0) {
        return(list(spectrum = rep(0.0, n), iterations = 0L,
                    converged = TRUE))
    }

    U <- matrix(0.0, nrow = m, ncol = 1L)
    U[, 1L] <- b / beta
    V <- matrix(0.0, nrow = n, ncol = 0L)
    alphas <- numeric(0L)
    betas  <- numeric(0L)

    best_x <- rep(0.0, n)
    iterations <- 0L
    converged <- FALSE

    for (k in seq_len(max_iterations)) {
        u <- U[, ncol(U)]
        if (k == 1L) {
            vv <- as.numeric(t(A) %*% u)
        } else {
            vv <- as.numeric(t(A) %*% u) - betas[length(betas)] * V[, ncol(V)]
        }
        alpha <- sqrt(sum(vv^2))
        if (alpha <= 1e-14) { converged <- TRUE; break }
        vv <- vv / alpha
        V <- cbind(V, vv)

        u2 <- as.numeric(A %*% vv) - alpha * u
        new_beta <- sqrt(sum(u2^2))
        if (new_beta <= 1e-14) {
            converged <- TRUE
        } else {
            U <- cbind(U, u2 / new_beta)
        }

        alphas <- c(alphas, alpha)
        betas  <- c(betas, new_beta)

        B <- matrix(0.0, nrow = k + 1L, ncol = k)
        for (i in seq_len(k)) B[i, i] <- alphas[i]
        if (k > 1L) {
            for (i in 2:k) B[i, i - 1L] <- betas[i - 1L]
        }
        bhat <- c(beta, rep(0.0, k))

        lam <- .lanczos_projected_gcv(B, bhat, m)
        if (!is.finite(lam) || lam <= 0) lam <- regularization

        # Reduced SVD of B (k+1 x k): U is (k+1, k), s is (k,), V is (k, k).
        svdB <- svd(B, nu = min(dim(B)), nv = min(dim(B)))
        Ub <- svdB$u
        s <- svdB$d
        Vb <- svdB$v
        c <- as.numeric(t(Ub) %*% bhat)
        # c has length min(k+1, k) = k; s also length k.
        s2 <- s^2
        y <- as.numeric(t(Vb) %*% (s * c / (s2 + lam)))
        best_x <- as.numeric(V %*% y)

        iterations <- k
        if (!is.null(noise_level)) {
            residual <- sqrt(sum((as.numeric(A %*% best_x) - b)^2))
            if (residual <= noise_level * sqrt(m)) {
                converged <- TRUE; break
            }
        }
        if (converged) break
    }
    list(spectrum = best_x, iterations = iterations, converged = converged)
}

.lanczos_projected_gcv <- function(B, bhat, m,
                                    n_lambdas = 200L,
                                    lambda_range = c(1e-12, 1e2)) {
    svdB <- svd(B, nu = min(dim(B)), nv = 0)
    Ub <- svdB$u
    s <- svdB$d
    c <- as.numeric(t(Ub) %*% bhat)
    orth_res <- sum(bhat^2) - sum(c^2)
    s2 <- s^2
    lambdas <- 10^seq(log10(lambda_range[1L]),
                       log10(lambda_range[2L]),
                       length.out = n_lambdas)
    gcv_values <- numeric(n_lambdas)
    for (i in seq_len(n_lambdas)) {
        lam <- lambdas[i]
        num <- sum((c * lam / (s2 + lam))^2) + orth_res
        den <- (m - sum(s2 / (s2 + lam)))^2
        gcv_values[i] <- if (den > 0) num / den else Inf
    }
    lambdas[which.min(gcv_values)]
}

#' Wrapper around \code{\link{solve_lanczos}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_lanczos
#' @param regularization_method Character; only \code{"gcv"} is supported.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_lanczos <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            regularization_method = "gcv",
                            max_iterations = NULL,
                            regularization = 1e-8, noise_level = NULL,
                            calculate_errors = FALSE,
                            n_montecarlo = 100L,
                            save_result = FALSE, random_state = NULL) {
    if (regularization_method != "gcv") {
        stop("Unsupported regularization method: ", regularization_method,
             ". The Lanczos hybrid method currently supports 'gcv'.")
    }
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_lanczos,
                                         max_iterations = max_iterations,
                                         regularization = regularization,
                                         noise_level = noise_level),
        solve_kwargs = list(),
        method_name = "Lanczos",
        extra_output = list(
            regularization_method = regularization_method,
            regularization = regularization
        ),
        calculate_errors = calculate_errors,
        noise_level = if (is.null(noise_level)) 0.01 else noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
