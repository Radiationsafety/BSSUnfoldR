#' Generalized Krylov Subspace (GKS) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_gks.py}. Builds a
#' Krylov subspace by Golub-Kahan bidiagonalization of \code{A} and projects
#' both \code{A} and the regularization operator \code{L} onto that subspace.
#' At each iteration a small projected Tikhonov problem is solved, with the
#' regularization parameter selected automatically by Generalized Cross
#' Validation (GCV), the Discrepancy Principle (DP) or the L-curve.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused; accepted for API compatibility.
#' @param smoothness_order Integer; 0 (identity), 1, or 2. Default 0.
#' @param regularization_method Character; \code{"gcv"} (default), \code{"dp"},
#'   \code{"lcurve"}, or \code{"manual"}.
#' @param max_iterations Optional integer; Krylov dimension. Default \code{NULL}
#'   = \code{min(nrow(A), ncol(A))}.
#' @param regularization Numeric fallback regularization parameter. Default 1e-8.
#' @param noise_level Optional numeric; relative noise level for DP. Default
#'   \code{NULL}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_gks(A, b, NULL, max_iterations = 3L)
solve_gks <- function(A, b, x0 = NULL, smoothness_order = 0L,
                       regularization_method = "gcv",
                       max_iterations = NULL,
                       regularization = 1e-8, noise_level = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A); n <- ncol(A)
    if (is.null(max_iterations)) max_iterations <- min(m, n)
    max_iterations <- max(1L, as.integer(max_iterations))

    L <- make_regularization_operator(n, smoothness_order,
                                       identity_for_zero = TRUE)
    beta <- sqrt(sum(b^2))
    if (beta == 0.0) {
        return(list(spectrum = rep(0.0, n), iterations = 0L,
                    converged = TRUE))
    }

    U <- matrix(0.0, nrow = m, ncol = 1L)
    U[, 1L] <- b / beta
    V <- matrix(0.0, nrow = n, ncol = 0L)
    alphas <- numeric(0L)
    betas <- numeric(0L)

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
        betas <- c(betas, new_beta)

        p <- ncol(U)
        B <- matrix(0.0, nrow = p, ncol = k)
        for (i in seq_len(min(p, k))) B[i, i] <- alphas[i]
        if (k > 1L && p > 1L) {
            sub <- min(k - 1L, p - 1L)
            if (sub >= 1L) {
                for (i in 1:sub) B[i + 1L, i] <- betas[i]
            }
        }
        # bhat has length p (= ncol(U)). When bidiag did not terminate,
        # p == k + 1; when it terminated, p == k, so we zero-pad to p.
        bhat <- c(beta, rep(0.0, p - 1L))

        RA <- as.numeric(t(U) %*% (A %*% V))
        dim(RA) <- c(p, k)
        RL <- as.matrix(L %*% V)
        dim(RL) <- c(nrow(L), k)

        lam <- switch(regularization_method,
            "gcv"     = .gks_projected_gcv(RA, RL, bhat),
            "dp"      = {
                nl <- if (is.null(noise_level)) 0.01 else noise_level
                .gks_projected_dp(RA, bhat, nl)
            },
            "lcurve"  = .gks_projected_lcurve(RA, bhat),
            "manual"  = regularization,
            stop("Unsupported regularization method: ", regularization_method)
        )
        if (!is.finite(lam) || lam <= 0) lam <- regularization

        lhs <- rbind(RA, sqrt(lam) * RL)
        rhs <- c(bhat, rep(0.0, nrow(RL)))
        y <- qr.solve(lhs, rhs)
        best_x <- as.numeric(V %*% y)
        iterations <- k
        if (converged) break
    }
    list(spectrum = pmax(best_x, 0), iterations = iterations,
         converged = converged)
}

.gks_projected_gcv <- function(RA, RL, bhat,
                                n_lambdas = 200L,
                                lambda_range = c(1e-12, 1e2)) {
    # Use reduced SVD: nu must not exceed min(dim(RA)) — R pads with zero
    # columns otherwise, breaking subsequent matrix-vector products.
    k_min <- min(dim(RA))
    sv <- svd(RA, nu = k_min, nv = 0)
    U <- sv$u
    s <- sv$d
    c <- as.numeric(t(U) %*% bhat)
    s2 <- s^2
    lambdas <- 10^seq(log10(lambda_range[1L]), log10(lambda_range[2L]),
                       length.out = n_lambdas)
    m_proj <- nrow(RA)
    gcv_values <- numeric(n_lambdas)
    for (i in seq_len(n_lambdas)) {
        lam <- lambdas[i]
        filt <- s2 / (s2 + lam)
        residual_coeff <- lam / (s2 + lam)
        residual_sq <- sum((residual_coeff * c)^2)
        trace_term <- sum(filt)
        denom <- (m_proj - trace_term)^2
        gcv_values[i] <- if (denom > 0) residual_sq / denom else Inf
    }
    lambdas[which.min(gcv_values)]
}

.gks_projected_dp <- function(RA, bhat, noise_level,
                               n_lambdas = 200L,
                               lambda_range = c(1e-12, 1e2)) {
    k_min <- min(dim(RA))
    sv <- svd(RA, nu = k_min, nv = k_min)
    U <- sv$u
    s <- sv$d
    Vt <- sv$v   # sv$v has dim (ncol(RA), k_min); use t(Vt) to get (k_min, ncol(RA))
    c <- as.numeric(t(U) %*% bhat)
    s2 <- s^2
    m_proj <- nrow(RA)
    target <- noise_level * sqrt(m_proj)
    lambdas <- 10^seq(log10(lambda_range[1L]), log10(lambda_range[2L]),
                       length.out = n_lambdas)
    residuals <- numeric(n_lambdas)
    for (i in seq_len(n_lambdas)) {
        lam <- lambdas[i]
        filt <- s / (s2 + lam)
        x <- as.numeric(t(Vt) %*% (filt * c))
        residuals[i] <- sqrt(sum((as.numeric(RA %*% x) - bhat)^2))
    }
    lambdas[which.min(abs(residuals - target))]
}

.gks_projected_lcurve <- function(RA, bhat,
                                    n_lambdas = 200L,
                                    lambda_range = c(1e-12, 1e2)) {
    k_min <- min(dim(RA))
    sv <- svd(RA, nu = k_min, nv = k_min)
    U <- sv$u
    s <- sv$d
    Vt <- sv$v
    c <- as.numeric(t(U) %*% bhat)
    s2 <- s^2
    lambdas <- 10^seq(log10(lambda_range[1L]), log10(lambda_range[2L]),
                       length.out = n_lambdas)
    residuals <- numeric(n_lambdas)
    norms <- numeric(n_lambdas)
    for (i in seq_len(n_lambdas)) {
        lam <- lambdas[i]
        filt <- s / (s2 + lam)
        x <- as.numeric(t(Vt) %*% (filt * c))
        residuals[i] <- sqrt(sum((as.numeric(RA %*% x) - bhat)^2))
        norms[i] <- sqrt(sum(x^2))
    }
    valid <- residuals > 0 & norms > 0
    if (sum(valid) < 3L) return(lambdas[floor(length(lambdas) / 2)])
    log_res <- log(pmax(residuals[valid], 1e-300))
    log_norm <- log(pmax(norms[valid], 1e-300))
    p1 <- c(log_res[1L], log_norm[1L])
    p2 <- c(log_res[length(log_res)], log_norm[length(log_norm)])
    d21 <- p2 - p1
    denom <- sqrt(sum(d21^2))
    if (denom < 1e-300) return(lambdas[floor(length(lambdas) / 2)])
    distances <- abs(d21[1L] * (p1[2L] - log_norm) -
                       d21[2L] * (p1[1L] - log_res)) / denom
    lambdas[valid][which.max(distances)]
}

#' Wrapper around \code{\link{solve_gks}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_gks
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_gks <- function(detector_names, n_energy_bins, E_MeV,
                        sensitivities, cc_icrp116, save_result_callback,
                        readings, initial_spectrum = NULL,
                        smoothness_order = 0L,
                        regularization_method = "gcv",
                        max_iterations = NULL,
                        regularization = 1e-8, noise_level = NULL,
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
        solve_func = make_solve_wrapper(solve_gks,
                                         smoothness_order = smoothness_order,
                                         regularization_method = regularization_method,
                                         max_iterations = max_iterations,
                                         regularization = regularization,
                                         noise_level = noise_level),
        solve_kwargs = list(),
        method_name = "GKS",
        extra_output = list(
            smoothness_order = as.integer(smoothness_order),
            regularization_method = regularization_method,
            regularization = regularization
        ),
        calculate_errors = calculate_errors,
        noise_level = if (is.null(noise_level)) 0.01 else noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
