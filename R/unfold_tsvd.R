#' TSVD (Truncated SVD) unfolding
#'
#' Solves \eqn{A x = b} via truncated SVD: \eqn{x = V_k S_k^{-1} U_k^T b}.
#' The truncation parameter k can be set explicitly, by threshold ratio, or
#' selected automatically by one of several criteria (discrepancy principle,
#' L-curve, GCV, etc.).
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Ignored; present for API compatibility.
#' @param method Character; k-selection method: 'discrepancy', 'l_curve',
#'   'gcv', 'energy', 'threshold_ratio', 'median_threshold', 'donoho',
#'   'mean_threshold'. Default 'discrepancy'.
#' @param k Optional integer; fixed truncation parameter. Overrides
#'   \code{method}. Default \code{NULL}.
#' @param threshold Optional numeric; threshold ratio for truncation.
#'   Default \code{NULL}.
#' @param noise_level Optional numeric noise level for discrepancy principle.
#'   Default \code{NULL}.
#' @return A list \code{list(spectrum, iterations = k, converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_tsvd(A, b, NULL, method = "gcv")
solve_tsvd <- function(A, b, x0 = NULL, method = "discrepancy",
                       k = NULL, threshold = NULL, noise_level = NULL) {
    v <- validate_system(A, b, x0 = x0)
    A <- v$A; b <- v$b
    m <- nrow(A); n <- ncol(A)
    sv <- svd(A, nu = min(m, n), nv = min(m, n))
    U <- sv$u; s <- sv$d; V <- sv$v
    Vt <- t(V)

    max_k <- min(m, n)
    if (!is.null(k)) {
        k <- min(as.integer(k), length(s))
    } else if (!is.null(threshold)) {
        k <- sum(s / s[1L] > threshold)
    } else {
        k <- .tsvd_select_k(s, A, b, method = method,
                            noise_level = noise_level)
    }
    k <- max(1L, min(k, m, n))
    s_k <- s[seq_len(k)]
    U_k <- U[, seq_len(k), drop = FALSE]
    V_k <- V[, seq_len(k), drop = FALSE]
    x <- as.numeric(V_k %*% ((1 / s_k) * (t(U_k) %*% b)))
    list(spectrum = pmax(x, 0), iterations = as.integer(k),
         converged = TRUE)
}

.tsvd_select_k <- function(s, A, b, method = "discrepancy",
                            noise_level = NULL) {
    m <- nrow(A); n <- ncol(A)
    max_k <- min(m, n)
    sv <- svd(A, nu = min(m, n), nv = min(m, n))
    U <- sv$u; s_full <- sv$d; V <- sv$v
    if (method == "discrepancy") {
        if (is.null(noise_level)) noise_level <- s_full[1L] * 1e-3
        for (i in seq_len(max_k)) {
            s_i <- s_full[seq_len(i)]
            U_i <- U[, seq_len(i), drop = FALSE]
            V_i <- V[, seq_len(i), drop = FALSE]
            x_i <- as.numeric(V_i %*% ((1 / s_i) * (t(U_i) %*% b)))
            res <- sqrt(sum((as.numeric(A %*% x_i) - b)^2))
            if (res <= noise_level * sqrt(max(m - i, 1L))) return(i)
        }
        return(max_k)
    }
    if (method == "energy") {
        cum_e <- cumsum(s_full^2) / sum(s_full^2)
        k <- which(cum_e >= 0.95)[1L]
        return(k)
    }
    if (method == "l_curve") {
        res_n <- numeric(max_k); sol_n <- numeric(max_k)
        for (i in seq_len(max_k)) {
            s_i <- s_full[seq_len(i)]
            U_i <- U[, seq_len(i), drop = FALSE]
            V_i <- V[, seq_len(i), drop = FALSE]
            x_i <- as.numeric(V_i %*% ((1 / s_i) * (t(U_i) %*% b)))
            res_n[i] <- sqrt(sum((as.numeric(A %*% x_i) - b)^2))
            sol_n[i] <- sqrt(sum(x_i^2))
        }
        lr <- log(pmax(res_n, 1e-300))
        ls <- log(pmax(sol_n, 1e-300))
        if (length(lr) >= 3L) {
            curv <- numeric(length(lr) - 2L)
            for (i in seq_len(length(curv))) {
                j <- i + 1L  # index into lr/ls, ranges 2..(n-1)
                dx1 <- lr[j] - lr[j - 1L]
                dy1 <- ls[j] - ls[j - 1L]
                dx2 <- lr[j + 1L] - lr[j]
                dy2 <- ls[j + 1L] - ls[j]
                curv[i] <- abs(dx1 * dy2 - dx2 * dy1) /
                           ((dx1^2 + dy1^2)^1.5 + 1e-10)
            }
            k_idx <- which.max(curv)
            return(min(k_idx + 1L, length(s_full)))
        }
        return(length(s_full) %/% 2L)
    }
    if (method == "gcv") {
        beta <- as.numeric(t(U) %*% b)
        ks <- seq_len(max_k)
        gcv <- sapply(ks, function(i) {
            resid <- sum(beta[(i + 1L):length(beta)]^2)
            eff_p <- m - i
            if (eff_p <= 0) Inf else resid / (eff_p^2)
        })
        return(ks[which.min(gcv)])
    }
    if (method == "threshold_ratio") {
        return(sum(s_full / s_full[1L] > 1e-2))
    }
    if (method == "median_threshold") {
        return(sum(s_full >= median(s_full)))
    }
    if (method == "donoho") {
        sigma_donoho <- 0.05
        donoho_rcond <- 4 / sqrt(3) * sqrt(n) * sigma_donoho
        return(sum(s_full > donoho_rcond))
    }
    # mean_threshold
    sum(s_full >= mean(s_full))
}

#' Wrapper around \code{\link{solve_tsvd}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_tsvd
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_tsvd <- function(detector_names, n_energy_bins, E_MeV,
                         sensitivities, cc_icrp116, save_result_callback,
                         readings, initial_spectrum = NULL,
                         method = "discrepancy", k = NULL,
                         threshold = NULL, noise_level = NULL,
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
        solve_func = make_solve_wrapper(solve_tsvd,
                                         method = method,
                                         k = k,
                                         threshold = threshold,
                                         noise_level = noise_level),
        solve_kwargs = list(),
        method_name = "TSVD",
        extra_output = list(k = k, k_method = method),
        calculate_errors = calculate_errors,
        noise_level = if (is.null(noise_level)) 0.01 else noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
