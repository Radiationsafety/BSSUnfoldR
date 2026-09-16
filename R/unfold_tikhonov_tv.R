#' Tikhonov-TV (Total Variation) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_tikhonov_tv.py}.
#' Noise-constrained ADMM solver for:
#' \deqn{\min f(m) \text{ subject to } \|A m - b\|^2 = \epsilon}
#' where the regularizer \eqn{f} is one of:
#' \describe{
#'   \item{TT}{\eqn{\|D_1 m\|_1 + \beta/2 \|D_1^{bar} g_2\|_2^2} (TV + Tikhonov, default)}
#'   \item{TV}{\eqn{\|D_1 m\|_1} (pure total variation)}
#'   \item{T}{\eqn{\beta/2 \|D_1^{bar} g_2\|_2^2} (pure Tikhonov)}
#' }
#' Here \eqn{D_1} is the first-derivative (TV) operator and \eqn{D_1^{bar}}
#' the second-derivative operator. The noise constraint is enforced through
#' an augmented Lagrangian / ADMM scheme.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused; accepted for API compatibility.
#' @param epsilon Optional numeric; squared 2-norm of the noise. If \code{NULL},
#'   derived from \code{noise_level} or from unregularized least-squares
#'   residuals.
#' @param mu Numeric vector of length 3; penalty parameters. Default
#'   \code{c(1.0, 1.0, 1.0)}.
#' @param max_iterations Positive integer; ADMM iteration count. Default 100.
#' @param type_ Character; \code{"TT"}, \code{"TV"}, or \code{"T"}. Default \code{"TT"}.
#' @param beta Numeric balancing parameter between TV and Tikhonov terms.
#'   Default 1.0.
#' @param zthr Numeric z-score threshold for adaptive beta (not used when
#'   beta is numeric). Default 2.5.
#' @param tolerance Numeric; stabilization stopping tolerance. Default 1e-4.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_tikhonov_tv(A, b, NULL, max_iterations = 50)
solve_tikhonov_tv <- function(A, b, x0 = NULL, epsilon = NULL,
                               mu = c(1.0, 1.0, 1.0), max_iterations = 100L,
                               type_ = "TT", beta = 1.0,
                               zthr = 2.5, tolerance = 1e-4) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m_len <- nrow(A); n <- ncol(A)
    if (!(type_ %in% c("TT", "TV", "T"))) {
        stop("Unsupported type_: ", type_, ". Choose from 'TT', 'TV', 'T'.")
    }
    adapt_beta <- is.character(beta) && beta == "adapt"
    mu1 <- mu[1L]; mu2 <- mu[2L]; mu3 <- mu[3L]
    D1 <- as.matrix(create_derivative_matrix(n, 1L))         # (n-1, n)
    if (n - 1L >= 3L) {
        D1_bar <- as.matrix(create_derivative_matrix(n - 1L, 2L))  # (n-3, n-1)
    } else {
        D1_bar <- matrix(0.0, nrow = 0L, ncol = n - 1L)
    }
    if (is.null(epsilon)) {
        x_ls <- qr.solve(A, b, tol = 1e-10)
        epsilon <- sum((b - as.numeric(A %*% x_ls))^2)
    }
    epsilon <- as.numeric(epsilon)
    if (epsilon <= 0) epsilon <- 1e-12

    B <- mu1 * crossprod(D1) + mu2 * crossprod(A)
    B_inv <- tryCatch(solve(B), error = function(e) {
        # Pseudo-inverse fallback
        sv <- svd(B, nu = nrow(B), nv = ncol(B))
        s_inv <- ifelse(sv$d > 1e-10 * max(sv$d), 1 / sv$d, 0)
        sv$u %*% (s_inv * (t(sv$u) %*% sv$v))
    })

    m_var <- rep(0.0, n)
    g1 <- rep(0.0, n - 1L)
    g2 <- rep(0.0, n - 1L)
    e <- rep(0.0, m_len)
    lambda_1 <- rep(0.0, n - 1L)
    lambda_2 <- rep(0.0, m_len)
    lambda_3 <- 0.0
    beta_k <- if (!adapt_beta) as.numeric(beta) else 1.0
    if (type_ == "TV") beta_k <- 0.0
    converged <- FALSE
    m_prev <- NULL
    stopit <- max_iterations

    for (k in seq_len(max_iterations)) {
        rhs <- mu1 * as.numeric(t(D1) %*% (g1 + g2 + lambda_1)) +
               mu2 * as.numeric(t(A) %*% (b + e + lambda_2))
        m_var <- as.numeric(B_inv %*% rhs)
        if (k == 1L) {
            m_prev <- m_var
        } else {
            norm_prev <- sqrt(sum(m_prev^2))
            diffm <- if (norm_prev > 0) {
                sqrt(sum((m_var - m_prev)^2)) / norm_prev
            } else {
                sqrt(sum(m_var^2))
            }
            m_prev <- m_var
            if (diffm < tolerance && !converged) {
                stopit <- k; converged <- TRUE; break
            }
        }
        # g1-subproblem
        if (type_ %in% c("TT", "TV")) {
            y1 <- as.numeric(D1 %*% m_var) - g2 - lambda_1
            g1 <- sign(y1) * pmax(abs(y1) - 1.0 / mu1, 0.0)
        }
        # g2-subproblem
        if (type_ %in% c("TT", "T")) {
            y2 <- as.numeric(D1 %*% m_var) - g1 - lambda_1
            if (beta_k > 0 && mu1 > 0 && nrow(D1_bar) > 0L) {
                lhs <- diag(n - 1L) + (beta_k / mu1) * crossprod(D1_bar)
                g2 <- qr.solve(lhs, y2)
            } else {
                g2 <- y2
            }
        }
        # e-subproblem
        y <- as.numeric(A %*% m_var) - b - lambda_2
        E <- sum(y * y)
        if (E > 0) {
            pp <- (mu2 - 2 * mu3 * (epsilon + lambda_3)) / (2 * mu3 * E)
            qq <- -mu2 / (2 * mu3 * E)
            gamma <- .ttv_gamma_from_cubic(pp, qq)
        } else {
            gamma <- 0.0
        }
        e <- gamma * y
        # Dual updates
        lambda_1 <- lambda_1 + g1 + g2 - as.numeric(D1 %*% m_var)
        lambda_2 <- lambda_2 + b + e - as.numeric(A %*% m_var)
        lambda_3 <- lambda_3 + epsilon - sum(e * e)
        # beta-update (only adaptive + TT)
        if (adapt_beta && type_ == "TT") {
            g <- as.numeric(D1 %*% m_var)
            target <- .ttv_zscore_max(g, zthr)
            value <- if (length(g2) > 0L) max(abs(g2)) else 0.0
            denom <- value + target
            if (denom > 0) {
                beta_k <- 2.0 * value / denom * beta_k
            }
        } else if (type_ == "T" && !adapt_beta) {
            beta_k <- as.numeric(beta)
        }
    }
    spectrum <- pmax(m_var, 0)
    list(spectrum = as.numeric(spectrum), iterations = stopit,
         converged = converged)
}

.ttv_zscore_max <- function(p, a = 2.5) {
    p <- sort(abs(p))
    if (length(p) == 0L) return(0.0)
    p <- p[ceiling(length(p) / 2):length(p)]
    mad <- 1.4826 * (mean(abs(p - mean(p))) + .Machine$double.eps)
    med <- median(p)
    z <- (p - med) / mad
    idx <- which(abs(z) < a)
    if (length(idx) == 0L) return(0.0)
    max(p[idx])
}

.ttv_gamma_from_cubic <- function(pp, qq) {
    # Largest real root of gamma^3 + pp*gamma + qq = 0
    roots <- polyroot(c(qq, pp, 0, 1))
    real_roots <- Re(roots)[abs(Im(roots)) < 1e-8]
    if (length(real_roots) == 0L) return(0.0)
    max(real_roots)
}

#' Wrapper around \code{\link{solve_tikhonov_tv}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_tikhonov_tv
#' @param noise_level Optional numeric; relative noise level used to derive
#'   \code{epsilon} when the latter is not given.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_tikhonov_tv <- function(detector_names, n_energy_bins, E_MeV,
                                 sensitivities, cc_icrp116, save_result_callback,
                                 readings, initial_spectrum = NULL,
                                 epsilon = NULL,
                                 mu = c(1.0, 1.0, 1.0), max_iterations = 100L,
                                 type_ = "TT", beta = 1.0, zthr = 2.5,
                                 tolerance = 1e-4, noise_level = NULL,
                                 calculate_errors = FALSE,
                                 n_montecarlo = 100L,
                                 save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    if (is.null(epsilon) && !is.null(noise_level)) {
        selected <- detector_names[detector_names %in% names(readings)]
        b_norm <- sqrt(sum(as.numeric(readings[selected])^2))
        epsilon <- (noise_level * b_norm)^2
    }
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_tikhonov_tv,
                                         epsilon = epsilon,
                                         mu = mu,
                                         max_iterations = max_iterations,
                                         type_ = type_,
                                         beta = beta,
                                         zthr = zthr,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "TikhonovTV",
        extra_output = list(
            type_ = type_,
            epsilon = if (!is.null(epsilon)) as.numeric(epsilon) else NULL,
            beta = if (!is.character(beta)) as.numeric(beta) else beta
        ),
        calculate_errors = calculate_errors,
        noise_level = if (is.null(noise_level)) 0.01 else noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
