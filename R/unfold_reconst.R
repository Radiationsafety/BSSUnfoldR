#' RECONST unfolding (Turchin/Vapnik statistical regularization)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_reconst.py}. Implements
#' the STREG1 algorithm from RECONST.FOR — solves
#' \eqn{(B \beta + \Omega \alpha) f = A_{vec} \beta} where \eqn{\Omega} is a
#' 5-diagonal smoothing matrix and \eqn{B = (A/S)^T (A/S)},
#' \eqn{A_{vec} = A^T (F/S^2)}, with \eqn{S} the measurement uncertainties
#' (\eqn{\sqrt{b}} when not supplied) normalised to geometric mean one.
#'
#' The regularisation weight \eqn{\alpha} and the data-fidelity weight
#' \eqn{\beta} are selected automatically (Vapnik bounds of the original
#' Fortran code) when \code{alpha < 0} and/or \code{beta <= 0}; see
#' \code{\link{solve_reconst}}.
#'
#' @note The automatic selection brackets \code{omega} and \code{delta}, each a
#'   difference of terms of order \eqn{10^{10}} against a result near zero, on a
#'   system matrix conditioned to \eqn{10^{11}}--\eqn{10^{13}}.  The branch the
#'   bisection takes is therefore decided below the precision of the evaluation,
#'   so it is not reproducible between BLAS implementations.  Both branches fit
#'   the readings equally well; pass \code{alpha} and \code{beta} explicitly for
#'   a reproducible spectrum.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused; accepted for API compatibility.
#' @param E_MeV Unused; accepted for API compatibility.
#' @param pp Numeric; smoothing parameter added to the diagonal of
#'   \eqn{\Omega}. Default 1e-3 (Python's value).
#' @param alpha Numeric regularization parameter for the smoothing matrix
#'   \eqn{\Omega}. \code{> 0} fixes it, \code{< 0} selects \code{|alpha|}
#'   automatically. Default -1 (auto), as in \code{bssunfold}.
#' @param beta Numeric data-fidelity weight. \code{> 0} fixes it,
#'   \code{<= 0} selects it automatically. Default 0 (auto).
#' @param sigma_b Optional measurement uncertainties (length m). Default
#'   \code{NULL} = \code{sqrt(b)}.
#' @return A list \code{list(spectrum, iterations = 1L, converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_reconst(A, b, NULL, alpha = 1.0)
solve_reconst <- function(A, b, x0 = NULL, E_MeV = NULL, pp = 1e-3,
                          alpha = -1.0, beta = 0.0, sigma_b = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A)
    n <- ncol(A)
    F <- b
    if (!is.null(sigma_b)) {
        S <- pmax(as.numeric(sigma_b), 1e-300)
    } else {
        S <- sqrt(pmax(F, 1e-10))
    }
    res <- .reconst_streg1(A, F, S, n, m, as.numeric(alpha), as.numeric(beta),
                           as.numeric(pp))
    list(spectrum = pmax(as.numeric(res$FI), 0.0), iterations = 1L,
         converged = TRUE)
}

## ainf[] of RECONST.FOR: relative stop criteria for the alpha search, the
## beta search and the alpha/beta alternation.
.reconst_ainf <- c(1.01, 1.01, 0.01, 0.0)

## Build the 5-diagonal smoothing matrix in band form: OMO rows 1..3 hold the
## second sub-diagonal, the first sub-diagonal and the diagonal (Python
## _build_omo_matrix, whose arrays are 0-based).
.reconst_build_omo <- function(n, pp) {
    XX <- seq(1.0, n + 1.0)          # XX[k] (0-based k) = k + 1
    AA <- numeric(n + 2L)
    BB <- numeric(n + 3L)
    CC <- numeric(n + 3L)
    ## Python: for i in range(2, n)  ->  1-based indices 3..n
    if (n >= 3L) {
        for (i in 3:(n + 1L - 1L)) {
            AA[i] <- 1.0 / (XX[i] - XX[i - 1L])
            CC[i] <- 1.0 / (XX[i - 1L] - XX[i - 2L])
            BB[i] <- -(AA[i] + CC[i])
        }
    }
    OMO <- matrix(0.0, nrow = 3L, ncol = n)
    for (i in seq_len(n)) {
        OMO[1L, i] <- AA[i] * CC[i]
        OMO[2L, i] <- AA[i] * BB[i] + BB[i + 1L] * CC[i + 1L]
        OMO[3L, i] <- AA[i]^2 + BB[i + 1L]^2 + CC[i + 2L]^2 +
                      pp * (XX[i + 1L] - XX[i])
    }
    OMO
}

## Band (3, n) -> full symmetric n x n matrix (Python _omo_to_full).
.reconst_omo_to_full <- function(OMO, n) {
    Omega <- matrix(0.0, n, n)
    for (i in seq_len(n)) {
        Omega[i, i] <- OMO[3L, i]
        if (i > 1L) Omega[i, i - 1L] <- OMO[2L, i]
        if (i > 2L) Omega[i, i - 2L] <- OMO[1L, i]
        if (i < n) Omega[i, i + 1L] <- OMO[2L, i + 1L]
        if (i < n - 1L) Omega[i, i + 2L] <- OMO[1L, i + 2L]
    }
    Omega
}

## Invert with the fallbacks of Python _invert_system: a Tikhonov shift when
## the condition number exceeds 1e12, escalating shifts when the factorisation
## fails, pseudo-inverse as the last resort.
.reconst_invert <- function(D) {
    n <- nrow(D)
    sv <- suppressWarnings(svd(D, nu = 0L, nv = 0L)$d)
    cond <- if (length(sv) == 0L || sv[n] == 0) Inf else sv[1L] / sv[n]
    ## numpy's `cond > 1e12` also catches an infinite condition number
    ## (exactly singular system matrix), so no is.finite() guard here.
    if (isTRUE(cond > 1e12)) {
        tr <- sum(diag(D))
        reg <- if (tr > 0) 1e-6 * tr / n else 1e-6
        return(solve(D + diag(n) * reg))
    }
    inv <- tryCatch(solve(D), error = function(e) NULL)
    if (!is.null(inv)) return(inv)
    for (reg in c(1e-6, 1e-4, 1e-2)) {
        inv <- tryCatch(solve(D + diag(n) * reg), error = function(e) NULL)
        if (!is.null(inv) && all(is.finite(inv))) return(inv)
    }
    ## numpy.linalg.pinv default rcond = 1e-15 relative to the largest
    ## singular value.
    sv <- svd(D)
    s_inv <- ifelse(sv$d > 1e-15 * max(sv$d), 1 / sv$d, 0)
    sv$v %*% (s_inv * t(sv$u))
}

## Python _reg1: build D = beta * B + alpha * Omega, invert it and form the
## regularised solution FI = D^-1 A_vec beta (ich > 0).
.reconst_reg1 <- function(B, OMO, A_vec, n, alpha, beta, ich = 2L, cache) {
    Omega <- .reconst_omo_to_full(OMO, n)
    D <- beta * B + alpha * Omega
    if (ich <= 0L) return(list(D = D, FI = NULL, SIGMA = NULL))
    D_inv <- .reconst_invert(D)
    FI <- as.numeric(D_inv %*% A_vec) * beta
    SIGMA <- sqrt(abs(diag(D_inv)))
    cache$Omega <- Omega
    cache$D_inv <- D_inv
    list(D = D_inv, FI = FI, SIGMA = SIGMA)
}

## omega(alpha) functional (Python _compute_omega), with the boundary handling
## of RECONST: rows 0, 1 only upper triangle + diagonal, row 2 one
## sub-diagonal, rows n-2 / n-1 trimmed.
.reconst_compute_omega <- function(Omega, D_inv, FI, n, alpha, cache) {
    code_trace <- 0.0
    for (i in seq_len(n)) {
        i0 <- i - 1L
        j_start <- if (i0 >= 3L) 0L else if (i0 == 2L) 1L else i0
        j_end <- if (i0 <= n - 3L) i0 + 2L else if (i0 == n - 2L) i0 + 1L else i0
        if (j_end < j_start) next
        for (j0 in j_start:j_end) {
            if (abs(i0 - j0) <= 2L) {
                code_trace <- code_trace +
                    Omega[i0 + 1L, j0 + 1L] * D_inv[j0 + 1L, i0 + 1L]
            }
        }
    }
    fof <- sum(as.numeric(crossprod(FI, Omega)) * FI)
    n / alpha - (code_trace + fof)
}

## delta(beta) discrepancy functional (Python _compute_delta).
.reconst_compute_delta <- function(B, D_inv, FI, A_vec, F, S, n, m, beta) {
    d1 <- sum(diag(B %*% D_inv))
    d2 <- sum(as.numeric(crossprod(FI, B)) * FI)
    d3 <- sum(A_vec * FI)
    d4 <- sum((F / S)^2)
    m / beta - (d1 + d2 - 2 * d3 + d4)
}

## Find alpha with omega(alpha) = 0 (Python _def_alpha).
.reconst_def_alpha <- function(B, OMO, A_vec, F, S, n, m, alpha, beta,
                               omega_init, cache) {
    alm <- if (omega_init >= 0) 4.0 else 1 / 4.0
    als <- omega_init
    for (rep in seq_len(50L)) {
        alpha <- alpha * alm
        r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
        omega <- .reconst_compute_omega(cache$Omega, r$D, r$FI, n, alpha, cache)
        if (omega * als <= 0) break
    }
    aln <- (alpha + alpha / alm) / 5.0
    alk <- 4.0 * aln
    for (rep in seq_len(100L)) {
        alpha <- (aln + alk) / 2.0
        r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
        omega <- .reconst_compute_omega(cache$Omega, r$D, r$FI, n, alpha, cache)
        if (omega < 0) alk <- alpha else aln <- alpha
        if (alk <= aln * .reconst_ainf[1L]) break
    }
    alpha
}

## Find beta with delta(beta) = 0 (Python _def_beta).
.reconst_def_beta <- function(B, OMO, A_vec, F, S, n, m, alpha, beta,
                              delta_init, cache) {
    betm <- if (delta_init >= 0) 4.0 else 1 / 4.0
    bets <- delta_init
    for (rep in seq_len(50L)) {
        beta <- beta * betm
        r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
        delta <- .reconst_compute_delta(B, r$D, r$FI, A_vec, F, S, n, m, beta)
        if (delta * bets <= 0) break
    }
    betn <- (beta + beta / betm) / 5.0
    betk <- 4.0 * betn
    for (rep in seq_len(100L)) {
        beta <- (betn + betk) / 2.0
        r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
        delta <- .reconst_compute_delta(B, r$D, r$FI, A_vec, F, S, n, m, beta)
        if (delta < 0) betk <- beta else betn <- beta
        if (betk <= betn * .reconst_ainf[2L]) break
    }
    beta
}

## Core STREG1 (Python _streg1). Returns list(FI, SIGMA).
.reconst_streg1 <- function(AK, F, S, n, m, alpha, beta, pp) {
    sa <- exp(mean(log(pmax(S, 1e-300))))
    S_norm <- S / sa
    ## B = (A / S_norm)^T (A / S_norm); Python divides every ROW i of A by
    ## S_norm[i], so the division must run along the row index (sweep margin
    ## 1).  Recycling a vector over a column-major matrix would only be
    ## correct when ncol(A) is a multiple of the number of detectors.
    W <- sweep(AK, 1L, S_norm, "/")
    B <- t(W) %*% W
    A_vec <- as.numeric(t(AK) %*% (F / S_norm^2))
    OMO <- .reconst_build_omo(n, pp)
    cache <- new.env(parent = emptyenv())

    if (beta > 0.0) {
        beta <- beta / sa^2
        if (alpha >= 0.0) {
            r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
            D_inv <- r$D; FI <- r$FI; SIGMA <- r$SIGMA
        } else {
            alpha <- -alpha
            r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
            omega_init <- .reconst_compute_omega(cache$Omega, r$D, r$FI, n,
                                                 alpha, cache)
            alpha <- .reconst_def_alpha(B, OMO, A_vec, F, S_norm, n, m, alpha,
                                       beta, omega_init, cache)
            r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
            D_inv <- r$D; FI <- r$FI; SIGMA <- r$SIGMA
        }
    } else {
        beta <- 1.0 / sa^2
        if (alpha >= 0.0) {
            r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
            delta_init <- .reconst_compute_delta(B, r$D, r$FI, A_vec, F,
                                                 S_norm, n, m, beta)
            beta <- .reconst_def_beta(B, OMO, A_vec, F, S_norm, n, m, alpha,
                                     beta, delta_init, cache)
            r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
            D_inv <- r$D; FI <- r$FI; SIGMA <- r$SIGMA
        } else {
            alpha <- -alpha
            bet_saved <- beta
            r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
            D_inv <- r$D; FI <- r$FI
            for (rep in seq_len(30L)) {
                omega_init <- .reconst_compute_omega(cache$Omega, D_inv, FI, n,
                                                     alpha, cache)
                alpha <- .reconst_def_alpha(B, OMO, A_vec, F, S_norm, n, m,
                                           alpha, beta, omega_init, cache)
                r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
                D_inv <- r$D; FI <- r$FI
                delta_init <- .reconst_compute_delta(B, D_inv, FI, A_vec, F,
                                                     S_norm, n, m, beta)
                beta <- .reconst_def_beta(B, OMO, A_vec, F, S_norm, n, m,
                                         alpha, beta, delta_init, cache)
                if (abs(bet_saved - beta) <= beta * .reconst_ainf[3L]) break
                r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
                D_inv <- r$D; FI <- r$FI
                bet_saved <- beta
            }
            cors <- 1.0 / (sqrt(beta) * sa)
            sa <- sa * cors
            S_norm <- S_norm * cors
            r <- .reconst_reg1(B, OMO, A_vec, n, alpha, beta, 2L, cache)
            D_inv <- r$D; FI <- r$FI; SIGMA <- r$SIGMA
        }
    }
    FI <- pmax(FI, 0)
    list(FI = FI, SIGMA = SIGMA, alpha = alpha, beta = beta)
}

#' Wrapper around \code{\link{solve_reconst}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_reconst
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_reconst <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              alpha = -1.0, pp = 1e-3, beta = 0.0,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_reconst,
                                        E_MeV = E_MeV, pp = pp,
                                        alpha = alpha, beta = beta),
        solve_kwargs = list(),
        method_name = "Reconst",
        extra_output = list(alpha = alpha, pp = pp, beta = beta),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
