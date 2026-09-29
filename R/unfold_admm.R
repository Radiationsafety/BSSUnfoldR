#' ADMM unfolding (Alternating Direction Method of Multipliers)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_admm.py}.
#' Solves the constrained regularized unfolding problem
#' \deqn{\min_x \frac12 \|Ax - b\|^2 + \lambda_1 \|x\|_1 + \lambda_{tv} \|Dx\|_1
#' \quad \text{s.t. } x \ge 0}
#' by consensus ADMM: the x-update is an exact non-negative least squares on
#' the augmented system, the z-updates are element-wise soft-thresholding and
#' the u-updates are scaled dual ascent.  The penalty \code{rho} is adapted
#' every ten iterations following Boyd et al., section 3.4.1.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Initial spectrum (length n); projected onto the nonnegative orthant.
#' @param max_iterations Positive integer; default 500L.
#' @param tolerance Positive numeric relative change tolerance; default 1e-6.
#' @param l1_penalty Numeric L1 (sparsity) penalty weight; default 0.
#' @param tv_penalty Numeric total-variation penalty weight \eqn{\|Dx\|_1}; default 0.
#' @param rho Optional ADMM penalty parameter; \code{NULL} (default) initializes
#'   and adapts it automatically.
#' @param adaptive_rho Logical; adapt \code{rho} every 10 iterations; default \code{TRUE}.
#' @param abstol Numeric absolute residual tolerance; default 1e-10.
#' @param reltol Numeric relative residual tolerance; default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
solve_admm <- function(A, b, x0, max_iterations = 500L, tolerance = 1e-6,
                       l1_penalty = 0.0, tv_penalty = 0.0, rho = NULL,
                       adaptive_rho = TRUE, abstol = 1e-10, reltol = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A); n <- ncol(A)

    l1_penalty <- max(as.numeric(l1_penalty), 0.0)
    tv_penalty <- max(as.numeric(tv_penalty), 0.0)

    x <- pmax(as.numeric(x0), 0.0)

    if (l1_penalty == 0.0 && tv_penalty == 0.0) {
        # Degenerate ADMM: consensus with no nonsmooth terms reduces to a
        # single (augmented) NNLS solve; perform it directly.
        return(list(spectrum = .bss_nnls(A, b), iterations = 1L, converged = TRUE))
    }

    D <- .admm_difference_matrix(n)

    if (is.null(rho)) {
        scale <- max(sqrt(sum(b^2)) / max(.bss_spectral_norm(A), 1e-30), 1e-12)
        rho <- scale * max(mean(A^2), 1e-12)
    }
    rho <- max(as.numeric(rho), 1e-12)

    # Fixed augmented design; only the RHS changes across iterations.
    blocks <- list(A, sqrt(rho) * diag(n))
    if (tv_penalty > 0) blocks <- c(blocks, list(sqrt(rho) * D))
    M <- do.call(rbind, blocks)
    Mr <- nrow(M)

    z1 <- x
    z2 <- if (nrow(D) > 0) as.numeric(D %*% x) else numeric(0)
    u1 <- numeric(n)
    u2 <- numeric(nrow(D))

    converged <- FALSE
    iterations <- 0L
    for (k in seq_len(max_iterations)) {
        z1_prev <- z1
        z2_prev <- z2

        # ---- x-update: exact NNLS on the augmented system -----------------
        rhs <- c(b, sqrt(rho) * (z1 - u1))
        if (tv_penalty > 0) rhs <- c(rhs, sqrt(rho) * (z2 - u2))
        x <- .bss_nnls(M, rhs, maxiter = 10L * ncol(M))

        # ---- z-updates: exact proximal operators --------------------------
        z1 <- .admm_soft_threshold(x + u1, l1_penalty / rho)
        if (tv_penalty > 0) {
            Dx <- as.numeric(D %*% x)
            z2 <- .admm_soft_threshold(Dx + u2, tv_penalty / rho)
        }

        # ---- dual updates --------------------------------------------------
        u1 <- u1 + x - z1
        if (tv_penalty > 0) u2 <- u2 + as.numeric(D %*% x) - z2

        iterations <- k

        # ---- Boyd primal/dual residual stopping rule -----------------------
        p <- sqrt(sum((x - z1)^2))
        if (tv_penalty > 0) p <- sqrt(p^2 + sum((as.numeric(D %*% x) - z2)^2))
        dz1 <- z1 - z1_prev
        dz2 <- z2 - z2_prev
        s <- sqrt(sum(dz1^2))
        if (tv_penalty > 0) s <- sqrt(s^2 + sum(as.numeric(t(D) %*% dz2)^2))
        s <- s * rho

        n_pri <- max(sqrt(sum(x^2)), sqrt(sum(z1^2)),
                     if (tv_penalty > 0) sqrt(sum((as.numeric(D %*% x))^2)) else 0.0)
        n_dual <- sqrt(sum(u1^2)) +
            (if (tv_penalty > 0) sqrt(sum((as.numeric(t(D) %*% u2))^2)) else 0.0)
        eps_pri <- sqrt(Mr) * abstol + reltol * max(n_pri, 1e-30)
        eps_dual <- n * abstol + reltol * max(n_dual, 1e-30) * rho

        if (p <= eps_pri && s <= eps_dual) {
            converged <- TRUE
            break
        }

        # ---- adaptive rho (Boyd et al., sec. 3.4.1) ------------------------
        if (isTRUE(adaptive_rho) && (k %% 10) == 0) {
            if (p > 10.0 * s) rho_new <- 2.0 * rho
            else if (s > 10.0 * p) rho_new <- 0.5 * rho
            else rho_new <- rho
            if (rho_new != rho) {
                factor <- rho / rho_new
                u1 <- u1 * factor
                u2 <- u2 * factor
                rho <- rho_new
                blocks <- list(A, sqrt(rho) * diag(n))
                if (tv_penalty > 0) blocks <- c(blocks, list(sqrt(rho) * D))
                M <- do.call(rbind, blocks)
                Mr <- nrow(M)
            }
        }
    }

    list(spectrum = as.numeric(x), iterations = iterations, converged = converged)
}

# ---- internal helpers ------------------------------------------------------

# ---- Lawson-Hanson NNLS primitives (port of scipy.optimize.nnls) ----------
# These reproduce scipy's compiled ``nnls.c`` (Lawson-Hanson active set with
# in-place Householder QR accumulation) bit-for-bit on the tested response
# matrices.  Unlike ``lsei::nnls`` / ``quadprog``, the incremental QR keeps the
# *exact* column basis that scipy selects on degenerate (near-collinear)
# energy grids, which is required for numerical parity of the NNLS-based
# methods (ADMM x-subproblem, Frank-Wolfe / mirror-descent total fluence).

.bss_nnls_copysign <- function(x, y) {
    if (y < 0 || (y == 0 && (1 / y) < 0)) -abs(x) else abs(x)
}

.bss_nnls_dnrm2 <- function(v) {
    # LAPACK-style scaled Euclidean norm (avoids overflow / underflow).
    if (length(v) == 0L) return(0)
    mx <- max(abs(v))
    if (mx == 0) return(0)
    mx * sqrt(sum((v / mx)^2))
}

.bss_nnls_dlapy2 <- function(a, b) .bss_nnls_dnrm2(c(a, b))

.bss_nnls_ulp <- function(u) {
    # nextafter(u, 2u) - u, i.e. the ulp (machine spacing) of a positive double.
    if (u <= 0) return(0)
    2^(floor(log2(u)) - 52)
}

# LAPACK dlarfgp: elementary reflector H = I - tau [1; v][1; v]^T with beta >= 0.
.bss_nnls_dlarfgp <- function(n, alpha, x) {
    x <- as.numeric(x)
    EPSM <- .Machine$double.eps
    SMLNUM <- .Machine$double.xmin / EPSM
    if (n <= 0L) return(list(alpha = alpha, x = x, tau = 0))
    XNORM <- .bss_nnls_dnrm2(x)
    if (XNORM <= EPSM * abs(alpha)) {
        if (alpha >= 0) {
            tau <- 0
        } else {
            tau <- 2
            x <- numeric(length(x))
            alpha <- -alpha
        }
        return(list(alpha = alpha, x = x, tau = tau))
    }
    BETA <- .bss_nnls_copysign(.bss_nnls_dlapy2(alpha, XNORM), alpha)
    KNT <- 0L
    if (abs(BETA) < SMLNUM) {
        BIGNUM <- 1 / SMLNUM
        repeat {
            KNT <- KNT + 1L
            x <- x * BIGNUM; BETA <- BETA * BIGNUM; alpha <- alpha * BIGNUM
            if (!(abs(BETA) < SMLNUM && KNT < 20L)) break
        }
        XNORM <- .bss_nnls_dnrm2(x)
        BETA <- .bss_nnls_copysign(.bss_nnls_dlapy2(alpha, XNORM), alpha)
    }
    SAVEALPHA <- alpha
    alpha <- alpha + BETA
    if (BETA < 0) {
        BETA <- -BETA
        TAU <- -alpha / BETA
    } else {
        alpha <- XNORM * (XNORM / alpha)
        TAU <- alpha / BETA
        alpha <- -alpha
    }
    if (abs(TAU) <= SMLNUM) {
        if (SAVEALPHA >= 0) {
            TAU <- 0
        } else {
            TAU <- 2
            x <- numeric(length(x))
            BETA <- -SAVEALPHA
        }
    } else {
        if (abs(alpha) < SMLNUM) {
            x <- x * (1 / SMLNUM)
            x <- x * (SMLNUM / alpha)
        } else {
            x <- x * (1 / alpha)
        }
    }
    for (i in seq_len(KNT)) BETA <- BETA * SMLNUM
    list(alpha = BETA, x = x, tau = TAU)
}

# LAPACK dlartgp: plane rotation with R >= 0.
.bss_nnls_dlartgp <- function(f, g) {
    if (g == 0) {
        cs <- .bss_nnls_copysign(1, f); sn <- 0; r <- abs(f)
    } else if (f == 0) {
        cs <- 0; sn <- .bss_nnls_copysign(1, g); r <- abs(g)
    } else {
        r <- sqrt(f^2 + g^2); cs <- f / r; sn <- g / r
        if (r < 0) { cs <- -cs; sn <- -sn; r <- -r }
    }
    list(cs = cs, sn = sn, r = r)
}

# Apply H = I - tau * v v^T from the left to one column (v has leading 1).
.bss_nnls_dlarf_col <- function(v, tau, C) {
    if (tau == 0) return(C)
    w <- tau * as.numeric(crossprod(v, C))
    C - v * w
}

# scipy.optimize.nnls equivalent: returns the nonnegative solution vector.
.bss_nnls <- function(A, b, maxiter = NULL) {
    a <- as.matrix(A); storage.mode(a) <- "double"
    bb <- as.numeric(b)
    m <- nrow(a); n <- ncol(a)
    if (m <= 0L || n <= 0L) return(numeric(n))
    if (is.null(maxiter)) maxiter <- 3L * n
    x <- numeric(n); w <- numeric(n); zz <- numeric(m)
    indices <- as.integer(seq_len(n))
    indz <- 0L; iteration <- 0L

    # back-substitution through the accumulated R factor -> prospective solution
    tri <- function() {
        z <- numeric(indz)
        for (i in 0:(indz - 1L)) z[i + 1L] <- bb[i + 1L]
        jjcol <- 0L
        for (k in 0:(indz - 1L)) {
            ip <- indz - 1L - k
            if (k != 0L) {
                for (i in 0:ip) z[i + 1L] <- z[i + 1L] - a[i + 1L, jjcol] * z[ip + 2L]
            }
            jjcol <- indices[ip + 1L]
            z[ip + 1L] <- z[ip + 1L] / a[ip + 1L, jjcol]
        }
        z
    }

    while (indz < min(m, n)) {
        # dual vector for the Z set: w[j] = a[indz:m-1, j] . b[indz:m-1]
        for (i in indz:(n - 1L)) {
            j <- indices[i + 1L]
            rr <- (indz + 1L):m
            w[j] <- as.numeric(crossprod(a[rr, j], bb[rr]))
        }

        repeat {
            wmax <- 0; izmax <- -1L
            for (kk in indz:(n - 1L)) {
                j <- indices[kk + 1L]
                if (w[j] > wmax) { wmax <- w[j]; izmax <- kk }
            }
            if (wmax <= 0) return(as.numeric(x))  # KKT certificate
            iz <- izmax; j <- indices[iz + 1L]
            unorm <- if (indz > 0L) .bss_nnls_dnrm2(a[1:indz, j]) else 0
            spacing <- if (unorm > 0) .bss_nnls_ulp(unorm) else 0
            rowsl <- (indz + 1L):m
            if (.bss_nnls_dnrm2(a[rowsl, j]) > 100 * spacing) {
                nn <- m - indz
                sub <- a[rowsl, j]
                res <- .bss_nnls_dlarfgp(nn, sub[1L], sub[-1L])
                pivot <- res$alpha; tau <- res$tau; vv <- res$x
                b0 <- bb[indz + 1L]
                tauvtb <- tau * (b0 + if (nn > 1) as.numeric(crossprod(vv, bb[(indz + 2L):m])) else 0)
                if ((b0 - tauvtb) / pivot > 0) break
            }
            w[j] <- 0
        }

        # accept column j: store the reflector and apply it to a, b
        zz[rowsl] <- c(1, vv)
        a[rowsl, j] <- zz[rowsl]
        bb[indz + 1L] <- bb[indz + 1L] - tauvtb
        if (nn > 1) bb[(indz + 2L):m] <- bb[(indz + 2L):m] - tauvtb * vv
        indices[iz + 1L] <- indices[indz + 1L]; indices[indz + 1L] <- j
        indz <- indz + 1L
        if (indz < n) {
            vstart <- indz:m
            vh <- a[vstart, j]
            for (kk in indz:(n - 1L)) {
                jjc <- indices[kk + 1L]
                a[vstart, jjc] <- .bss_nnls_dlarf_col(vh, tau, a[vstart, jjc])
            }
        }
        a[indz, j] <- pivot
        if (indz < m) a[(indz + 1L):m, j] <- 0
        w[j] <- 0
        zz[seq_len(indz)] <- tri()

        repeat {
            iteration <- iteration + 1L
            if (iteration >= maxiter) return(as.numeric(x))
            alpha <- 2; jjpos <- -1L
            for (ip in 0:(indz - 1L)) {
                k <- indices[ip + 1L]
                if (zz[ip + 1L] <= 0) {
                    Tv <- -x[k] / (zz[ip + 1L] - x[k])
                    if (alpha > Tv) { alpha <- Tv; jjpos <- ip }
                }
            }
            if (alpha == 2) break
            for (ip in 0:(indz - 1L)) { k <- indices[ip + 1L]; x[k] <- x[k] + alpha * (zz[ip + 1L] - x[k]) }
            i <- indices[jjpos + 1L]
            repeat {
                x[i] <- 0
                if (jjpos != indz - 1L) {
                    jjp <- jjpos + 1L
                    while (jjp <= indz - 1L) {
                        jpos0 <- jjp
                        ii <- indices[jpos0 + 1L]
                        indices[jpos0] <- ii
                        gg <- .bss_nnls_dlartgp(a[jpos0, ii], a[jpos0 + 1L, ii])
                        a[jpos0, ii] <- gg$r; a[jpos0 + 1L, ii] <- 0
                        cc <- gg$cs; ss <- gg$sn
                        for (kc in seq_len(n)) {
                            if (kc != ii) {
                                tmp <- a[jpos0, kc]
                                a[jpos0, kc] <- cc * tmp + ss * a[jpos0 + 1L, kc]
                                a[jpos0 + 1L, kc] <- -ss * tmp + cc * a[jpos0 + 1L, kc]
                            }
                        }
                        tmp <- bb[jpos0]
                        bb[jpos0] <- cc * tmp + ss * bb[jpos0 + 1L]
                        bb[jpos0 + 1L] <- -ss * tmp + cc * bb[jpos0 + 1L]
                        jjp <- jjp + 1L
                    }
                }
                indz <- indz - 1L
                indices[indz + 1L] <- i
                nobreak <- 0L
                for (jj2 in 0:(indz - 1L)) {
                    i2 <- indices[jj2 + 1L]
                    if (x[i2] <= 0) { jjpos <- jj2; break }
                    if (jj2 == indz - 1L) { nobreak <- 1L; jjpos <- jj2 }
                }
                if (nobreak == 1L) break
            }
            zz <- numeric(indz)
            for (i in 0:(indz - 1L)) zz[i + 1L] <- bb[i + 1L]
            jjcol <- 0L
            for (k in 0:(indz - 1L)) {
                ip <- indz - 1L - k
                if (k != 0L) {
                    for (i in 0:ip) zz[i + 1L] <- zz[i + 1L] - a[i + 1L, jjcol] * zz[ip + 2L]
                }
                jjcol <- indices[ip + 1L]
                zz[ip + 1L] <- zz[ip + 1L] / a[ip + 1L, jjcol]
            }
        }

        for (kk in 0:(indz - 1L)) { i <- indices[kk + 1L]; x[i] <- zz[kk + 1L] }
    }
    as.numeric(x)
}

.admm_soft_threshold <- function(v, threshold) {
    v <- as.numeric(v)
    if (threshold <= 0) return(v)
    sign(v) * pmax(abs(v) - threshold, 0.0)
}

.admm_difference_matrix <- function(n) {
    # First-order difference operator D with (D x)_i = x_{i+1} - x_i.
    if (n < 2L) return(matrix(numeric(0), nrow = 0L, ncol = n))
    D <- matrix(0.0, nrow = n - 1L, ncol = n)
    idx <- seq_len(n - 1L)
    D[cbind(idx, idx)] <- -1.0
    D[cbind(idx, idx + 1L)] <- 1.0
    D
}

.bss_spectral_norm <- function(A) {
    # Largest singular value of A (||A||_2), i.e. the 2-norm used by Python.
    svd(as.matrix(A), nu = 0L, nv = 0L)$d[1L]
}

#' @rdname solve_admm
#' @inheritParams run_unfolding
#' @param initial_spectrum Optional initial spectrum guess.
#' @param rho Optional ADMM penalty parameter.
#' @param adaptive_rho Logical; adapt \code{rho} every 10 iterations.
#' @param calculate_errors Logical; run Monte-Carlo uncertainty.
#' @param noise_level Numeric; MC noise level (must be in (0, 1]).
#' @param n_montecarlo Integer; MC samples.
#' @param save_result Logical; call the save callback.
#' @param random_state Optional integer seed.
#' @param max_neutron_energy Optional energy cutoff in MeV.
#' @return A result list produced by \code{\link{run_unfolding}}.
#' @export
unfold_admm <- function(detector_names, n_energy_bins, E_MeV,
                        sensitivities, cc_icrp116, save_result_callback,
                        readings, initial_spectrum = NULL,
                        max_iterations = 500L, tolerance = 1e-6,
                        l1_penalty = 0.0, tv_penalty = 0.0, rho = NULL,
                        adaptive_rho = TRUE,
                        calculate_errors = FALSE, noise_level = 0.01,
                        n_montecarlo = 100L, save_result = FALSE,
                        random_state = NULL, max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_admm,
                                        max_iterations = max_iterations,
                                        tolerance = tolerance,
                                        l1_penalty = l1_penalty,
                                        tv_penalty = tv_penalty,
                                        rho = rho,
                                        adaptive_rho = adaptive_rho),
        solve_kwargs = list(),
        method_name = "ADMM",
        extra_output = list(l1_penalty = l1_penalty,
                            tv_penalty = tv_penalty,
                            adaptive_rho = adaptive_rho),
        calculate_errors = calculate_errors,
        noise_level = noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
