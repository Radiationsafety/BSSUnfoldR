#' FRUIT-like parametric unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_fruit_like.py}.
#' Implements a parametric unfolding method inspired by FRUIT (Fast Real-time
#' Unfolding of neutron spectra with Iterative parameterNization Technique,
#' Bedogni et al., NIM A 580 (2007)).
#'
#' The parametric model consists of:
#' \itemize{
#'   \item Maxwellian thermal component: \eqn{A_{th} \sqrt{E} \exp(-E/T_{th})}
#'   \item 1/E epithermal component: \eqn{A_{epi} / E}
#'   \item Evaporation fast component: \eqn{A_f \exp(-E/T_{ev})}
#' }
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid (length n).
#' @param log_steps Numeric log-energy bin widths (length n).
#' @param initial_params Optional named list with \code{A_th}, \code{T_th},
#'   \code{A_epi}, \code{A_f}, \code{T_ev}.
#' @param method Character; optimizer method passed to lmfit-equivalent
#'   solver. Default \code{"leastsq"} (Levenberg-Marquardt via MINPACK,
#'   matching Python's lmfit default).
#' @return A list \code{list(spectrum, success, nfev, converged)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_fruit_like(A, b, E, compute_log_steps(E))
solve_fruit_like <- function(A, b, E, log_steps,
                                initial_params = NULL, method = "leastsq") {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E)
    log_steps <- as.numeric(log_steps)

    # Python: param name -> c(default, min, max); max = Inf means "min only"
    fl_bounds <- list(A_th = c(1e-6, 0, Inf),
                      T_th = c(0.025e-6, 1e-9, 1e-3),
                      A_epi = c(1e-6, 0, Inf),
                      A_f = c(1e-6, 0, Inf),
                      T_ev = c(2.0, 0.1, 20))

    if (is.null(initial_params)) {
        names_used <- names(fl_bounds)
    } else {
        # exactly like Python: only the keys the caller supplied become fit
        # variables (in the order they were given), unknown keys are ignored
        ip <- as.list(initial_params)
        names_used <- names(ip)[names(ip) %in% names(fl_bounds)]
        if (length(names_used) == 0L) {
            names_used <- names(fl_bounds)
            ip <- NULL
        }
    }
    p0 <- vapply(names_used,
                 function(nm) if (is.null(initial_params)) {
                     fl_bounds[[nm]][1L]
                 } else as.numeric(initial_params[[nm]]), numeric(1L))
    lower <- vapply(names_used, function(nm) fl_bounds[[nm]][2L], numeric(1L))
    upper <- vapply(names_used, function(nm) fl_bounds[[nm]][3L], numeric(1L))
    names(p0) <- names_used

    .residuals <- function(p) {
        spectrum <- .parametric_model(E, p[["A_th"]], p[["T_th"]],
                                        p[["A_epi"]], p[["A_f"]], p[["T_ev"]])
        as.numeric(A %*% (spectrum * log_steps)) - b
    }
    .spectrum_from <- function(p) {
        as.numeric(.parametric_model(E, p[["A_th"]], p[["T_th"]],
                                    p[["A_epi"]], p[["A_f"]], p[["T_ev"]]) *
                   log_steps)
    }

    if (method %in% c("leastsq", "lm", "l-bfgs-b", "L-BFGS-B")) {
        # lmfit/scipy leastsq: Levenberg-Marquardt needs at least as many
        # residuals as variables
        n_par <- length(p0)
        n_dat <- length(b)
        if (n_par > n_dat) {
            stop(sprintf(paste0("Improper input: func input vector length ",
                                "N=%d must not exceed func output vector ",
                                "length M=%d"), n_par, n_dat), call. = FALSE)
        }
        fit <- .fl_leastsq(.residuals, p0, lower, upper)
        return(list(spectrum = .spectrum_from(fit$values),
                    success = fit$success,
                    nfev = fit$nfev,
                    converged = fit$success,
                    message = fit$message,
                    info = fit$info,
                    params = as.list(fit$values)))
    }

    # explicit non-Levenberg-Marquardt method: bounded/unbounded optim
    .objective <- function(p) sum(.residuals(p)^2)
    if (method %in% c("Brent", "CG", "BFGS", "Nelder-Mead", "SANN")) {
        result <- stats::optim(unname(p0), .objective, method = method,
                                control = list(maxit = 1000))
    } else {
        result <- stats::optim(unname(p0), .objective, method = method,
                                lower = unname(lower),
                                upper = pmax(unname(upper), .Machine$double.xmax),
                                control = list(maxit = 1000))
    }
    p_opt <- setNames(as.numeric(result$par), names_used)
    list(spectrum = .spectrum_from(p_opt),
         success = (result$convergence == 0),
         nfev = as.integer(result$counts[1L]),
         converged = (result$convergence == 0),
         message = result$message,
         params = as.list(p_opt))
}

#' FRUIT-like three-component parametric spectrum (port of
#' \code{parametric_model})
#' @keywords internal
#' @noRd
.parametric_model <- function(E, A_th, T_th, A_epi, A_f, T_ev,
                                epi_max = 0.1) {
    E <- as.numeric(E)
    spectrum <- numeric(length(E))
    thermal <- E < 0.4e-6
    epithermal <- (E >= 0.4e-6) & (E < epi_max)
    fast <- E >= epi_max
    spectrum[thermal] <- spectrum[thermal] +
        (A_th * sqrt(E[thermal]) * exp(-E[thermal] / T_th))
    spectrum[epithermal] <- spectrum[epithermal] +
        (A_epi / (E[epithermal] + 1e-15))
    spectrum[fast] <- spectrum[fast] +
        (A_f * exp(-E[fast] / T_ev))
    spectrum
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------- #
#  Faithful MINPACK ``lmdif`` engine (verbatim transcription of the C     #
#  sources scipy compiles for ``scipy.optimize.leastsq``, so that the     #
#  fit trajectory -- and therefore the local minimum reached -- matches   #
#  the Python lmfit reference bit-for-bit.  Every accumulation is a       #
#  sequential IEEE-double loop (never ``sum()``, which R promotes to long #
#  double), every power matches the C operator (``pow`` via ``^2`` where  #
#  C writes ``pow(...,2.0)``, plain ``*`` where C writes ``x*x``).        #
#  Duplicated inside this file because \code{R/unfold_parametric.R} is    #
#  owned by another agent.                                                #
# ---------------------------------------------------------------------- #

# dpmpar constants (np.finfo(float64)): epsmch = 2.220446049250313e-16,
# dwarf = 2.2250738585072014e-308.
.fl_epsmch <- 2.220446049250313e-16
.fl_dwarf <- 2.2250738585072014e-308

#' ENORM -- scaled Euclidean norm (transcription of the C ``enorm``).
#' @noRd
.fl_enorm <- function(x) {
    x <- as.numeric(x)
    n <- length(x)
    rgiant <- 2^63.5
    rdwarf <- 0.5 / rgiant
    s1 <- 0; s2 <- 0; s3 <- 0
    x1max <- 0; x3max <- 0
    agiant <- rgiant / n
    for (i in seq_len(n)) {
        xabs <- abs(x[i])
        if (xabs <= rdwarf || xabs >= agiant) {
            if (xabs > rdwarf) {
                if (xabs > x1max) {
                    s1 <- 1 + s1 * (x1max / xabs)^2
                    x1max <- xabs
                } else {
                    s1 <- s1 + (xabs / x1max)^2
                }
            } else {
                if (xabs > x3max) {
                    s3 <- 1 + s3 * (x3max / xabs)^2
                    x3max <- xabs
                } else if (xabs != 0) {
                    s3 <- s3 + (xabs / x3max)^2
                }
            }
        } else {
            s2 <- s2 + xabs * xabs
        }
    }
    if (s1 != 0) return(x1max * sqrt(s1 + (s2 / x1max) / x1max))
    if (s2 != 0) {
        if (s2 >= x3max) return(sqrt(s2 * (1 + (x3max / s2) * (x3max * s3))))
        return(sqrt(x3max * ((s2 / x3max) + (x3max * s3))))
    }
    x3max * sqrt(s3)
}

#' FDJAC2 -- forward-difference Jacobian (transcription of the C ``fdjac2``).
#' @noRd
.fl_fdjac2 <- function(fcn, x, fvec, epsfcn) {
    x <- as.numeric(x); fvec <- as.numeric(fvec)
    m <- length(fvec); n <- length(x)
    eps <- sqrt(max(epsfcn, .fl_epsmch))
    fjac <- matrix(0, nrow = m, ncol = n)
    for (j in seq_len(n)) {
        temp <- x[j]
        h <- eps * abs(temp)
        if (h == 0) h <- eps
        x[j] <- temp + h
        wa <- as.numeric(fcn(x))
        x[j] <- temp
        for (i in seq_len(m)) fjac[i, j] <- (wa[i] - fvec[i]) / h
    }
    fjac
}

#' QRFAC -- Householder QR with column pivoting (transcription of ``qrfac``).
#' @noRd
.fl_qrfac <- function(a, pivot = TRUE) {
    storage.mode(a) <- "double"
    m <- nrow(a); n <- ncol(a)
    acnorm <- numeric(n); rdiag <- numeric(n); wa <- numeric(n)
    ipvt <- seq_len(n)
    for (j in seq_len(n)) {
        acnorm[j] <- .fl_enorm(a[, j])
        rdiag[j] <- acnorm[j]
        wa[j] <- rdiag[j]
    }
    minmn <- min(m, n)
    for (j in seq_len(minmn)) {
        if (pivot) {
            kmax <- j
            if (n > j) for (k in (j + 1L):n) if (rdiag[k] > rdiag[kmax]) kmax <- k
            if (kmax != j) {
                for (i in seq_len(m)) {
                    temp <- a[i, j]; a[i, j] <- a[i, kmax]; a[i, kmax] <- temp
                }
                rdiag[kmax] <- rdiag[j]
                wa[kmax] <- wa[j]
                k <- ipvt[j]; ipvt[j] <- ipvt[kmax]; ipvt[kmax] <- k
            }
        }
        ajnorm <- .fl_enorm(a[j:m, j])
        if (ajnorm == 0) {
            rdiag[j] <- -ajnorm
            next
        }
        if (a[j, j] < 0) ajnorm <- -ajnorm
        for (i in j:m) a[i, j] <- a[i, j] / ajnorm
        a[j, j] <- a[j, j] + 1
        if (j == n) {
            rdiag[j] <- -ajnorm
            next
        }
        for (k in (j + 1L):n) {
            ssum <- 0
            for (i in j:m) ssum <- ssum + a[i, j] * a[i, k]
            temp <- ssum / a[j, j]
            for (i in j:m) a[i, k] <- a[i, k] - temp * a[i, j]
            if (pivot && rdiag[k] != 0) {
                temp <- a[j, k] / rdiag[k]
                rdiag[k] <- rdiag[k] * sqrt(max(0, 1 - temp * temp))
                if (0.05 * (rdiag[k] / wa[k])^2 <= .fl_epsmch) {
                    rdiag[k] <- .fl_enorm(a[(j + 1L):m, k])
                    wa[k] <- rdiag[k]
                }
            }
        }
        rdiag[j] <- -ajnorm
    }
    list(a = a, acnorm = acnorm, rdiag = rdiag, wa = wa, ipvt = ipvt)
}

#' QRSOLV -- Givens-rotation subspace solve (transcription of ``qrsolv``).
#' @noRd
.fl_qrsolv <- function(r, ipvt, dvec, qtb, n) {
    x <- numeric(n); wa <- numeric(n); sdiag <- numeric(n)
    for (j in seq_len(n)) {
        for (i in j:n) r[i, j] <- r[j, i]
        x[j] <- r[j, j]
        wa[j] <- qtb[j]
    }
    for (j in seq_len(n)) {
        l <- ipvt[j]
        if (dvec[l] != 0) {
            for (k in j:n) sdiag[k] <- 0
            sdiag[j] <- dvec[l]
            qtbpj <- 0
            for (k in j:n) {
                if (sdiag[k] != 0) {
                    if (abs(r[k, k]) < abs(sdiag[k])) {
                        cotan <- r[k, k] / sdiag[k]
                        ssin <- 0.5 / sqrt(0.25 + 0.25 * cotan^2)
                        ccos <- ssin * cotan
                    } else {
                        tan <- sdiag[k] / r[k, k]
                        ccos <- 0.5 / sqrt(0.25 + 0.25 * tan^2)
                        ssin <- ccos * tan
                    }
                    r[k, k] <- r[k, k] * ccos
                    r[k, k] <- r[k, k] + ssin * sdiag[k]
                    temp <- ccos * wa[k] + ssin * qtbpj
                    qtbpj <- -ssin * wa[k] + ccos * qtbpj
                    wa[k] <- temp
                    if (k != n) for (i in (k + 1L):n) {
                        temp <- ccos * r[i, k] + ssin * sdiag[i]
                        sdiag[i] <- -ssin * r[i, k] + ccos * sdiag[i]
                        r[i, k] <- temp
                    }
                }
            }
        }
        sdiag[j] <- r[j, j]
        r[j, j] <- x[j]
    }
    nsing <- n
    for (j in seq_len(n)) {
        if (sdiag[j] == 0 && nsing == n) nsing <- j - 1L
        if (nsing < n) wa[j] <- 0
    }
    if (nsing > 0) for (k in seq_len(nsing)) {
        j <- nsing - k + 1L
        ssum <- 0
        if (nsing >= j + 1L) for (i in (j + 1L):nsing) ssum <- ssum + r[i, j] * wa[i]
        wa[j] <- (wa[j] - ssum) / sdiag[j]
    }
    for (j in seq_len(n)) { l <- ipvt[j]; x[l] <- wa[j] }
    list(r = r, x = x, sdiag = sdiag)
}

#' LMPAR -- LM parameter and direction (transcription of ``lmpar``).
#' @noRd
.fl_lmpar <- function(r, ipvt, diag, qtb, delta, par, n) {
    x <- numeric(n); sdiag <- numeric(n)
    wa1 <- numeric(n); wa2 <- numeric(n)
    nsing <- n
    for (j in seq_len(n)) {
        wa1[j] <- qtb[j]
        if (r[j, j] == 0 && nsing == n) nsing <- j - 1L
        if (nsing < n) wa1[j] <- 0
    }
    if (nsing > 0) for (k in seq_len(nsing)) {
        j <- nsing - k + 1L
        wa1[j] <- wa1[j] / r[j, j]
        temp <- wa1[j]
        if (j >= 2L) for (i in seq_len(j - 1L)) wa1[i] <- wa1[i] - r[i, j] * temp
    }
    for (j in seq_len(n)) { l <- ipvt[j]; x[l] <- wa1[j] }
    iter <- 0L
    for (j in seq_len(n)) wa2[j] <- diag[j] * x[j]
    dxnorm <- .fl_enorm(wa2)
    fp <- dxnorm - delta
    if (fp <= 0.1 * delta) {
        if (iter == 0L) par <- 0
        return(list(r = r, x = x, sdiag = sdiag, par = par))
    }
    parl <- 0
    if (nsing == n) {
        for (j in seq_len(n)) {
            l <- ipvt[j]
            wa1[j] <- diag[l] * (wa2[l] / dxnorm)
        }
        for (j in seq_len(n)) {
            ssum <- 0
            if (j >= 2L) for (i in seq_len(j - 1L)) ssum <- ssum + r[i, j] * wa1[i]
            wa1[j] <- (wa1[j] - ssum) / r[j, j]
        }
        temp <- .fl_enorm(wa1)
        parl <- ((fp / delta) / temp) / temp
    }
    for (j in seq_len(n)) {
        ssum <- 0
        for (i in seq_len(j)) ssum <- ssum + r[i, j] * qtb[i]
        l <- ipvt[j]
        wa1[j] <- ssum / diag[l]
    }
    gnorm <- .fl_enorm(wa1)
    paru <- gnorm / delta
    if (paru == 0) paru <- .fl_dwarf / min(delta, 0.1)
    par <- min(max(par, parl), paru)
    if (par == 0) par <- gnorm / dxnorm
    repeat {
        iter <- iter + 1L
        if (par == 0) par <- max(.fl_dwarf, paru * 0.001)
        temp <- sqrt(par)
        for (j in seq_len(n)) wa1[j] <- temp * diag[j]
        qs <- .fl_qrsolv(r, ipvt, wa1, qtb, n)
        r <- qs$r; x <- qs$x; sdiag <- qs$sdiag
        for (j in seq_len(n)) wa2[j] <- diag[j] * x[j]
        dxnorm <- .fl_enorm(wa2)
        temp <- fp
        fp <- dxnorm - delta
        if (abs(fp) <= 0.1 * delta ||
            (parl == 0 && fp <= temp && temp < 0) || iter == 10L) break
        for (j in seq_len(n)) {
            l <- ipvt[j]
            wa1[j] <- diag[l] * (wa2[l] / dxnorm)
        }
        for (j in seq_len(n)) {
            wa1[j] <- wa1[j] / sdiag[j]
            temp <- wa1[j]
            if (j != n) for (i in (j + 1L):n) wa1[i] <- wa1[i] - r[i, j] * temp
        }
        temp <- .fl_enorm(wa1)
        parc <- ((fp / delta) / temp) / temp
        if (fp > 0) parl <- max(parl, par)
        if (fp < 0) paru <- min(paru, par)
        par <- max(parl, par + parc)
    }
    if (iter == 0L) par <- 0
    list(r = r, x = x, sdiag = sdiag, par = par)
}

#' LMDIF -- Levenberg-Marquardt with forward-difference Jacobian
#' (transcription of the C ``LMDIF``).
#' @noRd
.fl_lmdif <- function(fcn, x0, ftol = 1.5e-8, xtol = 1.5e-8, gtol = 0,
                      maxfev = 4000L, epsfcn = 1e-10, factor = 100,
                      mode = 1L) {
    x <- as.numeric(x0)
    n <- length(x)
    info <- 0L; nfev <- 0L
    delta <- 0; xnorm <- 0
    fvec <- as.numeric(fcn(x))
    nfev <- 1L
    m <- length(fvec)
    fnorm <- .fl_enorm(fvec)
    par <- 0
    iter <- 1L
    diag <- numeric(n)
    qtf <- numeric(n); wa1 <- numeric(n); wa2 <- numeric(n)
    wa3 <- numeric(n); wa4 <- numeric(m); ipvt <- seq_len(n)
    fjac <- matrix(0, nrow = m, ncol = n)
    fvec_out <- fvec
    repeat {
        fjac <- .fl_fdjac2(fcn, x, fvec, epsfcn)
        nfev <- nfev + n
        qr <- .fl_qrfac(fjac, TRUE)
        fjac <- qr$a; wa1 <- qr$rdiag; wa2 <- qr$acnorm
        wa3 <- qr$wa; ipvt <- qr$ipvt
        if (iter == 1L) {
            if (mode != 2L) for (j in seq_len(n)) {
                diag[j] <- wa2[j]
                if (wa2[j] == 0) diag[j] <- 1
            }
            for (j in seq_len(n)) wa3[j] <- diag[j] * x[j]
            xnorm <- .fl_enorm(wa3)
            delta <- factor * xnorm
            if (delta == 0) delta <- factor
        }
        # Form Q' fvec, keep first n in qtf, fold Householder into wa4.
        wa4 <- fvec
        for (j in seq_len(n)) {
            if (fjac[j, j] != 0) {
                ssum <- 0
                for (i in j:m) ssum <- ssum + fjac[i, j] * wa4[i]
                temp <- -ssum / fjac[j, j]
                for (i in j:m) wa4[i] <- wa4[i] + fjac[i, j] * temp
            }
            fjac[j, j] <- wa1[j]
            qtf[j] <- wa4[j]
        }
        # Norm of scaled gradient.
        gnorm <- 0
        if (fnorm != 0) for (j in seq_len(n)) {
            l <- ipvt[j]
            if (wa2[l] != 0) {
                ssum <- 0
                for (i in seq_len(j)) ssum <- ssum + fjac[i, j] * (qtf[i] / fnorm)
                gnorm <- max(gnorm, abs(ssum / wa2[l]))
            }
        }
        if (gnorm <= gtol) info <- 4L
        if (info != 0L) break
        if (mode != 2L) for (j in seq_len(n)) diag[j] <- max(diag[j], wa2[j])
        ## inner loop
        repeat {
            lm <- .fl_lmpar(fjac, ipvt, diag, qtf, delta, par, n)
            fjac <- lm$r; par <- lm$par
            for (j in seq_len(n)) {
                wa1[j] <- -lm$x[j]
                wa2[j] <- x[j] + wa1[j]
                wa3[j] <- diag[j] * wa1[j]
            }
            pnorm <- .fl_enorm(wa3)
            if (iter == 1L) delta <- min(delta, pnorm)
            wa4 <- as.numeric(fcn(wa2))
            nfev <- nfev + 1L
            fnorm1 <- .fl_enorm(wa4)
            actred <- if (0.1 * fnorm1 < fnorm) 1 - (fnorm1 / fnorm)^2 else -1
            for (j in seq_len(n)) {
                wa3[j] <- 0
                l <- ipvt[j]
                temp <- wa1[l]
                for (i in seq_len(j)) wa3[i] <- wa3[i] + fjac[i, j] * temp
            }
            temp1 <- .fl_enorm(wa3) / fnorm
            temp2 <- (sqrt(par) * pnorm) / fnorm
            prered <- temp1 * temp1 + temp2 * temp2 / 0.5
            dirder <- -(temp1 * temp1 + temp2 * temp2)
            ratio <- if (prered != 0) actred / prered else 0
            if (ratio <= 0.25) {
                temp <- if (actred < 0) (0.5 * dirder) / (dirder + 0.5 * actred) else 0.5
                if (0.1 * fnorm1 >= fnorm || temp < 0.1) temp <- 0.1
                delta <- temp * min(delta, pnorm / 0.1)
                par <- par / temp
            } else {
                if (par == 0 || ratio >= 0.75) {
                    delta <- pnorm / 0.5
                    par <- par * 0.5
                }
            }
            if (ratio >= 1e-4) {
                for (j in seq_len(n)) {
                    x[j] <- wa2[j]
                    wa2[j] <- diag[j] * x[j]
                }
                fvec_out <- wa4
                xnorm <- .fl_enorm(wa2)
                fnorm <- fnorm1
                iter <- iter + 1L
            }
            info <- if (abs(actred) <= ftol && prered <= ftol &&
                        0.5 * ratio <= 1) 1L else 0L
            if (delta <= xtol * xnorm) info <- 2L
            if (abs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1 &&
                info == 2L) info <- 3L
            if (info != 0L) break
            if (nfev >= maxfev) info <- 5L
            if (abs(actred) <= .fl_epsmch && prered <= .fl_epsmch &&
                0.5 * ratio <= 1) info <- 6L
            if (delta <= .fl_epsmch * xnorm) info <- 7L
            if (gnorm <= .fl_epsmch) info <- 8L
            if (info != 0L) break
            if (ratio >= 1e-4) break
        }
        if (info != 0L) break
    }
    list(x = x, fvec = fvec_out, info = info, nfev = nfev, iter = iter,
         xnorm = xnorm, fnorm = fnorm)
}

#' Bounded least-squares driver mirroring
#' \code{lmfit.minimize(..., method = "leastsq")} on top of \code{.fl_lmdif}.
#' @noRd
.fl_leastsq <- function(resid, values, lower, upper) {
    blk <- .lmfit_setup_bounds(values, lower, upper)
    fcn <- function(xi) as.numeric(resid(.lmfit_from_internal(blk, xi)))
    n <- length(values)
    out <- .fl_lmdif(fcn, blk$internal, ftol = 1.5e-8, xtol = 1.5e-8,
                     gtol = 0, maxfev = 4000L * (n + 1L), epsfcn = 1e-10,
                     factor = 100, mode = 1L)
    ext <- .lmfit_from_internal(blk, out$x)
    names(ext) <- names(values)
    info <- out$info
    success <- info %in% c(1L, 2L, 3L, 4L)
    message <- if (info %in% c(1L, 2L, 3L)) {
        "Fit succeeded."
    } else if (info == 0L) {
        paste("Invalid Input Parameters. I.e. more variables than data",
              "points given, tolerance < 0.0, or no data provided.")
    } else if (info == 4L) {
        "One or more variable did not affect the fit."
    } else if (info == 5L) {
        sprintf("the maximum number of calls (%d) to the function",
                4000L * (n + 1L))
    } else {
        "Tolerance seems to be too small."
    }
    nfev <- max(out$nfev - 3L, 0L)
    list(values = ext, info = info, nfev = nfev, success = success,
         message = message)
}

#' Wrapper around \code{\link{solve_fruit_like}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_fruit_like
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_fruit_like <- function(detector_names, n_energy_bins, E_MeV,
                                 sensitivities, cc_icrp116, save_result_callback,
                                 readings, initial_spectrum = NULL,
                                 initial_params = NULL, method = "leastsq",
                                 calculate_errors = FALSE,
                                 noise_level = 0.01, n_montecarlo = 100L,
                                 save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    log_steps <- compute_log_steps(E_MeV) * log(10.0)
    x0_default <- rep(mean(b) / max(mean(rowSums(A)), 1e-30), n_energy_bins)
    solver <- function(A, b, x0 = NULL, ...) {
        res <- solve_fruit_like(A, b, E_MeV, log_steps,
                                initial_params = initial_params, method = method)
        list(spectrum = res$spectrum, iterations = res$nfev,
             converged = res$success)
    }
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver,
        solve_kwargs = list(),
        method_name = "fruit_like",
        extra_output = list(initial_params = initial_params, method = method),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
