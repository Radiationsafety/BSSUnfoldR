#' FRUIT parametric unfolding (full 3-component model)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_parametric.py} and its
#' model layer \code{bssunfold/src/bssunfold/core/_fruit.py}.
#' The spectrum is a weighted superposition of three components with
#' constraint P_th + P_epi + P_f = 1:
#' \describe{
#'   \item{Thermal}{\eqn{(E/T_0^2) \exp(-E/T_0)} for E < 1e-7 MeV}
#'   \item{Epithermal}{\eqn{[1-\exp(-(E/E_d)^2)] E^{b-1} \exp(-E/\beta')}
#'     for 1e-7 <= E < 0.1 MeV}
#'   \item{Fast}{\eqn{E^\alpha \exp(-E/\beta)} for E >= 0.1 MeV}
#' }
#'
#' @description The internal helpers in this file (\code{.fruit_model},
#'   \code{.fruit_find_initial_params}, \code{.fruit_compute_jacobian}, ...)
#'   are the exact port of the Python \code{_fruit.py} layer and are shared
#'   with \code{\link{solve_hybrid_parametric}} and the QP engine variants.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid (length n).
#' @param log_steps Numeric \eqn{d(\ln E)} bin widths (length n).  As in the
#'   Python original the model spectrum is multiplied by these steps
#'   \emph{on top of} the already lethargy-weighted response matrix, so the
#'   caller must pass \code{compute_log_steps(E) * log(10)}.
#' @param initial_params Optional named list or vector overriding any of
#'   \code{b}, \code{beta_prime}, \code{alpha}, \code{beta}, \code{P_th},
#'   \code{P_epi}.  When \code{NULL} a grid scan supplies the starting points.
#' @param method Character; backend name kept for parity with the Python
#'   \code{lmfit} selector.  The pure-R port always runs the ported
#'   Levenberg-Marquardt least-squares solver (\code{lmfit} \code{leastsq}).
#' @param alpha Numeric; Tikhonov weight penalising deviation from the
#'   starting parameters (\code{sqrt(alpha) * (p - p0)} appended to the
#'   residual vector, exactly like Python \code{_residuals}). Default 0.
#' @param alpha_auto Logical; select \code{alpha} by GCV on the linearised
#'   model (Python \code{_gcv_select_alpha}). Default \code{FALSE}.
#' @param n_restarts Integer; number of multistart runs taken from the best
#'   grid-scan candidates. Default 5.
#' @param max_iterations,tolerance Accepted for API compatibility but not
#'   used: the ported solver must apply the fixed
#'   \code{lmfit}/{@scipy} \code{leastsq} tolerances
#'   (\code{ftol = xtol = 1.5e-8}) and evaluation budget to reproduce the
#'   original result.
#' @param x0 Ignored: the parametric model builds its own starting point.
#'   Accepted so \code{run_unfolding()} can always forward the initial
#'   spectrum.
#' @return A list \code{list(spectrum, iterations, converged, params,
#'   success, message)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_parametric(A, b, E, compute_log_steps(E) * log(10),
#'                        max_iterations = 50)
solve_parametric <- function(A, b, E, log_steps, initial_params = NULL,
                                method = "leastsq", alpha = 0,
                                alpha_auto = FALSE, n_restarts = 5L,
                                max_iterations = 200L, tolerance = 1e-6,
                                x0 = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E); log_steps <- as.numeric(log_steps)

    # ---- starting points: grid scan over (P_th, P_epi) unless given ----
    if (is.null(initial_params)) {
        starts <- .fruit_find_initial_params(A, b, E, log_steps,
                                            n_grid = 7L,
                                            return_top = as.integer(n_restarts))
    } else {
        starts <- list(.fruit_get_initial_params(initial_params))
    }

    # GCV-selected Tikhonov weight (uses the best starting point)
    if (isTRUE(alpha_auto)) {
        alpha <- .fruit_gcv_select_alpha(A, b, E, log_steps, starts[[1L]])
    }

    bounds <- .fruit_get_param_bounds()
    lower <- vapply(bounds, function(x) x[1L], numeric(1L))
    upper <- vapply(bounds, function(x) x[2L], numeric(1L))

    best_spectrum <- NULL
    best_residual <- Inf
    best_params <- NULL
    best_success <- FALSE
    best_message <- ""
    best_iter <- 0L
    total_nfev <- 0L

    for (start in starts) {
        p0 <- unlist(start[.fruit_PARAM_NAMES])
        # exactly Python's _residuals: data residual plus the optional
        # sqrt(alpha) * (p - p0) Tikhonov tail
        resid <- function(p) {
            sp <- .fruit_model(E, p[1L], p[2L], p[3L], p[4L], p[5L], p[6L]) *
                log_steps
            res <- as.numeric(A %*% sp) - b
            if (alpha > 0) {
                res <- c(res, sqrt(alpha) * (as.numeric(p) - p0))
            }
            res
        }
        fit <- .lmfit_leastsq(resid, p0, lower, upper)
        total_nfev <- total_nfev + fit$nfev
        p_opt <- fit$values
        spectrum <- .fruit_model_vec(E, p_opt) * log_steps
        res_norm <- sqrt(sum((as.numeric(A %*% spectrum) - b)^2))
        if (res_norm < best_residual) {
            best_residual <- res_norm
            best_spectrum <- spectrum
            best_params <- p_opt
            best_success <- fit$success
            best_message <- fit$message
            best_iter <- fit$nfev
        }
    }

    .fruit_check_fit_quality(best_residual, b, "parametric")
    list(spectrum = as.numeric(best_spectrum),
         iterations = total_nfev,
         converged = is.finite(best_residual) && best_success,
         success = best_success,
         message = best_message,
         params = best_params,
         nfev = total_nfev, fit_iterations = best_iter)
}

# ------------------------------------------------------------------ #
#  QP engines: ports of solve_parametric_cvxpy / _qpsolvers / _combined
#
#  The cvxpy and the qpsolvers SQP paths build the very same
#  bound-constrained convex QP for the parameter update,
#
#    min_delta  ||A_eff delta + r||^2 + alpha ||delta||^2,
#    lo - p <= delta <= hi - p,
#
#  cvxpy writes it as cp.sum_squares(...) and qpsolvers in the library
#  convention 0.5 d'P d + q'd with P = A_eff'A_eff + alpha I,
#  q = A_eff'r; the two objectives differ only by the constant factor 2,
#  so their minimisers coincide and one engine serves both paths.  That
#  engine is the package's certified active-set box-QP solver
#  \code{.commercial_box_qp} (R/unfold_qp_commercial.R): for a PSD
#  quadratic with bound constraints it returns the exact minimiser plus a
#  KKT certificate, i.e. a global-optimality proof, in place of the
#  ECOS/OSQP iterates Python obtains from its licensed-free backends.
#  \code{.commercial_box_qp} minimises 0.5 x'Hx + c'x, hence H = P and
#  c = q below.
# ------------------------------------------------------------------ #

#' Box QP with an exact Jacobi preconditioning
#'
#' \code{.commercial_box_qp} measures its working tolerance and its
#' "on a bound" snapping against \code{max(1, |c|, diag(H))}.  The FRUIT
#' Jacobian columns span more than fifteen orders of magnitude (the thermal
#' component alone reaches \eqn{10^{8}}), so a raw call declares the
#' stationary point of the badly scaled problem converged at
#' \eqn{\delta = 0}.  Every commercial and open-source QP backend
#' (SCS/CLARABEL/OSQP/HiGHS included) prescales its problem for the same
#' reason, so the substitution \eqn{\delta = D y}{delta = D y} with
#' \eqn{D = diag(1/\sqrt{H_{ii}})}{D = diag(1/sqrt(H_ii))} is applied
#' first; it is an exact change of variables, so the minimiser returned is
#' the minimiser of the original QP, not of an approximation.
#'
#' @param H Symmetric PSD matrix.
#' @param cvec Gradient vector.
#' @param lo,hi Lower/upper bounds in the original variables.
#' @param x0 Optional feasible starting point in the original variables.
#' @return List \code{list(x, iterations, converged, kkt)}, with \code{x} and
#'   \code{kkt} expressed in the original variables.
#' @keywords internal
#' @noRd
.parametric_box_qp <- function(H, cvec, lo, hi, x0 = NULL) {
    dgl <- as.numeric(diag(H))
    s <- ifelse(is.finite(dgl) & dgl > 0, 1 / sqrt(pmax(dgl, .Machine$double.xmin)),
                1)
    Hs <- as.matrix(H * outer(s, s))
    cs <- as.numeric(cvec * s)
    scale <- max(1, max(abs(Hs)), max(abs(cs)))
    if (!is.finite(scale) || scale <= 0) scale <- 1
    sol <- .commercial_box_qp(Hs / scale, cs / scale, lo / s, hi / s,
                              x0 = if (is.null(x0)) NULL else as.numeric(x0) / s)
    x <- as.numeric(sol$x) * s
    grad <- as.numeric(H %*% x) + cvec
    at_lo <- is.finite(lo) & (x <= lo)
    at_hi <- is.finite(hi) & (x >= hi)
    kkt <- max(c(0, abs(grad[!(at_lo | at_hi)])),
               c(0, if (any(at_lo)) -grad[at_lo] else NULL),
               c(0, if (any(at_hi)) grad[at_hi] else NULL))
    list(x = x, iterations = sol$iterations, converged = sol$converged,
         kkt = kkt)
}

#' One SQP pass shared by the cvxpy and the qpsolvers optimizers
#' @keywords internal
#' @noRd
.parametric_solve_sqp <- function(A, b, E, log_steps, initial_params = NULL,
                                  alpha = 1e-4, max_iter = 50L, tol = 1e-6,
                                  method_name = "parametric_cvxpy") {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E)
    log_steps <- as.numeric(log_steps)
    alpha <- as.numeric(alpha)

    # Python: _find_initial_params(...) with its defaults n_grid=5, return_top=1
    params <- .fruit_find_initial_params(A, b, E, log_steps,
                                         n_grid = 5L, return_top = 1L)
    if (!is.null(initial_params)) {
        ip <- as.list(initial_params)
        for (nm in names(ip)) {
            if (nm %in% names(params)) params[[nm]] <- as.numeric(ip[[nm]])
        }
    }
    bounds <- .fruit_get_param_bounds()
    params <- .fruit_clamp_params(params, bounds)

    names_ <- .fruit_PARAM_NAMES
    n_params <- length(names_)
    lo <- vapply(names_, function(nm) bounds[[nm]][1L], numeric(1L))
    hi <- vapply(names_, function(nm) bounds[[nm]][2L], numeric(1L))
    pvec <- as.numeric(params[names_])

    message <- ""
    nfev <- 0L
    converged <- FALSE
    spectrum <- NULL

    for (k in seq_len(as.integer(max_iter))) {
        spectrum_k <- .fruit_model_vec(E, params) * log_steps
        residual <- as.numeric(A %*% spectrum_k) - b
        nfev <- nfev + 1L
        res_norm <- sqrt(sum(residual^2))
        if (res_norm < tol) {
            .fruit_check_fit_quality(res_norm, b, method_name)
            return(list(spectrum = as.numeric(spectrum_k),
                        iterations = nfev, converged = TRUE,
                        params = params,
                        message = sprintf("Converged in %d iterations",
                                          k - 1L)))
        }

        J <- .fruit_compute_jacobian(E, log_steps, params)
        A_eff <- as.matrix(A %*% J)
        H <- crossprod(A_eff) + alpha * diag(n_params)
        storage.mode(H) <- "double"
        cvec <- as.vector(crossprod(A_eff, residual))
        sol <- .parametric_box_qp(as.matrix(H), cvec, lo - pvec, hi - pvec)
        delta <- as.numeric(sol$x)
        if (!all(is.finite(delta))) {
            message <- sprintf("QP subproblem failed at iteration %d", k - 1L)
            break
        }

        pvec <- pmax(lo, pmin(hi, pvec + delta))
        names(pvec) <- names_
        params <- .fruit_clamp_params(as.list(pvec), bounds)

        if (sqrt(sum(delta^2)) < tol) {
            spectrum <- .fruit_model_vec(E, params) * log_steps
            converged <- TRUE
            message <- sprintf("Converged in %d iterations", k + 1L)
            return(list(spectrum = as.numeric(spectrum),
                        iterations = nfev, converged = TRUE,
                        params = params, message = message))
        }
    }

    spectrum <- .fruit_model_vec(E, params) * log_steps
    .fruit_check_fit_quality(sqrt(sum((as.numeric(A %*% spectrum) - b)^2)),
                             b, method_name)
    if (message == "") {
        message <- sprintf("Max iterations (%d) reached",
                           as.integer(max_iter))
    }
    list(spectrum = as.numeric(spectrum), iterations = nfev,
         converged = converged, params = params, message = message)
}

#' Parametric FRUIT unfolding by sequential QP (cvxpy optimizer)
#'
#' Port of \code{solve_parametric_cvxpy}.  The linearised QP is solved by
#' the package's internal active-set box-QP engine, which returns the same
#' global minimiser cvxpy obtains from ECOS/SCS/Clarabel but exactly and
#' with a KKT certificate.
#'
#' @inheritParams solve_parametric
#' @param log_steps Numeric \eqn{d(\ln E)} bin widths (length n).
#' @param solver_backend Character; backend label, kept for parity (the
#'   pure-R engine does not dispatch to an external solver).
#' @param max_iter Integer; maximum SQP iterations.
#' @param tol Numeric; SQP convergence tolerance.
#' @return A list \code{list(spectrum, iterations, converged, params,
#'   message)}.
#' @keywords internal
#' @noRd
solve_parametric_cvxpy <- function(A, b, E, log_steps,
                                   initial_params = NULL,
                                   method = "leastsq", alpha = 1e-4,
                                   solver_backend = "auto",
                                   max_iter = 50L, tol = 1e-6) {
    .parametric_solve_sqp(A, b, E, log_steps, initial_params = initial_params,
                          alpha = alpha, max_iter = max_iter, tol = tol,
                          method_name = "parametric_cvxpy")
}

#' Parametric FRUIT unfolding by sequential QP (qpsolvers optimizer)
#'
#' Port of \code{solve_parametric_qpsolvers}; identical QP formulation to
#' \code{\link{solve_parametric_cvxpy}} (see the note above the shared
#' engine), solved by the internal active-set box-QP solver.
#' @inheritParams solve_parametric_cvxpy
#' @keywords internal
#' @noRd
solve_parametric_qpsolvers <- function(A, b, E, log_steps,
                                       initial_params = NULL,
                                       method = "leastsq", alpha = 1e-4,
                                       solver_backend = "auto",
                                       max_iter = 50L, tol = 1e-6) {
    .parametric_solve_sqp(A, b, E, log_steps, initial_params = initial_params,
                          alpha = alpha, max_iter = max_iter, tol = tol,
                          method_name = "parametric_qpsolvers")
}

#' Parse a Python \code{solver_backend} string into its library name
#' (port of \code{_parse_solver_backend})
#' @keywords internal
#' @noRd
.parametric_solver_library <- function(solver_backend = "auto") {
    if (is.null(solver_backend) || !nzchar(solver_backend)) return("auto")
    parts <- strsplit(solver_backend, ":", fixed = TRUE)[[1L]]
    lib <- parts[1L]
    if (!lib %in% c("auto", "cvxpy", "qpsolvers")) {
        stop("Unknown solver library: '", lib,
             "'. Use 'cvxpy' or 'qpsolvers'.", call. = FALSE)
    }
    lib
}

#' FRUIT lmfit fit followed by a non-negative QP refinement
#'
#' Port of \code{solve_parametric_combined}.  Step 1 is the lmfit
#' (Levenberg-Marquardt) fit \code{\link{solve_parametric}}; step 2 solves
#'
#' \deqn{\min_{x \ge 0} \|Ax - b\|^2 + \alpha\|x - x_{\text{lmfit}}\|^2}
#'
#' which both the cvxpy and the qpsolvers refinement branch of the Python
#' original express as the same box QP
#' \eqn{0.5 x'(A'A + \alpha I)x - (A'b + \alpha x_{\text{lmfit}})'x}{0.5 x'(A'A + alpha*I)x - (A'b + alpha*x_lmfit)'x}.
#' Like Python, the refined lethargy-density vector is multiplied by the
#' log steps once more before it is returned.
#' @inheritParams solve_parametric_cvxpy
#' @keywords internal
#' @noRd
solve_parametric_combined <- function(A, b, E, log_steps,
                                      initial_params = NULL,
                                      method = "leastsq", alpha = 1e-4,
                                      solver_backend = "auto",
                                      max_iter = 50L, tol = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); log_steps <- as.numeric(log_steps)
    lib <- .parametric_solver_library(solver_backend)
    if (lib == "auto") lib <- "cvxpy"

    fit <- solve_parametric(A, b, E, log_steps,
                            initial_params = initial_params, method = method)
    spectrum_lmfit <- as.numeric(fit$spectrum)
    .fruit_check_fit_quality(
        sqrt(sum((as.numeric(A %*% spectrum_lmfit) - b)^2)), b,
        "parametric_combined(lmfit)")

    n <- ncol(A)
    qp <- .commercial_qp_data(A, b, alpha, norm = 2L, smoothness_order = 0L,
                              smoothness_weight = 1)
    cvec <- qp$q - alpha * spectrum_lmfit
    sol <- .parametric_box_qp(qp$P, cvec, rep(0, n), rep(Inf, n),
                              x0 = spectrum_lmfit)
    if (!all(is.finite(sol$x))) {
        return(list(spectrum = spectrum_lmfit, iterations = fit$iterations,
                    converged = isTRUE(fit$converged),
                    params = fit$params, message = "QP refinement failed"))
    }
    list(spectrum = as.numeric(pmax(sol$x, 0) * log_steps),
         iterations = fit$iterations,
         converged = isTRUE(fit$converged),
         params = fit$params,
         message = sprintf("lmfit (%s) + QP refinement OK", fit$message))
}


# ------------------------------------------------------------------ #
#  MINPACK LMDIF (Levenberg-Marquardt) -- pure-R port
#
#  Python's parametric family fits its models with
#  ``lmfit.minimize(..., method = "leastsq")``, which is a thin wrapper
#  around ``scipy.optimize.leastsq`` and therefore around the NETLIB
#  double-precision ``lmdif``/``fdjac2``/``qrfac``/``qrsolv``/``lmpar``
#  sources.  There is no CRAN equivalent that may be used here
#  (``minpack.lm`` is not available), so the algorithm is ported
#  statement by statement: the parametric models have several
#  near-degenerate local minima whose residual norms differ only in the
#  7th digit, so only the exact algorithm reproduces which minimum the
#  original reports.
# ------------------------------------------------------------------ #

#' ENORM: scaled Euclidean norm (port of the MINPACK ``enorm``)
#' @keywords internal
#' @noRd
.mp_enorm <- function(x) {
    x <- as.numeric(x)
    n <- length(x)
    if (n == 0L) return(0)
    one <- 1; zero <- 0
    rdwarf <- 3.834e-20; rgiant <- 1.304e19
    s1 <- zero; s2 <- zero; s3 <- zero
    x1max <- zero; x3max <- zero
    agiant <- rgiant / as.numeric(n)
    for (xa in abs(x)) {
        if (xa > rdwarf && xa < agiant) {
            s2 <- s2 + xa^2
        } else if (xa <= rdwarf) {
            if (xa > x3max) {
                s3 <- one + s3 * (x3max / xa)^2
                x3max <- xa
            } else if (xa != zero) {
                s3 <- s3 + (xa / x3max)^2
            }
        } else {
            if (xa > x1max) {
                s1 <- one + s1 * (x1max / xa)^2
                x1max <- xa
            } else {
                s1 <- s1 + (xa / x1max)^2
            }
        }
    }
    if (s1 != zero) return(x1max * sqrt(s1 + (s2 / x1max) / x1max))
    if (s2 == zero) return(x3max * sqrt(s3))
    if (s2 >= x3max) return(sqrt(s2 * (one + (x3max / s2) * (x3max * s3))))
    sqrt(x3max * ((s2 / x3max) + (x3max * s3)))
}

#' FDJAC2: forward-difference Jacobian (port of the MINPACK ``fdjac2``)
#' @keywords internal
#' @noRd
.mp_fdjac2 <- function(fcn, x, fvec, epsfcn, epsmch) {
    x <- as.numeric(x)
    m <- length(fvec); n <- length(x)
    eps <- sqrt(max(epsfcn, epsmch))
    fjac <- matrix(0, nrow = m, ncol = n)
    for (j in seq_len(n)) {
        temp <- x[j]
        h <- eps * abs(temp)
        if (h == 0) h <- eps
        x[j] <- temp + h
        wa <- fcn(x)
        x[j] <- temp
        if (is.null(wa)) return(list(fjac = fjac, iflag = -1L))
        fjac[, j] <- (as.numeric(wa) - fvec) / h
    }
    list(fjac = fjac, iflag = 1L)
}

#' QRFAC: Householder QR with optional column pivoting
#' (port of the MINPACK ``qrfac``)
#' @keywords internal
#' @noRd
.mp_qrfac <- function(a, pivot = TRUE) {
    a <- as.matrix(a); storage.mode(a) <- "double"
    m <- nrow(a); n <- ncol(a)
    epsmch <- .Machine$double.eps
    acnorm <- vapply(seq_len(n), function(j) .mp_enorm(a[, j]), numeric(1L))
    rdiag <- acnorm; wa <- acnorm
    ipvt <- seq_len(n)
    minmn <- min(m, n)
    for (j in seq_len(minmn)) {
        if (pivot) {
            kmax <- j
            if (n > j) for (k in (j + 1L):n) if (rdiag[k] > rdiag[kmax]) kmax <- k
            if (kmax != j) {
                tmp <- a[, j]; a[, j] <- a[, kmax]; a[, kmax] <- tmp
                rdiag[kmax] <- rdiag[j]
                wa[kmax] <- wa[j]
                kk <- ipvt[j]; ipvt[j] <- ipvt[kmax]; ipvt[kmax] <- kk
            }
        }
        ajnorm <- .mp_enorm(a[j:m, j])
        if (ajnorm != 0) {
            if (a[j, j] < 0) ajnorm <- -ajnorm
            a[j:m, j] <- a[j:m, j] / ajnorm
            a[j, j] <- a[j, j] + 1
            if (n >= j + 1L) for (k in (j + 1L):n) {
                s <- sum(a[j:m, j] * a[j:m, k])
                temp <- s / a[j, j]
                a[j:m, k] <- a[j:m, k] - temp * a[j:m, j]
                if (pivot && rdiag[k] != 0) {
                    temp <- a[j, k] / rdiag[k]
                    rdiag[k] <- rdiag[k] * sqrt(max(0, 1 - temp^2))
                    if (0.05 * (rdiag[k] / wa[k])^2 > epsmch) next
                    rdiag[k] <- .mp_enorm(a[(j + 1L):m, k])
                    wa[k] <- rdiag[k]
                }
            }
        }
        rdiag[j] <- -ajnorm
    }
    list(a = a, acnorm = acnorm, rdiag = rdiag, wa = wa, ipvt = ipvt)
}

#' QRSOLV: subspace solution of (R'R + D'D) x = Q'b by Givens rotations
#' (port of the MINPACK ``qrsolv``); ``r`` holds R in its upper triangle
#' @keywords internal
#' @noRd
.mp_qrsolv <- function(r, ipvt, dvec, qtb, n) {
    x <- numeric(n); wa <- numeric(n); sdiag <- numeric(n)
    for (j in seq_len(n)) {
        r[j:n, j] <- r[j, j:n]
        x[j] <- r[j, j]
        wa[j] <- qtb[j]
    }
    for (j in seq_len(n)) {
        l <- ipvt[j]
        if (dvec[l] != 0) {
            sdiag[j:n] <- 0
            sdiag[j] <- dvec[l]
            qtbpj <- 0
            for (k in j:n) {
                if (sdiag[k] != 0) {
                    if (abs(r[k, k]) >= abs(sdiag[k])) {
                        tan <- sdiag[k] / r[k, k]
                        cos <- 0.5 / sqrt(0.25 + 0.25 * tan^2)
                        sin <- cos * tan
                    } else {
                        cotan <- r[k, k] / sdiag[k]
                        sin <- 0.5 / sqrt(0.25 + 0.25 * cotan^2)
                        cos <- sin * cotan
                    }
                    r[k, k] <- cos * r[k, k] + sin * sdiag[k]
                    temp <- cos * wa[k] + sin * qtbpj
                    qtbpj <- -sin * wa[k] + cos * qtbpj
                    wa[k] <- temp
                    if (n >= k + 1L) for (i in (k + 1L):n) {
                        temp <- cos * r[i, k] + sin * sdiag[i]
                        sdiag[i] <- -sin * r[i, k] + cos * sdiag[i]
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
    if (nsing >= 1L) for (k in seq_len(nsing)) {
        j <- nsing - k + 1L
        s <- 0
        if (nsing >= j + 1L) for (i in (j + 1L):nsing) s <- s + r[i, j] * wa[i]
        wa[j] <- (wa[j] - s) / sdiag[j]
    }
    for (j in seq_len(n)) x[ipvt[j]] <- wa[j]
    list(r = r, x = x, sdiag = sdiag)
}

#' LMPAR: Levenberg-Marquardt parameter and direction
#' (port of the MINPACK ``lmpar``)
#' @keywords internal
#' @noRd
.mp_lmpar <- function(r, ipvt, diag, qtb, delta, par, n) {
    p1 <- 1e-1; p001 <- 1e-3; zero <- 0
    dwarf <- .Machine$double.xmin
    x <- numeric(n); sdiag <- numeric(n); wa1 <- numeric(n); wa2 <- numeric(n)
    nsing <- n
    for (j in seq_len(n)) {
        wa1[j] <- qtb[j]
        if (r[j, j] == 0 && nsing == n) nsing <- j - 1L
        if (nsing < n) wa1[j] <- 0
    }
    if (nsing >= 1L) for (k in seq_len(nsing)) {
        j <- nsing - k + 1L
        wa1[j] <- wa1[j] / r[j, j]
        temp <- wa1[j]
        if (j > 1L) for (i in seq_len(j - 1L)) wa1[i] <- wa1[i] - r[i, j] * temp
    }
    for (j in seq_len(n)) x[ipvt[j]] <- wa1[j]
    iter <- 0L
    for (j in seq_len(n)) wa2[j] <- diag[j] * x[j]
    dxnorm <- .mp_enorm(wa2)
    fp <- dxnorm - delta
    if (fp > p1 * delta) {
        parl <- zero
        if (nsing >= n) {
            for (j in seq_len(n)) {
                l <- ipvt[j]
                wa1[j] <- diag[l] * (wa2[l] / dxnorm)
            }
            for (j in seq_len(n)) {
                s <- zero
                if (j > 1L) for (i in seq_len(j - 1L)) s <- s + r[i, j] * wa1[i]
                wa1[j] <- (wa1[j] - s) / r[j, j]
            }
            temp <- .mp_enorm(wa1)
            parl <- ((fp / delta) / temp) / temp
        }
        for (j in seq_len(n)) {
            s <- zero
            for (i in seq_len(j)) s <- s + r[i, j] * qtb[i]
            wa1[j] <- s / diag[ipvt[j]]
        }
        gnorm <- .mp_enorm(wa1)
        paru <- gnorm / delta
        if (paru == 0) paru <- dwarf / min(delta, p1)
        par <- max(par, parl)
        par <- min(par, paru)
        if (par == 0) par <- gnorm / dxnorm
        repeat {
            iter <- iter + 1L
            if (par == 0) par <- max(dwarf, p001 * paru)
            temp <- sqrt(par)
            for (j in seq_len(n)) wa1[j] <- temp * diag[j]
            qs <- .mp_qrsolv(r, ipvt, wa1, qtb, n)
            r <- qs$r; x <- qs$x; sdiag <- qs$sdiag
            for (j in seq_len(n)) wa2[j] <- diag[j] * x[j]
            dxnorm <- .mp_enorm(wa2)
            temp <- fp
            fp <- dxnorm - delta
            if (abs(fp) <= p1 * delta ||
                (parl == 0 && fp <= temp && temp < 0) || iter == 10L) break
            for (j in seq_len(n)) {
                l <- ipvt[j]
                wa1[j] <- diag[l] * (wa2[l] / dxnorm)
            }
            for (j in seq_len(n)) {
                wa1[j] <- wa1[j] / sdiag[j]
                temp <- wa1[j]
                if (n >= j + 1L) for (i in (j + 1L):n) wa1[i] <- wa1[i] - r[i, j] * temp
            }
            temp <- .mp_enorm(wa1)
            parc <- ((fp / delta) / temp) / temp
            if (fp > 0) parl <- max(parl, par)
            if (fp < 0) paru <- min(paru, par)
            par <- max(parl, par + parc)
        }
    }
    if (iter == 0L) par <- zero
    list(r = r, x = x, sdiag = sdiag, par = par)
}

#' LMDIF: Levenberg-Marquardt least-squares minimiser
#' (port of the MINPACK ``lmdif`` with forward-difference Jacobian)
#'
#' @param fcn Function mapping the parameter vector to the residual vector;
#'   returning \code{NULL} aborts the fit (Fortran \code{iflag < 0}).
#' @param x0 Numeric starting point.
#' @param ftol,xtol,gtol Numeric convergence tolerances.
#' @param maxfev Integer; maximum number of function evaluations.
#' @param epsfcn Numeric; relative error assumed for the finite-difference
#'   Jacobian (the step is \code{sqrt(max(epsfcn, epsmch))}).
#' @param factor Numeric; initial step bound multiplier.
#' @param mode Integer; 1 = internal column-norm scaling, 2 = user \code{diag}.
#' @param diag Numeric scaling vector (only used when \code{mode = 2}).
#' @return List \code{list(x, fvec, info, nfev, xnorm, fnorm)}.
#' @keywords internal
#' @noRd
.mp_lmdif <- function(fcn, x0, ftol = 1.5e-8, xtol = 1.5e-8, gtol = 0,
                      maxfev = 4000L, epsfcn = 1e-10, factor = 100,
                      mode = 1L, diag = NULL) {
    x <- as.numeric(x0)
    n <- length(x)
    one <- 1; p1 <- 1e-1; p5 <- 5e-1; p25 <- 2.5e-1; p75 <- 7.5e-1
    p0001 <- 1e-4; zero <- 0
    epsmch <- .Machine$double.eps
    info <- 0L; nfev <- 0L
    fvec <- numeric(0)
    if (mode == 2L) {
        if (is.null(diag) || length(diag) != n || any(diag <= 0)) {
            return(list(x = x, fvec = fvec, info = 0L, nfev = 0L))
        }
        diag <- as.numeric(diag)
    } else {
        mode <- 1L
        diag <- rep(one, n)
    }
    fvec <- fcn(x)
    nfev <- 1L
    if (is.null(fvec) || !all(is.finite(fvec))) {
        return(list(x = x, fvec = fvec, info = -1L, nfev = nfev))
    }
    fvec <- as.numeric(fvec)
    m <- length(fvec)
    if (n <= 0L || m < n || ftol < zero || xtol < zero || gtol < zero ||
        maxfev <= 0L || factor <= zero) {
        return(list(x = x, fvec = fvec, info = 0L, nfev = nfev))
    }
    fnorm <- .mp_enorm(fvec)
    par <- zero
    iter <- 1L
    qtf <- numeric(n); wa1 <- numeric(n); wa2 <- numeric(n)
    wa3 <- numeric(n); wa4 <- numeric(m); ipvt <- seq_len(n)
    fjac <- matrix(0, nrow = m, ncol = n)
    xnorm <- 0; delta <- 0; gnorm <- 0
    repeat {
        fj <- .mp_fdjac2(fcn, x, fvec, epsfcn, epsmch)
        nfev <- nfev + n
        if (fj$iflag < 0L) { info <- -1L; break }
        fjac <- fj$fjac
        qr <- .mp_qrfac(fjac, TRUE)
        # Fortran: qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3) passes
        # rdiag -> wa1 and acnorm -> wa2 (the order of its own arguments).
        fjac <- qr$a; wa1 <- qr$rdiag; wa2 <- qr$acnorm; wa3 <- qr$wa
        ipvt <- qr$ipvt
        if (iter == 1L) {
            if (mode != 2L) for (j in seq_len(n)) {
                diag[j] <- wa2[j]
                if (wa2[j] == 0) diag[j] <- one
            }
            for (j in seq_len(n)) wa3[j] <- diag[j] * x[j]
            xnorm <- .mp_enorm(wa3)
            delta <- factor * xnorm
            if (delta == 0) delta <- factor
        }
        wa4 <- fvec
        for (j in seq_len(n)) {
            if (fjac[j, j] != 0) {
                s <- sum(fjac[j:m, j] * wa4[j:m])
                temp <- -s / fjac[j, j]
                wa4[j:m] <- wa4[j:m] + fjac[j:m, j] * temp
            }
            fjac[j, j] <- wa1[j]
            qtf[j] <- wa4[j]
        }
        gnorm <- zero
        if (fnorm != 0) for (j in seq_len(n)) {
            l <- ipvt[j]
            if (wa2[l] != 0) {
                s <- zero
                for (i in seq_len(j)) s <- s + fjac[i, j] * (qtf[i] / fnorm)
                gnorm <- max(gnorm, abs(s / wa2[l]))
            }
        }
        if (gnorm <= gtol) info <- 4L
        if (info != 0L) break
        if (mode != 2L) for (j in seq_len(n)) diag[j] <- max(diag[j], wa2[j])
        ## inner loop
        repeat {
            lm <- .mp_lmpar(fjac, ipvt, diag, qtf, delta, par, n)
            fjac <- lm$r; par <- lm$par
            for (j in seq_len(n)) {
                wa1[j] <- -lm$x[j]
                wa2[j] <- x[j] + wa1[j]
                wa3[j] <- diag[j] * wa1[j]
            }
            pnorm <- .mp_enorm(wa3)
            if (iter == 1L) delta <- min(delta, pnorm)
            wa4 <- fcn(wa2)
            nfev <- nfev + 1L
            if (is.null(wa4) || !all(is.finite(wa4))) { info <- -1L; break }
            wa4 <- as.numeric(wa4)
            fnorm1 <- .mp_enorm(wa4)
            actred <- -one
            if (p1 * fnorm1 < fnorm) actred <- one - (fnorm1 / fnorm)^2
            for (j in seq_len(n)) {
                wa3[j] <- zero
                temp <- wa1[ipvt[j]]
                for (i in seq_len(j)) wa3[i] <- wa3[i] + fjac[i, j] * temp
            }
            temp1 <- .mp_enorm(wa3) / fnorm
            temp2 <- (sqrt(par) * pnorm) / fnorm
            prered <- temp1^2 + temp2^2 / p5
            dirder <- -(temp1^2 + temp2^2)
            ratio <- zero
            if (prered != 0) ratio <- actred / prered
            if (ratio > p25) {
                if (par == 0 || ratio >= p75) {
                    delta <- pnorm / p5
                    par <- p5 * par
                }
            } else {
                if (actred >= 0) temp <- p5 else temp <- p5 * dirder / (dirder + p5 * actred)
                if (p1 * fnorm1 >= fnorm || temp < p1) temp <- p1
                delta <- temp * min(delta, pnorm / p1)
                par <- par / temp
            }
            if (ratio >= p0001) {
                for (j in seq_len(n)) {
                    x[j] <- wa2[j]
                    wa2[j] <- diag[j] * x[j]
                }
                fvec <- wa4
                xnorm <- .mp_enorm(wa2)
                fnorm <- fnorm1
                iter <- iter + 1L
            }
            info <- if (abs(actred) <= ftol && prered <= ftol &&
                        p5 * ratio <= one) 1L else 0L
            if (delta <= xtol * xnorm) info <- 2L
            if (abs(actred) <= ftol && prered <= ftol && p5 * ratio <= one &&
                info == 2L) info <- 3L
            if (info != 0L) break
            if (nfev >= maxfev) info <- 5L
            if (abs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= one) info <- 6L
            if (delta <= epsmch * xnorm) info <- 7L
            if (gnorm <= epsmch) info <- 8L
            if (info != 0L) break
            if (ratio >= p0001) break
        }
        if (info != 0L) break
    }
    list(x = x, fvec = fvec, info = info, nfev = nfev, iter = iter,
         xnorm = xnorm, fnorm = fnorm)
}

#' lmfit (Minuit-style) internal/external parameter transformation
#'
#' \code{lmfit} fits bounded parameters with the unbounded MINPACK
#' \code{leastsq} by mapping them through \code{sin}/\code{sqrt}
#' transformations (see \code{lmfit/parameter.py::setup_bounds}); the
#' residual function is evaluated at the transformed values.
#' @keywords internal
#' @noRd
.lmfit_setup_bounds <- function(values, lower, upper) {
    vnames <- names(values)
    if (is.null(vnames)) vnames <- rep("", length(values))
    values <- as.numeric(values); lower <- as.numeric(lower); upper <- as.numeric(upper)
    lower[is.na(lower)] <- -Inf
    upper[is.na(upper)] <- Inf
    type <- ifelse(is.infinite(lower) & is.infinite(upper), 0L,
            ifelse(is.infinite(upper), 1L,
            ifelse(is.infinite(lower), 2L, 3L)))
    internal <- numeric(length(values))
    for (i in seq_along(values)) {
        lo <- lower[i]; hi <- upper[i]
        val <- values[i]
        # lmfit does not clamp, but an out-of-bounds start would give NaN in
        # asin(); clamp defensively (a no-op for in-range starts).
        if (is.finite(lo) && val < lo) val <- lo
        if (is.finite(hi) && val > hi) val <- hi
        internal[i] <- switch(type[i] + 1L,
            val,
            sqrt((val - lo + 1)^2 - 1),
            sqrt((hi - val + 1)^2 - 1),
            asin(2 * (val - lo) / (hi - lo) - 1))
        if (abs(internal[i]) < .Machine$double.xmin) internal[i] <- 0
    }
    list(type = type, lower = lower, upper = upper, internal = internal,
         vname = vnames)
}

#' Inverse of \code{.lmfit_setup_bounds}: external values from internal ones
#' @keywords internal
#' @noRd
.lmfit_from_internal <- function(blk, internal) {
    internal <- as.numeric(internal)
    out <- vapply(seq_along(internal), function(i) {
        lo <- blk$lower[i]; hi <- blk$upper[i]
        val <- internal[i]
        switch(blk$type[i] + 1L,
            val,
            lo - 1 + sqrt(val * val + 1),
            hi + 1 - sqrt(val * val + 1),
            lo + (sin(val) + 1) * (hi - lo) / 2)
    }, numeric(1L))
    if (!is.null(blk$vname)) names(out) <- blk$vname
    out
}

#' Run one bounded least-squares fit exactly like
#' \code{lmfit.minimize(..., method = "leastsq")}
#'
#' @param resid Function of the external parameter vector returning the
#'   residual vector.
#' @param values Named numeric starting values (order = fit variables).
#' @param lower,upper Numeric bounds (may contain \code{Inf}/\code{-Inf}).
#' @return List \code{list(values, info, nfev, success, message)}.
#' @keywords internal
#' @noRd
.lmfit_leastsq <- function(resid, values, lower, upper) {
    blk <- .lmfit_setup_bounds(values, lower, upper)
    fcn <- function(xi) {
        r <- resid(.lmfit_from_internal(blk, xi))
        if (is.null(r)) return(NULL)
        as.numeric(r)
    }
    n <- length(values)
    out <- .mp_lmdif(fcn, blk$internal, ftol = 1.5e-8, xtol = 1.5e-8,
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
    # lmfit reports nfev corrected for the pre-fit initialisation checks
    nfev <- max(out$nfev - 3L, 0L)
    list(values = ext, info = info, nfev = nfev, success = success,
         message = message)
}

# Fixed constants from the FRUIT papers
.fruit_T0 <- 2.53e-8        # thermal peak energy (MeV)
.fruit_Ed <- 7.07e-8        # epithermal lower boundary parameter (MeV)
.fruit_THERMAL_MAX <- 1e-7  # MeV
.fruit_FAST_MIN <- 0.1      # MeV

.fruit_PARAM_NAMES <- c("b", "beta_prime", "alpha", "beta", "P_th", "P_epi")

# name -> c(default, lower, upper)
.fruit_PARAM_DEFAULTS <- list(
    b = c(1.0, 0.5, 2.0),
    beta_prime = c(0.01, 1e-4, 1.0),
    alpha = c(0.5, 0.0, 5.0),
    beta = c(2.0, 0.1, 20.0),
    P_th = c(1.0, 0.0, 1.0),
    P_epi = c(1.0, 0.0, 1.0)
)

#' Thermal / epithermal / fast components of the FRUIT model
#' @keywords internal
#' @noRd
.fruit_thermal <- function(E) (E / (.fruit_T0^2)) * exp(-E / .fruit_T0)

.fruit_epithermal <- function(E, b, beta_prime) {
    (1 - exp(-(E / .fruit_Ed)^2)) * E^(b - 1) * exp(-E / beta_prime)
}

.fruit_fast <- function(E, alpha, beta) E^alpha * exp(-E / beta)

#' FRUIT three-component parametric spectrum (port of parametric_model)
#' @keywords internal
#' @noRd
.fruit_model <- function(E, b, beta_prime, alpha, beta, P_th, P_epi) {
    E <- as.numeric(E)
    P_f <- max(0, 1 - P_th - P_epi)
    thermal <- numeric(length(E))
    epithermal <- numeric(length(E))
    fast <- numeric(length(E))
    m_th <- E < .fruit_THERMAL_MAX
    m_epi <- (E >= .fruit_THERMAL_MAX) & (E < .fruit_FAST_MIN)
    m_f <- E >= .fruit_FAST_MIN
    if (any(m_th)) thermal[m_th] <- .fruit_thermal(E[m_th])
    if (any(m_epi)) epithermal[m_epi] <- .fruit_epithermal(E[m_epi], b,
                                                           beta_prime)
    if (any(m_f)) fast[m_f] <- .fruit_fast(E[m_f], alpha, beta)
    P_th * thermal + P_epi * epithermal + P_f * fast
}

#' Same model evaluated from a named parameter vector/list
#' @keywords internal
#' @noRd
.fruit_model_vec <- function(E, params) {
    .fruit_model(E, params[["b"]], params[["beta_prime"]], params[["alpha"]],
                 params[["beta"]], params[["P_th"]], params[["P_epi"]])
}

#' Default parameter values with user overrides applied
#' @keywords internal
#' @noRd
.fruit_get_initial_params <- function(initial_params = NULL) {
    params <- vapply(.fruit_PARAM_DEFAULTS, function(x) x[1L], numeric(1L))
    if (!is.null(initial_params)) {
        ip <- as.list(initial_params)
        for (nm in names(params)) {
            if (nm %in% names(ip)) params[[nm]] <- as.numeric(ip[[nm]])
        }
    }
    params
}

#' Parameter bounds as list(name = c(lower, upper))
#' @keywords internal
#' @noRd
.fruit_get_param_bounds <- function() {
    lapply(.fruit_PARAM_DEFAULTS, function(x) c(x[2L], x[3L]))
}

#' Clamp parameters into their bounds
#' @keywords internal
#' @noRd
.fruit_clamp_params <- function(params, bounds = .fruit_get_param_bounds()) {
    for (nm in names(bounds)) {
        if (!is.null(params[[nm]])) {
            params[[nm]] <- max(bounds[[nm]][1L],
                                min(bounds[[nm]][2L], params[[nm]]))
        }
    }
    params
}

#' Forward finite-difference Jacobian of (model * log_steps) w.r.t. the
#' parameters, with bound-aware backward fallback (port of _compute_jacobian)
#' @keywords internal
#' @noRd
.fruit_compute_jacobian <- function(E, log_steps, params, delta = 1e-8) {
    bounds <- .fruit_get_param_bounds()
    n_params <- length(.fruit_PARAM_NAMES)
    J <- matrix(0, nrow = length(E), ncol = n_params)
    s0 <- .fruit_model_vec(E, params) * log_steps

    for (i in seq_len(n_params)) {
        nm <- .fruit_PARAM_NAMES[i]
        lo <- bounds[[nm]][1L]
        hi <- bounds[[nm]][2L]
        p_val <- params[[nm]]
        d <- delta
        if (!is.na(hi) && p_val + d > hi) d <- max(0, hi - p_val) * 0.5
        if (!is.na(lo) && p_val + d < lo) d <- 0
        if (d < 1e-15) {
            d <- delta
            if (!is.na(lo) && p_val - d >= lo) {
                pert <- params
                pert[[nm]] <- p_val - d
                s_pert <- .fruit_model_vec(E, pert) * log_steps
                J[, i] <- (s0 - s_pert) / d
            } else {
                J[, i] <- 0
            }
            next
        }
        pert <- params
        pert[[nm]] <- p_val + d
        s_plus <- .fruit_model_vec(E, pert) * log_steps
        J[, i] <- (s_plus - s0) / d
    }
    J
}

#' Brute-force scan over (P_th, P_epi): best (or top N) starting points
#' @keywords internal
#' @noRd
.fruit_find_initial_params <- function(A, b, E, log_steps, n_grid = 5L,
                                       return_top = 1L) {
    n_grid <- as.integer(n_grid)
    res_norms <- numeric(0L)
    cand_params <- list()
    p_th_vals <- seq(0, 1, length.out = n_grid)
    p_epi_vals <- seq(0, 1, length.out = n_grid)
    for (p_th in p_th_vals) {
        for (p_epi in p_epi_vals) {
            if (p_th + p_epi > 1) next
            params <- .fruit_get_initial_params(NULL)
            params[["P_th"]] <- p_th
            params[["P_epi"]] <- p_epi
            spectrum <- .fruit_model_vec(E, params) * log_steps
            residual <- as.numeric(A %*% spectrum) - b
            res_norms <- c(res_norms, sqrt(sum(residual^2)))
            cand_params[[length(cand_params) + 1L]] <- params
        }
    }
    if (length(cand_params) == 0L) {
        base <- .fruit_get_initial_params(NULL)
        return(if (return_top == 1L) base else list(base))
    }
    # Stable sort by residual, like Python's list.sort
    ord <- order(res_norms, method = "radix")
    if (return_top == 1L) return(cand_params[[ord[1L]]])
    lapply(ord[seq_len(min(return_top, length(ord)))],
           function(k) cand_params[[k]])
}

#' GCV selection of the Tikhonov weight on the linearised model
#' (port of _gcv_select_alpha)
#' @keywords internal
#' @noRd
.fruit_gcv_select_alpha <- function(A, b, E, log_steps, initial_params,
                                    n_coarse = 50L, n_refine = 20L) {
    J <- .fruit_compute_jacobian(E, log_steps, initial_params)
    A_eff <- as.matrix(A) %*% J
    m <- nrow(A_eff); n <- ncol(A_eff)
    if (m < 2L || n < 2L) return(1e-4)
    svd_parts <- compute_svd_components(A_eff)
    U <- svd_parts$U; s_sq <- svd_parts$s_sq
    UTb <- as.vector(crossprod(U, b))
    gcv_value <- function(alpha_) {
        filt <- s_sq / (s_sq + alpha_)
        residual_coeff <- alpha_ / (s_sq + alpha_)
        residual_sq <- sum((residual_coeff * UTb)^2)
        denom <- (m - sum(filt))^2
        if (denom < 1e-30) return(Inf)
        residual_sq / denom
    }
    alphas_coarse <- 10^seq(-8, 2, length.out = as.integer(n_coarse))
    gcv_coarse <- vapply(alphas_coarse, gcv_value, numeric(1L))
    alpha_best <- alphas_coarse[which.min(gcv_coarse)]
    lo <- max(alpha_best / 10, 1e-10)
    hi <- alpha_best * 10
    alphas_refine <- seq(lo, hi, length.out = as.integer(n_refine))
    gcv_refine <- vapply(alphas_refine, gcv_value, numeric(1L))
    alphas_refine[which.min(gcv_refine)]
}

#' Warn on a poor parametric fit (port of _check_fit_quality)
#' @keywords internal
#' @noRd
.fruit_check_fit_quality <- function(residual_norm, b, method_name = "parametric",
                                     model_descr = "3-component parametric model") {
    b_norm <- sqrt(sum(as.numeric(b)^2))
    if (b_norm > 0) {
        relative_residual <- residual_norm / b_norm
        if (relative_residual > 10) {
            warning(sprintf(paste0("%s: large residual (%.2e / %.2e = %.1fx). ",
                                   "The %s may not represent this spectrum well."),
                            method_name, residual_norm, b_norm,
                            relative_residual, model_descr),
                    call. = FALSE)
        }
    }
    invisible(NULL)
}

#' Wrapper around \code{\link{solve_parametric}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_parametric
#' @param optimizer Character; backend selector, as in Python:
#'   \code{"lmfit"} (default, bounded least squares in R), \code{"cvxpy"},
#'   \code{"qpsolvers"} or \code{"combined"} (SQP variants implemented in
#'   \code{R/unfold_parametric_engine.R}).
#' @param solver_backend Character; QP backend string, kept for parity with
#'   Python (\code{"auto"}).
#' @param max_iter Integer; maximum SQP iterations for the QP optimizers.
#' @param tol Numeric; SQP convergence tolerance.
#' @param calculate_errors Logical; run Monte-Carlo uncertainty estimation.
#' @param noise_level Numeric; relative noise for Monte-Carlo.
#' @param n_montecarlo Integer; Monte-Carlo sample count.
#' @param save_result Logical; hand the result to \code{save_result_callback}.
#' @param random_state Optional integer seed for the Monte-Carlo part.
#' @param max_neutron_energy Optional energy cutoff in MeV.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_parametric <- function(detector_names, n_energy_bins, E_MeV,
                                 sensitivities, cc_icrp116, save_result_callback,
                                 readings, initial_spectrum = NULL,
                                 initial_params = NULL,
                                 method = "leastsq", optimizer = "lmfit",
                                 alpha = 1e-4, alpha_auto = FALSE,
                                 solver_backend = "auto",
                                 max_iter = 50L, tol = 1e-6,
                                 max_iterations = 200L, tolerance = 1e-6,
                                 calculate_errors = FALSE,
                                 noise_level = 0.01, n_montecarlo = 100L,
                                 save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b
    # Python: log_steps = _matrix_utils.compute_log_steps(E, n)  (dex)
    # The package-level helper is the exact port of that function (log10 of
    # E + 1e-15, edge differences at the ends, central differences inside), so
    # this method always follows the core (+1e-15) convention.
    log_steps <- compute_log_steps(E_MeV)
    ln_steps <- log_steps * log(10)

    optimizer <- match.arg(tolower(optimizer),
                           c("lmfit", "cvxpy", "qpsolvers", "combined"))
    if (optimizer == "lmfit") {
        # lmfit path: tiny Tikhonov term for numerical stability
        lmfit_alpha <- if (isTRUE(alpha_auto)) alpha else 1e-8
        solver <- make_solve_wrapper(solve_parametric, E = E_MeV,
                                     log_steps = ln_steps,
                                     initial_params = initial_params,
                                     method = method, alpha = lmfit_alpha,
                                     alpha_auto = alpha_auto,
                                     max_iterations = max_iterations,
                                     tolerance = tolerance)
        method_name <- "parametric"
        extra <- list(initial_params = initial_params, lmfit_method = method,
                      alpha_auto = alpha_auto, T0 = .fruit_T0, Ed = .fruit_Ed)
    } else {
        engine_fun <- switch(optimizer,
            cvxpy = solve_parametric_cvxpy,
            qpsolvers = solve_parametric_qpsolvers,
            combined = solve_parametric_combined)
        solver <- function(A, b, x0 = NULL, ...) {
            if (optimizer == "combined") {
                res <- engine_fun(A, b, E_MeV, ln_steps,
                                  initial_params = initial_params,
                                  method = method, alpha = alpha,
                                  solver_backend = solver_backend)
            } else {
                res <- engine_fun(A, b, E_MeV, ln_steps,
                                  initial_params = initial_params,
                                  alpha = alpha,
                                  solver_backend = solver_backend,
                                  max_iter = as.integer(max_iter),
                                  tol = tol)
            }
            list(spectrum = res$spectrum, iterations = res$iterations,
                 converged = res$converged)
        }
        method_name <- paste0("parametric_", optimizer)
        extra <- list(initial_params = initial_params, optimizer = optimizer,
                      alpha = alpha, solver_backend = solver_backend,
                      max_iter = max_iter, tol = tol,
                      T0 = .fruit_T0, Ed = .fruit_Ed)
    }

    x0_default <- rep(mean(b) / max(mean(rowSums(A)), .Machine$double.xmin),
                      n_energy_bins)

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver, solve_kwargs = list(),
        method_name = method_name, extra_output = extra,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}

# ------------------------------------------------------------------ #
#  BON95 model layer (port of _bon95.py / _parametric_shared.py)
# ------------------------------------------------------------------ #

# Thermal peak temperature of BON95 (MeV), = 0.035 eV
.bon95_Tth <- 3.5e-8
.bon95_SHAPE_NAMES <- c("b", "Tf", "c")
.bon95_SHAPE_BOUNDS <- list(b = c(0.5, 2.0), Tf = c(0.5, 10.0), c = c(0.5, 3.0))
.bon95_B_RANGE <- c(0.5, 2.0, 5)
.bon95_TF_RANGE <- c(0.5, 10.0, 5)
.bon95_C_RANGE <- c(0.5, 3.0, 4)

.bon95_Fth <- function(E, Tth = .bon95_Tth) {
    Xth <- E / Tth
    Xth^1.5 * exp(-Xth)
}

.bon95_Fepi <- function(E, b, Tth = .bon95_Tth) {
    Xth <- E / Tth
    E^(-b) * (1 - exp(-Xth))
}

.bon95_Fint <- function(E, Tth = .bon95_Tth) {
    Xth <- E / Tth
    1 - exp(-Xth)
}

.bon95_Ff <- function(E, Tf, c) {
    Xf <- (E / Tf)^c
    Xf^1.5 * exp(-Xf)
}

#' BON95 lethargy spectrum E*Phi(E)
#' @keywords internal
#' @noRd
.bon95_model <- function(E, b, Tf, c, a1, a2, a3, a4) {
    E <- as.numeric(E)
    a1 * .bon95_Fth(E) + a2 * .bon95_Fepi(E, b) + a3 * .bon95_Fint(E) +
        a4 * .bon95_Ff(E, Tf, c)
}

#' BON95 fluence spectrum Phi(E) = E*Phi(E) / E
#' @keywords internal
#' @noRd
.bon95_spectrum <- function(E, b, Tf, c, a1, a2, a3, a4) {
    E <- as.numeric(E)
    lethargy <- .bon95_model(E, b, Tf, c, a1, a2, a3, a4)
    ifelse(E > 0, lethargy / E, 0)
}

#' numpy.linalg.lstsq(rcond=None) equivalent: SVD pseudo-inverse with the
#' LAPACK gelsd rank cut-off eps * max(M, N) * s_max.
#' @keywords internal
#' @noRd
.numpy_lstsq <- function(X, y) {
    X <- as.matrix(X); y <- as.numeric(y)
    sv <- svd(X, nu = min(dim(X)), nv = min(dim(X)))
    if (length(sv$d) == 0L) return(rep(0, ncol(X)))
    tol <- .Machine$double.eps * max(dim(X)) * sv$d[1L]
    keep <- which(sv$d > tol)
    if (length(keep) == 0L) return(rep(0, ncol(X)))
    z <- as.vector(crossprod(sv$u[, keep, drop = FALSE], y)) / sv$d[keep]
    as.vector(sv$v[, keep, drop = FALSE] %*% z)
}

#' Weighted NLS for the linear coefficients a1..a4 given shape parameters
#' (port of _solve_linear_coefficients)
#' @keywords internal
#' @noRd
.bon95_solve_linear_coefficients <- function(A, b, E, ln_steps, b_param, Tf,
                                             c_param, weights = NULL) {
    A <- as.matrix(A)
    n_det <- nrow(A)
    if (is.null(weights)) weights <- rep(1, n_det)
    E_safe <- ifelse(E > 0, E, 1)
    F_cols <- cbind(.bon95_Fth(E), .bon95_Fepi(E, b_param),
                    .bon95_Fint(E), .bon95_Ff(E, Tf, c_param))
    weighted_F <- (F_cols / E_safe) * ln_steps
    B <- A %*% weighted_F
    sw <- sqrt(weights)
    a <- pmax(.numpy_lstsq(B * sw, as.numeric(b) * sw), 0)
    residual <- as.vector(B %*% a) - as.numeric(b)
    chi2 <- if (all(weights > 0)) {
        mean(residual^2 * weights)
    } else {
        mean(residual^2)
    }
    list(a = a, chi2 = chi2)
}

#' Best-fit spectrum at fixed shape parameters (port of _solve_shape_nls)
#' @keywords internal
#' @noRd
.bon95_solve_shape_nls <- function(A, b, E, ln_steps, b_param, Tf, c_param,
                                   weights) {
    fit <- .bon95_solve_linear_coefficients(A, b, E, ln_steps, b_param, Tf,
                                            c_param, weights)
    phi <- .bon95_spectrum(E, b_param, Tf, c_param, fit$a[1L], fit$a[2L],
                           fit$a[3L], fit$a[4L])
    phi <- .clean_edge_bins(pmax(phi, 0))
    list(spectrum = phi * ln_steps, chi2 = fit$chi2, a = fit$a)
}

#' Grid search + NLS over (b, Tf, c) (port of solve_bon95_parametric)
#' @keywords internal
#' @noRd
.bon95_grid_search <- function(A, b, E, ln_steps, b_range = .bon95_B_RANGE,
                               Tf_range = .bon95_TF_RANGE,
                               c_range = .bon95_C_RANGE, b_meas = NULL,
                               top_n = 5L) {
    A <- as.matrix(A)
    b_vals <- seq(b_range[1L], b_range[2L],
                  length.out = as.integer(b_range[3L]))
    Tf_vals <- seq(Tf_range[1L], Tf_range[2L],
                   length.out = as.integer(Tf_range[3L]))
    c_vals <- seq(c_range[1L], c_range[2L],
                  length.out = as.integer(c_range[3L]))
    weights <- if (!is.null(b_meas)) {
        ifelse(b_meas > 0, 1 / (b_meas^2), 1)
    } else {
        rep(1, nrow(A))
    }
    chi2s <- numeric(0L)
    candidates <- list()
    for (bv in b_vals) {
        for (Tfv in Tf_vals) {
            for (cv in c_vals) {
                fit <- .bon95_solve_linear_coefficients(A, b, E, ln_steps, bv,
                                                        Tfv, cv, weights)
                chi2s <- c(chi2s, fit$chi2)
                candidates[[length(candidates) + 1L]] <- list(
                    b = bv, Tf = Tfv, c = cv,
                    a1 = fit$a[1L], a2 = fit$a[2L], a3 = fit$a[3L],
                    a4 = fit$a[4L], chi2 = fit$chi2)
            }
        }
    }
    ord <- order(chi2s, method = "radix")
    best <- candidates[[ord[1L]]]
    top <- candidates[ord[seq_len(min(as.integer(top_n), length(ord)))]]
    list(best = best, best_chi2 = best$chi2, top = top,
         n_candidates = length(ord))
}

#' Jacobian of the best-fit spectrum w.r.t. the shape parameters
#' (port of _compute_bon95_shape_jacobian)
#' @keywords internal
#' @noRd
.bon95_shape_jacobian <- function(A, b, E, ln_steps, params, weights,
                                  delta = 1e-6) {
    A <- as.matrix(A)
    s0 <- .bon95_solve_shape_nls(A, b, E, ln_steps, params[["b"]],
                                params[["Tf"]], params[["c"]],
                                weights)$spectrum
    residual <- as.numeric(A %*% s0) - as.numeric(b)
    J <- matrix(0, nrow = length(E), ncol = 3L)
    for (i in seq_along(.bon95_SHAPE_NAMES)) {
        nm <- .bon95_SHAPE_NAMES[i]
        lo <- .bon95_SHAPE_BOUNDS[[nm]][1L]
        hi <- .bon95_SHAPE_BOUNDS[[nm]][2L]
        p_val <- params[[nm]]
        d <- delta
        if (p_val + d > hi) d <- max(0, hi - p_val) * 0.5
        if (d < 1e-15) {
            d <- delta
            if (p_val - d >= lo) {
                pert <- params
                pert[[nm]] <- p_val - d
                s_pert <- .bon95_solve_shape_nls(A, b, E, ln_steps, pert[["b"]],
                                                 pert[["Tf"]], pert[["c"]],
                                                 weights)$spectrum
                J[, i] <- (s0 - s_pert) / d
            } else {
                J[, i] <- 0
            }
            next
        }
        pert <- params
        pert[[nm]] <- p_val + d
        s_pert <- .bon95_solve_shape_nls(A, b, E, ln_steps, pert[["b"]],
                                         pert[["Tf"]], pert[["c"]],
                                         weights)$spectrum
        J[, i] <- (s_pert - s0) / d
    }
    list(J = J, residual = residual)
}

#' Clamp BON95 shape parameters to their bounds
#' @keywords internal
#' @noRd
.bon95_clamp_shape <- function(params, bounds = .bon95_SHAPE_BOUNDS) {
    for (nm in names(bounds)) {
        params[[nm]] <- max(bounds[[nm]][1L], min(bounds[[nm]][2L],
                                                 params[[nm]]))
    }
    params
}

#' Zero anomalously large edge bins (port of _clean_edge_bins)
#' @keywords internal
#' @noRd
.clean_edge_bins <- function(phi, factor = 10) {
    phi <- as.numeric(phi)
    n <- length(phi)
    if (n < 3L) return(phi)
    neighbor_mean <- mean(phi[2L:3L])
    if (neighbor_mean > 0 && phi[1L] > factor * neighbor_mean) phi[1L] <- 0
    neighbor_mean <- mean(phi[(n - 2L):(n - 1L)])
    if (neighbor_mean > 0 && phi[n] > factor * neighbor_mean) phi[n] <- 0
    phi
}

#' Estimate measurement uncertainties from the readings
#' (port of _build_measurement_uncertainties)
#' @keywords internal
#' @noRd
.build_measurement_uncertainties <- function(b, noise_level = 0.05) {
    abs(as.numeric(b)) * noise_level + 1e-30
}

#' Directed-divergence (I-divergence) refinement of a BON95 spectrum
#'
#' Port of \code{directed_divergence_iteration} from
#' \code{unfold_parametric2.py}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid (length n).
#' @param ln_steps Numeric \eqn{d(\ln E)} bin widths (length n).
#' @param phi0 Numeric initial fluence-spectrum guess (length n).
#' @param b_meas Optional measurement uncertainties (weights are
#'   \code{1/b_meas^2} where \code{b_meas > 0}).
#' @param max_iter Integer; maximum iterations. Default 200.
#' @param tol_chi2 Numeric; stop when the weighted mean chi-square drops
#'   below this value. Default 1.
#' @param tol_rel Numeric; stop when the relative change of the spectrum
#'   falls below this value. Default 1e-6.
#' @return A list \code{list(spectrum, iterations, chi2, converged)}.
#' @export
directed_divergence_iteration <- function(A, b, E, ln_steps, phi0,
                                          b_meas = NULL, max_iter = 200L,
                                          tol_chi2 = 1.0, tol_rel = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); ln_steps <- as.numeric(ln_steps)
    phi <- pmax(as.numeric(phi0), 1e-30)
    n_det <- nrow(A)
    weights <- if (!is.null(b_meas)) {
        ifelse(b_meas > 0, 1 / (b_meas^2), 1)
    } else {
        rep(1, n_det)
    }
    denom <- pmax(as.numeric(colSums(A)), 1e-30)

    for (iteration in seq_len(as.integer(max_iter))) {
        M_p <- as.numeric(A %*% (phi * ln_steps))
        M_p_safe <- pmax(M_p, 1e-30)
        residual <- M_p - b
        chi2 <- mean(residual^2 * weights)
        if (chi2 < tol_chi2) {
            return(list(spectrum = phi, iterations = iteration,
                        chi2 = chi2, converged = TRUE))
        }
        ratios <- b / M_p_safe
        numerator <- as.numeric(crossprod(A, ratios))
        phi_new <- pmax(phi * numerator / denom, 1e-30)
        rel_change <- max(abs(phi_new - phi)) / (max(phi) + 1e-30)
        phi <- phi_new
        if (rel_change < tol_rel) {
            phi <- .clean_edge_bins(phi)
            M_p_final <- as.numeric(A %*% (phi * ln_steps))
            chi2_final <- mean((M_p_final - b)^2 * weights)
            return(list(spectrum = phi, iterations = iteration,
                        chi2 = chi2_final, converged = TRUE))
        }
    }
    phi <- .clean_edge_bins(phi)
    M_p_final <- as.numeric(A %*% (phi * ln_steps))
    chi2_final <- mean((M_p_final - b)^2 * weights)
    list(spectrum = phi, iterations = as.integer(max_iter),
         chi2 = chi2_final, converged = chi2_final < tol_chi2)
}

#' BON95 parametric unfolding (4-component model)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_parametric2.py} plus
#' its model layer \code{_bon95.py}.
#' The spectrum E*Phi(E) is a linear combination of four components
#' (thermal, epithermal, intermediate, fast) whose shape parameters
#' \code{(b, Tf, c)} are found by grid search, the linear coefficients
#' \code{a1..a4} by weighted NLS, and the result is refined by
#' directed-divergence iterations.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid (length n).
#' @param ln_steps Numeric \eqn{d(\ln E)} bin widths (length n). As in Python
#'   the model spectrum is weighted by these steps once more, i.e. the caller
#'   passes \code{compute_log_steps(E) * log(10)}.
#' @param b_meas Optional measurement uncertainties used as NLS weights.
#' @param optimizer Character: \code{"grid"} (default), \code{"cvxpy"},
#'   \code{"qpsolvers"} or \code{"combined"}.
#' @param b_range,Tf_range,c_range Length-3 numeric vectors
#'   \code{c(min, max, n_points)} for the shape grid search.
#' @param alpha Numeric; Tikhonov weight for the SQP optimizers. Default 1e-4.
#' @param solver_backend Character; QP backend label (accepted for parity;
#'   the pure-R engine always uses the internal box-QP active set).
#' @param max_iter_qp Integer; maximum SQP iterations. Default 50.
#' @param tol_qp Numeric; SQP convergence tolerance. Default 1e-6.
#' @param max_iter Integer; maximum directed-divergence iterations. Default 200.
#' @param tol_chi2 Numeric; chi-square stop threshold for the DD phase.
#'   Default 1.
#' @return A list \code{list(spectrum, iterations, converged, chi2, params,
#'   message)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_parametric2(A, b, E, compute_log_steps(E) * log(10))
solve_parametric2 <- function(A, b, E, ln_steps, b_meas = NULL,
                              optimizer = "grid",
                              b_range = .bon95_B_RANGE,
                              Tf_range = .bon95_TF_RANGE,
                              c_range = .bon95_C_RANGE,
                              alpha = 1e-4, solver_backend = "auto",
                              max_iter_qp = 50L, tol_qp = 1e-6,
                              max_iter = 200L, tol_chi2 = 1.0) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E); ln_steps <- as.numeric(ln_steps)
    optimizer <- match.arg(tolower(optimizer),
                           c("grid", "cvxpy", "qpsolvers", "combined"))

    if (optimizer == "grid") {
        grid <- .bon95_grid_search(A, b, E, ln_steps, b_range = b_range,
                                  Tf_range = Tf_range, c_range = c_range,
                                  b_meas = b_meas, top_n = 5L)
        best <- grid$best
        phi_param <- .bon95_spectrum(E, best$b, best$Tf, best$c,
                                     best$a1, best$a2, best$a3, best$a4)
        nfev <- length(grid$top)
        best_chi2 <- grid$best_chi2
    } else {
        qp <- .bon95_solve_sqp(A, b, E, ln_steps, b_meas = b_meas,
                               optimizer = optimizer, alpha = alpha,
                               max_iter = as.integer(max_iter_qp),
                               tol = tol_qp)
        phi_param <- qp$spectrum / pmax(ln_steps, 1e-30)
        nfev <- qp$nfev
        best_chi2 <- 0
    }

    phi_param <- .clean_edge_bins(pmax(as.numeric(phi_param), 0))

    dd <- directed_divergence_iteration(A, b, E, ln_steps, phi_param,
                                        b_meas = b_meas,
                                        max_iter = as.integer(max_iter),
                                        tol_chi2 = tol_chi2)
    phi_refined <- .clean_edge_bins(dd$spectrum)
    spectrum <- phi_refined * ln_steps

    computed <- as.numeric(A %*% spectrum)
    .fruit_check_fit_quality(sqrt(sum((computed - b)^2)), b, "parametric2",
                             model_descr = "4-component BON95 model")
    list(spectrum = spectrum, iterations = nfev, converged = dd$converged,
         chi2 = dd$chi2, dd_iterations = dd$iterations,
         best_chi2 = best_chi2,
         message = sprintf("BON95 %s fit (chi2=%.4f) + DD (%d iters, chi2=%.4f)",
                           optimizer, best_chi2, dd$iterations, dd$chi2))
}

#' BON95 SQP step: box-constrained QP on the shape-parameter update.
#'
#' Pure-R replacement for the cvxpy / qpsolvers backends of
#' \code{_bon95.py}.  Solves
#' \eqn{\min_\delta \|A_{eff}\delta + r\|^2 + \alpha\|\delta\|^2}
#' subject to the per-parameter box, with an exact active-set method.
#' @keywords internal
#' @noRd
.bon95_solve_sqp <- function(A, b, E, ln_steps, b_meas = NULL,
                             optimizer = "cvxpy", alpha = 1e-4,
                             max_iter = 50L, tol = 1e-6) {
    A <- as.matrix(A)
    weights <- if (!is.null(b_meas)) {
        ifelse(b_meas > 0, 1 / (b_meas^2), 1)
    } else {
        rep(1, nrow(A))
    }
    grid <- .bon95_grid_search(A, b, E, ln_steps, b_meas = b_meas, top_n = 1L)
    params <- list(b = grid$best$b, Tf = grid$best$Tf, c = grid$best$c)
    params <- .bon95_clamp_shape(params)
    nfev <- 0L
    message <- ""
    for (k in seq_len(as.integer(max_iter))) {
        spec_k <- .bon95_solve_shape_nls(A, b, E, ln_steps, params[["b"]],
                                        params[["Tf"]], params[["c"]],
                                        weights)$spectrum
        nfev <- nfev + 1L
        residual <- as.numeric(A %*% spec_k) - b
        if (sqrt(sum(residual^2)) < tol) {
            return(list(spectrum = spec_k, converged = TRUE, nfev = nfev,
                        message = sprintf("Converged in %d iterations", k)))
        }
        jac <- .bon95_shape_jacobian(A, b, E, ln_steps, params, weights)
        nfev <- nfev + 2L * 3L
        A_eff <- A %*% jac$J
        H <- crossprod(A_eff) + alpha * diag(3L)
        g <- as.vector(crossprod(A_eff, residual))
        lo <- vapply(.bon95_SHAPE_NAMES,
                     function(nm) .bon95_SHAPE_BOUNDS[[nm]][1L] - params[[nm]],
                     numeric(1L))
        hi <- vapply(.bon95_SHAPE_NAMES,
                     function(nm) .bon95_SHAPE_BOUNDS[[nm]][2L] - params[[nm]],
                     numeric(1L))
        delta <- .solve_box_qp(as.matrix(H), g, lo, hi)
        for (i in seq_along(.bon95_SHAPE_NAMES)) {
            params[[.bon95_SHAPE_NAMES[i]]] <-
                params[[.bon95_SHAPE_NAMES[i]]] + delta[i]
        }
        params <- .bon95_clamp_shape(params)
        if (sqrt(sum(delta^2)) < tol) {
            spec_f <- .bon95_solve_shape_nls(A, b, E, ln_steps, params[["b"]],
                                             params[["Tf"]], params[["c"]],
                                             weights)$spectrum
            return(list(spectrum = spec_f, converged = TRUE, nfev = nfev,
                        message = sprintf("Converged in %d iterations", k)))
        }
    }
    spec_f <- .bon95_solve_shape_nls(A, b, E, ln_steps, params[["b"]],
                                     params[["Tf"]], params[["c"]],
                                     weights)$spectrum
    list(spectrum = spec_f, converged = FALSE, nfev = nfev,
         message = sprintf("Max iterations (%d) reached", max_iter))
}

#' Exact box-constrained quadratic program
#' \eqn{\min 0.5 x'Hx + g'x}{min 0.5 x'Hx + g'x} with \eqn{lo <= x <= hi}.
#'
#' Lawson-Hanson style active set for a strictly convex QP with only bound
#' constraints; used in place of the cvxpy / qpsolvers backends.
#'
#' @param H Symmetric positive definite matrix.
#' @param g Gradient vector.
#' @param lo,hi Bound vectors.
#' @param max_active_sets Integer; safety limit on active-set swaps.
#' @return Numeric solution vector.
#' @keywords internal
#' @noRd
.solve_box_qp <- function(H, g, lo, hi, max_active_sets = 200L) {
    n <- length(g)
    x <- pmin(pmax(-as.vector(solve(H, g)), lo), hi)
    if (!all(is.finite(x))) x <- rep(0, n)
    # state: 0 free, -1 at lo, +1 at hi
    state <- ifelse(x <= lo + 0, -1L, ifelse(x >= hi - 0, 1L, 0L))
    state <- as.integer(state)
    for (it in seq_len(as.integer(max_active_sets))) {
        free <- which(state == 0L)
        fixed <- which(state != 0L)
        if (length(free) > 0L) {
        rhs <- -g[free]
        if (length(fixed) > 0L) rhs <- rhs - H[free, fixed, drop = FALSE] %*% x[fixed]
        sol <- tryCatch(solve(H[free, free, drop = FALSE], rhs),
                        error = function(e) NULL)
        if (is.null(sol)) sol <- .numpy_lstsq(H[free, free, drop = FALSE], rhs)
        x_new <- x
        x_new[free] <- as.vector(sol)
        } else {
            x_new <- x
        }
        viol <- which(state == 0L & (x_new < lo | x_new > hi))
        if (length(viol) > 0L) {
            # fix the most violated variable at its nearest bound
            over <- pmax(x_new - hi, 0)
            under <- pmax(lo - x_new, 0)
            amount <- over[viol] + under[viol]
            j <- viol[which.max(amount)]
            state[j] <- if (x_new[j] > hi[j]) 1L else -1L
            x[j] <- if (state[j] == 1L) hi[j] else lo[j]
            next
        }
        x <- pmin(pmax(x_new, lo), hi)
        # KKT check on the fixed variables
        w <- as.vector(H %*% x) + g
        bad_lo <- which(state == -1L & w < 0)
        bad_hi <- which(state == 1L & w > 0)
        if (length(bad_lo) == 0L && length(bad_hi) == 0L) return(x)
        cand <- c(bad_lo, bad_hi)
        scores <- abs(w[cand])
        j <- cand[which.max(scores)]
        state[j] <- 0L
    }
    x
}

#' Wrapper around \code{\link{solve_parametric2}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_parametric2
#' @param noise_level Numeric; relative measurement uncertainty feeding the
#'   BON95 weights (\code{b_meas = |b| * noise_level + 1e-30}). Default 0.05.
#' @param calculate_errors Logical; run Monte-Carlo uncertainty estimation.
#' @param n_montecarlo Integer; Monte-Carlo sample count.
#' @param save_result Logical; hand the result to \code{save_result_callback}.
#' @param random_state Optional integer seed for the Monte-Carlo part.
#' @param max_neutron_energy Optional energy cutoff in MeV.
#' @param initial_params,max_iterations,tolerance Legacy arguments kept for API
#'   compatibility with the Python workflow; unused by the parametric fit.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_parametric2 <- function(detector_names, n_energy_bins, E_MeV,
                                  sensitivities, cc_icrp116, save_result_callback,
                                  readings, initial_spectrum = NULL,
                                  optimizer = "grid",
                                  b_range = .bon95_B_RANGE,
                                  Tf_range = .bon95_TF_RANGE,
                                  c_range = .bon95_C_RANGE,
                                  alpha = 1e-4, solver_backend = "auto",
                                  max_iter_qp = 50L, tol_qp = 1e-6,
                                  noise_level = 0.05,
                                  max_iter = 200L, tol_chi2 = 1.0,
                                  initial_params = NULL,
                                  max_iterations = 200L, tolerance = 1e-6,
                                  calculate_errors = FALSE,
                                  n_montecarlo = 100L,
                                  save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b
    log_steps <- compute_log_steps(E_MeV)
    ln_steps <- log_steps * log(10)

    b_meas_local <- .build_measurement_uncertainties(b, noise_level)

    solver <- function(A, b, x0 = NULL, ...) {
        bm <- .build_measurement_uncertainties(b, noise_level)
        res <- solve_parametric2(A, b, E_MeV, ln_steps, b_meas = bm,
                                 optimizer = optimizer,
                                 b_range = b_range, Tf_range = Tf_range,
                                 c_range = c_range, alpha = alpha,
                                 solver_backend = solver_backend,
                                 max_iter_qp = max_iter_qp, tol_qp = tol_qp,
                                 max_iter = max_iter, tol_chi2 = tol_chi2)
        list(spectrum = res$spectrum, iterations = res$iterations,
             converged = res$converged)
    }

    x0_default <- rep(mean(b) / max(mean(rowSums(A)), .Machine$double.xmin),
                      n_energy_bins)

    out <- run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver, solve_kwargs = list(),
        method_name = "parametric2",
        extra_output = list(optimizer = optimizer, b_range = b_range,
                            Tf_range = Tf_range, c_range = c_range,
                            alpha = alpha, solver_backend = solver_backend,
                            noise_level = noise_level, bon95_Tth = .bon95_Tth,
                            b_meas = b_meas_local),
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
    out
}
