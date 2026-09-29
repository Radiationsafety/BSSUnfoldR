#' P-spline REML unfolding (mixed-model smoothing)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_pspline_reml.py}.
#' The unfolded spectrum is represented in a clamped cubic B-spline basis on
#' the (energy-dependent) knot grid and the smoothness is selected by
#' restricted maximum likelihood (REML) in the linear-mixed-model
#' formulation of \code{LMMsolver} (Eilers & Marx 1996; Wand & Ormerod 2008;
#' Currie et al. 2004): the P-spline penalty \code{G = D^(d)T D^(d)} is
#' diagonalised, its null space becomes the unpenalised fixed effect and its
#' range space the penalised random effect, Henderson's mixed-model
#' equations are solved for the coefficients and the variance ratio
#' \code{lam} is found by maximising the REML profile along
#' \code{log10(lam/lam_ref)}.
#'
#' @name pspline-reml-methods
NULL

# Numerical guards (mirrors of the Python module constants) --------------
.pspline_TINY <- 1e-300
# Search interval for the *relative* smoothing parameter (log10 units)
.pspline_LAMBDA_REL_BOUNDS <- c(1e-6, 1e6)
# Eigenvalues of G below this relative threshold form the null space.
.pspline_NULLSPACE_RTOL <- 1e9 * .Machine$double.eps
# scipy.optimize.minimize_scalar(method = "bounded") defaults: the bounded
# Brent search uses xatol = 1e-5 on log10(lam/lam_ref) and maxfun = 500.
.pspline_XATOL <- 1e-5
.pspline_MAXFUN <- 500L

# ---- internal helpers ------------------------------------------------------

# Cox-de Boor design matrix for the knot vector `t` of degree `degree`,
# evaluated at `x`. Mirrors scipy's BSpline.design_matrix(x, t, degree)
# (clamped knot vector, values clipped to [t[1], t[n]], right end of the
# domain assigned to the last non-empty knot interval).
.pspline_design_matrix <- function(x, t, degree) {
    nt <- length(t)
    n_basis <- nt - degree - 1L
    x <- pmin(pmax(as.numeric(x), t[1]), t[nt])
    # degree-0 basis functions N_{i,0} = 1 on [t_i, t_{i+1})
    N <- matrix(0.0, nrow = length(x), ncol = nt - 1L)
    last_int <- nt - degree - 1L
    for (i in seq_len(nt - 1L)) {
        sel <- (x >= t[i]) & (x < t[i + 1L])
        if (i == last_int) sel <- sel | (x == t[nt])
        N[, i] <- as.numeric(sel)
    }
    for (j in seq_len(degree)) {
        ncols <- nt - j - 1L
        Nn <- matrix(0.0, nrow = length(x), ncol = ncols)
        for (i in seq_len(ncols)) {
            d1 <- t[i + j] - t[i]
            d2 <- t[i + j + 1L] - t[i + 1L]
            v <- numeric(length(x))
            if (d1 > 0) v <- v + (x - t[i]) / d1 * N[, i]
            if (d2 > 0) v <- v + (t[i + j + 1L] - x) / d2 * N[, i + 1L]
            Nn[, i] <- v
        }
        N <- Nn
    }
    N[, seq_len(n_basis), drop = FALSE]
}

# knot_spacing == "auto" resolution (Python _resolve_knot_spacing): log-spaced
# interior knots when the grid spans a ratio > 100, uniform knots otherwise.
.pspline_resolve_knot_spacing <- function(knot_spacing, E) {
    if (knot_spacing != "auto") return(knot_spacing)
    if (length(E) == 0 || min(E) <= 0) return("log")
    if ((max(E) / min(E)) > 100) "log" else "uniform"
}

# Clamped B-spline basis, mirror of Python's
# core/unfold_mlem_bs.py:build_bspline_basis(E, n_basis, spline_order, spacing)
.pspline_bspline_basis <- function(E, n_basis, spline_order = 4L,
                                   knot_spacing = "auto") {
    E <- as.numeric(E)
    n_basis <- as.integer(n_basis); spline_order <- as.integer(spline_order)
    spacing <- .pspline_resolve_knot_spacing(knot_spacing, E)
    emin <- E[1]; emax <- E[length(E)]
    n_interior <- n_basis - spline_order
    if (n_interior > 0) {
        grid <- if (spacing == "log") {
            10^seq(log10(emin), log10(emax), length.out = n_interior + 2L)
        } else {
            seq(emin, emax, length.out = n_interior + 2L)
        }
        interior <- grid[-c(1L, n_interior + 2L)]
    } else {
        interior <- numeric(0)
    }
    t <- c(rep(emin, spline_order), interior, rep(emax, spline_order))
    .pspline_design_matrix(E, t, spline_order - 1L)
}

# Python's difference_matrix(n, order) = np.diff(np.eye(n), order, axis=0)
.pspline_difference_matrix <- function(n, order = 2L) {
    n <- as.integer(n); order <- as.integer(order)
    if (n < 1) stop("n must be a positive integer, got ", n)
    if (order < 1 || order > 4) stop("order must be in [1, 4], got ", order)
    if (n <= order) stop("n (", n, ") must be greater than the difference order (", order, ")")
    D <- diag(as.integer(n))
    as.matrix(diff(D, differences = order, margin = 1L))
}

# eigen-based fixed/random split of the P-spline penalty (mixed_model_split)
.pspline_mixed_model_split <- function(n, order = 2L) {
    G <- as.matrix(crossprod(.pspline_difference_matrix(n, order)))
    ev <- eigen(G, symmetric = TRUE)
    o <- order(ev$values)                # ascending, like numpy.linalg.eigh
    g <- ev$values[o]
    U <- ev$vectors[, o, drop = FALSE]
    g_max <- max(g[length(g)], .pspline_TINY)
    is_fixed <- g <= .pspline_NULLSPACE_RTOL * g_max
    if (!any(is_fixed)) {
        is_fixed[which.min(g)] <- TRUE
    }
    list(U_fixed = U[, is_fixed, drop = FALSE],
         U_random = U[, !is_fixed, drop = FALSE],
         g_random = pmax(g[!is_fixed], .pspline_TINY))
}

# REML profile log-likelihood for one smoothing value (reml_profile).
# Returns c(loglik, sigma2_hat).
.pspline_reml_profile <- function(y_w, X_w, Z_w, g_random, lam) {
    m <- length(y_w)
    p_f <- ncol(X_w)
    df <- m - p_f
    if (df < 1) return(c(-Inf, NaN))
    L_inv <- 1 / pmax(g_random, .pspline_TINY)
    V <- diag(m) + lam * as.matrix(sweep(Z_w, 2L, L_inv, "*") %*% t(Z_w))
    ok <- tryCatch({chol(V); TRUE}, error = function(e) FALSE,
                   warning = function(e) FALSE)
    if (!isTRUE(ok)) return(c(-Inf, NaN))
    Vinv_X <- solve(V, X_w)
    XtVinvX <- as.matrix(crossprod(X_w, Vinv_X))
    dX <- determinant(XtVinvX, logarithm = TRUE)
    if (dX$sign <= 0) return(c(-Inf, NaN))
    logdet_XtVX <- as.numeric(dX$modulus)
    beta <- tryCatch(solve(XtVinvX, as.numeric(crossprod(X_w, solve(V, y_w)))),
                     error = function(e) NULL)
    if (is.null(beta)) return(c(-Inf, NaN))
    resid <- as.numeric(y_w) - as.numeric(Vinv_X %*% beta)
    ss <- as.numeric(crossprod(resid, solve(V, resid)))
    if (!is.finite(ss) || ss <= 0) ss <- max(ss, .pspline_TINY)
    sigma2 <- ss / df
    logdet_V <- as.numeric(determinant(V, logarithm = TRUE)$modulus)
    loglik <- -0.5 * (df * log(sigma2) + logdet_V + logdet_XtVX)
    c(loglik, sigma2)
}

# Exact port of scipy.optimize's bounded Brent scalar minimiser
# (_minimize_scalar_bounded, xatol default 1e-5, maxfun = maxiter = 500).
# `f` must be vectorised over a scalar; NaN comparisons follow R's isTRUE
# semantics so that they evaluate FALSE, as in numpy.
.pspline_min_bounded <- function(f, x1, x2, xatol = 1e-5, maxfun = 500L) {
    sqrt_eps <- sqrt(2.2e-16)
    golden_mean <- 0.5 * (3 - sqrt(5))
    a <- x1; b <- x2
    fulc <- a + golden_mean * (b - a)
    nfc <- fulc; xf <- fulc
    rat <- 0; e <- 0
    x <- xf
    fx <- f(x)
    num <- 1L
    fu <- Inf
    ffulc <- fx; fnfc <- fx
    xm <- 0.5 * (a + b)
    tol1 <- sqrt_eps * abs(xf) + xatol / 3
    tol2 <- 2 * tol1
    while (isTRUE(abs(xf - xm) > (tol2 - 0.5 * (b - a)))) {
        golden <- 1L
        if (isTRUE(abs(e) > tol1)) {
            golden <- 0L
            r <- (xf - nfc) * (fx - ffulc)
            q <- (xf - fulc) * (fx - fnfc)
            p <- (xf - fulc) * q - (xf - nfc) * r
            q <- 2 * (q - r)
            if (isTRUE(q > 0)) p <- -p
            q <- abs(q)
            r <- e
            e <- rat
            if (isTRUE(abs(p) < abs(0.5 * q * r)) &&
                isTRUE(p > q * (a - xf)) &&
                isTRUE(p < q * (b - xf))) {
                rat <- (p + 0) / q
                x <- xf + rat
                if (isTRUE((x - a) < tol2) || isTRUE((b - x) < tol2)) {
                    si <- sign(xm - xf) + as.numeric(xm == xf)
                    rat <- tol1 * si
                }
            } else {
                golden <- 1L
            }
        }
        if (golden == 1L) {
            e <- if (isTRUE(xf >= xm)) a - xf else b - xf
            rat <- golden_mean * e
        }
        si <- sign(rat) + as.numeric(rat == 0)
        stp1 <- si * pmax(abs(rat), tol1)
        if (is.na(stp1)) stp1 <- NA_real_
        x <- xf + stp1
        fu <- f(x)
        num <- num + 1L
        if (isTRUE(fu <= fx)) {
            if (isTRUE(x >= xf)) a <- xf else b <- xf
            fulc <- nfc; ffulc <- fnfc
            nfc <- xf; fnfc <- fx
            xf <- x; fx <- fu
        } else {
            if (isTRUE(x < xf)) a <- x else b <- x
            if (isTRUE(fu <= fnfc) || isTRUE(nfc == xf)) {
                fulc <- nfc; ffulc <- fnfc
                nfc <- x; fnfc <- fu
            } else if (isTRUE(fu <= ffulc) || isTRUE(fulc == xf) ||
                       isTRUE(fulc == nfc)) {
                fulc <- x; ffulc <- fu
            }
        }
        xm <- 0.5 * (a + b)
        tol1 <- sqrt_eps * abs(xf) + xatol / 3
        tol2 <- 2 * tol1
        if (num >= maxfun) break
    }
    list(x = xf, fun = fx, nfev = num)
}

# select_lambda_reml: maximise the REML profile over the relative smoothing
# parameter (cached neg-loglik, bounded Brent on log10 units).
.pspline_select_lambda_reml <- function(y_w, X_w, Z_w, g_random, lam_ref,
                                        lam_bounds = .pspline_LAMBDA_REL_BOUNDS,
                                        xatol = 1e-5, maxfun = 500L) {
    lo <- log10(max(lam_bounds[1], .pspline_TINY))
    hi <- log10(lam_bounds[2])
    n_fev <- 0L
    cache <- new.env(parent = emptyenv())
    neg_loglik <- function(t) {
        n_fev <<- n_fev + 1L
        key <- format(t, digits = 17)
        hit <- cache[[key]]
        if (!is.null(hit)) return(hit)
        lam <- lam_ref * (10^t)
        ll <- .pspline_reml_profile(y_w, X_w, Z_w, g_random, lam)[1]
        val <- -ll
        cache[[key]] <- val
        val
    }
    res <- .pspline_min_bounded(neg_loglik, lo, hi, xatol = xatol,
                                maxfun = maxfun)
    converged <- is.finite(res$fun)
    t_opt <- res$x
    lam_rel <- 10^t_opt
    lam <- lam_ref * lam_rel
    prof <- .pspline_reml_profile(y_w, X_w, Z_w, g_random, lam)
    list(lam = lam, lam_relative = lam_rel, sigma2 = prof[2],
         reml_loglik = prof[1], converged = converged,
         n_iterations = if (n_fev > 0L) n_fev else res$nfev)
}

# Core solver, mirroring Python's solve_pspline_reml_full + solve_pspline_reml
.pspline_solve <- function(A, b, E_MeV, n_basis = NULL, spline_order = 4L,
                           diff_order = 2L, knot_spacing = "auto",
                           weights = "uniform", lam_relative = NULL,
                           lam_bounds = .pspline_LAMBDA_REL_BOUNDS,
                           nonneg = TRUE, xatol = 1e-5, maxfun = 500L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A); n <- ncol(A)
    E <- as.numeric(E_MeV)
    if (is.null(n_basis)) {
        n_basis <- max(diff_order + 2L, min(n %/% 2L, 30L))
    }
    n_basis <- as.integer(n_basis)
    if (n_basis < diff_order + 2L) {
        stop(sprintf("n_basis (%d) must be >= diff_order + 2 (%d)",
                     n_basis, diff_order + 2L), call. = FALSE)
    }
    if (n_basis > n) {
        stop(sprintf(paste0("n_basis (%d) cannot exceed the number of energy ",
                            "bins (%d)"), n_basis, n), call. = FALSE)
    }
    B <- .pspline_bspline_basis(E, n_basis, spline_order = spline_order,
                                knot_spacing = knot_spacing)
    split <- .pspline_mixed_model_split(n_basis, diff_order)
    X <- A %*% (B %*% split$U_fixed)
    Z <- A %*% (B %*% split$U_random)
    if (is.character(weights) && identical(weights, "poisson")) {
        floor_b <- 1e-3 * max(max(b), .pspline_TINY)
        w <- 1 / pmax(b, floor_b)
    } else if (is.numeric(weights) && length(weights) == m) {
        w <- as.numeric(weights)
    } else {
        w <- rep(1, m)
    }
    sw <- sqrt(w)
    y_w <- sw * b
    X_w <- X * sw
    Z_w <- Z * sw
    data_scale <- sum(Z_w^2) / max(ncol(Z_w), 1L)
    pen_scale <- mean(split$g_random)
    lam_ref <- max(data_scale / max(pen_scale, .pspline_TINY), .pspline_TINY)
    if (is.null(lam_relative)) {
        selection <- .pspline_select_lambda_reml(y_w, X_w, Z_w, split$g_random,
                                                 lam_ref, lam_bounds = lam_bounds,
                                                 xatol = xatol, maxfun = maxfun)
    } else {
        lam_rel <- as.numeric(lam_relative)
        lam <- lam_ref * lam_rel
        prof <- .pspline_reml_profile(y_w, X_w, Z_w, split$g_random, lam)
        selection <- list(lam = lam, lam_relative = lam_rel, sigma2 = prof[2],
                          reml_loglik = prof[1], converged = is.finite(prof[1]),
                          n_iterations = 0L)
    }
    lam <- selection$lam
    # Henderson mixed-model equations
    XtWX <- as.matrix(crossprod(X, X * w))
    XtWZ <- as.matrix(crossprod(X, Z * w))
    ZtWZ <- as.matrix(crossprod(Z, Z * w))
    rhs <- c(as.numeric(crossprod(X, w * b)), as.numeric(crossprod(Z, w * b)))
    p_f <- ncol(split$U_fixed)
    saddle <- rbind(cbind(XtWX, XtWZ),
                    cbind(t(XtWZ), ZtWZ + lam * diag(split$g_random)))
    coef <- tryCatch(as.numeric(solve(saddle, rhs)),
                     error = function(e) as.numeric(qr.solve(saddle, rhs)))
    beta_hat <- coef[seq_len(p_f)]
    b_random_hat <- coef[p_f + seq_len(ncol(split$U_random))]
    spectrum <- as.numeric(B %*% (split$U_fixed %*% beta_hat +
                                      split$U_random %*% b_random_hat))
    ed <- tryCatch({
        ZtZ_lam <- ZtWZ + lam * diag(split$g_random)
        p_f + sum(diag(solve(ZtZ_lam, ZtWZ)))
    }, error = function(e) NaN)
    out <- list(spectrum = spectrum,
                iterations = as.integer(selection$n_iterations),
                converged = isTRUE(selection$converged) &&
                    all(is.finite(pmax(spectrum, 0))),
                lambda = lam,
                coefficients = coef,
                n_basis = n_basis,
                spline_order = spline_order,
                diff_order = diff_order,
                knot_spacing = .pspline_resolve_knot_spacing(knot_spacing, E),
                lam_relative = selection$lam_relative,
                lam_ref = lam_ref,
                sigma2 = selection$sigma2,
                reml_loglik = selection$reml_loglik,
                ed = ed,
                ed_norm = ed / n_basis)
    if (isTRUE(nonneg)) out$spectrum <- pmax(out$spectrum, 0)
    out
}

#' Solve by P-spline REML
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n); reserved for API
#'   compatibility (the fit is linear in the basis).
#' @param E_MeV Numeric energy grid (length n), used for knot placement.
#'   Required, as in Python; knot placement on a \code{NULL} grid is not
#'   meaningful.
#' @param n_basis Integer; dimension of the B-spline space. Default
#'   \code{NULL} = \code{max(diff_order + 2, min(n \%/\% 2, 30))}, the same
#'   rule Python uses.
#' @param spline_order Integer; B-spline order. Default 4 (cubic).
#' @param diff_order Integer; difference order of the penalty. Default 2.
#' @param knot_spacing One of \code{"auto"}, \code{"uniform"} or
#'   \code{"log"}. Default \code{"auto"}.
#' @param weights \code{"uniform"}, \code{"poisson"} or a numeric vector of
#'   length \code{m}. Default \code{"uniform"}.
#' @param lam_relative Optional fixed relative smoothing parameter; when set,
#'   REML selection is skipped. Default \code{NULL}.
#' @param lam_bounds Length-2 search interval for the relative smoothing
#'   parameter.
#' @param max_iterations Integer; maximum number of REML profile evaluations
#'   (function evaluations of the bounded Brent search). Default 50.
#' @param tolerance Numeric; absolute tolerance of the bounded Brent search on
#'   \code{log10(lam/lam_ref)}. Default 1e-3; when left at its default the
#'   search instead uses the reference tolerance 1e-5, which matches the
#'   Python implementation.
#' @return A list \code{list(spectrum, iterations, converged, lambda,
#'   coefficients, ...)} with the REML diagnostics \code{lam_relative},
#'   \code{lam_ref}, \code{sigma2}, \code{reml_loglik}, \code{ed} and
#'   \code{ed_norm}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 12)
#' A <- matrix(runif(3 * 12), nrow = 3)
#' b <- as.numeric(A %*% rep(0.5, 12))
#' r <- solve_pspline_reml(A, b, NULL, E_MeV = E, n_basis = 5,
#'                         max_iterations = 5L)
solve_pspline_reml <- function(A, b, x0 = NULL, E_MeV = NULL, n_basis = NULL,
                               spline_order = 4L, diff_order = 2L,
                               knot_spacing = "auto", weights = "uniform",
                               lam_relative = NULL,
                               lam_bounds = .pspline_LAMBDA_REL_BOUNDS,
                               max_iterations = 50L, tolerance = 1e-3) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (is.null(E_MeV)) {
        stop("E_MeV is required by solve_pspline_reml", call. = FALSE)
    }
    E <- as.numeric(E_MeV)
    if (length(E) != n) {
        stop("E_MeV must have length ncol(A)", call. = FALSE)
    }
    nb <- if (is.null(n_basis)) NULL else as.integer(n_basis)[1L]
    xatol <- if (missing(tolerance)) .pspline_XATOL else
        max(as.numeric(tolerance)[1], .Machine$double.eps)
    maxfun <- if (missing(max_iterations)) .pspline_MAXFUN else
        max(as.integer(max_iterations)[1], 1L)
    .pspline_solve(A, b, E, n_basis = nb, spline_order = spline_order,
                   diff_order = diff_order, knot_spacing = knot_spacing,
                   weights = weights, lam_relative = lam_relative,
                   lam_bounds = lam_bounds,
                   nonneg = TRUE, xatol = xatol, maxfun = maxfun)
}

#' Wrapper around \code{\link{solve_pspline_reml}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_pspline_reml
#' @export
unfold_pspline_reml <- function(detector_names, n_energy_bins, E_MeV,
                                sensitivities, cc_icrp116,
                                save_result_callback, readings,
                                initial_spectrum = NULL, n_basis = NULL,
                                spline_order = 4L, diff_order = 2L,
                                knot_spacing = "auto", weights = "uniform",
                                lam_relative = NULL,
                                lam_bounds = .pspline_LAMBDA_REL_BOUNDS,
                                max_iterations = 50L, tolerance = 1e-3,
                                calculate_errors = FALSE,
                                noise_level = 0.01, n_montecarlo = 100L,
                                save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    nb <- if (is.null(n_basis)) NULL else as.integer(n_basis)[1L]
    xatol <- if (missing(tolerance)) .pspline_XATOL else
        max(as.numeric(tolerance)[1], .Machine$double.eps)
    maxfun <- if (missing(max_iterations)) .pspline_MAXFUN else
        max(as.integer(max_iterations)[1], 1L)
    E <- as.numeric(E_MeV)
    ps <- function(A, b, nonneg) {
        .pspline_solve(A, b, E, n_basis = nb, spline_order = spline_order,
                       diff_order = diff_order, knot_spacing = knot_spacing,
                       weights = weights, lam_relative = lam_relative,
                       lam_bounds = lam_bounds,
                       nonneg = nonneg, xatol = xatol, maxfun = maxfun)
    }
    solver <- function(A, b, x0 = NULL, ...) ps(A, b, nonneg = TRUE)
    keep <- detector_names[detector_names %in% names(readings)]
    diag <- tryCatch(ps(do.call(rbind, lapply(keep,
                                              function(n) {
                                                  as.numeric(sensitivities[[n]])
                                              })),
                        as.numeric(readings[keep]), nonneg = FALSE),
                     error = function(e) NULL)
    extra_output <- if (is.null(diag)) NULL else list(
        n_basis = diag$n_basis,
        lam = diag$lambda,
        lam_relative = diag$lam_relative,
        lam_ref = diag$lam_ref,
        sigma2 = diag$sigma2,
        reml_loglik = diag$reml_loglik,
        ed = diag$ed,
        ed_norm = diag$ed_norm,
        reml_converged = diag$converged)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(0, n_energy_bins),
        solve_func = solver,
        solve_kwargs = list(),
        method_name = "P-spline REML",
        extra_output = extra_output,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
