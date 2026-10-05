#' LMfit / QP / CVXPY-equivalent unfolding methods
#'
#' R ports of \code{bssunfold/src/bssunfold/core/unfold_lmfit.py},
#' \code{unfold_qpsolvers.py}, \code{unfold_cvxpy.py}.
#'
#' @name optimization-methods
NULL

#' Levenberg-Marquardt unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_lmfit.py}.
#' Minimises \eqn{\|A x - b\|^2 + \alpha \|L x\|^2} using
#' \code{\link[stats]{optim}} with \code{L-BFGS-B} method and box constraints.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial guess. Default \code{NULL} = flat.
#' @param alpha Numeric; regularization parameter. Default 0.01.
#' @param smoothness_order Integer; 0 (identity), 1, or 2. Default 0.
#' @param max_iterations Integer; default 1000.
#' @param tolerance Numeric; default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_lmfit(A, b, rep(1, 3), alpha = 0.01, max_iterations = 50)
solve_lmfit <- function(A, b, x0 = NULL, alpha = 0.01,
                          smoothness_order = 0L, max_iterations = 1000L,
                          tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    L <- make_regularization_operator(n, smoothness_order,
                                       identity_for_zero = FALSE)
    .objective <- function(x) {
        x <- pmax(x, 0)
        resid <- as.numeric(A %*% x) - b
        data <- sum(resid^2)
        reg <- if (!is.null(L)) alpha * sum(as.numeric(L %*% x)^2) else 0
        data + reg
    }
    .gradient <- function(x) {
        x <- pmax(x, 0)
        resid <- as.numeric(A %*% x) - b
        grad <- as.numeric(t(A) %*% resid)
        if (!is.null(L)) grad <- grad + alpha * as.numeric(t(L) %*% (L %*% x))
        grad
    }
    if (is.null(x0)) x0 <- rep(mean(b) / max(mean(rowSums(A)), 1e-10), n)
    result <- stats::optim(x0, .objective, .gradient, method = "L-BFGS-B",
                            lower = rep(0, n), upper = rep(Inf, n),
                            control = list(maxit = max_iterations))
    list(spectrum = pmax(as.numeric(result$par), 0),
         iterations = as.integer(result$counts[1]),
         converged = (result$convergence == 0))
}

#' Wrapper around \code{\link{solve_lmfit}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_lmfit
#' @export
unfold_lmfit <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            alpha = 0.01, smoothness_order = 0L,
                            max_iterations = 1000L, tolerance = 1e-6,
                            calculate_errors = FALSE,
                            noise_level = 0.01, n_montecarlo = 100L,
                            save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.5, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_lmfit,
                                         alpha = alpha,
                                         smoothness_order = smoothness_order,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "LMfit",
        extra_output = list(alpha = alpha,
                            smoothness_order = as.integer(smoothness_order)),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}

#' QP solver unfolding (Tikhonov-regularized NNLS)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_qpsolvers.py}.
#' Solves \eqn{\min \|A x - b\|^2 + \alpha \|L x\|^2} with \eqn{x \geq 0}
#' via the augmented-matrix NNLS approach.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused.
#' @param alpha Numeric; regularization. Default 0.01.
#' @param smoothness_order Integer; 0, 1, or 2. Default 0.
#' @return A list \code{list(spectrum, iterations = 1, converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_qpsolvers(A, b, NULL, alpha = 0.01)
solve_qpsolvers <- function(A, b, x0 = NULL, alpha = 0.01,
                              smoothness_order = 0L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    L <- make_regularization_operator(n, smoothness_order,
                                       identity_for_zero = FALSE)
    if (!is.null(L) && alpha > 0) {
        Aw <- rbind(A, sqrt(alpha) * L)
        bw <- c(b, rep(0, nrow(L)))
    } else {
        Aw <- A; bw <- b
    }
    x <- tryCatch(as.numeric(lsei::nnls(Aw, bw)$x),
                  error = function(e) as.numeric(qr.solve(A, b)))
    list(spectrum = pmax(x, 0), iterations = 1L, converged = TRUE)
}

#' Wrapper around \code{\link{solve_qpsolvers}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_qpsolvers
#' @export
unfold_qpsolvers <- function(detector_names, n_energy_bins, E_MeV,
                                 sensitivities, cc_icrp116, save_result_callback,
                                 readings, initial_spectrum = NULL,
                                 alpha = 0.01, smoothness_order = 0L,
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
        solve_func = make_solve_wrapper(solve_qpsolvers,
                                         alpha = alpha,
                                         smoothness_order = smoothness_order),
        solve_kwargs = list(),
        method_name = "QPsolvers",
        extra_output = list(alpha = alpha,
                            smoothness_order = as.integer(smoothness_order)),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}

#' CVXPY-equivalent convex optimization unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_cvxpy.py}.
#' The Python code hands the conic program
#' \deqn{\min_x \; \|A x - b\|_2 + \alpha \, \|x\|_p \quad \text{s.t. } x \ge 0}
#' with \eqn{p} (\code{norm}) equal to 1 or 2, plus per-bin upper bounds that
#' encode the \code{max_neutron_energy} cutoff, to CVXPY, which solves it with
#' its first installed conic solver.  Note that the data term is the
#' \emph{unsquared} 2-norm, so this is a second-order cone program, not a QP:
#' \code{regularization} multiplies a norm, never a squared norm, and there is
#' no smoothness operator (\eqn{L} is the identity) and no sum or dose
#' normalisation constraint.
#'
#' The program is reproduced here with the package's own active-set engine,
#' solved to (near) exact optimality.  For \code{norm = 2} the KKT conditions
#' of the cone program are exactly the KKT conditions of the Tikhonov problem
#' \eqn{\min_{x \ge 0} \|A x - b\|^2 + \lambda \|x\|^2}{min ||Ax-b||^2 + lambda ||x||^2}
#' with \eqn{\lambda = \alpha \|A x - b\|_2 / \|x\|_2}{lambda = alpha*||Ax-b||/||x||};
#' that scalar equation has a single root (\eqn{\lambda \mapsto \alpha g/s}
#' crosses the diagonal once), so it is solved by bracketing plus bisection
#' followed by a Newton polish on the fixed support.  When an exact
#' non-negative fit exists the multiplier vanishes and the cone optimum is the
#' minimum-norm exact fit, which the engine selects from a multiplier ladder
#' scored by the true cone objective.  For \code{norm = 1} the non-negativity
#' makes \eqn{\|x\|_1}{||x||_1} the linear form \eqn{\mathbf{1}'x}{1'x}, and the
#' analogous stationarity identity reads
#' \eqn{\mu = 2 \alpha \|A x(\mu) - b\|_2}{mu = 2*alpha*||Ax(mu)-b||} with
#' \eqn{x(\mu)}{x(mu)} solving \eqn{\min_{x \ge 0} \|A x - b\|^2 - \mu \mathbf{1}'x}{%
#' min ||Ax-b||^2 - mu 1'x}; \eqn{\mu \mapsto 2 \alpha g(\mu)} is non-decreasing
#' and bounded by \eqn{2 \alpha \|b\|}{2*alpha*||b||}, so the monotone iteration
#' from \eqn{\mu = 0}{mu = 0} converges to the same point.
#'
#' Because the reference environment has no ECOS and CVXPY therefore falls back
#' to SCS (a first-order splitting solver) at its default tolerances, the
#' reference spectra can be measurably short of optimality on badly scaled
#' grids; this engine returns the exact minimiser of the same program, so on
#' those grids it differs from the reference by the reference's own solver
#' error rather than by a different formulation.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused (accepted for the unified solver API and Python parity).
#' @param alpha Numeric; regularization parameter multiplying the norm of the
#'   spectrum.  Default 1e-4, the reference default.
#' @param norm Integer; 1 (L1) or 2 (L2). Default 2.
#' @param solver Character; CVXPY solver name.  Accepted for API parity: the
#'   port uses its own internal conic engine, so the name only selects which
#'   label is reported and any recognised CVXPY solver name is accepted.
#'   Default \code{"ECOS"} (the reference default).
#' @param ub Optional numeric per-bin upper bound of length \code{n}.
#'   Bins with \code{ub == 0} are forced to zero (the
#'   \code{max_neutron_energy} cutoff), finite positive entries are honoured as
#'   upper bounds and the returned spectrum is clipped to them, exactly like
#'   \code{_solve_cvxpy_problem}.  Default \code{NULL} = unbounded.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_cvxpy(A, b, NULL, alpha = 1e-4, norm = 2)
solve_cvxpy <- function(A, b, x0 = NULL, alpha = 1e-4, norm = 2L,
                        solver = "ECOS", ub = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    p <- if (is.null(norm)) 2L else as.integer(norm)
    if (!p %in% c(1L, 2L)) {
        stop("solve_cvxpy: `norm` must be 1 or 2, got ", norm, call. = FALSE)
    }
    alpha <- if (is.null(alpha) || !is.finite(alpha)) 0 else as.numeric(alpha)
    if (alpha < 0) alpha <- 0
    if (!is.null(solver) && !identical(solver, "default")) {
        known <- c("ECOS", "SCS", "CLARABEL", "OSQP", "MOSEK", "GUROBI",
                   "SCIPY", "HiGHS", "diffcp", "glpk", "glpk_mi", "cplex")
        if (!solver %in% known) {
            warning("solve_cvxpy: unknown CVXPY solver '", solver,
                    "'; the internal conic engine is used.", call. = FALSE)
        }
    }
    ## Upper bounds: only finite entries constrain (Python honours
    ## `x[finite] <= ub[finite]` and post-clips `ub == 0` to exact zero).
    box <- NULL
    if (!is.null(ub)) {
        ubv <- as.numeric(ub)
        if (length(ubv) != n) {
            stop("solve_cvxpy: `ub` must have length ", n,
                 ", got ", length(ubv), call. = FALSE)
        }
        box <- ifelse(is.finite(ubv), pmax(ubv, 0), NA_real_)
    }
    spec <- if (p == 1L) .cvxpy_socp_l1(A, b, alpha, box)
            else .cvxpy_socp_l2(A, b, alpha, box)
    spec <- pmax(as.numeric(spec), 0)
    if (!is.null(box)) {
        hit <- which(!is.na(box))
        spec[hit] <- pmin(spec[hit], box[hit])
    }
    list(spectrum = spec, iterations = 1L, converged = TRUE)
}

## ---- internal engines for the CVXPY cone program (unfold_cvxpy) -----------

## Plain non-negative least squares on the augmented (Tikhonov) system:
##   argmin_{0 <= x <= box} ||A x - b||^2 + lambda ||x||^2
.cvx_nnls_tikh <- function(A, b, lambda, box = NULL) {
    n <- ncol(A)
    G <- crossprod(A)
    h <- as.numeric(crossprod(A, as.numeric(b)))
    x <- tryCatch(as.numeric(lsei::nnls(rbind(A, if (lambda > 0)
                                               sqrt(lambda) * diag(n)
                                               else matrix(0, 0, n)),
                                        c(b, numeric(if (lambda > 0) n else 0L))
                                        )$x),
                  error = function(e) NULL)
    if (is.null(x) || !all(is.finite(x)) || !is.null(box)) {
        x <- .cvx_nnls_gram(G, h, box = box, lambda = lambda)
    }
    x
}

## Symmetric positive (semi-)definite solve.  With a strictly positive definite
## block (any Tikhonov multiplier lambda > 0) the Cholesky factorisation is
## used, which keeps the small eigen-directions of A'A that the minimum-norm
## exact fit lives in; the truncated eigen solver is the fallback for the
## singular lambda = 0 systems (m < n_energy_bins makes A'A singular there).
.psolve <- function(M, v) {
    decomp <- eigen(M, symmetric = TRUE)
    d <- decomp$values
    if (!length(d) || all(d <= 0)) return(rep(0, length(v)))
    tol <- max(dim(M)) * max(d) * .Machine$double.eps * 8
    if (min(d) > tol) {
        cf <- tryCatch(chol(M), error = function(e) NULL)
        if (!is.null(cf)) {
            sol <- tryCatch(as.numeric(backsolve(cf, forwardsolve(t(cf), v))),
                            error = function(e) NULL)
            if (!is.null(sol) && all(is.finite(sol))) return(sol)
        }
    }
    keep <- which(d > tol)
    if (!length(keep)) return(rep(0, length(v)))
    V <- decomp$vectors[, keep, drop = FALSE]
    as.numeric(V %*% ((t(V) %*% v) / d[keep]))
}

## Lawson-Hanson active set for
##   argmin_{0 <= x <= box} x'(G + lambda I)x - 2 h' x   (G psd).
## Works straight from the Gram pair, which is what a linear penalty needs.
## `box = NULL` (the reference case) keeps the classical lower-bound-only walk.
.cvx_nnls_gram <- function(G, h, tol = 1e-13, box = NULL, lambda = 0) {
    n <- length(h)
    H <- if (lambda > 0) G + lambda * diag(n) else G
    ub <- if (is.null(box)) rep(Inf, n)
          else ifelse(is.na(box), Inf, pmax(as.numeric(box), 0))
    has_ub <- any(is.finite(ub))
    ## state: -1 at the lower bound, 0 free, 1 at a finite upper bound
    state <- rep(-1L, n)
    free <- rep(FALSE, n)
    x <- numeric(n)
    hs <- max(abs(h))
    if (!is.finite(hs) || hs == 0) return(x)
    scale <- hs + 1e-300
    for (it in seq_len(6L * n + 60L)) {
        w <- h - as.numeric(H %*% x)
        w[free] <- NA_real_
        j <- which.max(w)                      # best improving lower-bound var
        gain <- if (length(j) == 0L || is.na(w[j])) -Inf else w[j]
        if (has_ub) {
            up <- which(state == 1L)
            if (length(up)) {
                k <- up[which.min(w[up])]
                if (w[k] < -tol * scale && -w[k] > gain) {
                    free[k] <- TRUE; state[k] <- 0L; gain <- -w[k]
                }
            }
        }
        if (gain <= tol * scale) break
        if (length(j) == 1L && !is.na(w[j]) && w[j] > tol * scale) {
            free[j] <- TRUE; state[j] <- -1L
        }
        for (inner in seq_len(6L * n + 60L)) {
            idx <- which(free)
            if (!length(idx)) break
            s <- .psolve(H[idx, idx, drop = FALSE], h[idx])
            lo_ok <- s > 0
            hi_ok <- s <= ub[idx]
            if (all(lo_ok & hi_ok)) {
                x[idx] <- s
                break
            }
            blocked <- which(!lo_ok | !hi_ok)
            down <- !lo_ok[blocked]                 # heading below zero
            ratio <- ifelse(down,
                            ifelse(x[idx[blocked]] > 0,
                                   x[idx[blocked]] /
                                       (x[idx[blocked]] - s[blocked]), Inf),
                            ifelse(ub[idx[blocked]] > x[idx[blocked]],
                                   (ub[idx[blocked]] - x[idx[blocked]]) /
                                       (s[blocked] - x[idx[blocked]]), Inf))
            a <- min(ratio)
            if (!is.finite(a) || a <= 0) {
                x[idx] <- pmin(pmax(s, 0), ub[idx])
                break
            }
            x[idx] <- x[idx] + a * (s - x[idx])
            gone_lo <- idx[down & abs(x[idx][down]) <= tol * scale]
            gone_hi <- idx[!down & has_ub &
                               abs(x[idx][!down] - ub[idx][!down]) <=
                               tol * scale]
            if (length(gone_lo)) {
                x[gone_lo] <- 0; free[gone_lo] <- FALSE; state[gone_lo] <- -1L
            }
            if (length(gone_hi)) {
                x[gone_hi] <- ub[gone_hi]
                free[gone_hi] <- FALSE; state[gone_hi] <- 1L
            }
            if (!any(free)) break
        }
        if (has_ub) {
            fin <- which(is.finite(ub) & !free)
            x[fin] <- ifelse(state[fin] == 1L, ub[fin], 0)
        }
    }
    idx <- which(free)
    if (length(idx)) x[idx] <- pmin(pmax(x[idx], 0), ub[idx])
    pmax(x, 0)
}

## Minimum-norm exact non-negative fit:
##   argmin ||x||_2  s.t.  A x = b,  0 <= x <= box
## Primal-dual active set (the QP is strictly convex, so the KKT system below
## characterises the unique optimum): with multiplier y,
##   x_j = (A_j' y) > 0 on the support,  (A_j' y) <= 0 off it,  A_S x_S = b.
## This is the cone optimum whenever an exact non-negative fit exists (the
## multi-dimensional grids of this package almost always do), because then the
## data term vanishes and only the penalty is left to minimise.
.cvx_min_norm_exact <- function(A, b, box = NULL) {
    n <- ncol(A)
    m <- nrow(A)
    keep <- rep(TRUE, n)
    x <- numeric(n)
    fixed <- rep(FALSE, n)          # variables pinned at a positive upper bound
    bb <- b
    if (!is.null(box)) {
        zero <- which(is.finite(box) & box <= 0)
        if (length(zero)) keep[zero] <- FALSE
    }
    for (it in seq_len(8L * n + 60L)) {
        idx <- which(keep)
        if (!length(idx)) return(NULL)
        Aa <- A[, idx, drop = FALSE]
        y <- tryCatch(as.numeric(.psolve(Aa %*% t(Aa), bb)),
                      error = function(e) NULL)
        if (is.null(y) || !all(is.finite(y))) return(NULL)
        xs <- as.numeric(crossprod(Aa, y))
        neg <- which(xs <= 0)
        if (length(neg)) {                    # drop the most negative variable
            keep[idx[neg[which.min(xs[neg])]]] <- FALSE
            next
        }
        w <- as.numeric(crossprod(A, y))       # value an excluded column would take
        add <- which(!keep & !fixed & w > 0)
        if (length(add)) {
            keep[add[which.max(w[add])]] <- TRUE
            next
        }
        cand <- numeric(n)
        cand[idx] <- xs
        if (!is.null(box)) {
            over <- which(is.finite(box) & cand > box + 1e-12 * max(1, max(box)))
            if (length(over)) {                # pin at the bound and re-fit
                j <- over[which.max(cand[over] - box[over])]
                cand <- numeric(n); cand[j] <- box[j]
                fixed[j] <- TRUE; keep[j] <- FALSE
                bb <- bb - as.numeric(A[, j]) * box[j]
                next
            }
        }
        if (.cvx_g(A, b, cand) > 1e-6 * norm(b, "2")) return(NULL)
        return(cand)
    }
    NULL
}

## ||A x - b||_2 and ||x||_2 helpers
.cvx_g <- function(A, b, x) norm(as.numeric(A %*% x) - b, "2")

.cvxpy_socp_l2 <- function(A, b, alpha, box = NULL) {
    n <- ncol(A)
    b <- as.numeric(b)
    nb <- norm(b, "2")
    if (nb == 0) return(numeric(n))
    if (alpha == 0) return(.cvx_nnls_tikh(A, b, 0, box))
    ## x = 0 is optimal once alpha reaches the directional derivative bound
    ## ||(A'b)_+||_2 / ||b||_2 (a converged ECOS/SCS/CLARABEL see the same cut).
    grad0 <- as.numeric(crossprod(A, b))
    if (alpha >= norm(pmax(grad0, 0), "2") / nb) return(numeric(n))

    G <- crossprod(A)
    h <- as.numeric(crossprod(A, b))
    gG <- max(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
    sol_of <- function(lambda) .cvx_nnls_gram(G, h, box = box, lambda = lambda)
    obj_of <- function(v) .cvx_g(A, b, v) + alpha * norm(v, "2")
    ## sign of  lambda*s(lambda) - alpha*g(lambda);  positive above the root
    resid_of <- function(lambda) {
        x <- sol_of(lambda)
        s <- norm(x, "2")
        lambda * s - alpha * .cvx_g(A, b, x)
    }
    ## Degenerate case: a non-negative exact fit exists, the multiplier vanishes
    ## and the cone optimum is the minimum-norm exact fit.
    g0 <- .cvx_g(A, b, sol_of(0))
    degenerate <- g0 <= 1e-10 * nb

    lo <- 0
    hi <- max(alpha * nb, .Machine$double.xmin^0.25)
    guard <- 0L
    while (resid_of(hi) < 0 && guard < 200L) { hi <- hi * 10; guard <- guard + 1L }
    if (resid_of(hi) < 0) return(numeric(n))
    lam_root <- hi
    if (!degenerate) {
        for (it in seq_len(120L)) {
            mid <- 0.5 * (lo + hi)
            if (mid <= lo || mid >= hi) break
            if (resid_of(mid) < 0) lo <- mid else hi <- mid
            if (hi - lo <= 1e-15 * max(1, hi)) break
        }
        lam_root <- 0.5 * (lo + hi)
        ## Newton polish on the fixed support: with S = supp(x) the map
        ## F(lam) = lam*s - alpha*g is smooth and
        ##   dF/dlam = s - lam*q*(1/s + alpha/g),  q = x'(G_SS + lam I)^-1 x_S
        ## (both derivatives follow from dx/dlam = -(G_SS + lam I)^-1 x_S and
        ## the stationarity identity G_SS x_S - h_S = -lam x_S).
        x <- sol_of(lam_root)
        sup <- which(x > 0)
        if (length(sup) > 0 && length(sup) < n) {
            Gss <- G[sup, sup, drop = FALSE]
            hs <- h[sup]
            for (k in seq_len(12L)) {
                z <- .psolve(Gss + lam_root * diag(length(sup)), hs)
                if (any(z <= 0)) break
                q <- as.numeric(crossprod(z, .psolve(Gss + lam_root *
                                                     diag(length(sup)), z)))
                xtmp <- numeric(n); xtmp[sup] <- z
                s <- norm(xtmp, "2"); g <- .cvx_g(A, b, xtmp)
                F <- lam_root * s - alpha * g
                dF <- s - lam_root * q * (1 / s + alpha / max(g, 1e-300))
                if (!is.finite(F) || !is.finite(dF) || dF <= 0) break
                step <- F / dF
                cand <- lam_root - step
                if (!is.finite(cand) || cand <= 0) break
                if (abs(step) > 0.5 * lam_root) cand <- 0.5 * lam_root
                if (abs(cand - lam_root) <= 1e-16 * lam_root) {
                    lam_root <- cand; break
                }
                lam_root <- cand
            }
        }
    }
    x <- sol_of(lam_root)
    ## Score the vanishing-multiplier family too: on an exact fit only the
    ## penalty breaks the tie, so a ladder of absolute multipliers (anchored on
    ## the scale of A'A, which is where the degeneracy bites) is cheap insurance
    ## and can only lower the reported objective.
    cands <- list(x, sol_of(0))
    for (f in 10^(-seq(-2L, 16L))) cands[[length(cands) + 1L]] <- sol_of(lam_root * f)
    for (f in 10^(-seq_len(16L))) cands[[length(cands) + 1L]] <- sol_of(hi * f)
    if (is.finite(gG) && gG > 0) {
        for (k in seq_len(14L)) {
            cands[[length(cands) + 1L]] <- sol_of(gG * 10^(-k))
        }
    }
    if (!is.null(box)) {
        hit <- which(is.finite(ifelse(is.na(box), Inf, box)))
        for (i in seq_along(cands)) {
            cands[[i]][hit] <- pmin(cands[[i]][hit],
                                    ifelse(is.na(box[hit]), Inf, box[hit]))
        }
    }
    objs <- vapply(cands, obj_of, numeric(1))
    cands[[which.min(objs)]]
}

.cvxpy_socp_l1 <- function(A, b, alpha, box = NULL) {
    n <- ncol(A)
    b <- as.numeric(b)
    nb <- norm(b, "2")
    if (nb == 0) return(numeric(n))
    if (alpha == 0) return(.cvx_nnls_tikh(A, b, 0, box))
    ## x = 0 optimal iff alpha >= max_i (A'b)_i / ||b||_2 (unit direction e_i).
    grad0 <- as.numeric(crossprod(A, b))
    if (alpha >= max(grad0) / nb) return(numeric(n))
    G <- crossprod(A)
    h0 <- as.numeric(crossprod(A, b))
    gG <- max(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
    scale <- max(abs(h0))
    if (!is.finite(scale) || scale == 0) scale <- 1
    ## In the flat part of the trade-off curve the linear penalty is the only
    ## thing that breaks the tie between equivalent fits, so the active-set
    ## tolerance must stay well below its magnitude.
    gram_tol <- function(mu) max(4 * .Machine$double.eps,
                                min(1e-13, 1e-8 * 0.5 * mu / scale))
    xof <- function(mu) .cvx_nnls_gram(G, h0 - 0.5 * mu, tol = gram_tol(mu),
                                       box = box)
    muf <- function(mu) 2 * alpha * .cvx_g(A, b, xof(mu))
    obj_of <- function(v) .cvx_g(A, b, v) + alpha * sum(v)
    ## mu -> 2*alpha*g(mu) is non-decreasing while the left hand side grows
    ## linearly, so mu - 2*alpha*g(mu) changes sign exactly once.
    lo <- 0
    hi <- max(2 * alpha * nb, .Machine$double.xmin^0.25)
    guard <- 0L
    while (muf(hi) > hi && guard < 200L) { hi <- hi * 10; guard <- guard + 1L }
    if (muf(hi) > hi) hi <- lo
    for (it in seq_len(120L)) {
        mid <- 0.5 * (lo + hi)
        if (mid <= lo || mid >= hi) break
        if (muf(mid) > mid) lo <- mid else hi <- mid
        if (hi - lo <= 1e-15 * max(1, hi)) break
    }
    mu_root <- 0.5 * (lo + hi)
    x <- xof(mu_root)
    sup <- which(x > 0)
    if (length(sup) > 0) {
        Gss <- G[sup, sup, drop = FALSE]
        mu <- mu_root
        for (it in seq_len(60L)) {
            xs <- .psolve(Gss, h0[sup] - 0.5 * mu)
            if (any(xs <= 0)) break
            xtmp <- numeric(n); xtmp[sup] <- xs
            mu_new <- 2 * alpha * .cvx_g(A, b, xtmp)
            if (!is.finite(mu_new) || mu_new <= 0) break
            if (abs(mu_new - mu) <= 1e-16 * max(1, mu_new)) { mu <- mu_new; break }
            mu <- 0.5 * mu + 0.5 * mu_new   # damped, keeps the iteration stable
        }
        xs <- .psolve(Gss, h0[sup] - 0.5 * mu)
        if (all(xs > 0)) {
            x_new <- numeric(n); x_new[sup] <- xs
            if (obj_of(x_new) < obj_of(x)) x <- x_new
        }
    }
    ## Same vanishing-multiplier guard as the L2 branch: when an exact
    ## non-negative fit exists the L1 root collapses to zero and any small
    ## positive penalty selects the minimum-1'x fit, so score a wide ladder of
    ## multipliers with the true cone objective and keep the best.
    cands <- list(x, xof(0))
    for (f in 10^(-seq(-2L, 16L))) cands[[length(cands) + 1L]] <- xof(mu_root * f)
    for (f in 10^(-seq_len(16L))) cands[[length(cands) + 1L]] <- xof(hi * f)
    if (is.finite(gG) && gG > 0) {
        for (k in seq_len(14L)) cands[[length(cands) + 1L]] <- xof(gG * 10^(-k))
    }
    if (!is.null(box)) {
        hit <- which(is.finite(ifelse(is.na(box), Inf, box)))
        for (i in seq_along(cands)) {
            cands[[i]][hit] <- pmin(cands[[i]][hit],
                                    ifelse(is.na(box[hit]), Inf, box[hit]))
        }
    }
    objs <- vapply(cands, obj_of, numeric(1))
    cands[[which.min(objs)]]
}

#' Resolve the CVXPY regularization parameter
#'
#' Port of the \code{regularization_method} branch of
#' \code{bssunfold.core.unfold_cvxpy}: \code{"manual"} keeps the supplied
#' \code{regularization}, \code{"cosine"} requires \code{initial_spectrum}, and
#' the automatic selectors evaluate the documented fallback criteria of
#' \code{bssunfold.core.regularization} (\code{pytikhonov}, which the reference
#' environment lacks, is not required here).
#' @keywords internal
.cvxpy_select_regularization <- function(A, b, method, regularization,
                                         noise_var = NULL,
                                         initial_spectrum = NULL) {
    if (identical(method, "manual")) return(as.numeric(regularization))
    if (identical(method, "cosine")) {
        if (is.null(initial_spectrum)) {
            stop("For 'cosine' method, initial_spectrum must be provided.",
                 call. = FALSE)
        }
        return(.cvxpy_cosine_selection(A, b, initial_spectrum))
    }
    if (identical(method, "lcurve")) return(.cvxpy_lcurve_selection(A, b))
    if (identical(method, "gcv")) return(.cvxpy_gcv_selection(A, b))
    if (identical(method, "dp")) {
        if (is.null(noise_var)) noise_var <- .cvxpy_noise_variance(A, b)
        return(.cvxpy_dp_selection(A, b, noise_var))
    }
    stop("Unknown regularization selection method: ", method,
         ". Choose from 'manual', 'cosine', 'lcurve', 'gcv', 'dp'.",
         call. = FALSE)
}

## Shared alpha grid of the fallback selectors: np.logspace(-9, 2, n_alphas).
.cvxpy_alpha_grid <- function(n_alphas = 50L,
                              alpha_range = c(1e-9, 1e2)) {
    10^seq(log10(alpha_range[1]), log10(alpha_range[2]),
           length.out = as.integer(n_alphas))
}

## Unconstrained Tikhonov solve with a non-negativity clip, the oracle used by
## every fallback selector (np.linalg.solve then np.maximum(x, 0)).
.cvxpy_tikh_clip <- function(A, b, alpha) {
    n <- ncol(A)
    x <- tryCatch(as.numeric(.psolve(crossprod(A) + alpha * diag(n),
                                     as.numeric(crossprod(A, as.numeric(b))))),
                  error = function(e) NULL)
    if (is.null(x) || !all(is.finite(x))) return(NULL)
    pmax(x, 0)
}

#' @keywords internal
.cvxpy_lcurve_selection <- function(A, b, n_alphas = 50L,
                                    alpha_range = c(1e-9, 1e2)) {
    alphas <- .cvxpy_alpha_grid(n_alphas, alpha_range)
    res <- numeric(0); nrm <- numeric(0)
    for (alpha in alphas) {
        x <- .cvxpy_tikh_clip(A, b, alpha)
        if (is.null(x)) next
        res <- c(res, .cvx_g(A, b, x))
        nrm <- c(nrm, norm(x, "2"))           # L = identity
    }
    if (length(res) < 3) return(1)
    lr <- log(res); ln <- log(nrm)
    dvec <- c(lr[length(lr)] - lr[1], ln[length(ln)] - ln[1])
    den <- sqrt(sum(dvec^2))
    if (!is.finite(den) || den == 0) return(as.numeric(alphas[1]))
    ## Distance of each log-log point from the chord joining the endpoints; the
    ## reference uses np.cross on 2-vectors, which NumPy >= 2 rejects, so the
    ## scalar 2-D cross product is evaluated directly.
    dist <- abs(dvec[1] * (lr[1] - lr) - dvec[2] * (ln[1] - ln)) / den
    as.numeric(alphas[which.max(dist)])
}

#' @keywords internal
.cvxpy_gcv_selection <- function(A, b, n_alphas = 50L,
                                 alpha_range = c(1e-9, 1e2)) {
    alphas <- .cvxpy_alpha_grid(n_alphas, alpha_range)
    m <- nrow(A)
    sv <- compute_svd_components(A)
    s_sq <- sv$s_sq
    UTb <- as.numeric(t(sv$U) %*% as.numeric(b))
    gcv <- vapply(alphas, function(alpha) {
        rc <- alpha / (s_sq + alpha)
        den <- (m - sum(s_sq / (s_sq + alpha)))^2
        if (den == 0) Inf else sum((rc * UTb)^2) / den
    }, numeric(1))
    if (!length(gcv) || !any(is.finite(gcv))) return(1)
    as.numeric(alphas[which.min(gcv)])
}

#' @keywords internal
.cvxpy_dp_selection <- function(A, b, noise_var, n_alphas = 50L,
                                alpha_range = c(1e-9, 1e2)) {
    alphas <- .cvxpy_alpha_grid(n_alphas, alpha_range)
    target <- sqrt(noise_var) * sqrt(length(b))
    res <- vapply(alphas, function(alpha) {
        x <- .cvxpy_tikh_clip(A, b, alpha)
        if (is.null(x)) Inf else .cvx_g(A, b, x)
    }, numeric(1))
    as.numeric(alphas[which.min(abs(res - target))])
}

#' @keywords internal
.cvxpy_cosine_selection <- function(A, b, initial_spectrum, n_alphas = 100L,
                                    alpha_range = c(-9, 2)) {
    ref <- pmax(as.numeric(initial_spectrum), 0)
    nrm <- norm(ref, "2")
    if (nrm == 0) stop("Initial spectrum has zero norm.", call. = FALSE)
    ref <- ref / nrm
    alphas <- 10^seq(alpha_range[1], alpha_range[2],
                     length.out = as.integer(n_alphas))
    sv <- compute_svd_components(A)
    s_sq <- sv$s_sq
    UTb <- as.numeric(t(sv$U) %*% as.numeric(b))
    sim <- vapply(alphas, function(alpha) {
        filt <- sv$s / (s_sq + alpha)
        x <- pmax(as.numeric(t(sv$Vt) %*% (filt * UTb)), 0)
        nx <- norm(x, "2")
        if (nx == 0) 0 else as.numeric(crossprod(x, ref)) / nx
    }, numeric(1))
    as.numeric(alphas[which.max(sim)])
}

#' @keywords internal
.cvxpy_noise_variance <- function(A, b) {
    x <- as.numeric(qr.solve(as.matrix(A), as.numeric(b)))
    r <- as.numeric(b) - as.numeric(A %*% x)
    stats::var(r)
}

#' Wrapper around \code{\link{solve_cvxpy}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_cvxpy
#' @param regularization Numeric; the CVXPY regularization parameter passed to
#'   \code{\link{solve_cvxpy}} as \code{alpha}.  Default 1e-4.
#' @param regularization_method Character; \code{"manual"} (default) uses
#'   \code{regularization}, \code{"cosine"} needs \code{initial_spectrum}, and
#'   \code{"lcurve"}, \code{"gcv"} or \code{"dp"} select the parameter from the
#'   response system.
#' @param noise_var Optional noise variance for \code{regularization_method =
#'   "dp"}; estimated from the least-squares residual when \code{NULL}.
#' @param ln_steps Optional numeric \eqn{d(\ln E)} bin widths.  Accepted for
#'   Python API parity; the sensitivities held by \code{Detector} are already
#'   lethargy-weighted, so this method does not consume it.
#' @param reading_uncertainties Optional per-detector reading uncertainties.
#'   Accepted for Python API parity.
#' @param reading_covariance Optional reading covariance matrix.  Accepted for
#'   Python API parity.
#' @param noise_model Character noise model name.  Accepted for Python API
#'   parity.  Default \code{"gaussian"}.
#' @param measurement_time Optional measurement time in seconds.  Accepted for
#'   Python API parity.
#' @param solver Character; CVXPY solver name, \code{"default"} resolves to the
#'   reference preference order.  The port always uses its own internal conic
#'   engine, so this only labels the reported field.
#' @export
unfold_cvxpy <- function(detector_names, n_energy_bins, E_MeV,
                         sensitivities, cc_icrp116, save_result_callback,
                         readings, ln_steps = NULL, initial_spectrum = NULL,
                         regularization = 1e-4, norm = 2L,
                         solver = "default",
                         calculate_errors = FALSE,
                         noise_level = 0.01, n_montecarlo = 100L,
                         save_result = FALSE,
                         regularization_method = "manual",
                         noise_var = NULL, random_state = NULL,
                         max_neutron_energy = NULL,
                         reading_uncertainties = NULL,
                         reading_covariance = NULL,
                         noise_model = "gaussian",
                         measurement_time = NULL) {
    ## Python: solver == "default" picks the first installed conic solver of
    ## ECOS, SCS, CLARABEL.  None of them ships with this package, so the
    ## internal conic engine stands in for the resolved name.
    if (identical(solver, "default")) solver <- "CLARABEL"
    if (!regularization_method %in% c("manual", "cosine", "lcurve", "gcv",
                                      "dp")) {
        stop("Unknown regularization selection method: ", regularization_method,
             ". Choose from 'manual', 'cosine', 'lcurve', 'gcv', 'dp'.",
             call. = FALSE)
    }
    sys <- .build_system(readings, detector_names, sensitivities)
    alpha <- .cvxpy_select_regularization(
        sys$A, sys$b, regularization_method, regularization,
        noise_var = noise_var, initial_spectrum = initial_spectrum)
    x0_default <- numeric(as.integer(n_energy_bins))
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_cvxpy,
                                         alpha = alpha, norm = norm,
                                         solver = solver),
        solve_kwargs = list(),
        method_name = "cvxpy",
        extra_output = list(norm = as.integer(norm), solver = solver,
                            regularization_method = regularization_method,
                            selected_regularization = as.numeric(alpha)),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
