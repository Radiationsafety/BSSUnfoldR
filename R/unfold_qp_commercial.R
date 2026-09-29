#' Shared numerical backend for the commercial-license QP solvers
#'
#' R port of \code{bssunfold/src/bssunfold/core/_commercial_qp.py} plus the
#' solver-independent half of \code{core/unfold_commercial.py} and
#' \code{core/regularization.py}.
#'
#' @details
#' The Python package exposes five parallel unfolding methods
#' (\code{unfold_gurobi}, \code{unfold_mosek}, \code{unfold_cplex},
#' \code{unfold_copt}, \code{unfold_xpress}) that are \code{functools.partial}
#' bindings of one single implementation, \code{unfold_commercial}, differing
#' only in the \code{solver} alias.  Each of them hands the \emph{same}
#' convex quadratic program to the cvxpy interface of the licensed engine,
#' with \strong{no open-source fallback}: if the engine or its license is
#' missing the solve returns \code{None} and the wrapper warns and yields a
#' zero spectrum.  In the reference environment all five therefore return an
#' all-zero spectrum and there is nothing to compare against.
#'
#' This module \emph{ports the optimisation problem itself} and solves it with
#' a pure-R active-set quadratic-programming engine, so the five methods
#' produce a genuine solution instead of the zero-spectrum placeholder.  The
#' problem is written exactly as in \code{_commercial_qp.py}:
#'
#' \deqn{\min_x \; \tfrac12 x^T P x + q^T x
#'        \quad \text{s.t.}\; x \ge 0 \;(\text{when } \code{nonneg}),\;
#'        x \le ub}
#'
#' with the library's canonical QP convention
#' \itemize{
#'   \item \eqn{P = A^T A}{P = A'A}, \eqn{q = -A^T b}{q = -A'b};
#'   \item \code{norm = 2}: \eqn{P \mathrel{+}= \alpha I} for
#'     \code{smoothness_order = 0}, else
#'     \eqn{P \mathrel{+}= \alpha\,w\,L^T L}{P += alpha*w*L'L} with \eqn{L} the
#'     \code{create_derivative_matrix} operator
#'     (\code{build_smoothness_penalty});
#'   \item \code{norm = 1}: \eqn{q \mathrel{+}= \alpha \mathbf{1}}{q += alpha*1}
#'     and, for \code{smoothness_order} 1 or 2, the same
#'     \eqn{\alpha\,w\,L^T L}{alpha*w*L'L} term on \eqn{P}.  The L1 penalty
#'     equals \eqn{\alpha\sum_j x_j} only under non-negativity, so
#'     \code{norm = 1} with \code{nonneg = FALSE} is rejected exactly as
#'     Python rejects it;
#'   \item \eqn{P} is symmetrised as \eqn{\tfrac12(P + P^T)} before use.
#' }
#'
#' Because \eqn{P} is positive semidefinite by construction, any point that
#' satisfies the Karush-Kuhn-Tucker conditions is a \emph{global} minimiser:
#' stationarity on the free set, a non-negative gradient on variables pinned
#' at their lower bound and a non-positive gradient on variables pinned at an
#' upper bound.  The engine below therefore terminates on a KKT certificate
#' rather than on an iteration count, and \code{timeout} is only an advisory
#' wall-clock guard (the licensed engines use it as a true time limit).
#'
#' @name commercial-qp-backend
NULL

# ---------------------------------------------------------------------------
# backend metadata (COMMERCIAL_SOLVER_ALIASES / _PY_MODULES)
# ---------------------------------------------------------------------------

.COMMERCIAL_SOLVER_ALIASES <- c(
    gurobi = "GUROBI", mosek = "MOSEK", cplex = "CPLEX",
    copt = "COPT", xpress = "XPRESS"
)

.COMMERCIAL_SOLVER_PY_MODULES <- c(
    gurobi = "gurobipy", mosek = "mosek", cplex = "cplex",
    copt = "coptpy", xpress = "xpress"
)

#' License / backend metadata for one commercial solver alias
#'
#' Internal helper mirroring \code{commercial_solver_info}.  \code{available}
#' is \code{FALSE} for all five engines here: no proprietary engine or license
#' is present, and the R port never probes for one -- the problem is solved by
#' the internal QP engine instead.
#'
#' @param alias Character; one of \code{"gurobi"}, \code{"mosek"},
#'   \code{"cplex"}, \code{"copt"}, \code{"xpress"}.
#' @return A named list with \code{alias}, \code{cvxpy_solver},
#'   \code{pip_package}, \code{license_required} and \code{available}.
#' @keywords internal
.commercial_solver_info <- function(alias) {
    if (!is.character(alias) || length(alias) != 1L ||
        !alias %in% names(.COMMERCIAL_SOLVER_ALIASES)) {
        stop("Unknown commercial solver '", alias, "'. Supported: ",
             paste(sort(names(.COMMERCIAL_SOLVER_ALIASES)), collapse = ", "),
             call. = FALSE)
    }
    list(
        alias = alias,
        cvxpy_solver = .COMMERCIAL_SOLVER_ALIASES[[alias]],
        pip_package = .COMMERCIAL_SOLVER_PY_MODULES[[alias]],
        license_required = TRUE,
        available = FALSE
    )
}

#' Per-bin upper bounds for a maximum-neutron-energy cutoff
#'
#' Mirrors \code{bssunfold.core._max_energy.upper_bounds}: bins with
#' \code{E_MeV > max_neutron_energy} get \code{0} (which, together with
#' \code{lb = 0}, forces them to zero), all other bins get \code{Inf}.
#'
#' @param E_MeV Numeric energy grid.
#' @param max_neutron_energy Numeric cutoff in MeV, or \code{NULL} to disable.
#' @return Numeric vector the same length as \code{E_MeV}.
#' @keywords internal
.commercial_upper_bounds <- function(E_MeV, max_neutron_energy) {
    E_MeV <- as.numeric(E_MeV)
    ub <- rep(Inf, length(E_MeV))
    if (!is.null(max_neutron_energy)) {
        cutoff <- as.numeric(max_neutron_energy)
        if (length(cutoff) == 1L && is.finite(cutoff)) {
            ub[E_MeV > cutoff] <- 0
        }
    }
    ub
}

# ---------------------------------------------------------------------------
# QP construction
# ---------------------------------------------------------------------------

#' Build the canonical QP data \eqn{P, q}{P, q}
#'
#' Mirrors the \code{P}/\code{q} block of \code{solve_commercial_qp},
#' including the \code{build_smoothness_penalty} semantics (order 0 adds the
#' identity term only for \code{norm = 2}; orders 1 and 2 add
#' \eqn{\alpha w L^T L}{alpha*w*L'L} for both norms).
#' @keywords internal
.commercial_qp_data <- function(A, b, alpha, norm, smoothness_order,
                                smoothness_weight) {
    n <- ncol(A)
    alpha <- as.numeric(alpha)
    if (is.finite(alpha) && alpha < 0) {
        warning("negative regularization alpha=", alpha,
                " is not a convex QP; clamped to 0.", call. = FALSE)
        alpha <- 0
    } else if (!is.finite(alpha)) {
        alpha <- 0
    }
    weight <- as.numeric(smoothness_weight)
    order <- as.integer(smoothness_order)
    # build_smoothness_penalty(): alpha * w * L'L for orders 1 and 2, else
    # NULL (the identity term is added separately and only for norm = 2).
    deriv <- if (order %in% c(1L, 2L)) {
        D <- if (n > order) {
            as.matrix(create_derivative_matrix(n, order))
        } else {
            matrix(0.0, nrow = 0L, ncol = n)   # empty operator, as in scipy
        }
        (alpha * weight) * crossprod(D)
    } else {
        NULL
    }
    P <- crossprod(A)
    q <- -as.numeric(crossprod(A, b))
    if (norm == 2L) {
        P <- P + (if (is.null(deriv)) alpha * diag(n) else deriv)
    } else {
        q <- q + alpha * rep(1, n)
        if (!is.null(deriv)) P <- P + deriv
    }
    P <- 0.5 * (P + t(P))
    storage.mode(P) <- "double"
    list(P = P, q = q)
}

# ---------------------------------------------------------------------------
# linear algebra helpers
# ---------------------------------------------------------------------------

# Symmetric (possibly semidefinite) linear solve.  Uses a Cholesky factor when
# it is numerically trustworthy and otherwise the eigen pseudo-inverse, which
# returns the minimum-norm solution.
.commercial_sym_solve <- function(M, rhs, rcond = 1e-11) {
    rhs <- as.numeric(rhs)
    if (length(rhs) == 0L) return(numeric(0L))
    out <- tryCatch({
        ct <- chol(M)
        dg <- abs(diag(ct))
        if (min(dg) <= 0 || max(dg) / min(dg) > 1e7) {
            NULL
        } else {
            as.numeric(chol2inv(ct) %*% rhs)
        }
    }, error = function(e) NULL)
    if (!is.null(out) && all(is.finite(out))) {
        res <- as.numeric(M %*% out) - rhs
        nrm <- max(sqrt(sum(res^2)) / max(sqrt(sum(rhs^2)), 1e-300), 0)
        if (nrm < 1e-7) return(out)
    }
    ev <- eigen(M, symmetric = TRUE, only.values = FALSE)
    d <- ev$values
    smax <- if (length(d)) max(abs(d)) else 0
    keep <- abs(d) > rcond * max(smax, 1e-300)
    z <- as.numeric(t(ev$vectors) %*% rhs)
    s <- numeric(length(d))
    s[keep] <- z[keep] / d[keep]
    as.numeric(ev$vectors %*% s)
}

# Projection onto the box [lo, hi].
.commercial_box_project <- function(x, lo, hi) pmin(pmax(x, lo), hi)

#' FISTA-style projected-gradient polish
#'
#' Safety net for the active-set method: accelerated projected gradient on the
#' box.  Only used when the active set cannot certify KKT optimality (for
#' example a badly conditioned semidefinite case).
#' @keywords internal
.commercial_pg_polish <- function(H, g, lo, hi, x, max_iterations = 20000L,
                                  deadline = Inf) {
    n <- length(x)
    L <- tryCatch(max(eigen(H, symmetric = TRUE, only.values = TRUE)$values),
                  error = function(e) max(abs(diag(H))) * n)
    L <- max(as.numeric(L), 1e-300)
    step <- 1 / L
    y <- x; tk <- 1
    for (it in seq_len(max_iterations)) {
        grad <- as.numeric(H %*% y) + g
        z <- .commercial_box_project(y - step * grad, lo, hi)
        tn <- (1 + sqrt(1 + 4 * tk^2)) / 2
        y <- z + ((tk - 1) / tn) * (z - x)
        x <- z
        tk <- tn
        if (it %% 500L == 0L && proc.time()[3] > deadline) break
    }
    .commercial_box_project(x, lo, hi)
}

#' Solve the box-regularised convex QP by a primal active-set method
#'
#' \deqn{\min \tfrac12 x^T H x + c^T x \quad \text{s.t. } lo \le x \le hi}
#' with \eqn{H} positive semidefinite.  Variables whose two bounds coincide
#' are eliminated exactly (this is how \code{ub = 0} above a
#' \code{max_neutron_energy} cutoff and the \code{nonneg} lower bound are
#' handled), the remaining ones are solved by the classical Lawson-Hanson /
#' Ibnd-("bounds") iteration: release the working-set member with the most
#' violating multiplier, otherwise take the exact Newton step on the free
#' block and shorten it to the first blocking bound.  Termination is a KKT
#' certificate, which is a global optimality certificate because \eqn{H} is
#' PSD.
#'
#' @return A list \code{list(x, iterations, converged, kkt)}.
#' @keywords internal
.commercial_box_qp <- function(H, cvec, lo, hi, x0 = NULL, tol = 1e-9,
                               max_iter = 4000L, deadline = Inf) {
    n <- length(cvec)
    lo <- as.numeric(lo); hi <- as.numeric(hi)
    storage.mode(H) <- "double"

    # --- eliminate exactly-fixed variables (lo == hi, both finite) ----------
    fixed <- is.finite(lo) & is.finite(hi) & (lo == hi)
    freev <- which(!fixed)
    xfull <- numeric(n)
    xfull[fixed] <- hi[fixed]
    if (length(freev) == 0L) {
        return(list(x = xfull, iterations = 0L, converged = TRUE,
                    kkt = 0))
    }
    Hf <- H[freev, freev, drop = FALSE]
    cf <- as.numeric(H[freev, fixed, drop = FALSE] %*% xfull[fixed]) +
        cvec[freev]
    lof <- lo[freev]; hif <- hi[freev]
    nr <- length(freev)

    # --- feasible start ------------------------------------------------------
    if (!is.null(x0) && length(x0) == n) {
        xs <- .commercial_box_project(as.numeric(x0)[freev], lof, hif)
    } else {
        xs <- .commercial_box_project(rep(0, nr), lof, hif)
    }
    hsnap <- max(1, max(abs(cvec[freev]), 0), max(abs(diag(Hf)), 0))
    atol <- tol * hsnap
    gtol <- tol * hsnap

    at_lo <- is.finite(lof) & (xs <= lof + atol)
    at_hi <- is.finite(hif) & (xs >= hif - atol)
    at_lo[at_hi] <- FALSE
    st <- ifelse(at_lo, -1L, ifelse(at_hi, 1L, 0L))
    xs[st == -1L] <- lof[st == -1L]
    xs[st ==  1L] <- hif[st ==  1L]

    converged <- FALSE
    it <- 0L
    for (it in seq_len(max_iter)) {
        if (proc.time()[3] > deadline) break
        grad <- as.numeric(Hf %*% xs) + cf

        # (1) multiplier / release test on the working set
        cand <- integer(0)
        cand_lo <- which(st == -1L & grad < -gtol)
        cand_hi <- which(st ==  1L & grad >  gtol)
        mlo <- if (length(cand_lo)) -min(grad[cand_lo]) else -Inf
        mhi <- if (length(cand_hi))  max(grad[cand_hi]) else -Inf
        if (is.finite(mlo) || is.finite(mhi)) {
            if (mlo >= mhi) {
                j <- cand_lo[which.min(grad[cand_lo])]
            } else {
                j <- cand_hi[which.max(grad[cand_hi])]
            }
            st[j] <- 0L
            next
        }

        # (2) exact Newton step on the free block
        F <- which(st == 0L)
        if (length(F) == 0L) { converged <- TRUE; break }
        d <- -.commercial_sym_solve(Hf[F, F, drop = FALSE], grad[F])
        if (!all(is.finite(d))) break
        nd <- sqrt(sum(d^2))
        if (nd <= gtol) { converged <- TRUE; break }

        # (3) longest feasible step
        tt <- rep(Inf, length(F))
        for (k in seq_along(F)) {
            j <- F[k]
            if (d[k] < 0 && is.finite(lof[j])) {
                tt[k] <- (lof[j] - xs[j]) / (-d[k])
            } else if (d[k] > 0 && is.finite(hif[j])) {
                tt[k] <- (hif[j] - xs[j]) / d[k]
            }
        }
        tmin <- min(tt[is.finite(tt)], 1)
        tmin <- max(tmin, 0)
        alpha <- min(1, tmin)
        xs[F] <- xs[F] + alpha * d
        if (alpha < 1) {
            blockers <- which(is.finite(tt) & tt <= alpha + atol * 1e-2)
            for (k in blockers) {
                j <- F[k]
                if (d[k] < 0) {
                    st[j] <- -1L; xs[j] <- lof[j]
                } else {
                    st[j] <-  1L; xs[j] <- hif[j]
                }
            }
        } else {
            # no bound touched: snap any value that landed on a bound
            hit_lo <- which(st == 0L & is.finite(lof) & xs <= lof + atol)
            hit_hi <- which(st == 0L & is.finite(hif) & xs >= hif - atol)
            if (length(hit_lo) + length(hit_hi) > 0L) {
                st[hit_lo] <- -1L; xs[hit_lo] <- lof[hit_lo]
                st[hit_hi] <-  1L; xs[hit_hi] <- hif[hit_hi]
            }
        }
    }

    # --- KKT certificate ----------------------------------------------------
    xfull[freev] <- xs
    grad_full <- as.numeric(H %*% xfull) + cvec
    bad_lo <- is.finite(lo) & (xfull <= lo + atol) & (grad_full < -gtol)
    bad_hi <- is.finite(hi) & (xfull >= hi - atol) & (grad_full >  gtol)
    interior <- !(bad_lo | bad_hi) &
        !((is.finite(lo) & xfull <= lo + atol) |
          (is.finite(hi) & xfull >= hi - atol))
    kkt_stat <- if (any(interior)) max(abs(grad_full[interior])) else 0
    kkt <- max(kkt_stat,
               if (any(bad_lo)) max(-grad_full[bad_lo]) else 0,
               if (any(bad_hi)) max(grad_full[bad_hi] ) else 0)
    if (!(converged && kkt <= 2 * gtol)) {
        xs2 <- .commercial_pg_polish(Hf, cf, lof, hif, xs,
                                     deadline = deadline)
        xfull[freev] <- xs2
        grad_full <- as.numeric(H %*% xfull) + cvec
        bad_lo <- is.finite(lo) & (xfull <= lo + atol) & (grad_full < -gtol)
        bad_hi <- is.finite(hi) & (xfull >= hi - atol) & (grad_full >  gtol)
        interior <- !(bad_lo | bad_hi) &
            !((is.finite(lo) & xfull <= lo + atol) |
              (is.finite(hi) & xfull >= hi - atol))
        kkt <- max(if (any(interior)) max(abs(grad_full[interior])) else 0,
                   if (any(bad_lo)) max(-grad_full[bad_lo]) else 0,
                   if (any(bad_hi)) max(grad_full[bad_hi]) else 0)
        converged <- kkt <= 2 * gtol
    }
    # mirror Python's `result[ub == 0] = 0` clean-up
    if (any(fixed)) xfull[fixed & (hi == 0)] <- 0
    list(x = xfull, iterations = as.integer(it),
         converged = converged, kkt = kkt)
}

# ---------------------------------------------------------------------------
# regularization-parameter selection (core/regularization.py)
# ---------------------------------------------------------------------------

#' Population noise variance of the least-squares residual
#' @keywords internal
.commercial_estimate_noise_variance <- function(A, b) {
    r <- .commercial_lstsq_residual(A, b)
    sum(r^2) / length(r)                      # numpy.var -> denominator n
}

# residual b - A x with x the minimum-norm least-squares solution (the R
# equivalent of np.linalg.lstsq(..., rcond = NULL)).
.commercial_lstsq_residual <- function(A, b) {
    sv <- svd(as.matrix(A), nu = 0L, nv = 0L)
    b <- as.numeric(b)
    rcond <- max(dim(A)) * .Machine$double.eps * max(sv$d, 1e-300)
    keep <- sv$d > rcond
    # x = V diag(1/s) U' b ; residual via the orthogonal projector, which is
    # numerically stabler than forming x explicitly.
    U <- sv$u
    z <- as.numeric(t(U) %*% b)
    rz <- numeric(length(z))
    rz[!keep] <- z[!keep]
    as.numeric(U %*% rz)
}

#' Cosine-similarity alpha selection (\code{cosine_similarity_selection})
#' @keywords internal
.commercial_cosine_similarity_selection <- function(A, b, initial_spectrum,
                                                    n_alphas = 100L,
                                                    alpha_range = c(-9, 2),
                                                    norm = 2L) {
    alphas <- 10^seq(alpha_range[1], alpha_range[2], length.out = n_alphas)
    norm_init <- sqrt(sum(initial_spectrum^2))
    if (norm_init == 0) stop("Initial spectrum has zero norm.", call. = FALSE)
    init_n <- as.numeric(initial_spectrum) / norm_init
    sv <- compute_svd_components(A)
    U <- sv$U; s <- sv$s; Vt <- sv$Vt; s_sq <- sv$s_sq
    UTb <- as.numeric(t(U) %*% as.numeric(b))
    similarities <- numeric(0)
    for (alpha in alphas) {
        filt <- s / (s_sq + alpha)
        x <- as.numeric(t(Vt) %*% (filt * UTb))
        x <- pmax(x, 0)
        norm_x <- sqrt(sum(x^2))
        if (norm_x == 0) similarities <- c(similarities, 0)
        similarities <- c(similarities, sum(x * init_n) / norm_x)
    }
    idx <- which.max(similarities)             # like np.argmax: first maximum
    if (is.na(idx) || idx > length(alphas)) {
        stop("list index out of range", call. = FALSE)
    }
    as.numeric(alphas[idx])
}

#' L-curve corner alpha, fallback branch (\code{_lcurve_fallback})
#'
#' The Python fallback uses \code{np.cross} on two 2-component vectors, which
#' NumPy >= 2 rejects outright ("Both input arrays must be (arrays of)
#' 3-dimensional vectors"), so \code{regularization_method = "lcurve"} raises
#' \code{ValueError("Regularization selection failed: ...")} in the reference
#' environment.  This port evaluates the same documented chord-distance corner
#' detector with the scalar 2-D cross product \eqn{a_1 b_2 - a_2 b_1}{a1*b2 -
#' a2*b1}, which is what the code intends.
#' @keywords internal
.commercial_lcurve_selection <- function(A, b, n_alphas = 50L,
                                         alpha_range = c(1e-9, 1e2)) {
    alphas <- 10^seq(log10(alpha_range[1]), log10(alpha_range[2]),
                     length.out = n_alphas)
    n <- ncol(A)
    ATA <- crossprod(A); ATb <- as.numeric(crossprod(A, as.numeric(b)))
    residuals <- numeric(0); norms <- numeric(0)
    for (alpha in alphas) {
        x <- tryCatch(as.numeric(.commercial_sym_solve(ATA + alpha * diag(n),
                                                       ATb)),
                      error = function(e) NULL)
        if (is.null(x) || !all(is.finite(x))) next
        x <- pmax(x, 0)
        residuals <- c(residuals, sqrt(sum((as.numeric(A %*% x) -
                                            as.numeric(b))^2)))
        norms <- c(norms, sqrt(sum(x^2)))       # L = identity -> ||L x||
    }
    if (length(residuals) < 3) return(1)
    log_res <- log(residuals); log_norm <- log(norms)
    p1 <- c(log_res[1], log_norm[1])
    p2 <- c(log_res[length(log_res)], log_norm[length(log_norm)])
    dvec <- p2 - p1
    den <- sqrt(sum(dvec^2))
    if (!is.finite(den) || den == 0) return(as.numeric(alphas[1]))
    dist <- mapply(function(r, nn) {
        p <- c(r, nn); v1 <- c(p1[1] - p[1], p1[2] - p[2])
        abs(dvec[1] * v1[2] - dvec[2] * v1[1]) / den
    }, log_res, log_norm)
    as.numeric(alphas[which.max(dist)])
}

#' GCV alpha selection, fallback branch (\code{_gcv_fallback})
#' @keywords internal
.commercial_gcv_selection <- function(A, b, n_alphas = 50L,
                                      alpha_range = c(1e-9, 1e2)) {
    alphas <- 10^seq(log10(alpha_range[1]), log10(alpha_range[2]),
                     length.out = n_alphas)
    m <- nrow(A)
    sv <- compute_svd_components(A)
    s_sq <- sv$s_sq
    UTb <- as.numeric(t(sv$U) %*% as.numeric(b))
    gcv <- vapply(alphas, function(alpha) {
        filt <- s_sq / (s_sq + alpha)
        rc <- alpha / (s_sq + alpha)
        num <- sum((rc * UTb)^2)
        den <- (m - sum(filt))^2
        if (den == 0) Inf else num / den
    }, numeric(1))
    if (!length(gcv) || !any(is.finite(gcv))) return(1)
    as.numeric(alphas[which.min(gcv)])
}

#' Discrepancy-principle alpha selection, fallback branch (\code{_dp_fallback})
#' @keywords internal
.commercial_dp_selection <- function(A, b, noise_var, n_alphas = 50L,
                                     alpha_range = c(1e-9, 1e2)) {
    alphas <- 10^seq(log10(alpha_range[1]), log10(alpha_range[2]),
                     length.out = n_alphas)
    delta <- sqrt(noise_var)
    m <- length(b)
    target <- delta * sqrt(m)
    n <- ncol(A)
    ATA <- crossprod(A); ATb <- as.numeric(crossprod(A, as.numeric(b)))
    residuals <- vapply(alphas, function(alpha) {
        x <- tryCatch(as.numeric(.commercial_sym_solve(ATA + alpha * diag(n),
                                                       ATb)),
                      error = function(e) NULL)
        if (is.null(x) || !all(is.finite(x))) return(Inf)
        x <- pmax(x, 0)
        sqrt(sum((as.numeric(A %*% x) - as.numeric(b))^2))
    }, numeric(1))
    as.numeric(alphas[which.min(abs(residuals - target))])
}

#' Dispatch for the automatic selection methods
#' @keywords internal
.commercial_select_regularization_parameter <- function(A, b, method,
                                                        noise_var = NULL,
                                                        initial_spectrum = NULL) {
    switch(method,
        lcurve = .commercial_lcurve_selection(A, b),
        gcv    = .commercial_gcv_selection(A, b),
        dp     = {
            if (is.null(noise_var)) {
                noise_var <- .commercial_estimate_noise_variance(A, b)
            }
            .commercial_dp_selection(A, b, noise_var)
        },
        cosine = .commercial_cosine_similarity_selection(
            A, b, initial_spectrum),
        stop("Unknown regularization selection method: ", method,
             ". Choose from 'lcurve', 'gcv', 'dp', 'cosine'.", call. = FALSE)
    )
}

#' Resolve the regularization parameter from the requested method
#'
#' Port of \code{bssunfold.core.regularization.resolve_regularization_parameter}
#' restricted to the five methods the commercial wrappers expose.  The
#' \code{pytikhonov} branch of \code{lcurve}/\code{gcv}/\code{dp} is absent
#' from the reference environment, so the documented fallback formulas are
#' always used (see \code{\link{.commercial_lcurve_selection}} for the one
#' numeric divergence that had to be repaired).
#'
#' @keywords internal
.commercial_resolve_regularization_parameter <- function(A, b,
                                                         regularization_method,
                                                         regularization,
                                                         n_energy_bins,
                                                         initial_spectrum = NULL,
                                                         norm = 2L,
                                                         noise_var = NULL) {
    if (identical(regularization_method, "manual")) {
        return(as.numeric(regularization))
    }
    if (identical(regularization_method, "cosine")) {
        if (is.null(initial_spectrum)) {
            stop("For 'cosine' regularization method, initial_spectrum must ",
                 "be provided.", call. = FALSE)
        }
        if (norm != 2) {
            warning("Cosine regularization selection method assumes L2 norm, ",
                    "but norm=", norm, " was requested. Using L2 for ",
                    "selection.", call. = FALSE)
        }
        init <- pmax(as.numeric(initial_spectrum), 0)
        if (length(init) != n_energy_bins) {
            stop("Initial spectrum length (", length(initial_spectrum),
                 ") must match number of energy bins (", n_energy_bins, ")",
                 call. = FALSE)
        }
        return(.commercial_cosine_similarity_selection(A, b, init, norm = 2L))
    }
    if (norm != 2) {
        warning("Automatic regularization selection methods assume L2 norm, ",
                "but norm=", norm, " was requested. Using L2 for selection.",
                call. = FALSE)
    }
    sel <- tryCatch(
        .commercial_select_regularization_parameter(
            A, b, regularization_method, noise_var = noise_var,
            initial_spectrum = initial_spectrum),
        error = function(e) stop("Regularization selection failed: ",
                                 conditionMessage(e),
                                 ". Consider using manual regularization.",
                                 call. = FALSE))
    as.numeric(sel)
}

# ---------------------------------------------------------------------------
# public shared entry points
# ---------------------------------------------------------------------------

#' Solve the unfolding QP with one commercial-solver formulation
#'
#' R port of \code{bssunfold.core.unfold_commercial.solve_commercial} and its
#' backend \code{bssunfold.core._commercial_qp.solve_commercial_qp}: the QP is
#' built exactly as in the Python backend and solved by the package's own
#' active-set engine (see \code{\link{commercial-qp-backend}}).  \code{solver}
#' selects the reported backend only -- the problem is identical for all five
#' aliases, mirroring the \code{functools.partial} bindings of the Python
#' module.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional numeric warm start (length n).  Default \code{NULL}.
#' @param alpha Numeric regularization parameter.  Default \code{1e-4}.
#' @param norm Integer norm type, 1 (L1) or 2 (L2).  Default \code{2}.
#' @param solver Character commercial alias: \code{"gurobi"}, \code{"mosek"},
#'   \code{"cplex"}, \code{"copt"} or \code{"xpress"}.  Default
#'   \code{"gurobi"}.
#' @param timeout Numeric time limit in seconds.  Default \code{10}.  In the
#'   licensed engines this is a true solver time limit; here it is an advisory
#'   wall-clock guard on the active-set loop, so the KKT certificate is
#'   normally reached well inside it.
#' @param smoothness_order Integer smoothness order, 0, 1 or 2.  Default
#'   \code{0}.
#' @param smoothness_weight Numeric weight on the derivative penalty.  Default
#'   \code{1}.
#' @param nonneg Logical; constrain \eqn{x \ge 0}.  Default \code{TRUE}.
#' @param random_state Optional integer seed.  The QP itself is deterministic;
#'   the seed is only forwarded to \code{\link{set.seed}} for parity with the
#'   Python \code{seed=} keyword.
#' @param ub Optional numeric per-bin upper bounds; non-finite entries mean
#'   "unbounded".  Default \code{NULL}.
#' @return A list \code{list(spectrum, iterations, converged, kkt)}.
#'   \code{converged} is the KKT certificate, \code{kkt} its residual.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_commercial(A, b, alpha = 1e-4)
#' round(r$spectrum, 6)
solve_commercial <- function(A, b, x0 = NULL, alpha = 1e-4, norm = 2L,
                             solver = "gurobi", timeout = 10.0,
                             smoothness_order = 0L, smoothness_weight = 1.0,
                             nonneg = TRUE, random_state = NULL, ub = NULL) {
    .commercial_qp_solve(A = A, b = b, solver = solver, x0 = x0, alpha = alpha,
                         norm = norm, timeout = timeout,
                         smoothness_order = smoothness_order,
                         smoothness_weight = smoothness_weight,
                         nonneg = nonneg, random_state = random_state, ub = ub)
}

# The single shared workhorse behind all six entry points.
.commercial_qp_solve <- function(A, b, solver = "gurobi", x0 = NULL,
                                 alpha = 1e-4, norm = 2L, timeout = 10,
                                 smoothness_order = 0L,
                                 smoothness_weight = 1.0, nonneg = TRUE,
                                 random_state = NULL, ub = NULL) {
    info <- .commercial_solver_info(solver)      # validates the alias
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    if (length(dim(A)) != 2L || !is.null(dim(b)) || length(b) != nrow(A)) {
        stop(solver, " solver: received ill-formed input.", call. = FALSE)
    }
    norm <- as.integer(norm)
    if (!norm %in% c(1L, 2L)) {
        stop("Unsupported norm type: ", norm, call. = FALSE)
    }
    if (norm == 1L && !isTRUE(nonneg)) {
        stop(solver, " solver: L1 penalty equals alpha * sum(x) only under ",
             "the non-negativity constraint (pass nonneg=TRUE or use ",
             "norm=2).", call. = FALSE)
    }
    n <- ncol(A)
    if (!is.null(random_state)) {
        suppressWarnings(set.seed(as.integer(random_state)))
    }
    qp <- .commercial_qp_data(A, b, alpha, norm, smoothness_order,
                              smoothness_weight)
    lo <- if (isTRUE(nonneg)) rep(0, n) else rep(-Inf, n)
    hi <- if (is.null(ub)) rep(Inf, n) else as.numeric(ub)
    if (length(hi) != n) {
        stop(solver, " solver: received ill-formed input.", call. = FALSE)
    }
    hi[is.na(hi)] <- Inf
    if (any(is.finite(hi) & hi < lo)) {
        warning("Commercial solver '", info$cvxpy_solver, "' did not find a ",
                "solution (infeasible upper bounds). Returning zero spectrum.",
                call. = FALSE)
        return(list(spectrum = numeric(n), iterations = 0L, converged = FALSE,
                    kkt = Inf))
    }
    deadline <- proc.time()[3] + max(as.numeric(timeout), 1e-3)
    sol <- .commercial_box_qp(qp$P, qp$q, lo, hi, x0 = x0,
                              deadline = deadline)
    x <- as.numeric(sol$x)
    # mirror Python's `result[ub == 0.0] = 0.0` clean-up
    if (!is.null(ub)) x[as.numeric(ub) == 0] <- 0
    if (isTRUE(nonneg)) x <- pmax(x, 0)
    list(spectrum = x, iterations = sol$iterations,
         converged = isTRUE(sol$converged), kkt = sol$kkt)
}

# Shared workflow behind the five unfold_<name> wrappers.  Mirrors Python
# `unfold_commercial`: the alias is validated, alpha is resolved from the
# FULL (untrimmed) system, the max-neutron-energy cutoff is expressed as QP
# upper bounds (`upper_bounds(E_MeV, max_neutron_energy)`) exactly like the
# Python wrapper rather than by dropping columns, and the default initial
# spectrum is all zeros.
.unfold_commercial_common <- function(solver, detector_names, n_energy_bins,
                                      E_MeV, sensitivities, cc_icrp116,
                                      save_result_callback, readings, ln_steps,
                                      reading_uncertainties,
                                      reading_covariance, noise_model,
                                      measurement_time, initial_spectrum,
                                      regularization, norm, timeout,
                                      smoothness_order, smoothness_weight,
                                      nonneg, calculate_errors, noise_level,
                                      n_montecarlo, save_result,
                                      regularization_method, noise_var,
                                      random_state, max_neutron_energy) {
    info <- .commercial_solver_info(solver)      # validates the alias
    sys <- .build_system(readings, detector_names, sensitivities)
    alpha <- .commercial_resolve_regularization_parameter(
        sys$A, sys$b, regularization_method, regularization,
        as.integer(n_energy_bins), initial_spectrum = initial_spectrum,
        norm = norm, noise_var = noise_var)
    ub <- .commercial_upper_bounds(E_MeV, max_neutron_energy)
    x0_default <- rep(0, as.integer(n_energy_bins))

    solver_fun <- function(A, b, x0 = NULL, ...) {
        .commercial_qp_solve(A, b, solver = solver, x0 = x0, alpha = alpha,
                             norm = norm, timeout = timeout,
                             smoothness_order = smoothness_order,
                             smoothness_weight = smoothness_weight,
                             nonneg = nonneg, random_state = random_state,
                             ub = ub)
    }

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver_fun, solve_kwargs = list(),
        method_name = solver,
        extra_output = list(
            norm = norm,
            solver = solver,
            cvxpy_solver = info$cvxpy_solver,
            pip_package = info$pip_package,
            license_required = TRUE,
            regularization = regularization,
            regularization_method = regularization_method,
            selected_regularization = as.numeric(alpha),
            smoothness_order = smoothness_order,
            smoothness_weight = smoothness_weight,
            timeout = timeout,
            nonneg = nonneg
        ),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result)
}

# ---- documented-but-unused Python-only arguments ---------------------------
# Shared roxygen fragment; `@inheritParams` cannot reach a NULL placeholder,
# so the five entry points document them through `commercial-extra-args`.
#' Arguments accepted for Python API parity
#'
#' @param ln_steps Optional numeric \eqn{d(\\ln E)} bin widths.  Accepted for
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
#' @param initial_spectrum Optional numeric initial spectrum guess (length
#'   \code{n_energy_bins}).  Also the reference spectrum required by
#'   \code{regularization_method = "cosine"}.
#' @param regularization Numeric regularization parameter, used when
#'   \code{regularization_method = "manual"}.  Default \code{1e-4}.
#' @param regularization_method Character; one of \code{"manual"},
#'   \code{"cosine"}, \code{"lcurve"}, \code{"gcv"}, \code{"dp"}.  Default
#'   \code{"manual"}.
#' @param noise_var Optional noise variance for the discrepancy principle
#'   (\code{regularization_method = "dp"}); estimated from the least-squares
#'   residual when \code{NULL}.
#' @name commercial-extra-args
NULL

# ---------------------------------------------------------------------------
# The five public entry points
# ---------------------------------------------------------------------------

#' Unfold a neutron spectrum with the Gurobi QP formulation
#'
#' R port of \code{bssunfold.core.unfold_commercial(solver = "gurobi")}.
#' Gurobi itself is proprietary and is not required here: the identical
#' convex quadratic program is solved by the package's internal active-set
#' QP engine (see \code{\link{commercial-qp-backend}}), so this method
#' returns a real spectrum in an environment with no Gurobi license, where
#' the Python original warns and returns zeros.
#'
#' @inheritParams unfold_fista
#' @inheritParams commercial-extra-args
#' @param timeout Advisory wall-clock limit in seconds.  Default \code{10}.
#' @param norm Norm type, \code{1} or \code{2}.  Default \code{2}.
#' @param smoothness_order Derivative order of the smoothness penalty, 0, 1
#'   or 2.  Default \code{0}.
#' @param smoothness_weight Weight of the smoothness penalty.  Default
#'   \code{1}.
#' @param nonneg Constrain the spectrum to be non-negative.  Default
#'   \code{TRUE}.
#' @return A list as produced by \code{run_unfolding}, additionally carrying
#'   \code{license_required = TRUE} and the backend metadata.
#' @seealso \code{\link{commercial-qp-backend}} for the problem formulation.
#' @family commercial solvers
#' @export
unfold_gurobi <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, initial_spectrum = NULL,
                          regularization = 1e-4, norm = 2L, timeout = 10.0,
                          smoothness_order = 0L, smoothness_weight = 1.0,
                          nonneg = TRUE, calculate_errors = FALSE,
                          noise_level = 0.01, n_montecarlo = 100L,
                          save_result = FALSE,
                          regularization_method = "manual", noise_var = NULL,
                          random_state = NULL, max_neutron_energy = NULL,
                          ln_steps = NULL, reading_uncertainties = NULL,
                          reading_covariance = NULL,
                          noise_model = "gaussian",
                          measurement_time = NULL) {
    .unfold_commercial_common(
        solver = "gurobi", detector_names = detector_names,
        n_energy_bins = n_energy_bins, E_MeV = E_MeV,
        sensitivities = sensitivities, cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback, readings = readings,
        ln_steps = ln_steps,
        reading_uncertainties = reading_uncertainties,
        reading_covariance = reading_covariance, noise_model = noise_model,
        measurement_time = measurement_time,
        initial_spectrum = initial_spectrum,
        regularization = regularization, norm = norm, timeout = timeout,
        smoothness_order = smoothness_order,
        smoothness_weight = smoothness_weight, nonneg = nonneg,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, save_result = save_result,
        regularization_method = regularization_method, noise_var = noise_var,
        random_state = random_state,
        max_neutron_energy = max_neutron_energy)
}

#' Unfold a neutron spectrum with the MOSEK QP formulation
#'
#' R port of \code{bssunfold.core.unfold_commercial(solver = "mosek")}; the
#' problem is identical to \code{\link{unfold_gurobi}} and is solved by the
#' package's internal QP engine.
#'
#' @inheritParams unfold_gurobi
#' @family commercial solvers
#' @export
unfold_mosek <- function(detector_names, n_energy_bins, E_MeV,
                         sensitivities, cc_icrp116, save_result_callback,
                         readings, initial_spectrum = NULL,
                         regularization = 1e-4, norm = 2L, timeout = 10.0,
                         smoothness_order = 0L, smoothness_weight = 1.0,
                         nonneg = TRUE, calculate_errors = FALSE,
                         noise_level = 0.01, n_montecarlo = 100L,
                         save_result = FALSE,
                         regularization_method = "manual", noise_var = NULL,
                         random_state = NULL, max_neutron_energy = NULL,
                         ln_steps = NULL, reading_uncertainties = NULL,
                         reading_covariance = NULL,
                         noise_model = "gaussian",
                         measurement_time = NULL) {
    .unfold_commercial_common(
        solver = "mosek", detector_names = detector_names,
        n_energy_bins = n_energy_bins, E_MeV = E_MeV,
        sensitivities = sensitivities, cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback, readings = readings,
        ln_steps = ln_steps,
        reading_uncertainties = reading_uncertainties,
        reading_covariance = reading_covariance, noise_model = noise_model,
        measurement_time = measurement_time,
        initial_spectrum = initial_spectrum,
        regularization = regularization, norm = norm, timeout = timeout,
        smoothness_order = smoothness_order,
        smoothness_weight = smoothness_weight, nonneg = nonneg,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, save_result = save_result,
        regularization_method = regularization_method, noise_var = noise_var,
        random_state = random_state,
        max_neutron_energy = max_neutron_energy)
}

#' Unfold a neutron spectrum with the CPLEX QP formulation
#'
#' R port of \code{bssunfold.core.unfold_commercial(solver = "cplex")}; the
#' problem is identical to \code{\link{unfold_gurobi}} and is solved by the
#' package's internal QP engine.
#'
#' @inheritParams unfold_gurobi
#' @family commercial solvers
#' @export
unfold_cplex <- function(detector_names, n_energy_bins, E_MeV,
                         sensitivities, cc_icrp116, save_result_callback,
                         readings, initial_spectrum = NULL,
                         regularization = 1e-4, norm = 2L, timeout = 10.0,
                         smoothness_order = 0L, smoothness_weight = 1.0,
                         nonneg = TRUE, calculate_errors = FALSE,
                         noise_level = 0.01, n_montecarlo = 100L,
                         save_result = FALSE,
                         regularization_method = "manual", noise_var = NULL,
                         random_state = NULL, max_neutron_energy = NULL,
                         ln_steps = NULL, reading_uncertainties = NULL,
                         reading_covariance = NULL,
                         noise_model = "gaussian",
                         measurement_time = NULL) {
    .unfold_commercial_common(
        solver = "cplex", detector_names = detector_names,
        n_energy_bins = n_energy_bins, E_MeV = E_MeV,
        sensitivities = sensitivities, cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback, readings = readings,
        ln_steps = ln_steps,
        reading_uncertainties = reading_uncertainties,
        reading_covariance = reading_covariance, noise_model = noise_model,
        measurement_time = measurement_time,
        initial_spectrum = initial_spectrum,
        regularization = regularization, norm = norm, timeout = timeout,
        smoothness_order = smoothness_order,
        smoothness_weight = smoothness_weight, nonneg = nonneg,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, save_result = save_result,
        regularization_method = regularization_method, noise_var = noise_var,
        random_state = random_state,
        max_neutron_energy = max_neutron_energy)
}

#' Unfold a neutron spectrum with the COPT QP formulation
#'
#' R port of \code{bssunfold.core.unfold_commercial(solver = "copt")}; the
#' problem is identical to \code{\link{unfold_gurobi}} and is solved by the
#' package's internal QP engine.
#'
#' @inheritParams unfold_gurobi
#' @family commercial solvers
#' @export
unfold_copt <- function(detector_names, n_energy_bins, E_MeV,
                        sensitivities, cc_icrp116, save_result_callback,
                        readings, initial_spectrum = NULL,
                        regularization = 1e-4, norm = 2L, timeout = 10.0,
                        smoothness_order = 0L, smoothness_weight = 1.0,
                        nonneg = TRUE, calculate_errors = FALSE,
                        noise_level = 0.01, n_montecarlo = 100L,
                        save_result = FALSE,
                        regularization_method = "manual", noise_var = NULL,
                        random_state = NULL, max_neutron_energy = NULL,
                        ln_steps = NULL, reading_uncertainties = NULL,
                        reading_covariance = NULL,
                        noise_model = "gaussian",
                        measurement_time = NULL) {
    .unfold_commercial_common(
        solver = "copt", detector_names = detector_names,
        n_energy_bins = n_energy_bins, E_MeV = E_MeV,
        sensitivities = sensitivities, cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback, readings = readings,
        ln_steps = ln_steps,
        reading_uncertainties = reading_uncertainties,
        reading_covariance = reading_covariance, noise_model = noise_model,
        measurement_time = measurement_time,
        initial_spectrum = initial_spectrum,
        regularization = regularization, norm = norm, timeout = timeout,
        smoothness_order = smoothness_order,
        smoothness_weight = smoothness_weight, nonneg = nonneg,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, save_result = save_result,
        regularization_method = regularization_method, noise_var = noise_var,
        random_state = random_state,
        max_neutron_energy = max_neutron_energy)
}

#' Unfold a neutron spectrum with the XPRESS QP formulation
#'
#' R port of \code{bssunfold.core.unfold_commercial(solver = "xpress")}; the
#' problem is identical to \code{\link{unfold_gurobi}} and is solved by the
#' package's internal QP engine.
#'
#' @inheritParams unfold_gurobi
#' @family commercial solvers
#' @export
unfold_xpress <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, initial_spectrum = NULL,
                          regularization = 1e-4, norm = 2L, timeout = 10.0,
                          smoothness_order = 0L, smoothness_weight = 1.0,
                          nonneg = TRUE, calculate_errors = FALSE,
                          noise_level = 0.01, n_montecarlo = 100L,
                          save_result = FALSE,
                          regularization_method = "manual", noise_var = NULL,
                          random_state = NULL, max_neutron_energy = NULL,
                          ln_steps = NULL, reading_uncertainties = NULL,
                          reading_covariance = NULL,
                          noise_model = "gaussian",
                          measurement_time = NULL) {
    .unfold_commercial_common(
        solver = "xpress", detector_names = detector_names,
        n_energy_bins = n_energy_bins, E_MeV = E_MeV,
        sensitivities = sensitivities, cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback, readings = readings,
        ln_steps = ln_steps,
        reading_uncertainties = reading_uncertainties,
        reading_covariance = reading_covariance, noise_model = noise_model,
        measurement_time = measurement_time,
        initial_spectrum = initial_spectrum,
        regularization = regularization, norm = norm, timeout = timeout,
        smoothness_order = smoothness_order,
        smoothness_weight = smoothness_weight, nonneg = nonneg,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, save_result = save_result,
        regularization_method = regularization_method, noise_var = noise_var,
        random_state = random_state,
        max_neutron_energy = max_neutron_energy)
}
