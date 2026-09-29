#' qpmad-based unfolding method for neutron spectrum reconstruction
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_qpmad.py}.
#' Method built on the qpmad solver of Alexander Sherikov
#' (https://github.com/asherikov/qpmad), a header-only C++14 implementation
#' of the Goldfarb-Idnani dual active-set algorithm for (strictly convex)
#' quadratic programming:
#' \deqn{\minimize \frac{1}{2} x' H x + g' x}{min 1/2 x'Hx + g'x}
#' subject to simple bounds, general inequality and equality constraints.
#'
#' For BSS unfolding the strictly convex regularised least-squares problem
#' \deqn{\minimize \frac{1}{2} ||A x - b||^2 + \frac{\alpha}{2} ||L x||^2
#'   + \frac{\alpha_0}{2} ||x||^2}{min 1/2 ||Ax-b||^2 + alpha/2 ||Lx||^2 +
#'   alpha0/2 ||x||^2}
#' subject to \eqn{x \ge 0} (default) or \eqn{lb \le x \le ub}{lb <= x <= ub}
#' is solved, which has Hessian \eqn{H = A'A + \alpha L'L + \alpha_0 I}{H =
#' A'A + alpha L'L + alpha0 I} (symmetric positive definite, as required by
#' Goldfarb-Idnani) and linear term \eqn{g = -A'b}{g = -A'b}.
#'
#' The original qpmad is a C++ library (built on Eigen).  Like the Python
#' module, this port provides two interchangeable backends:
#' \itemize{
#'   \item \code{backend = "python"} --- a self-contained base-R port of the
#'     active-set QP solver (Nocedal & Wright, \emph{Numerical Optimization},
#'     ch. 16.4) that the Python module implements.  The algorithmic spirit
#'     is the same as Goldfarb-Idnani (it solves the same strictly convex QP
#'     with inequality constraints); the search path through constraint
#'     activations/deactivations is the one of the Python implementation and
#'     is reproduced here step by step, including its termination behaviour.
#'     This is the default and needs no external dependency.
#'   \item \code{backend = "qpmad"} (\code{"cpp"}) --- would call the
#'     upstream C++ library through its bindings.  No such engine exists in
#'     this R package, so, exactly like the Python module when the binding is
#'     missing, a warning is emitted and the \code{"python"} backend is used.
#' }
#'
#' @name qpmad
NULL

# Least-squares solve with the minimum-norm solution (numpy
# linalg.lstsq(rcond=None)): SVD based, cutoff = max(M, N) * eps.
.qp_lstsq <- function(M, rhs) {
    s <- svd(M)
    cutoff <- .Machine$double.eps * max(nrow(M), ncol(M)) * max(s$d)
    keep <- s$d > cutoff
    n <- sum(keep)
    if (n == 0L) return(rep(0, ncol(M)))
    u <- s$u[, seq_len(n), drop = FALSE]
    v <- s$v[, seq_len(n), drop = FALSE]
    ds <- 1 / s$d[seq_len(n)]
    as.numeric(v %*% (ds * as.numeric(t(u) %*% rhs)))
}

# Complete QR factorisation of t(M): returns list(Q = n x n orthonormal,
# rank = number of rows of M assumed independent) mirroring
# np.linalg.qr(Cw.T, mode = "complete").
.qp_qr_complete <- function(M) {
    n <- ncol(M)
    mt <- matrix(as.numeric(t(M)), nrow = n)
    qr <- qr(mt, pivot = FALSE)
    list(Q = as.matrix(qr.Q(qr, complete = TRUE)), rank = qr$rank)
}

# Solve a strictly convex QP with an active-set method.  Minimises
# 0.5 x'Hx + g'x subject to lb <= x <= ub and lb_A <= A x <= ub_A.
# H must be symmetric positive-definite.  Port of the Python
# _solve_qp_goldfarb_idnani (primal active-set framework, Nocedal & Wright
# ch. 16.4): algorithmically simpler than the dual Goldfarb-Idnani scheme but
# solving the same QP; the search path through the working set is the one of
# the Python implementation.
#
# Returns list(x, status) with status "OK", "INFEASIBLE" or "MAX_ITER".
.qp_solve_goldfarb_idnani <- function(H, g, lb = NULL, ub = NULL, A = NULL,
                                      lb_A = NULL, ub_A = NULL,
                                      tol = 1e-9, max_iterations = 10000L) {
    H <- as.matrix(H); storage.mode(H) <- "double"
    g <- as.numeric(g)
    n <- nrow(H)
    if (ncol(H) != n) {
        stop("H must be square, got ", nrow(H), "x", ncol(H))
    }
    H <- 0.5 * (H + t(H))  # symmetrise defensively

    # ---- Cholesky factorisation of H (with jitter retries) --------------- #
    jitter <- 0
    L_chol <- NULL
    for (trial in 1:5) {
        L_chol <- tryCatch({
            if (jitter > 0) {
                ch <- chol(H + jitter * diag(n))
            } else {
                ch <- chol(H)
            }
            t(ch)   # lower-triangular factor, H = L L'
        }, error = function(e) NULL)
        if (!is.null(L_chol)) break
        jitter <- max(jitter * 10, 1e-10)
        L_chol <- NULL
    }
    if (is.null(L_chol)) {
        # Final fallback: pseudo-inverse based unconstrained minimum
        x <- tryCatch(.qp_lstsq(H, -g), error = function(e) rep(0, n))
        # Apply projection to bounds as best-effort
        if (!is.null(lb) && !is.null(ub)) {
            x <- pmin(pmax(x, as.numeric(lb)), as.numeric(ub))
        }
        return(list(x = x, status = "INFEASIBLE"))
    }

    # Helper: solve H x = rhs using the cached Cholesky factor.
    solve_H <- function(rhs) {
        y <- forwardsolve(L_chol, as.numeric(rhs))
        as.numeric(backsolve(t(L_chol), y))
    }

    # ---- Collect inequality constraints into stacked form C x >= d ------- #
    rows_C <- list()
    vals_d <- numeric(0)
    if (!is.null(lb) && !is.null(ub)) {
        lb_arr <- as.numeric(lb); ub_arr <- as.numeric(ub)
        for (i in seq_len(n)) {
            if (is.finite(lb_arr[i])) {
                r <- numeric(n); r[i] <- 1.0
                rows_C[[length(rows_C) + 1L]] <- r
                vals_d <- c(vals_d, lb_arr[i])
            }
            if (is.finite(ub_arr[i])) {
                r <- numeric(n); r[i] <- -1.0
                rows_C[[length(rows_C) + 1L]] <- r
                vals_d <- c(vals_d, -ub_arr[i])
            }
        }
    }
    if (!is.null(A) && !is.null(lb_A) && !is.null(ub_A)) {
        A_arr <- as.matrix(A); storage.mode(A_arr) <- "double"
        lbA_arr <- as.numeric(lb_A); ubA_arr <- as.numeric(ub_A)
        for (i in seq_len(nrow(A_arr))) {
            if (is.finite(lbA_arr[i])) {
                rows_C[[length(rows_C) + 1L]] <- A_arr[i, , drop = TRUE]
                vals_d <- c(vals_d, lbA_arr[i])
            }
            if (is.finite(ubA_arr[i])) {
                rows_C[[length(rows_C) + 1L]] <- -A_arr[i, , drop = TRUE]
                vals_d <- c(vals_d, -ubA_arr[i])
            }
        }
    }
    if (length(rows_C) == 0L) {
        # No inequality constraints -- return the unconstrained minimum.
        return(list(x = solve_H(-g), status = "OK"))
    }
    C <- do.call(rbind, rows_C)
    storage.mode(C) <- "double"
    d <- as.numeric(vals_d)
    m <- nrow(C)

    # Each row's single non-zero entry, when the row is a scaled unit vector
    # (the simple-bound case).  Rows that are not unit vectors get index 0.
    row_idx <- integer(m); row_sgn <- numeric(m)
    for (i in seq_len(m)) {
        nz <- which(C[i, ] != 0)
        if (length(nz) == 1L) {
            row_idx[i] <- nz
            row_sgn[i] <- C[i, nz]
        }
    }
    bounds_only <- all(row_idx > 0L) && all(abs(row_sgn) == 1)

    # ---- Primal active-set algorithm ------------------------------------- #
    # Start from a feasible point: the projection of the unconstrained
    # minimum onto the constraint box.
    x <- solve_H(-g)
    if (!is.null(lb) && !is.null(ub)) {
        x <- pmin(pmax(x, as.numeric(lb)), as.numeric(ub))
    }
    # If general constraints are present, make sure they are feasible too:
    # a simple alternating projection scheme (50 sweeps).
    for (proj in 1:50) {
        violations <- as.numeric(C %*% x) - d
        if (all(violations >= -tol)) break
        j <- which.min(violations)
        n_j <- solve_H(C[j, , drop = TRUE])
        denom <- as.numeric(C[j, , drop = TRUE] %*% n_j)
        if (abs(denom) < 1e-15) break
        x <- x + max(0, -violations[j]) / denom * n_j
        if (!is.null(lb) && !is.null(ub)) {
            x <- pmin(pmax(x, as.numeric(lb)), as.numeric(ub))
        }
    }
    if (any(as.numeric(C %*% x) - d < -tol * 100)) {
        return(list(x = x, status = "INFEASIBLE"))
    }

    # Active set: W = { i : C_i x = d_i }  (logical over the m rows)
    W <- logical(m)
    for (i in seq_len(m)) {
        if (abs(as.numeric(C[i, , drop = TRUE] %*% x) - d[i]) <= tol * 10) {
            W[i] <- TRUE
        }
    }

    # Helper: solve the equality-constrained QP on the working set,
    #   minimise 0.5 x'Hx + g'x  s.t.  C_W x = d_W
    # by the null-space method x = x_part + Z p.
    solve_eqp <- function(x_cur, W) {
        W_list <- which(W)
        if (length(W_list) == 0L) return(solve_H(-g))
        Cw <- C[W_list, , drop = FALSE]
        dw <- d[W_list]
        if (bounds_only && length(unique(row_idx[W_list])) == length(W_list)) {
            # Closed form for unit constraint rows: x_part is exact on the
            # active coordinates and zero elsewhere, Z is the identity on the
            # free coordinates.
            x_part <- numeric(n)
            x_part[row_idx[W_list]] <- dw / row_sgn[W_list]
            free <- setdiff(seq_len(n), row_idx[W_list])
            if (length(free) == 0L) return(x_part)
            Hred <- H[free, free, drop = FALSE]
            grad_red <- as.numeric(H %*% x_part + g)[free]
            p <- tryCatch(as.numeric(solve(Hred, -grad_red)),
                          error = function(e) .qp_lstsq(Hred, -grad_red))
            x_new <- x_part
            x_new[free] <- p
            return(x_new)
        }
        # Compute a particular solution by solving C_W x = d_W via least
        # squares (gives the closest feasible point to the current iterate).
        x_part <- tryCatch(.qp_lstsq(Cw, dw), error = function(e) x_cur)
        # Null-space basis Z via QR of C_W'
        q <- tryCatch(.qp_qr_complete(Cw), error = function(e) NULL)
        if (is.null(q)) return(x_cur)
        k <- nrow(Cw)
        Z <- q$Q[, (k + 1L):n, drop = FALSE]
        # Reduced problem: min 0.5 p'(Z'HZ) p + (H x_part + g)' Z p
        Hz <- H %*% Z
        Hred <- crossprod(Z, Hz)
        grad_red <- as.numeric(t(Z) %*% (as.numeric(H %*% x_part) + g))
        p <- tryCatch(as.numeric(solve(Hred, -grad_red)),
                      error = function(e) .qp_lstsq(Hred, -grad_red))
        as.numeric(x_part + Z %*% p)
    }

    # Multipliers of the working-set constraints from
    #   H x + g = -sum_i mu_i C_i'   ->   C_W' mu_W = -(H x + g)
    # and the KKT sign test (constraints are c_i'x >= d_i, so mu >= 0).
    multiplier_drop <- function(x, W_list) {
        Cw <- C[W_list, , drop = FALSE]
        rhs <- -(as.numeric(H %*% x) + g)
        mu <- if (bounds_only &&
                  length(unique(row_idx[W_list])) == length(W_list)) {
            row_sgn[W_list] * rhs[row_idx[W_list]]
        } else {
            tryCatch(.qp_lstsq(t(Cw), rhs), error = function(e) rep(0, length(W_list)))
        }
        list(mu = mu, drop = which.min(mu))
    }

    # Main active-set loop.
    for (iter in seq_len(as.integer(max_iterations))) {
        x_new <- solve_eqp(x, W)

        # Maximum feasible step length alpha in [0, 1] such that
        # x + alpha (x_new - x) stays feasible for the constraints that are
        # not in the working set (min-of-ratios blocking test).
        direction <- x_new - x
        if (sqrt(sum(direction^2)) <= tol) {
            # KKT point on the current working set.  Check multipliers.
            W_list <- which(W)
            if (length(W_list) > 0L) {
                md <- multiplier_drop(x, W_list)
                if (all(md$mu >= -tol)) {
                    return(list(x = x, status = "OK"))
                }
                # Drop the most-negative multiplier.
                W[W_list[md$drop]] <- FALSE
                next
            }
            return(list(x = x, status = "OK"))
        }

        alpha <- 1.0
        blocking_idx <- -1L
        for (i in seq_len(m)) {
            if (W[i]) next
            c_i <- C[i, , drop = TRUE]
            directional <- as.numeric(c_i %*% direction)
            slack <- as.numeric(c_i %*% x) - d[i]  # >= 0 in the feasible set
            if (directional < -tol && slack < -directional * alpha) {
                # Constraint i becomes active at alpha = -slack / directional
                a_i <- -slack / directional
                if (a_i < alpha) {
                    alpha <- a_i
                    blocking_idx <- i
                }
            }
        }
        alpha <- max(0, min(1, alpha))
        x <- x + alpha * direction

        if (blocking_idx >= 0L && alpha < 1) {
            W[blocking_idx] <- TRUE
        }

        if (alpha >= 1) {
            # Full step taken; check KKT conditions on the new working set.
            W_list <- which(W)
            if (length(W_list) > 0L) {
                md <- multiplier_drop(x, W_list)
                if (all(md$mu >= -tol)) {
                    return(list(x = x, status = "OK"))
                }
                W[W_list[md$drop]] <- FALSE
            } else {
                return(list(x = x, status = "OK"))
            }
        }
    }

    # Iteration cap reached: if the solution is feasible and the KKT
    # residual is tiny, declare success anyway.  This guards against
    # pathological cycling between working sets on numerically-degenerate
    # problems where the active-set oscillation does not affect the
    # practical quality of the solution.
    W_list <- which(W)
    grad_lag <- as.numeric(H %*% x) + g
    if (length(W_list) > 0L) {
        Cw <- C[W_list, , drop = FALSE]
        q <- tryCatch(.qp_qr_complete(Cw), error = function(e) NULL)
        grad_proj <- if (is.null(q)) {
            grad_lag
        } else {
            Zk <- q$Q[, (length(W_list) + 1L):n, drop = FALSE]
            as.numeric(t(Zk) %*% grad_lag)
        }
    } else {
        grad_proj <- grad_lag
    }

    kkt_residual <- sqrt(sum(grad_proj^2))
    feasibility <- min(as.numeric(C %*% x) - d)  # should be >= -tol
    if (kkt_residual <= 1e-6 && feasibility >= -1e-6) {
        return(list(x = x, status = "OK"))
    }
    list(x = x, status = "MAX_ITER")
}

# Try the (hypothetical) compiled qpmad engine.  Mirrors
# _try_import_qpmad(): nothing is available here, so NULL is returned and the
# caller falls back to the R port with a warning.
.qp_try_qpmad_engine <- function() NULL

#' Solve the BSS unfolding problem with qpmad (active-set QP)
#'
#' Core solver mirroring \code{solve_qpmad} in
#' \code{bssunfold/src/bssunfold/core/unfold_qpmad.py}.  Recasts the
#' regularised non-negative least-squares problem
#' \deqn{\minimize \frac{1}{2} ||A x - b||^2 + \frac{\alpha}{2} ||L x||^2
#'   + \frac{\alpha_0}{2} ||x||^2}{min 1/2 ||Ax-b||^2 + alpha/2 ||Lx||^2 +
#'   alpha0/2 ||x||^2}
#' subject to \eqn{lb \le x \le ub}{lb <= x <= ub} (default \eqn{x \ge 0})
#' as the QP \eqn{\min 0.5 x'Hx + g'x} with \eqn{H = A'A + \alpha L'L +
#' \alpha_0 I}{H = A'A + alpha L'L + alpha0 I} and \eqn{g = -A'b}{g = -A'b},
#' and solves it with the qpmad algorithm of Sherikov.
#'
#' @rdname qpmad
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional spectrum accepted for API compatibility --- the
#'   Goldfarb-Idnani algorithm does not use a warm start (it always starts
#'   from the unconstrained minimum).  Default \code{NULL}.
#' @param regularization Numeric; Tikhonov / smoothness regularisation
#'   weight.  Default 1e-4.
#' @param smoothness_order Integer; smoothness penalty order (0, 1 or 2).
#'   Default 0.
#' @param smoothness_weight Numeric; weight for the smoothness term.
#'   Default 1.0.
#' @param floor Numeric; diagonal regularisation floor added to \code{H} to
#'   guarantee strict positive-definiteness.  Default 1e-6.
#' @param lb,ub Optional numeric simple bounds (length n).  If both are
#'   \code{NULL} (default) the method enforces \eqn{x \ge 0}{x >= 0}.
#' @param backend Character; \code{"python"} (default) uses the base-R port
#'   of the active-set solver; \code{"qpmad"} / \code{"cpp"} request the
#'   upstream C++ library and fall back to \code{"python"} with a warning
#'   when it is unavailable.
#' @param tol Numeric; numerical tolerance.  Default 1e-9.
#' @param max_iterations Positive integer; iteration cap.  Default 10000.
#' @return A list \code{list(spectrum, iterations, converged)} where
#'   \code{iterations} is the status code: 0 (OK), 1 (infeasible, treated as
#'   not converged) or 2 (max iterations hit, treated as not converged).
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_qpmad(A, b, rep(0, 3))
solve_qpmad <- function(A, b, x0 = NULL, regularization = 1e-4,
                        smoothness_order = 0L, smoothness_weight = 1.0,
                        floor = 1e-6, lb = NULL, ub = NULL,
                        backend = "python", tol = 1e-9,
                        max_iterations = 10000L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)

    if (!smoothness_order %in% c(0L, 1L, 2L)) {
        stop("Unsupported smoothness order: ", smoothness_order,
             ". Use 0, 1 or 2.")
    }
    if (!backend %in% c("python", "qpmad", "cpp")) {
        stop("Unsupported backend: ", deparse(backend),
             ". Use 'python' or 'qpmad'.")
    }

    # Build H = A'A + alpha * L'L + alpha0 * I (symmetric PD).
    H <- t(A) %*% A + floor * diag(n)
    smoothness_order <- as.integer(smoothness_order)
    if (smoothness_order %in% c(1L, 2L) && regularization > 0) {
        L <- as.matrix(create_derivative_matrix(n, smoothness_order))
        H <- H + as.numeric(regularization) * as.numeric(smoothness_weight) *
            (t(L) %*% L)
    }
    H <- 0.5 * (H + t(H))

    g <- -as.numeric(t(A) %*% b)

    # Default bounds: x >= 0.
    if (is.null(lb) && is.null(ub)) {
        lb_arr <- rep(0, n)
        ub_arr <- rep(Inf, n)
    } else {
        lb_arr <- if (!is.null(lb)) as.numeric(lb) else rep(-Inf, n)
        ub_arr <- if (!is.null(ub)) as.numeric(ub) else rep(Inf, n)
    }

    # Try the C++ backend first if requested.
    if (backend %in% c("qpmad", "cpp")) {
        engine <- .qp_try_qpmad_engine()
        if (is.null(engine)) {
            warning("qpmad C++ bindings not available; using the python backend. ",
                    "Install qpmad with bindings to use the C++ backend.",
                    call. = FALSE)
        }
    }

    # Pure-R implementation.
    res <- .qp_solve_goldfarb_idnani(H, g, lb = lb_arr, ub = ub_arr,
                                     tol = tol, max_iterations = max_iterations)
    converged <- identical(res$status, "OK")
    code <- if (identical(res$status, "OK")) 0L else
        if (identical(res$status, "INFEASIBLE")) 1L else 2L
    list(spectrum = pmax(as.numeric(res$x), 0), iterations = code,
         converged = converged)
}

#' qpmad unfolding (unified workflow wrapper)
#'
#' Thin wrapper around \code{\link{solve_qpmad}} for the unified workflow,
#' mirroring \code{unfold_qpmad} in
#' \code{bssunfold/src/bssunfold/core/unfold_qpmad.py}.
#'
#' @rdname qpmad
#' @inheritParams run_unfolding
#' @param calculate_errors Logical; if \code{TRUE}, run Monte-Carlo
#'   uncertainty estimation.  Default \code{FALSE}.
#' @param noise_level Numeric; relative Gaussian noise level for Monte-Carlo.
#'   Default 0.01.
#' @param n_montecarlo Integer; number of Monte-Carlo samples.
#'   Default 100.
#' @param random_state Optional integer seed for Monte-Carlo.
#' @param max_neutron_energy Optional numeric energy cutoff in MeV.
#'   Default \code{NULL} = no cutoff.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_qpmad <- function(detector_names, n_energy_bins, E_MeV,
                         sensitivities, cc_icrp116, save_result_callback,
                         readings, initial_spectrum = NULL,
                         regularization = 1e-4, smoothness_order = 0L,
                         smoothness_weight = 1.0, floor = 1e-6,
                         lb = NULL, ub = NULL, backend = "python",
                         tol = 1e-9, max_iterations = 10000L,
                         calculate_errors = FALSE,
                         noise_level = 0.01,
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
        solve_func = make_solve_wrapper(solve_qpmad,
                                        regularization = regularization,
                                        smoothness_order = smoothness_order,
                                        smoothness_weight = smoothness_weight,
                                        floor = floor,
                                        lb = lb, ub = ub,
                                        backend = backend,
                                        tol = tol,
                                        max_iterations = max_iterations),
        solve_kwargs = list(),
        method_name = "qpmad",
        extra_output = list(regularization = regularization,
                            smoothness_order = as.integer(smoothness_order),
                            smoothness_weight = smoothness_weight,
                            floor = floor,
                            backend = backend,
                            tol = tol,
                            max_iterations = as.integer(max_iterations)),
        calculate_errors = calculate_errors,
        noise_level = noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
