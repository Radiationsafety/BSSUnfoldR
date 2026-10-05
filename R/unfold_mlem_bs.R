#' MLEM-BS (B-spline MLEM with sieve regularization)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_mlem_bs.py}.
#' Implements the regularized MLEM iteration on B-spline coefficients of
#' Mazankova et al., CNDGS'2026, \url{https://doi.org/10.47459/cndcgs.2026.61}.
#'
#' The spectrum is represented as \eqn{x(E) = \sum_s b_s B_s(E)} so the
#' effective system matrix is \eqn{RB = R B}. The B-spline coefficients
#' are updated by the multiplicative MLEM iteration with a second-derivative
#' penalty \eqn{P(b) = \|D^{(2)} b\|^2_2} and a sieve restriction to
#' non-negative coefficients:
#'
#' \deqn{b_s^{(k+1)} = \frac{b_s^{(k)}}{\sum_i (RB)_{is} + \beta\, dP/db_s}
#'        \sum_i (RB)_{is} \frac{n_i}{\sum_{s'} (RB)_{is'} b_{s'}^{(k)}}.}
#'
#' @section Scaling note:
#' The absolute penalty \code{beta} of Eq. 4 is problem-scale dependent (the
#' paper uses \code{beta = 1e-17}). \code{beta_relative} defines the effective
#' penalty as \code{beta = beta_relative * mean(colSums(RB))}; exactly one of
#' the two may be non-zero, and with neither the penalty is zero (pure sieve
#' MLEM on the B-spline parameterization).
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n), projected onto the
#'   non-negative B-spline sieve to obtain the starting coefficients.
#'   \code{NULL} uses a flat 0.5 spectrum.
#' @param E_MeV Numeric energy grid (length n).
#' @param n_basis Integer; number of B-spline basis functions \eqn{N_s}.
#'   Default \code{NULL} = \code{max(spline_order + 1, min(n \%/\% 2, 40))}.
#' @param beta Numeric; absolute penalty strength. Exactly one of
#'   \code{beta} / \code{beta_relative} may be non-zero. Default 0.
#' @param beta_relative Numeric; penalty strength relative to the mean
#'   column sum of \eqn{RB}. Default 0.
#' @param max_iterations Positive integer; default 1000.
#' @param tolerance Positive numeric; default 1e-6.
#' @param knot_spacing Character; \code{"auto"}, \code{"uniform"}, or
#'   \code{"log"}. Default \code{"auto"} = log when the grid spans more
#'   than two decades.
#' @param spline_order Integer B-spline order \eqn{p} (degree \eqn{p-1});
#'   the paper uses \code{4L} (cubic). Allowed range 2..8.
#' @param auto_params Logical; select \eqn{N_s}, the penalty strength and
#'   the iteration count by minimising the \eqn{K_S} statistic (Eq. 6).
#' @param ks_patience Patience (iterations without \eqn{K_S} improvement)
#'   used in the \eqn{K_S}-minimising mode.
#' @return A list \code{list(spectrum, iterations, converged, coefficients,
#'   ks_history, ks_final, chi2_pearson, n_basis, spline_order, knot_spacing,
#'   interior_knots, beta_effective, beta_relative)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_mlem_bs(A, b, rep(1, 60), E, n_basis = 10L,
#'                    beta_relative = 1e-3, max_iterations = 50L)
solve_mlem_bs <- function(A, b, x0 = NULL, E_MeV = NULL, n_basis = NULL,
                            beta = 0.0, beta_relative = 0.0,
                            max_iterations = 1000L, tolerance = 1e-6,
                            knot_spacing = "auto", spline_order = 4L,
                            auto_params = FALSE, ks_patience = 50L) {
    out <- .mlem_bs_solve(A, b, x0 = x0, E_MeV = E_MeV, n_basis = n_basis,
                          beta = beta, beta_relative = beta_relative,
                          max_iterations = max_iterations,
                          tolerance = tolerance, knot_spacing = knot_spacing,
                          spline_order = spline_order,
                          auto_params = auto_params,
                          ks_patience = ks_patience)
    out$spectrum <- pmax(as.numeric(out$spectrum), 0.0)
    out
}

# Numerical floors / guards of the Python module ------------------------------
.MB_TINY <- 1e-300           # absolute floor for denominators
.MB_COUNT_FLOOR <- 1e-10     # floor for forward-model counts (div-by-zero)
.MB_AUTO_BETA_GRID <- c(0.0, 1e-4, 1e-3, 1e-2, 1e-1)
.MB_VALID_KNOT_SPACING <- c("auto", "uniform", "log")

# core/unfold_mlem_bs.py:_resolve_knot_spacing
.mlem_bs_resolve_spacing <- function(knot_spacing, E) {
    knot_spacing <- tolower(as.character(knot_spacing)[1L])
    if (!knot_spacing %in% .MB_VALID_KNOT_SPACING) {
        stop("knot_spacing must be one of ",
             paste(.MB_VALID_KNOT_SPACING, collapse = ", "),
             ", got '", knot_spacing, "'", call. = FALSE)
    }
    if (knot_spacing != "auto") return(knot_spacing)
    ratio <- if (length(E) && min(E) > 0) max(E) / min(E) else Inf
    if (isTRUE(ratio > 100)) "log" else "uniform"
}

# core/unfold_mlem_bs.py:build_bspline_basis
.mlem_bs_basis <- function(E_MeV, n_basis, spline_order = 4L,
                           knot_spacing = "auto") {
    E <- as.numeric(E_MeV)
    if (length(E) < 2L) {
        stop("E_MeV must have >= 2 points, got ", length(E), call. = FALSE)
    }
    if (any(!is.finite(E))) stop("E_MeV contains non-finite values", call. = FALSE)
    if (any(diff(E) <= 0)) stop("E_MeV must be strictly increasing", call. = FALSE)
    if (any(E <= 0)) stop("E_MeV must contain positive energies", call. = FALSE)

    n_basis <- as.integer(n_basis)
    spline_order <- as.integer(spline_order)
    if (spline_order < 2L || spline_order > 8L) {
        stop("spline_order must be in [2, 8], got ", spline_order, call. = FALSE)
    }
    if (n_basis < spline_order) {
        stop("n_basis (N_s = ", n_basis, ") must be >= spline_order (p = ",
             spline_order, ") for a clamped B-spline basis", call. = FALSE)
    }
    degree <- spline_order - 1L
    emin <- E[1L]; emax <- E[length(E)]
    knot_spacing <- .mlem_bs_resolve_spacing(knot_spacing, E)

    n_interior <- n_basis - spline_order
    if (n_interior > 0L) {
        keep <- seq_len(n_interior + 2L)[-c(1L, n_interior + 2L)]
        if (knot_spacing == "log") {
            # numpy.geomspace(emin, emax, n_interior + 2)
            t_int <- exp(seq(log(emin), log(emax),
                             length.out = n_interior + 2L))[keep]
        } else {
            t_int <- seq(emin, emax, length.out = n_interior + 2L)[keep]
        }
    } else {
        t_int <- numeric(0L)
    }
    knots <- c(rep(emin, spline_order), t_int, rep(emax, spline_order))

    n <- length(E)
    B <- matrix(0.0, nrow = n, ncol = n_basis)
    x <- pmin(pmax(E, emin), emax)          # numpy clip, guards float fuzz
    for (i in seq_len(n)) {
        B[i, ] <- .mlem_bs_basis_eval(x[i], knots, degree, n_basis)
    }
    B
}

# core/unfold_mlem_bs.py:_interior_knots
.mlem_bs_interior_knots <- function(E, n_basis, spline_order, knot_spacing) {
    n_interior <- as.integer(n_basis) - as.integer(spline_order)
    if (n_interior <= 0L) return(numeric(0L))
    emin <- E[1L]; emax <- E[length(E)]
    keep <- seq_len(n_interior + 2L)[-c(1L, n_interior + 2L)]
    if (knot_spacing == "log") {
        exp(seq(log(emin), log(emax), length.out = n_interior + 2L))[keep]
    } else {
        seq(emin, emax, length.out = n_interior + 2L)[keep]
    }
}

.mlem_bs_basis_eval <- function(x, knots, degree, n_basis) {
    # Clean Cox-de Boor recursion. Returns the n_basis-long vector of basis
    # function values N_i(x) at point x. Partition of unity holds for
    # clamped B-splines at every interior point.
    n_knots <- length(knots)
    # Find span index i (1-based) such that knots[i] <= x < knots[i+1],
    # and i in [degree+1, n_basis] (interior span).
    if (x >= knots[n_knots]) {
        i <- n_basis
    } else if (x <= knots[1L]) {
        i <- degree + 1L
    } else {
        # Binary search could be used; linear is fine for small knot vectors
        i <- degree + 1L
        for (k in (degree + 1L):(n_basis)) {
            if (knots[k] <= x && x < knots[k + 1L]) { i <- k; break }
        }
    }
    # Degree-0 basis functions (only the one containing x is non-zero)
    N0 <- numeric(degree + 1L)
    for (j in 0:degree) {
        kk <- i - degree + j  # 1-based knot index
        if (x >= knots[kk] && x < knots[kk + 1L]) N0[j + 1L] <- 1.0
    }
    # Special case: x at the right endpoint
    if (x == knots[n_knots]) N0[degree + 1L] <- 1.0
    # Cox-de Boor recursion. At degree 0, only N0[degree+1] (= N_{i, 0}) is
    # 1. At degree d, the non-zero offsets are j = (degree-d)..degree (the
    # LAST d+1 entries of N). We use a separate N_prev snapshot to avoid
    # overwriting d-1 values we still need.
    N <- N0
    for (d in 1:degree) {
        N_prev <- N  # snapshot of d-1 values
        for (j in (degree - d):degree) {
            kk <- i - degree + j  # 1-based knot index
            denom_a <- knots[kk + d] - knots[kk]
            denom_b <- knots[kk + d + 1L] - knots[kk + 1L]
            term_a <- if (denom_a > 0) (x - knots[kk]) / denom_a * N_prev[j + 1L] else 0
            right_N <- if (j + 2L <= degree + 1L) N_prev[j + 2L] else 0
            term_b <- if (denom_b > 0) (knots[kk + d + 1L] - x) / denom_b * right_N else 0
            N[j + 1L] <- term_a + term_b
        }
        # Zero out the FIRST (degree-d) entries -- they fall outside the
        # active support at degree d.
        if (degree - d >= 1L) {
            for (j in 0:(degree - d - 1L)) {
                N[j + 1L] <- 0
            }
        }
    }
    # Map N[1..degree+1] to basis vector at positions (i-degree)..(i)
    basis <- numeric(n_basis)
    for (j in 0:degree) {
        idx <- i - degree + j  # 1-based position in basis vector
        # idx is in [i-degree, i] which should be in [1, n_basis]
        if (idx >= 1L && idx <= n_basis) basis[idx] <- N[j + 1L]
    }
    basis
}

#' Build B-spline basis matrix (clamped)
#'
#' @param E_MeV Numeric energy grid.
#' @param n_basis Integer; number of basis functions \eqn{N_s}.
#' @param knot_spacing Character; \code{"auto"}, \code{"uniform"} or
#'   \code{"log"}.
#' @param spline_order Integer B-spline order \eqn{p} (degree \eqn{p-1});
#'   default \code{4L} (cubic, as in the paper).
#' @return Numeric matrix of shape \code{(length(E_MeV), n_basis)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' B <- build_bspline_basis(E, 10L, "log")
#' dim(B)
build_bspline_basis <- function(E_MeV, n_basis, knot_spacing = "auto",
                                spline_order = 4L) {
    .mlem_bs_basis(E_MeV, n_basis, spline_order = spline_order,
                   knot_spacing = knot_spacing)
}

#' Second-difference matrix (for MLEM-BS penalty)
#'
#' Returns the second-difference operator of shape \eqn{(n-2, n)} with
#' stencil \code{[1, -2, 1]}.
#'
#' @param n Integer; length of the coefficient vector.
#' @return Numeric matrix of shape \code{(n-2, n)}.
#' @export
#' @examples
#' D2 <- second_difference_matrix(6L)
#' dim(D2)
second_difference_matrix <- function(n) {
    n <- as.integer(n)
    if (n < 3L) stop("second_difference_matrix requires n >= 3, got ", n,
                     call. = FALSE)
    as.matrix(create_derivative_matrix(n, 2L))
}

#' K_S goodness-of-fit statistic
#'
#' \eqn{K_S = | sum_i (n_i - model_i)^2 / sum_i model_i - 1 |}
#'
#' @param measurements Numeric vector.
#' @param estimate Numeric vector (model prediction).
#' @return Numeric scalar.
#' @export
#' @examples
#' ks_statistic(c(1, 2, 3), c(1.1, 1.9, 3.1))
ks_statistic <- function(measurements, estimate) {
    measurements <- as.numeric(measurements)
    estimate <- as.numeric(estimate)
    denom <- sum(estimate)
    if (denom <= .MB_TINY) return(Inf)
    resid <- measurements - estimate
    abs(sum(resid * resid) / denom - 1.0)
}

# Pearson chi-square diagnostic (core/unfold_mlem_bs.py:_chi2_pearson)
.mlem_bs_chi2_pearson <- function(measurements, model) {
    measurements <- as.numeric(measurements)
    model <- as.numeric(model)
    mask <- model > .MB_COUNT_FLOOR
    if (!any(mask)) return(Inf)
    resid <- measurements[mask] - model[mask]
    sum(resid * resid / model[mask])
}

# core/unfold_mlem_bs.py:_sieve_projection
.mlem_bs_sieve_projection <- function(B, x0) {
    ns <- ncol(B)
    mean_x0 <- if (length(x0)) mean(as.numeric(x0)) else 0.0
    fallback <- rep(max(mean_x0, .MB_COUNT_FLOOR), ns)
    b0 <- tryCatch(as.numeric(lsei::nnls(B, as.numeric(x0))$x),
                   error = function(e) {
                       # unconstrained least squares clipped at zero
                       v <- tryCatch(
                           as.numeric(solve(crossprod(B),
                                            crossprod(B, as.numeric(x0)))),
                           error = function(e2) rep(NA_real_, ns))
                       pmax(v, 0.0)
                   })
    if (length(b0) != ns || any(!is.finite(b0))) return(fallback)
    if (!any(b0 > 0)) return(fallback)
    b0
}

# core/unfold_mlem_bs.py:_resolve_beta.  The R API keeps the documented
# convention that a zero argument means "not supplied" (Python uses None),
# so both arguments may not be non-zero at the same time.
.mlem_bs_resolve_beta <- function(beta, beta_relative, RB) {
    beta <- as.numeric(beta)[1L]
    beta_relative <- as.numeric(beta_relative)[1L]
    if (beta != 0 && beta_relative != 0) {
        stop("Provide either 'beta' (absolute, paper Eq. 4) or ",
             "'beta_relative' (scaled by mean sensitivity), not both",
             call. = FALSE)
    }
    if (beta_relative != 0) {
        if (beta_relative < 0) {
            stop("beta_relative must be >= 0, got ", beta_relative,
                 call. = FALSE)
        }
        s_bar <- if (length(RB)) mean(colSums(RB)) else 0.0
        return(list(beta_eff = beta_relative * s_bar,
                    beta_relative = beta_relative))
    }
    beta_eff <- beta
    if (beta_eff < 0) stop("beta must be >= 0, got ", beta_eff, call. = FALSE)
    list(beta_eff = beta_eff, beta_relative = NULL)
}

# core/unfold_mlem_bs.py:_argmin_index (0-based index of the first minimum,
# non-finite values never minimal)
.mlem_bs_argmin <- function(values) {
    best_i <- 1L
    best_v <- Inf
    for (i in seq_along(values)) {
        v <- values[i]
        if (is.finite(v) && v < best_v) { best_v <- v; best_i <- i }
    }
    best_i
}

# core/unfold_mlem_bs.py:_mlem_bs_iterate
.mlem_bs_iterate <- function(RB, b, b_coef, pen_mat, beta, max_iterations,
                             tolerance, ks_min_mode, ks_patience) {
    s_s <- colSums(RB)                       # MLEM sensitivity sum_i (RB)_is
    s_s_safe <- pmax(s_s, .MB_TINY)

    best_coef <- b_coef
    best_ks <- ks_statistic(b, as.numeric(RB %*% b_coef))
    ks_history <- best_ks
    argmin <- 1L                              # index of min in ks_history
    converged <- FALSE
    iterations <- 0L

    for (k in seq_len(max_iterations)) {
        # forward projection: sum_s' (RB)_is' b_s'^(k)
        fwd <- pmax(as.numeric(RB %*% b_coef), .MB_COUNT_FLOOR)
        ratio <- b / fwd
        # multiplicative ML correction
        correction <- as.numeric(crossprod(RB, ratio))
        # penalized denominator: sum_i (RB)_is + beta * dP/db_s  (Eq. 4)
        if (beta != 0 && !is.null(pen_mat)) {
            grad_P <- as.numeric(pen_mat %*% b_coef)
            denom <- s_s + beta * grad_P
            # one-step-late guard: fall back to the unpenalized sensitivity
            # where the penalty gradient would flip the denominator
            # non-positive (standard OSL safeguard)
            denom <- ifelse(denom > .MB_TINY, denom, s_s_safe)
        } else {
            denom <- s_s_safe
        }
        denom <- pmax(denom, .MB_TINY)

        b_new <- pmax(b_coef * (correction / denom), 0.0)  # sieve: b_s >= 0
        iterations <- k

        # K_S statistic of the new iterate (Eq. 6)
        ks <- ks_statistic(b, as.numeric(RB %*% b_new))
        ks_history <- c(ks_history, ks)

        if (ks_min_mode) {
            if (is.finite(ks) && ks < best_ks) {
                best_ks <- ks
                best_coef <- b_new
                argmin <- k + 1L
            } else if ((k - (argmin - 1L)) >= ks_patience) {
                b_coef <- b_new
                break
            }
        } else {
            best_coef <- b_new
        }

        # convergence on relative coefficient change
        norm_old <- sqrt(sum(b_coef^2))
        diff <- sqrt(sum((b_new - b_coef)^2)) / (norm_old + .MB_TINY)
        b_coef <- b_new
        if (diff < tolerance) { converged <- TRUE; break }
    }

    list(best_coef = best_coef, iterations = iterations,
         converged = converged, ks_history = ks_history, argmin = argmin)
}

# core/unfold_mlem_bs.py:_auto_ns_grid
.mlem_bs_auto_ns_grid <- function(n_energy_bins, spline_order, n_points = 8L) {
    lo <- max(spline_order + 1L, 8L)
    hi <- max(min(as.integer(n_energy_bins), 150L), lo)
    if (hi <= lo) return(lo)
    grid <- exp(seq(log(lo), log(hi), length.out = as.integer(n_points)))
    sizes <- sort(unique(c(round(grid), lo, hi)))
    sizes[as.integer(sizes) >= spline_order + 1L]
}

# core/unfold_mlem_bs.py:solve_mlem_bs_full
.mlem_bs_solve <- function(A, b, x0 = NULL, E_MeV = NULL, n_basis = NULL,
                           beta = 0.0, beta_relative = 0.0,
                           max_iterations = 1000L, tolerance = 1e-6,
                           knot_spacing = "auto", spline_order = 4L,
                           auto_params = FALSE, ks_patience = 50L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    spline_order <- as.integer(spline_order)
    max_iterations <- as.integer(max_iterations)
    validate_system(A, b, x0 = x0, max_iterations = max_iterations,
                    tolerance = tolerance)

    n_bins <- ncol(A)
    if (is.null(E_MeV)) {
        stop("E_MeV is required to build the B-spline basis", call. = FALSE)
    }
    E <- as.numeric(E_MeV)
    if (length(E) != n_bins) {
        stop("Length of E_MeV (", length(E),
             ") must match the number of columns of A (", n_bins, ")",
             call. = FALSE)
    }
    if (is.null(x0)) x0 <- rep(0.5, n_bins)
    x0 <- as.numeric(x0)
    knot_spacing <- .mlem_bs_resolve_spacing(knot_spacing, E)

    # --- single configuration (Python's _run_single closure) -------------
    run_single <- function(ns, beta_eff, beta_rel, max_iter, ks_mode) {
        B <- .mlem_bs_basis(E, ns, spline_order = spline_order,
                            knot_spacing = knot_spacing)
        RB <- A %*% B
        pen_mat <- NULL
        if (beta_eff != 0) {
            D2 <- second_difference_matrix(ns)
            pen_mat <- 2.0 * as.matrix(crossprod(D2))   # 2 (D2)^T D2
        }
        b0 <- .mlem_bs_sieve_projection(B, x0)
        it <- .mlem_bs_iterate(RB, b, b0, pen_mat, beta_eff, max_iter,
                               tolerance, ks_mode, ks_patience)
        spectrum <- as.numeric(B %*% it$best_coef)
        model <- as.numeric(RB %*% it$best_coef)
        list(spectrum = spectrum, coefficients = it$best_coef,
             iterations = as.integer(it$iterations),
             converged = it$converged, ks_history = it$ks_history,
             ks_final = if (length(it$ks_history))
                            as.numeric(it$ks_history[length(it$ks_history)])
                        else Inf,
             chi2_pearson = .mlem_bs_chi2_pearson(b, model),
             n_basis = as.integer(ns), spline_order = spline_order,
             knot_spacing = knot_spacing,
             interior_knots = .mlem_bs_interior_knots(E, ns, spline_order,
                                                     knot_spacing),
             beta_effective = as.numeric(beta_eff),
             beta_relative = beta_rel, method = "MLEM-BS")
    }

    # --- auto selection of (N_s, beta, iterations) by minimizing K_S -----
    if (isTRUE(auto_params)) {
        ns_grid <- .mlem_bs_auto_ns_grid(n_bins, spline_order)
        scan_iter <- as.integer(min(max_iterations, 400L))
        candidates <- list()
        best <- NULL
        for (ns in ns_grid) {
            for (rho in .MB_AUTO_BETA_GRID) {
                beta_eff <- 0.0
                if (rho > 0) {
                    Btmp <- .mlem_bs_basis(E, ns, spline_order = spline_order,
                                           knot_spacing = knot_spacing)
                    beta_eff <- rho * mean(colSums(A %*% Btmp))
                }
                res <- run_single(ns, beta_eff, if (rho > 0) rho else 0.0,
                                  scan_iter, TRUE)
                entry <- list(n_basis = as.integer(ns),
                              beta_relative = as.numeric(rho),
                              ks = min(unlist(res$ks_history)),
                              ks_iteration = .mlem_bs_argmin(res$ks_history) - 1L,
                              iterations = res$iterations)
                candidates[[length(candidates) + 1L]] <- entry
                if (is.null(best) || entry$ks < best$ks) best <- entry
            }
        }
        chosen_ns <- as.integer(best$n_basis)
        chosen_rho <- as.numeric(best$beta_relative)
        beta_eff <- 0.0
        if (chosen_rho > 0) {
            Bc <- .mlem_bs_basis(E, chosen_ns, spline_order = spline_order,
                                 knot_spacing = knot_spacing)
            beta_eff <- chosen_rho * mean(colSums(A %*% Bc))
        }
        result <- run_single(chosen_ns, beta_eff,
                             if (chosen_rho > 0) chosen_rho else 0.0,
                             max_iterations, TRUE)
        result$auto_selection <- list(candidates = candidates, chosen = best,
                                      ns_grid = ns_grid,
                                      beta_relative_grid = .MB_AUTO_BETA_GRID)
        result$beta_relative <- if (chosen_rho > 0) chosen_rho else NULL
        result$beta_effective <- as.numeric(beta_eff)
        return(result)
    }

    # --- fixed parameters -------------------------------------------------
    if (is.null(n_basis)) {
        n_basis <- max(spline_order + 1L, min(n_bins %/% 2L, 40L))
    }
    n_basis <- as.integer(n_basis)
    rb_placeholder <- matrix(0.0, nrow = length(b), ncol = 1L)
    resolved <- .mlem_bs_resolve_beta(beta, beta_relative, rb_placeholder)
    beta_eff <- resolved$beta_eff
    beta_rel <- resolved$beta_relative
    # resolve beta against the actual RB (s_bar depends on the basis)
    if (!is.null(beta_rel)) {
        B <- .mlem_bs_basis(E, n_basis, spline_order = spline_order,
                            knot_spacing = knot_spacing)
        beta_eff <- beta_rel * mean(colSums(A %*% B))
    }
    run_single(n_basis, beta_eff, beta_rel, max_iterations, FALSE)
}

#' Poisson bootstrap confidence intervals for MLEM-BS
#'
#' Port of \code{bssunfold} core \code{unfold_mlem_bs.py} Eqs. 7-9: the
#' backward-reconstructed counts \eqn{n^{(0)} = A x^{(0)}} define Poisson
#' means, every bootstrap replicate \eqn{n^{*}} is unfolded with the same
#' solver configuration and the interval is formed from the
#' \eqn{\alpha/2}/\eqn{1-\alpha/2} quantiles.
#'
#' @param A Numeric response matrix (m x n).
#' @param x0 Numeric initial spectrum (length n) or \code{NULL}.
#' @param E_MeV Numeric energy grid (length n).
#' @param spectrum Numeric unfolded spectrum \eqn{x^{(0)}} (length n).
#' @param n_bootstrap Integer; number of replicates.
#' @param ci_alpha Numeric in \eqn{(0, 1)}; 0.05 gives a 95 percent interval.
#' @param random_state Integer seed or \code{NULL}.
#' @param ... Solver configuration forwarded to \code{\link{solve_mlem_bs}}.
#' @return A list with \code{ci_low}, \code{ci_high}, \code{ci_level},
#'   \code{bootstrap_samples}, \code{bootstrap_alpha}, \code{bootstrap_mean},
#'   \code{bootstrap_std}.
#' @keywords internal
.mlem_bs_bootstrap_ci <- function(A, x0, E_MeV, spectrum, n_bootstrap,
                                  ci_alpha, random_state, ...) {
    n_bootstrap <- as.integer(n_bootstrap)
    if (n_bootstrap < 1L) {
        stop("n_bootstrap must be >= 1, got ", n_bootstrap, call. = FALSE)
    }
    ci_alpha <- as.numeric(ci_alpha)
    if (!(ci_alpha > 0 && ci_alpha < 1)) {
        stop("ci_alpha must be in (0, 1), got ", ci_alpha, call. = FALSE)
    }
    if (!is.null(random_state)) set.seed(as.integer(random_state))
    # Eq. 8: backward reconstruction of the initial count spectrum
    n0 <- pmax(as.numeric(A %*% as.numeric(spectrum)), 0.0)
    n_bins <- length(spectrum)
    replicates <- matrix(0.0, nrow = n_bootstrap, ncol = n_bins)
    for (i in seq_len(n_bootstrap)) {
        # Eq. 9: Poisson resampling preserving counting statistics
        n_star <- as.numeric(stats::rpois(n = length(n0), lambda = n0))
        replicates[i, ] <- .mlem_bs_solve(A, n_star, x0 = x0, E_MeV = E_MeV,
                                         ...)$spectrum
    }
    list(ci_low = stats::quantile(replicates, probs = ci_alpha / 2,
                                  type = 7, names = FALSE),
         ci_high = stats::quantile(replicates,
                                   probs = 1 - ci_alpha / 2,
                                   type = 7, names = FALSE),
         ci_level = 1 - ci_alpha,
         bootstrap_samples = n_bootstrap,
         bootstrap_alpha = ci_alpha,
         bootstrap_mean = colMeans(replicates),
         bootstrap_std = if (n_bootstrap > 1L)
                             apply(replicates, 2, stats::sd)
                         else rep(0.0, n_bins))
}

#' Wrapper around \code{\link{solve_mlem_bs}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_mlem_bs
#' @param spline_order Integer; B-spline order (degree \code{spline_order -
#'   1}); the paper uses \code{4L}.
#' @param auto_params Logical; select \eqn{N_s}, penalty strength and
#'   iteration count by minimising \eqn{K_S}.
#' @param ks_patience Integer; patience of the \eqn{K_S}-based early stop.
#' @param bootstrap_ci Logical; add Poisson-bootstrap confidence intervals
#'   (Eqs. 7-9) to the output.
#' @param n_bootstrap Integer; number of bootstrap replicates.
#' @param ci_alpha Numeric; CI significance level.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_mlem_bs <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            n_basis = NULL, beta = 0.0,
                            beta_relative = 0.0,
                            max_iterations = 1000L, tolerance = 1e-6,
                            knot_spacing = "auto",
                            spline_order = 4L, auto_params = FALSE,
                            ks_patience = 50L, bootstrap_ci = FALSE,
                            n_bootstrap = 100L, ci_alpha = 0.05,
                            calculate_errors = FALSE,
                            noise_level = 0.01, n_montecarlo = 100L,
                            save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.5, n_energy_bins)
    diag_env <- new.env(parent = emptyenv())
    # Use closure to forward E_MeV and the basis/penalty configuration.
    # \code{run_unfolding()} drops the bins above \code{max_neutron_energy}
    # from A (they are the trailing ones, the grid is increasing) and pads the
    # returned spectrum, so the solver must see the matching energy prefix --
    # in Python the Detector already hands the core a truncated grid.
    solver_with_E <- function(A, b, x0 = NULL, ...) {
        nb <- ncol(A)
        E_use <- as.numeric(E_MeV)[seq_len(nb)]
        x_in <- if (is.null(x0)) x0_default else as.numeric(x0)
        x_in <- x_in[seq_len(nb)]
        out <- .mlem_bs_solve(A = A, b = b,
                              x0 = x_in,
                              E_MeV = E_use, n_basis = n_basis,
                              beta = beta, beta_relative = beta_relative,
                              max_iterations = max_iterations,
                              tolerance = tolerance,
                              knot_spacing = knot_spacing,
                              spline_order = spline_order,
                              auto_params = auto_params,
                              ks_patience = ks_patience)
        diag_env$out <- out
        out
    }
    x0_main <- if (is.null(initial_spectrum)) {
        x0_default
    } else {
        .normalize_initial(initial_spectrum, x0_default, n_energy_bins)
    }
    out <- run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver_with_E,
        solve_kwargs = list(),
        method_name = "MLEM-BS",
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)

    used <- diag_env$out
    if (!is.null(used)) {
        out <- c(out, list(ks_history = list(used$ks_history),
                           ks_final = used$ks_final,
                           chi2_pearson = used$chi2_pearson,
                           n_basis = used$n_basis,
                           spline_order = used$spline_order,
                           knot_spacing = used$knot_spacing,
                           interior_knots = list(used$interior_knots),
                           beta_effective = used$beta_effective,
                           beta_relative = used$beta_relative,
                           coefficients = list(used$coefficients)))
        if (!is.null(used$auto_selection)) {
            out$auto_selection <- used$auto_selection
        }
    }

    if (isTRUE(bootstrap_ci)) {
        sys <- .build_system(readings, detector_names, sensitivities)
        A_sys <- sys$A
        spec <- as.numeric(out$spectrum)
        E_boot <- as.numeric(E_MeV)
        x0_boot <- x0_main
        if (!is.null(max_neutron_energy) &&
            is.finite(as.numeric(max_neutron_energy)) &&
            as.numeric(max_neutron_energy) > 0) {
            keep <- E_boot <= as.numeric(max_neutron_energy)
            A_sys <- A_sys[, keep, drop = FALSE]
            E_boot <- E_boot[keep]
            spec <- spec[keep]
            x0_boot <- x0_boot[keep]
        }
        ns <- if (!is.null(used)) used$n_basis else
            max(spline_order + 1L, min(as.integer(n_energy_bins) %/% 2L, 40L))
        boot <- .mlem_bs_bootstrap_ci(
            A_sys, x0 = x0_boot, E_MeV = E_boot,
            spectrum = spec, n_bootstrap = n_bootstrap,
            ci_alpha = ci_alpha, random_state = random_state,
            n_basis = ns, beta = if (!is.null(used)) used$beta_effective else beta,
            beta_relative = if (is.null(used)) beta_relative else 0,
            max_iterations = max_iterations, tolerance = tolerance,
            knot_spacing = if (!is.null(used)) used$knot_spacing else knot_spacing,
            spline_order = spline_order)
        out <- c(out, boot)
    }
    out
}
