#' EPIC Tikhonov regularization unfolding
#'
#' R port of \code{bssunfold/core/unfold_epic.py} (the EPIC_LS package of
#' Ortega-Culaciati et al., JGR Solid Earth, 2021).
#' The Equal Posterior Information Condition (EPIC) selects the prior
#' variances of a regularization operator \eqn{H} such that the a posteriori
#' variances of the model parameters match the user supplied target
#' sigmas. The betas (natural logarithms of the reciprocal prior variances)
#' are obtained by solving the nonlinear system
#' \deqn{(\mathrm{diag}((P + H' \exp(\beta) H)^{-1}) - s^2)/s^2 = 0}{(diag((P +
#' H' exp(beta) H)^{-1}) - s^2)/s^2 = 0}
#' in the least-squares sense, with \eqn{P = A'A} the precision matrix and
#' \eqn{s} the target sigmas. Once the weights are known the (optionally
#' non-negative) regularized least squares problem is solved by stacking the
#' misfit and regularization blocks.
#'
#' The returned spectrum is a lethargy density (fluence per \eqn{d \ln E}),
#' like every other solver in this package: \code{A} and \code{b} arrive
#' already weighted by the per-bin log steps.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused (kept for API compatibility with the other solvers).
#' @param sigma_frac Numeric; fraction of max LS solution as target sigma. Default 0.1.
#' @param max_iterations Integer; max EPIC iterations. Each unit allows 100
#'   function evaluations of the EPIC weight search. Default 30.
#' @param tolerance Numeric; convergence threshold on the EPIC betas. Default 1e-4.
#' @param nonneg Logical; enforce non-negativity via NNLS. Default TRUE.
#' @return A list \code{list(spectrum, iterations, converged, betas,
#'   posterior_vars)}. \code{betas} holds the EPIC log-weights and
#'   \code{posterior_vars} the a posteriori variances they produce.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_epic(A, b, NULL, max_iterations = 10)
solve_epic <- function(A, b, x0 = NULL, sigma_frac = 0.1,
                        max_iterations = 30L, tolerance = 1e-4,
                        nonneg = TRUE) {
    A <- .epic_as_matrix(A, b)
    b <- as.numeric(b)
    weights <- .epic_weights(A, b, sigma_frac = sigma_frac,
                             max_iterations = max_iterations,
                             tolerance = tolerance)
    x <- .epic_final_solve(A, b, weights$Wx, weights$H,
                           rep(0, nrow(weights$H)), weights$Wh,
                           non_neg = nonneg)
    list(spectrum = x,
         iterations = weights$meta$iterations,
         converged = weights$meta$epic_converged,
         betas = weights$betas,
         posterior_vars = weights$posterior_vars)
}

# ---- internal helpers, mirroring the Python module ----

#' Dense double matrix used by the EPIC solver.
#' @noRd
.epic_as_matrix <- function(A, b) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    if (nrow(A) != length(b)) {
        stop("b length (", length(b), ") must match the number of A rows (",
             nrow(A), ")")
    }
    A
}

#' Port of beta_bounds.compute_bounds from EPIC_LS: the largest k such that
#' exp(k_center - k) + exp(k_center + k) is still well represented, keeping
#' `distance` away from that limit.
#' @noRd
.epic_beta_bounds <- function(k_center = 0, distance = 2) {
    eps <- .Machine$double.eps
    k_test <- 0.0
    for (i in seq_len(999999L)) {
        k_test <- k_test + 0.01
        if (abs(exp(k_center - k_test) + exp(k_center + k_test) -
                exp(k_center + k_test)) < eps) break
    }
    k_test <- k_test - distance
    c(k_center - k_test, k_center + k_test)
}

#' min-norm least squares solution, matching numpy.linalg.lstsq(rcond=None).
#' @noRd
.epic_lstsq <- function(A, b) {
    sv <- svd(A)
    if (length(sv$d) == 0L) return(rep(0, ncol(A)))
    cutoff <- max(dim(A)) * .Machine$double.eps * max(sv$d)
    keep <- sv$d > cutoff
    if (!any(keep)) return(rep(0, ncol(A)))
    utb <- as.numeric(crossprod(sv$u, b))
    as.numeric(sv$v[, keep, drop = FALSE] %*% (utb[keep] / sv$d[keep]))
}

#' Port of _default_target_sigmas(): scale from the naive least-squares
#' solution, falling back to the measurement scale when it is degenerate.
#' @noRd
.epic_default_target_sigmas <- function(A, b, n, sigma_frac) {
    x_ls <- tryCatch(.epic_lstsq(A, b), error = function(e) numeric(0))
    scale <- if (length(x_ls)) max(abs(x_ls)) else 1
    if (!is.finite(scale) || scale <= 0) {
        scale <- if (length(b)) max(abs(b)) else 1
        if (!is.finite(scale) || scale <= 0) scale <- 1
    }
    rep(sigma_frac * scale, n)
}

#' Port of _build_precision(): P = A' inv(Cx) A and the misfit weight Wx.
#' @noRd
.epic_precision <- function(A, noise_var = NULL) {
    m <- nrow(A)
    if (is.null(noise_var)) {
        list(P = crossprod(A), Wx = diag(m))
    } else {
        if (!is.finite(noise_var) || noise_var <= 0) {
            stop("noise_var must be positive, got ", noise_var)
        }
        inv_cx <- (1 / noise_var) * diag(m)
        list(P = crossprod(A, inv_cx %*% A),
             Wx = as.matrix(t(chol(inv_cx))))
    }
}

#' Port of _build_regularization_matrix(): order 0 is the identity, order 1/2
#' the plain finite difference operators (grid independent, legacy EPIC_LS
#' behaviour).
#' @noRd
.epic_reg_matrix <- function(n, order) {
    n <- as.integer(n); order <- as.integer(order)
    if (!order %in% c(0L, 1L, 2L)) {
        stop("Unsupported regularization_order: ", order,
             ". Use 0 (identity), 1 or 2.")
    }
    if (order == 0L) return(diag(n))
    as.matrix(create_derivative_matrix(n, order))
}

#' Port of calc_F from _calc_epic_ch: the EPIC residual vector at betas beta.
#' @noRd
.epic_calc_F <- function(beta, P, H, target_var) {
    if (any(!is.finite(beta))) {
        beta <- pmin(pmax(beta, -700), 700)
        beta[!is.finite(beta)] <- 0
    }
    M <- P + crossprod(H, exp(beta) * H)
    inv_diag <- .epic_inv_diag(M)
    (inv_diag - target_var) / target_var
}

#' Diagonal of the inverse of a symmetric matrix, with a robust fallback so a
#' numerically singular precision never reaches LAPACK as NaN/Inf.
#' @noRd
.epic_inv_diag <- function(M) {
    if (any(!is.finite(M))) return(rep(Inf, nrow(M)))
    d <- tryCatch({
        cm <- chol(M)
        as.numeric(diag(chol2inv(cm)))
    }, error = function(e) {
        tryCatch(as.numeric(diag(solve(M))), error = function(e2) {
            as.numeric(diag(Matrix::solve(M, Matrix::Diagonal(nrow(M)),
                                         method = "LU")))
        })
    })
    if (any(!is.finite(d))) return(rep(Inf, nrow(M)))
    d
}

#' Port of calc_JF from _calc_epic_ch: d F_i / d beta_j.
#' @noRd
.epic_calc_JF <- function(beta, P, H, target_var) {
    eb <- exp(pmin(pmax(beta, -700), 700))
    eb[!is.finite(eb)] <- 0
    M <- P + crossprod(H, eb * H)
    invA <- tryCatch(solve(M), error = function(e) {
        Mm <- as.matrix(M)
        sv <- svd(Mm)
        tol <- .Machine$double.eps * max(dim(Mm)) * sv$d[1]
        di <- ifelse(sv$d > tol, 1 / pmax(sv$d, tol), 0)
        sv$v %*% (di * t(sv$u))
    })
    if (any(!is.finite(invA))) {
        return(matrix(0, nrow = length(target_var), ncol = length(beta)))
    }
    B <- H %*% invA
    # JF[i, j] = -(1/target_var[i]) * exp(beta[j]) * B[j, i]^2
    JF <- -t(eb * (B * B))
    JF * (1 / target_var)
}

#' scipy.optimize.least_squares(method='trf', tr_solver='exact') port.
#'
#' Shared by \code{\link{solve_epic}} (which uses the \code{box} argument for
#' the EPIC beta bounds) and \code{\link{solve_express}}. Both Python methods
#' call \code{scipy.optimize.least_squares}, so this reproduces that algorithm
#' rather than a generic Levenberg-Marquardt: the trust-region subproblem is
#' solved with the Moré-Sorenson secular equation on the SVD of the scaled
#' Jacobian (\code{solve_lsq_trust_region}), the radius is updated with
#' \code{update_tr_radius} and termination with scipy's \code{check_termination}
#' (which needs BOTH \code{dF < ftol * F} and \code{ratio > 0.25}). That matters
#' because the published Express answer is an ftol stopping point, not a
#' converged optimum, so the iteration path has to be replicated.
#'
#' \code{fun(x)} returns the residual vector of length \code{m}; the objective
#' is \code{0.5 * sum(fun^2)}. \code{jac = NULL} uses scipy's \code{jac =
#' '2-point'} differencing (\code{rel_step = sqrt(eps)} times
#' \code{sign(x) * max(1, |x|)}, forward difference, with the denominator taken
#' as \code{(x + h) - x}). \code{nfev} counts the initial evaluation plus every
#' trial evaluation of \code{fun}, as scipy's local counter does; Jacobian
#' evaluations are not counted.
#' @noRd
.bss_lsq_fd <- function(fun, x, f0) {
    rs <- .Machine$double.eps^0.5
    n <- length(x); m <- length(f0)
    J <- matrix(0, nrow = m, ncol = n)
    for (i in seq_len(n)) {
        ## np.sign() replacement: scipy uses (x0 >= 0) * 2 - 1, so 0 steps up.
        sg <- if (x[i] >= 0) 1 else -1
        h <- rs * sg * max(1, abs(x[i]))
        xi <- x[i] + h
        xp <- x
        xp[i] <- xi
        J[, i] <- (as.numeric(fun(xp)) - f0) / (xi - x[i])
    }
    J
}

#' evaluate_quadratic(): 0.5 * s' (J' J) s + g' s.
#' @noRd
.bss_lsq_quad <- function(J, g, s) {
    Js <- as.numeric(J %*% s)
    0.5 * sum(Js * Js) + sum(s * g)
}

#' compute_jac_scale(): column scale from the Jacobian norms.
#' @noRd
.bss_lsq_jac_scale <- function(J, scale_inv_old = NULL) {
    si <- sqrt(colSums(J * J))
    if (is.null(scale_inv_old)) {
        si[si == 0] <- 1
    } else {
        si <- pmax(si, scale_inv_old)
    }
    list(scale = 1 / si, scale_inv = si)
}

#' update_tr_radius(): new trust radius and step quality ratio.
#' @noRd
.bss_lsq_radius <- function(Delta, actual, predicted, step_norm, bound_hit) {
    ratio <- if (predicted > 0) actual / predicted
             else if (predicted == 0 && actual == 0) 1 else 0
    if (ratio < 0.25) {
        Delta <- 0.25 * step_norm
    } else if (ratio > 0.75 && bound_hit) {
        Delta <- Delta * 2.0
    }
    list(Delta = Delta, ratio = ratio)
}

#' check_termination(): scipy's combined ftol/xtol test (status 2, 3 or 4).
#' @noRd
.bss_lsq_term <- function(dF, F, dx_norm, x_norm, ratio, ftol, xtol) {
    ftol_ok <- dF < ftol * F && ratio > 0.25
    xtol_ok <- dx_norm < xtol * (xtol + x_norm)
    if (ftol_ok && xtol_ok) 4L else if (ftol_ok) 2L else if (xtol_ok) 3L else NULL
}

#' solve_lsq_trust_region(): the exact More-Sorenson subproblem solver.
#' @noRd
.bss_lsq_subproblem <- function(n, m, uf, s, V, Delta, initial_alpha = 0,
                                rtol = 0.01, max_iter = 10L) {
    eps <- .Machine$double.eps
    k <- length(s)
    suf <- s * uf
    full_rank <- if (m >= n) (s[k] > eps * m * s[1]) else FALSE
    if (full_rank) {
        p <- -V %*% (uf / s)
        if (sqrt(sum(p * p)) <= Delta) {
            return(list(p = as.numeric(p), alpha = 0, n_iter = 0L))
        }
    }
    alpha_upper <- sqrt(sum(suf * suf)) / Delta
    if (full_rank) {
        denom <- s * s
        p_norm <- sqrt(sum((suf / denom)^2))
        phi <- p_norm - Delta
        phi_prime <- -sum(suf^2 / denom^3) / p_norm
        alpha_lower <- -phi / phi_prime
    } else {
        alpha_lower <- 0
    }
    fallback <- max(0.001 * alpha_upper, sqrt(alpha_lower * alpha_upper))
    alpha <- if (initial_alpha == 0 && !full_rank) fallback else initial_alpha
    it <- 0L
    for (i in seq_len(max_iter)) {
        it <- i
        if (alpha < alpha_lower || alpha > alpha_upper) alpha <- fallback
        denom <- s * s + alpha
        p_norm <- sqrt(sum((suf / denom)^2))
        phi <- p_norm - Delta
        phi_prime <- -sum(suf^2 / denom^3) / p_norm
        if (phi < 0) alpha_upper <- alpha
        ratio <- phi / phi_prime
        alpha_lower <- max(alpha_lower, alpha - ratio)
        alpha <- alpha - (phi + Delta) * ratio / Delta
        if (abs(phi) < rtol * Delta) break
    }
    p <- -V %*% (suf / (s * s + alpha))
    ## Rescale onto the sphere so the step never leaves the trust region.
    p <- p * (Delta / sqrt(sum(p * p)))
    list(p = as.numeric(p), alpha = alpha, n_iter = it)
}

#' Main trust-region loop, mirroring scipy's \code{trf_no_bounds} plus the box
#' clipping that \code{trf_bounds} would apply for finite bounds.
#' @noRd
.bss_lsq_trf <- function(fun, x0, m, n, jac = NULL, ftol = 1e-8, xtol = 1e-8,
                         gtol = 1e-8, max_nfev = 100L, x_scale = 1,
                         box = NULL) {
    x <- as.numeric(x0)
    if (!is.null(box)) x <- pmin(pmax(x, box$lo), box$hi)
    f_of <- function(v) {
        r <- suppressWarnings(as.numeric(fun(v)))
        if (length(r) != m) return(rep(NaN, m))
        r
    }
    f <- f_of(x)
    nfev <- 1L
    cost <- 0.5 * sum(f * f)
    J <- if (is.null(jac)) .bss_lsq_fd(fun, x, f) else as.matrix(jac(x))
    g <- as.numeric(crossprod(J, f))
    jac_scale <- is.character(x_scale) && x_scale == "jac"
    if (jac_scale) {
        sc <- .bss_lsq_jac_scale(J)
        scale <- sc$scale; scale_inv <- sc$scale_inv
    } else {
        scale <- rep(as.numeric(x_scale), length.out = n)
        scale_inv <- rep(1 / as.numeric(x_scale), length.out = n)
    }
    Delta <- sqrt(sum((x * scale_inv)^2))
    if (Delta == 0) Delta <- 1
    alpha <- 0
    status <- NULL
    actual_reduction <- -1
    iteration <- 0L
    while (TRUE) {
        g_norm <- if (length(g)) max(abs(g)) else 0
        if (g_norm < gtol) status <- 1L
        if (!is.null(status) || nfev == max_nfev) break
        d <- scale
        g_h <- d * g
        J_h <- J * matrix(d, nrow = m, ncol = n, byrow = TRUE)
        sv <- svd(J_h)
        uf <- as.numeric(crossprod(sv$u, f))
        actual_reduction <- -1
        while (actual_reduction <= 0 && nfev < max_nfev) {
            step_h <- .bss_lsq_subproblem(n, m, uf, sv$d, sv$v, Delta,
                                          initial_alpha = alpha)$p
            predicted <- -.bss_lsq_quad(J_h, g_h, step_h)
            step <- d * step_h
            x_new <- x + step
            if (!is.null(box)) x_new <- pmin(pmax(x_new, box$lo), box$hi)
            f_new <- f_of(x_new)
            nfev <- nfev + 1L
            step_h_norm <- sqrt(sum(step_h * step_h))
            if (!all(is.finite(f_new))) {
                Delta <- 0.25 * step_h_norm
                next
            }
            cost_new <- 0.5 * sum(f_new * f_new)
            actual_reduction <- cost - cost_new
            ur <- .bss_lsq_radius(Delta, actual_reduction, predicted,
                                  step_h_norm, step_h_norm > 0.95 * Delta)
            Delta_new <- ur$Delta
            step_norm <- sqrt(sum(step * step))
            ts <- .bss_lsq_term(actual_reduction, cost, step_norm,
                                sqrt(sum(x * x)), ur$ratio, ftol, xtol)
            if (!is.null(ts)) {
                status <- ts
                break
            }
            alpha <- alpha * Delta / Delta_new
            Delta <- Delta_new
        }
        if (actual_reduction > 0) {
            x <- x_new
            f <- f_new
            cost <- cost_new
            J <- if (is.null(jac)) .bss_lsq_fd(fun, x, f) else as.matrix(jac(x))
            g <- as.numeric(crossprod(J, f))
            if (jac_scale) {
                sc <- .bss_lsq_jac_scale(J, scale_inv)
                scale <- sc$scale; scale_inv <- sc$scale_inv
            }
        } else {
            actual_reduction <- 0
        }
        iteration <- iteration + 1L
    }
    list(x = x, cost = cost, fun = f, nfev = nfev, nit = iteration,
         status = if (is.null(status)) 0L else status,
         success = !is.null(status) && status > 0)
}

#' Port of _calc_epic_ch(): solve the EPIC condition for the log prior
#' variances. Python runs scipy.optimize.least_squares twice (a homogeneous
#' 1-D search with TolX1/TolFun1/TolG1, then the full heterogeneous betas with
#' TolX2/TolFun2/TolG2, x_scale = 'jac' and tr_solver = 'exact'), both confined
#' to the representability bounds from _compute_bounds().
#' @noRd
.epic_calc_ch <- function(P, H, target_sigmas, homogeneous_step = TRUE,
                          beta_shift_k = 0, beta_distance = 2,
                          max_iterations = 30L, tolerance = 1e-4) {
    Nh <- nrow(H); Nm <- ncol(H)
    target_var <- as.numeric(target_sigmas)^2
    bounds <- .epic_beta_bounds(beta_shift_k, beta_distance)
    # Python: X0 = ones(Nh) * (bounds[0] + bounds[1]) / 2
    X0 <- rep((bounds[1] + bounds[2]) / 2, Nh)
    # scipy tolerances from _calc_epic_ch (Nh <= Nm for every unfolding grid
    # here, because the first-derivative operator has one row fewer than
    # columns); the R 'tolerance' argument may only tighten them.
    tol <- min(as.numeric(tolerance), 1e-4)
    if (Nh > Nm) {
        t1 <- c(1e-6, 1e-6, 1e-6); t2 <- c(1e-6, 1e-6, 1e-8)
    } else {
        t1 <- c(1e-6, 1e-6, 1e-6); t2 <- c(1e-8, 1e-8, 1e-10)
    }
    t1 <- pmin(t1, tol); t2 <- pmin(t2, tol)
    box <- list(lo = rep(bounds[1], Nh), hi = rep(bounds[2], Nh))
    fun <- function(beta) .epic_calc_F(beta, P, H, target_var)
    jac <- function(beta) .epic_calc_JF(beta, P, H, target_var)

    if (homogeneous_step) {
        sol0 <- .bss_lsq_trf(
            function(k) fun(as.numeric(k) + X0), 0, m = Nm, n = 1L,
            ftol = t1[2], xtol = t1[1], gtol = t1[3], max_nfev = 100L,
            x_scale = 1, box = list(lo = bounds[1], hi = bounds[2]))
        Xnext <- as.numeric(sol0$x) + X0
    } else {
        Xnext <- X0
    }

    sol <- .bss_lsq_trf(fun, Xnext, m = Nm, n = Nh, jac = jac,
                        ftol = t2[2], xtol = t2[1], gtol = t2[3],
                        max_nfev = max(1L, as.integer(max_iterations)) * 100L,
                        x_scale = "jac", box = box)
    list(beta = sol$x, cost = sol$cost, iterations = sol$nfev,
         converged = sol$success)
}

#' Port of _epic_weights(): EPIC log-weights plus the matrices of the final
#' regularized least squares problem.
#' @noRd
.epic_weights <- function(A, b, target_sigmas = NULL, sigma_frac = 0.1,
                          regularization_order = 1L, noise_var = NULL,
                          homogeneous_step = TRUE,
                          beta_shift_k = 0, beta_distance = 2,
                          max_iterations = 30L, tolerance = 1e-4) {
    n <- ncol(A)
    prec <- .epic_precision(A, noise_var)
    P <- prec$P; Wx <- prec$Wx
    H <- .epic_reg_matrix(n, regularization_order)

    if (is.null(target_sigmas)) {
        target_sigmas <- .epic_default_target_sigmas(A, b, n, sigma_frac)
    } else {
        target_sigmas <- as.numeric(target_sigmas)
    }
    if (length(target_sigmas) != n) {
        stop("target_sigmas length (", length(target_sigmas),
             ") must match the number of parameters (", n, ")")
    }
    if (any(!is.finite(target_sigmas)) || any(target_sigmas <= 0)) {
        stop("target_sigmas must be finite and strictly positive")
    }

    sol <- .epic_calc_ch(P, H, target_sigmas,
                         homogeneous_step = homogeneous_step,
                         beta_shift_k = beta_shift_k,
                         beta_distance = beta_distance,
                         max_iterations = max_iterations,
                         tolerance = tolerance)
    beta <- sol$beta
    posterior_vars <- .epic_inv_diag(P + crossprod(H, exp(beta) * H))
    list(Wx = Wx, H = H, Wh = diag(exp(beta / 2), nrow = nrow(H)),
         betas = beta,
         posterior_vars = posterior_vars,
         target_sigmas = target_sigmas,
         ncol_design = n,
         meta = list(epic_converged = all(is.finite(sol$cost)) &&
                         sol$cost < 1e300,
                     epic_cost = sol$cost,
                     epic_nfev = sol$iterations,
                     beta_min = min(beta), beta_max = max(beta),
                     target_sigmas = target_sigmas,
                     iterations = sol$iterations))
}

#' Port of _final_solve(): stack misfit and regularization blocks into one
#' least squares problem, optionally with a non-negativity constraint.
#' @noRd
.epic_final_solve <- function(A, b, Wx, H, ho, Wh, non_neg = TRUE) {
    WxG <- Wx %*% A
    Wxd <- as.numeric(Wx %*% b)
    WhH <- Wh %*% H
    Whho <- as.numeric(Wh %*% ho)
    F <- rbind(WxG, WhH)
    D <- c(Wxd, Whho)
    if (any(!is.finite(F)) || any(!is.finite(D))) {
        return(rep(0, ncol(A)))
    }
    x <- if (isTRUE(non_neg)) {
        tryCatch(as.numeric(lsei::nnls(F, D)$x),
                 error = function(e) as.numeric(.epic_lstsq(F, D)))
    } else {
        .epic_lstsq(F, D)
    }
    x <- as.numeric(x)
    if (length(x) != ncol(A)) x <- rep(0, ncol(A))
    pmax(ifelse(is.finite(x), x, 0), 0)
}

#' Wrapper around \code{\link{solve_epic}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_epic
#' @export
unfold_epic <- function(detector_names, n_energy_bins, E_MeV,
                           sensitivities, cc_icrp116, save_result_callback,
                           readings, initial_spectrum = NULL,
                           sigma_frac = 0.1, max_iterations = 30L,
                           tolerance = 1e-4, nonneg = TRUE,
                           calculate_errors = FALSE,
                           noise_level = 0.01, n_montecarlo = 100L,
                           save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    # Python computes the EPIC weights once from the (noise free) system and
    # reuses them for every Monte-Carlo sample.
    sys <- .build_system(readings, detector_names, sensitivities)
    weights <- .epic_weights(sys$A, sys$b, sigma_frac = sigma_frac,
                             max_iterations = max_iterations,
                             tolerance = tolerance)

    epic_solver <- function(A, b, x0 = NULL) {
        w <- if (ncol(A) == weights$ncol_design) {
            weights
        } else {
            # Energy cutoff trimmed the design: rebuild weights for it.
            .epic_weights(A, b, sigma_frac = sigma_frac,
                          max_iterations = max_iterations,
                          tolerance = tolerance)
        }
        ho <- rep(0, nrow(w$H))
        spectrum <- .epic_final_solve(A, b, w$Wx, w$H, ho, w$Wh,
                                      non_neg = nonneg)
        list(spectrum = spectrum,
             iterations = w$meta$iterations,
             converged = w$meta$epic_converged,
             betas = w$betas,
             posterior_vars = w$posterior_vars)
    }

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = epic_solver,
        solve_kwargs = list(),
        method_name = "EPIC",
        extra_output = list(sigma_frac = sigma_frac,
                            regularization_order = 1L,
                            non_neg = as.logical(nonneg),
                            epic_converged = weights$meta$epic_converged,
                            epic_cost = weights$meta$epic_cost,
                            beta_min = weights$meta$beta_min,
                            beta_max = weights$meta$beta_max,
                            target_sigmas = weights$meta$target_sigmas),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
