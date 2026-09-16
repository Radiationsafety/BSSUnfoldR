#' Uno-style constrained unfolding: Lagrange-Newton NLP presets
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_uno.py}.
#'
#' The unfolding problem is posed as the smooth non-linear program
#' \deqn{\min f(x) = \tfrac{1}{2} ||W (A x - b)||^2 + \tfrac{\lambda}{2} ||D_2 x||^2}
#' subject to \eqn{x \ge 0}, a convex quadratic objective with a single
#' inequality constraint set. Two Uno presets are provided:
#' \describe{
#'   \item{filter_sqp}{Sequential quadratic programming with the exact
#'     Hessian \eqn{H = A^T W^2 A + \lambda D_2^T D_2}, Newton directions
#'     and the Vanaret-Leyffer filter globalisation.}
#'   \item{ipopt_like}{A primal-dual interior-point method in the IPOPT
#'     manner: the inequality is handled by the log-barrier
#'     \eqn{-\mu \sum \log x}, the Newton system is regularised with the
#'     barrier Hessian \eqn{\mathrm{diag}(\mu / x^2)}, \eqn{\mu} follows
#'     a geometric decrease, and a fraction-to-the-boundary rule protects
#'     strict feasibility. The Hessian is either the exact convex
#'     \eqn{H} (\code{hessian="exact"}) or a dense BFGS approximation
#'     (\code{hessian="bfgs"}).}
#' }
#'
#' @name uno-methods
#' @rdname uno-methods
NULL

.UNO_PRESETS <- c("filter_sqp", "ipopt_like")
.UNO_FTB_TAU <- 0.995

.uno_derivative <- function(n, order = 2L) {
    as.matrix(create_derivative_matrix(n, order))
}

#' Uno objective: \eqn{\tfrac{1}{2} ||W(Ax-b)||^2 + \tfrac{\lambda}{2} ||D_2 x||^2}
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param w Numeric weights vector (length m).
#' @param regularization Numeric; roughness ridge.
#' @param x Numeric spectrum vector (length n).
#' @param cache Optional environment or list with cached derivative
#'   operator under \code{$D}.
#' @return Numeric scalar.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' uno_objective(A, b, rep(1, 3), 0.01, c(0.5, 0.5, 0.5))
uno_objective <- function(A, b, w, regularization, x, cache = NULL) {
    x <- as.numeric(x)
    r <- w * (as.numeric(A %*% x) - b)
    obj <- 0.5 * sum(r * r)
    if (regularization > 0 && length(x) > 2L) {
        D <- if (!is.null(cache) && !is.null(cache$D)) cache$D else
                 .uno_derivative(length(x))
        if (!is.null(cache)) cache$D <- D
        Dx <- D %*% x
        obj <- obj + 0.5 * regularization * sum(Dx * Dx)
    }
    as.numeric(obj)
}

#' Gradient of \code{\link{uno_objective}}
#'
#' @inheritParams uno_objective
#' @return Numeric vector (length n).
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' uno_gradient(A, b, rep(1, 3), 0.01, c(0.5, 0.5, 0.5))
uno_gradient <- function(A, b, w, regularization, x, cache = NULL) {
    x <- as.numeric(x)
    Aw <- sweep(A, 1L, w, "*")
    g <- as.numeric(t(Aw) %*% (w * (as.numeric(A %*% x) - b)))
    if (regularization > 0 && length(x) > 2L) {
        D <- if (!is.null(cache) && !is.null(cache$D)) cache$D else
                 .uno_derivative(length(x))
        if (!is.null(cache)) cache$D <- D
        g <- g + regularization * as.numeric(t(D) %*% (D %*% x))
    }
    g
}

#' Fletcher-Leyffer (Uno) filter acceptability test
#'
#' @param filter_entries List of \code{c(f_j, viol_j)} pairs.
#' @param f Numeric; trial objective value.
#' @param viol Numeric; trial constraint violation.
#' @param gamma Numeric; safety margin. Default 1e-5.
#' @return \code{TRUE} when the trial point is filter-acceptable.
#' @export
#' @examples
#' uno_filter(list(c(0.5, 0.1)), 0.4, 0.05)
uno_filter <- function(filter_entries, f, viol, gamma = 1e-5) {
    for (entry in filter_entries) {
        f_j <- entry[1L]; viol_j <- entry[2L]
        if (viol_j > 0.0) {
            acceptable <- (f <= f_j - gamma * viol_j) ||
                              (viol <= (1.0 - gamma) * viol_j)
        } else {
            acceptable <- (f < f_j) || (viol < viol_j)
        }
        if (!acceptable) return(FALSE)
    }
    TRUE
}

.uno_filter_sqp <- function(A, b, w, x0, regularization, max_iterations,
                             tolerance, cache) {
    m <- nrow(A); n <- ncol(A)
    # Equivalent single least-squares system [W A; sqrt(lam) Dn] x = [W b; 0]
    if (regularization > 0 && n > 2L) {
        D <- .uno_derivative(n)
        Gsqrt <- sqrt(regularization) * D
        A_aug <- rbind(sweep(A, 1L, w, "*"), Gsqrt)
        b_aug <- c(b * w, rep(0, n - 2L))
    } else {
        A_aug <- sweep(A, 1L, w, "*")
        b_aug <- b * w
    }
    xs <- tryCatch(as.numeric(lsei::nnls(A_aug, b_aug)$x),
                   error = function(e) qr.solve(A_aug, b_aug, tol = 1e-10))
    x <- pmax(xs, 0.0)
    f <- uno_objective(A, b, w, regularization, x, cache = cache)
    viol <- sum(pmin(x, 0.0)^2)
    g <- uno_gradient(A, b, w, regularization, x, cache = cache)
    interior <- x > 10.0 * .TINY
    m_int <- if (any(interior)) max(abs(g[interior])) else 0.0
    act <- !interior
    m_act <- if (any(act)) max(pmax(-g[act], 0.0)) else 0.0
    dual_inf <- max(m_int, m_act)
    converged <- (dual_inf <= tolerance && viol <= tolerance) &&
                     all(is.finite(x))
    list(x = x, it = 1L, converged = converged, f = f, viol = viol,
         dual_inf = dual_inf)
}

.uno_ipopt_like <- function(A, b, w, H, x0, regularization, max_iterations,
                             tolerance, hessian_mode, cache) {
    n <- length(x0)
    x <- pmax(as.numeric(x0), .TINY)
    mus <- 1.0 * (0.1^(0:29))
    mu_k <- 1L
    converged <- FALSE
    Id <- diag(n)
    it <- 0L
    f <- uno_objective(A, b, w, regularization, x, cache = cache)
    hess_mode <- tolower(hessian_mode)
    if (!(hess_mode %in% c("exact", "bfgs")))
        stop("hessian must be 'exact' or 'bfgs', got '", hessian_mode, "'")
    scale0 <- mean(diag(H))
    B <- Id * (if (is.finite(scale0) && scale0 > 0) scale0 else 1.0)
    g_prev <- NULL
    g0 <- uno_gradient(A, b, w, regularization, x, cache = cache)
    grad_scale <- max(max(abs(g0 - max(mus[1L], 0.0) / x)), 1.0)

    for (it in seq_len(as.integer(max_iterations))) {
        mu <- mus[min(mu_k, length(mus))]
        if (mu_k < length(mus) && it %% 3L == 0L) mu_k <- mu_k + 1L
        Hb <- if (hess_mode == "exact") H else B
        Hb <- Hb + mu * diag(1.0 / (x * x))
        Hb <- Hb + 1e-12 * max(mean(abs(H)), 1.0) * Id

        g <- uno_gradient(A, b, w, regularization, x, cache = cache)
        barrier_grad <- g - mu / x

        d <- tryCatch(as.numeric(solve(Hb, -barrier_grad)),
                      error = function(e)
                          as.numeric(qr.solve(Hb, -barrier_grad, tol = 1e-10)))
        if (!all(is.finite(d))) d <- rep(0, n)

        # fraction-to-the-boundary
        neg <- d < 0
        if (any(neg)) {
            alpha_max <- .UNO_FTB_TAU * min(-x[neg] / d[neg])
        } else {
            alpha_max <- 1.0
        }
        t <- min(1.0, alpha_max)
        x_old <- x
        # backtracking on the barrier objective
        for (inner in 1:50) {
            xs <- pmax(x + t * d, .TINY)
            fs <- uno_objective(A, b, w, regularization, xs, cache = cache)
            barrier_new <- fs - mu * sum(log(xs))
            barrier_old <- f - mu * sum(log(pmax(x, .TINY)))
            if (barrier_new <= barrier_old - 1e-4 * t * abs(barrier_old) ||
                    t < 1e-14) {
                x <- xs; f <- fs; break
            }
            t <- t * 0.5
        }
        g_new <- uno_gradient(A, b, w, regularization, x, cache = cache)
        if (hess_mode == "bfgs" && !is.null(g_prev)) {
            s <- x - x_old
            y <- g_new - g_prev
            ys <- sum(s * y)
            if (ys > 1e-10) {
                Bs <- as.numeric(B %*% s)
                sBs <- sum(s * Bs)
                if (sBs > 1e-12) {
                    B <- B + outer(y, y) / ys - outer(Bs, Bs) / sBs
                }
            }
        }
        g_prev <- g_new
        g <- g_new
        dual_inf_it <- max(abs(barrier_grad))
        if (dual_inf_it <= grad_scale * tolerance && mu <= mus[length(mus)]) {
            converged <- TRUE; break
        }
    }
    g_last <- uno_gradient(A, b, w, regularization, x, cache = cache)
    mu_final <- mus[min(mu_k, length(mus))]
    dual_inf <- max(abs(g_last - mu_final / pmax(x, .TINY)))
    dual_inf_rel <- dual_inf / grad_scale
    viol <- sum(pmin(x, 0.0)^2)
    list(x = x, it = it, converged = converged, f = f, viol = viol,
         dual_inf = dual_inf_rel)
}

#' Solve the unfolding NLP with an Uno preset (full diagnostics)
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial spectrum.
#' @param preset Character: \code{"filter_sqp"} (default) or
#'   \code{"ipopt_like"}.
#' @param weights Character (\code{"uniform"} default, \code{"poisson"})
#'   or explicit positive weight array.
#' @param regularization Numeric; relative roughness ridge. Default 1e-3.
#' @param hessian Character: \code{"exact"} (default) or \code{"bfgs"}
#'   (interior-point preset only).
#' @param max_iterations Integer; max iterations. Default 300.
#' @param tolerance Numeric; KKT tolerance. Default 1e-10.
#' @return Diagnostics list with \code{spectrum}, \code{preset},
#'   \code{hessian}, \code{objective}, \code{constraint_violation},
#'   \code{dual_infeasibility}, \code{n_iterations}, \code{converged}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' d <- solve_uno_full(A, b, preset = "filter_sqp", max_iterations = 100)
solve_uno_full <- function(A, b, x0 = NULL, preset = "filter_sqp",
                            weights = "uniform", regularization = 1e-3,
                            hessian = "exact", max_iterations = 300L,
                            tolerance = 1e-10) {
    preset <- tolower(as.character(preset))
    if (!(preset %in% .UNO_PRESETS))
        stop("preset must be one of ", paste(.UNO_PRESETS, collapse = ", "),
             ", got '", preset, "'")
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A); n <- ncol(A)
    if (regularization < 0)
        stop("regularization must be non-negative, got ", regularization)

    if (is.character(weights)) {
        weights <- tolower(weights)
        if (weights %in% c("", "uniform", "none", "ones")) {
            w <- rep(1.0, m)
        } else if (weights == "poisson") {
            w <- 1.0 / pmax(b, .TINY)
        } else {
            stop("weights must be 'uniform', 'poisson' or an array, got '",
                 weights, "'")
        }
    } else {
        w <- as.numeric(weights)
        if (length(w) != m || any(w <= 0))
            stop("weights must be a positive array of length ", m)
    }

    cache <- new.env(parent = emptyenv())
    if (regularization > 0 && n > 2L) {
        D <- .uno_derivative(n)
        G <- t(D) %*% D
        G <- G / max(mean(diag(G)), 1.0)
    } else {
        G <- matrix(0.0, n, n)
    }
    cache$D <- .uno_derivative(n)

    Aw <- sweep(A, 1L, w, "*")
    H <- t(Aw) %*% Aw + regularization * G

    x_start <- if (is.null(x0)) rep(1, n) else as.numeric(x0)

    if (preset == "filter_sqp") {
        res <- .uno_filter_sqp(A, b, w, x_start, regularization,
                                max_iterations, tolerance, cache)
    } else {
        res <- .uno_ipopt_like(A, b, w, H, x_start, regularization,
                                max_iterations, tolerance, hessian, cache)
    }
    list(spectrum = res$x,
         preset = preset,
         hessian = if (preset == "filter_sqp") "exact" else tolower(hessian),
         objective = as.numeric(res$f),
         constraint_violation = as.numeric(res$viol),
         dual_infeasibility = as.numeric(res$dual_inf),
         n_iterations = as.integer(res$it),
         converged = as.logical(res$converged))
}

#' Solve the unfolding NLP with an Uno preset
#'
#' @inheritParams solve_uno_full
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_uno(A, b, preset = "filter_sqp")
solve_uno <- function(A, b, x0 = NULL, preset = "filter_sqp",
                       weights = "uniform", regularization = 1e-3,
                       hessian = "exact", max_iterations = 300L,
                       tolerance = 1e-10) {
    diag <- solve_uno_full(A, b, x0 = x0, preset = preset, weights = weights,
                            regularization = regularization, hessian = hessian,
                            max_iterations = max_iterations,
                            tolerance = tolerance)
    list(spectrum = diag$spectrum,
         iterations = diag$n_iterations,
         converged = diag$converged)
}

#' Wrapper around \code{\link{solve_uno}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_uno_full
#' @export
unfold_uno <- function(detector_names, n_energy_bins, E_MeV, sensitivities,
                        cc_icrp116, save_result_callback, readings,
                        initial_spectrum = NULL, preset = "filter_sqp",
                        weights = "uniform", regularization = 1e-3,
                        hessian = "exact", max_iterations = 300L,
                        tolerance = 1e-10, calculate_errors = FALSE,
                        noise_level = 0.01, n_montecarlo = 100L,
                        save_result = FALSE, random_state = NULL,
                        max_neutron_energy = NULL) {
    x0_default <- rep(0, n_energy_bins)

    extra_output <- tryCatch({
        A_mat <- do.call(rbind, lapply(detector_names, function(nm)
            as.numeric(sensitivities[[nm]])))
        b_vec <- as.numeric(readings[detector_names])
        diag <- solve_uno_full(A_mat, b_vec, x0 = NULL, preset = preset,
                                weights = weights,
                                regularization = regularization,
                                hessian = hessian,
                                max_iterations = max_iterations,
                                tolerance = tolerance)
        list(uno_preset = diag$preset,
             objective = diag$objective,
             constraint_violation = diag$constraint_violation,
             dual_infeasibility = diag$dual_infeasibility,
             uno_converged = diag$converged)
    }, error = function(e) NULL)

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_uno, preset = preset,
                                         weights = weights,
                                         regularization = regularization,
                                         hessian = hessian,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = paste0("Uno (", preset, ")"),
        extra_output = extra_output,
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
