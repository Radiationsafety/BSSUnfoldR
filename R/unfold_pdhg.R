#' Primal-Dual Hybrid Gradient (Chambolle-Pock) unfolding with TV
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_pdhg.py}
#' (ODL-independent, pure R). Solves
#' \eqn{\min_x \tfrac{1}{2}\|Ax-b\|^2 + w_\mathrm{tv}\,\mathrm{TV}(x)}
#' (+ non-negativity indicator) with Chambolle & Pock's primal-dual hybrid
#' gradient on the one-dimensional energy axis.
#'
#' @name pdhg-methods
NULL

.tv_operator <- function(n) {
    D <- matrix(0.0, nrow = n, ncol = n)
    for (i in seq_len(max(n - 1L, 0L))) {
        D[i, i] <- -1
        D[i, i + 1L] <- 1
    }
    D
}

.shrink <- function(z, s) sign(z) * pmax(0, abs(z) - s)

#' Solve by PDHG (L2 + TV)
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param tau Numeric; primal step. Default \code{NULL} = auto.
#' @param sigma Numeric; dual step. Default \code{NULL} = auto (tau chosen
#'   so that \code{tau*sigma*||K||^2 < 1}).
#' @param use_tv Logical; include the TV dual block. Default TRUE.
#' @param tv_weight Numeric TV weight. Default 0.05.
#' @param nonnegativity Logical; clamp negative bins each primal update.
#'   Default TRUE.
#' @param tolerance Numeric relative tolerance. Default 1e-6.
#' @param max_iterations Integer iterations. Default 500.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_pdhg(A, b, rep(1, 3), max_iterations = 100L)
solve_pdhg <- function(A, b, x0 = NULL, tau = NULL, sigma = NULL,
                       use_tv = TRUE, tv_weight = 0.05, nonnegativity = TRUE,
                       tolerance = 1e-6, max_iterations = 500L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (is.null(x0)) x0 <- rep(mean(b) / max(sum(A), 1e-10) * n, n)
    use_tv <- isTRUE(use_tv) && tv_weight > 0
    D <- if (use_tv) .tv_operator(n) else matrix(0.0, nrow = 0L, ncol = n)
    # operator K = [A; D]; PDHG solves min 1/2||Ax-b||^2 + w TV(x)
    smax <- tryCatch(norm(A, "2"),
                     error = function(e)
                         sqrt(max(eigen(crossprod(A), symmetric = TRUE,
                                        only.values = TRUE)$values)))
    if (is.null(tau) || is.null(sigma)) {
        l2norm_D <- if (use_tv) sqrt(4) else 1
        normK <- sqrt(max(smax^2 + l2norm_D^2, 1e-30))
        if (is.null(tau)) tau <- 1.0 / normK      # primal step
        if (is.null(sigma)) sigma <- 1.0 / normK  # dual step; tau*sigma*||K||^2 < 1
    }
    p <- if (use_tv) rep(0, n) else numeric(0)
    x <- pmax(as.numeric(x0), 0)
    iterations <- 0L
    converged <- FALSE
    bnorm <- max(sqrt(sum(b^2)), 1e-30)
    creator <- NULL
    for (it in seq_len(max_iterations)) {
        iterations <- it
        # dual ascent on TV block: p^{k+1} = prox_{w/sigma}(p + sigma D x)
        if (use_tv) {
            p <- .shrink(p + sigma * as.numeric(D %*% x),
                         tv_weight / sigma)
        }
        # primal descent: x^{k+1} = prox_{tau*1/2||A.-b||^2}(x - tau D' p)
        xd <- if (use_tv) x - tau * as.numeric(t(D) %*% p) else x
        # prox of the quadratic: (I + tau A'A) z = x + tau A'b
        if (it == 1L) {
            lhs <- crossprod(A) * tau + diag(n)
            creator <- function(yy)
                as.numeric(solve(lhs, yy, tol = 1e-10))
        }
        rhs <- xd + tau * as.numeric(t(A) %*% b)
        x_new <- creator(rhs)
        if (isTRUE(nonnegativity)) x_new <- pmax(x_new, 0)
        x <- x_new
        resid <- sqrt(sum((as.numeric(A %*% x) - b)^2)) / bnorm
        if (resid < tolerance) { converged <- TRUE; break }
    }
    spectrum <- if (isTRUE(nonnegativity)) pmax(x, 0) else x
    list(spectrum = spectrum, iterations = as.integer(iterations),
         converged = converged)
}

#' Wrapper around \code{\link{solve_pdhg}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_pdhg
#' @export
unfold_pdhg <- function(detector_names, n_energy_bins, E_MeV, sensitivities,
                        cc_icrp116, save_result_callback, readings,
                        initial_spectrum = NULL, tau = NULL, sigma = NULL,
                        use_tv = TRUE, tv_weight = 0.05,
                        nonnegativity = TRUE, tolerance = 1e-6,
                        max_iterations = 500L, method_name = "PDHG",
                        calculate_errors = FALSE, noise_level = 0.01,
                        n_montecarlo = 100L, save_result = FALSE,
                        random_state = NULL,
                              max_neutron_energy = NULL) {
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1, n_energy_bins),
        solve_func = solve_pdhg,
        solve_kwargs = list(tau = tau, sigma = sigma, use_tv = use_tv,
                            tv_weight = tv_weight,
                            nonnegativity = nonnegativity,
                            tolerance = tolerance,
                            max_iterations = max_iterations),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}

#' Douglas-Rachford splitting for composite L2 + TV (+ non-negativity)
#'
#' R port of \code{bssunfold/core/unfold_douglas_rachford.py}: splits
#' \eqn{f = \tfrac12\|Ax-b\|^2} and \eqn{g = w\mathrm{TV} + \iota_{\{x\ge0\}}}
#' with the classical Douglas-Rachford iteration
#' \eqn{u \leftarrow (1-\gamma)u + \gamma\;\mathrm{prox}_f(
#' 2\,\mathrm{prox}_g(u) - u)}.
#'
#' @inheritParams solve_pdhg
#' @param relaxation Numeric relaxation \eqn{\gamma \in (0,2)}, default 1.0.
#' @param tv_weight Numeric weight of the TV term (default 0.05; use 0 for
#'   a pure non-negative L2 solve).
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_douglas_rachford(A, b, rep(1, 3), max_iterations = 200L)
solve_douglas_rachford <- function(A, b, x0 = NULL, use_tv = TRUE,
                                   tv_weight = 0.05, relaxation = 1.0,
                                   tolerance = 1e-6,
                                   max_iterations = 500L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (is.null(x0)) x0 <- rep(mean(b) / max(sum(A), 1e-10) * n, n)
    gamma <- pmin(pmax(as.numeric(relaxation), 0.01), 1.99)
    # prox of the quadratic f with unit strength:
    # prox_f(u) = (I + A'A)^{-1} (u + A'b)
    lhs <- crossprod(A) + diag(n)
    lstsq <- as.numeric(solve(lhs, as.numeric(t(A) %*% b), tol = 1e-10))
    D <- if (isTRUE(use_tv) && tv_weight > 0) .tv_operator(n) else
         matrix(0.0, nrow = 0L, ncol = n)
    # prox of TV + indicator g (TV prox by dual shrinkage, then projection)
    .prox_g <- function(u) {
        if (nrow(D) > 0) {
            dv <- as.numeric(D %*% u)
            sv <- as.numeric(.shrink(dv, tv_weight))
            v <- u - as.numeric(t(D) %*% (dv - sv))
        } else {
            v <- u
        }
        pmax(v, 0)
    }
    u <- as.numeric(x0)
    bnorm <- max(sqrt(sum(b^2)), 1e-30)
    iterations <- 0L
    converged <- FALSE
    .prox_f <- function(u) as.numeric(solve(lhs, u + as.numeric(t(A) %*% b)))
    for (it in seq_len(max_iterations)) {
        iterations <- it
        pf <- .prox_f(u)
        y2 <- .prox_g(2 * pf - u)
        u <- u + gamma * (y2 - pf)
        resid <- sqrt(sum((as.numeric(A %*% pmax(pf, 0)) - b)^2)) / bnorm
        if (resid < tolerance) { converged <- TRUE; break }
    }
    list(spectrum = pmax(pf, 0), iterations = as.integer(iterations),
         converged = converged)
}

#' Wrapper around \code{\link{solve_douglas_rachford}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_douglas_rachford
#' @export
unfold_douglas_rachford <- function(detector_names, n_energy_bins, E_MeV,
                                       sensitivities, cc_icrp116,
                                       save_result_callback, readings,
                                       initial_spectrum = NULL, use_tv = TRUE,
                                       tv_weight = 0.05, relaxation = 1.0,
                                       tolerance = 1e-6,
                                       max_iterations = 500L,
                                       method_name = "Douglas-Rachford",
                                       calculate_errors = FALSE,
                                       noise_level = 0.01,
                                       n_montecarlo = 100L,
                                       save_result = FALSE,
                                       random_state = NULL,
                                       max_neutron_energy = NULL) {
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1, n_energy_bins),
        solve_func = solve_douglas_rachford,
        solve_kwargs = list(use_tv = use_tv, tv_weight = tv_weight,
                            relaxation = relaxation, tolerance = tolerance,
                            max_iterations = max_iterations),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
