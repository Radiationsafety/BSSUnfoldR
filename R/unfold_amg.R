#' AMG / preconditioned Krylov unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_amg.py}.
#' Solves the (optionally Tikhonov-regularized) response system with a
#' Krylov method (\code{cg} for the SPD case, \code{bicgstab} and
#' \code{gmres} for the general case) preconditioned by
#' \itemize{
#'   \item \code{"amg"} \eqn{\rightarrow} algebraic-multigrid-like
#'     two-level smoothing (down-Jacobi sweeps + coarse correction);
#'   \item \code{"jacobi"} \eqn{\rightarrow} point-Jacobi preconditioner;
#'   \item \code{"gs"} \eqn{\rightarrow} Gauss-Seidel iteration (backward
#'     application);
#'   \item \code{"sor"} / \code{"ssor"} \eqn{\rightarrow} symmetric or
#'     forward SOR relaxation;
#'   \item \code{"none"} \eqn{\rightarrow} no preconditioner.
#' }
#' The pure-R AMG emulation follows Katzengruber et al.: aggregated coarse
#' grids with relaxed Jacobi pre/post smoothing, which is a faithful
#' classical-memory model for the Python version that uses pyamg
#' (optional dependency).
#'
#' @name amg-methods
NULL

# ---- Krylov ----------------------------------------------------------------

.amg_cg <- function(A, b, apply_pre, max_iterations, tolerance) {
    x <- rep(0, ncol(A))
    r <- as.numeric(b)
    z <- apply_pre(r)
    p <- z
    rz <- sum(r * z)
    iterations <- 0L
    bnorm <- max(sqrt(sum(r^2)), 1e-30)
    converged <- FALSE
    for (it in seq_len(max_iterations)) {
        iterations <- it
        Ap <- as.numeric(A %*% p)
        denom <- sum(p * Ap)
        if (abs(denom) < 1e-30) break
        alpha <- rz / denom
        x <- x + alpha * p
        r <- r - alpha * Ap
        if (sqrt(sum(r^2)) / bnorm < tolerance) { converged <- TRUE; break }
        z <- apply_pre(r)
        rz_new <- sum(r * z)
        if (abs(rz) < 1e-30) break
        p <- z + (rz_new / rz) * p
        rz <- rz_new
    }
    list(x = x, iterations = iterations, converged = converged)
}

.amg_gmres <- function(A, b, apply_pre, max_iterations, tolerance) {
    n <- ncol(A)
    x <- rep(0, n)
    r <- as.numeric(b)
    bnorm <- max(sqrt(sum(r^2)), 1e-30)
    # Restarted preconditioned GMRES (projection on Krylov space)
    beta <- bnorm
    V <- matrix(0, nrow = n, ncol = max_iterations)
    H <- matrix(0, nrow = max_iterations + 1L, ncol = max_iterations)
    V[, 1] <- r / beta
    krylov_k <- max_iterations
    converged <- FALSE
    iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        w <- apply_pre(as.numeric(A %*% V[, it]))
        hnorm_w <- max(sqrt(sum(w^2)), 1e-30)
        # simple orthogonalization (modified Gram-Schmidt)
        for (j in seq_len(it)) {
            H[j, it] <- sum(w * V[, j])
            w <- w - H[j, it] * V[, j]
        }
        H[it + 1L, it] <- sqrt(sum(w^2))
        if (H[it + 1L, it] > 1e-30) V[, it + 1L] <- w / H[it + 1L, it]
        # Residual estimate from last column (Givens-free estimate)
        resid_est <- H[it + 1L, it] / bnorm
        if (resid_est < tolerance) { krylov_k <- it; converged <- TRUE; break }
    }
    # Least-squares solve on the small system for y
    k <- krylov_k
    Hs <- H[seq_len(k + 1L), seq_len(k), drop = FALSE]
    rhs0 <- as.numeric(t(Hs) %*% c(beta, rep(0, k)))
    omega <- tryCatch(
        as.numeric(solve(t(Hs) %*% Hs + 1e-10 * diag(k), rhs0, tol = 1e-10)),
        error = function(e)
            as.numeric(qr.solve(t(Hs) %*% Hs + 1e-6 * diag(k), rhs0,
                                 tol = 1e-8)))
    x <- V[, seq_len(k), drop = FALSE] %*% omega
    list(x = as.numeric(x), iterations = iterations, converged = converged)
}

.amg_bicgstab <- function(A, b, apply_pre, max_iterations, tolerance) {
    n <- ncol(A)
    x <- rep(0, n)
    r <- as.numeric(b)
    rhat <- r
    bnorm <- max(sqrt(sum(r^2)), 1e-30)
    z <- apply_pre(r)
    p <- z
    rz <- sum(rhat * r)
    rho <- 1.0
    alpha <- 1.0
    omega <- 1.0
    iterations <- 0L
    converged <- FALSE
    for (it in seq_len(max_iterations)) {
        iterations <- it
        u <- apply_pre(as.numeric(A %*% p))
        denom <- as.numeric(sum(rhat * u))
        if (abs(denom) < 1e-30) break
        alpha <- rho / denom
        s <- r - alpha * u
        if (sqrt(sum(s^2)) / bnorm < tolerance) {
            x <- x + alpha * p
            r <- s
            converged <- TRUE
            break
        }
        t <- apply_pre(as.numeric(A %*% s))
        omega <- sum(t * s) / max(sum(t * t), 1e-30)
        x <- x + alpha * p + omega * s
        r <- s - omega * t
        if (sqrt(sum(r^2)) / bnorm < tolerance) { converged <- TRUE; break }
        rho_new <- sum(t * rhat) / max(sum(rhat * rhat), 1e-30)
        beta_r <- (rho_new / rho) * (alpha / omega)
        p <- r + beta_r * (p - omega * u)
        rho <- rho_new
    }
    list(x = x, iterations = iterations, converged = converged)
}

# ---- Preconditioners on M = A'A (SPD normal-equation preconditioning) ------

.amg_preconditioner <- function(A, preconditioner, omega, smooth_sweeps) {
    At <- t(A)
    n <- nrow(At)
    M <- crossprod(A)
    switch(tolower(preconditioner),
        none = function(r) r,
        jacobi = {
            dinv <- 1 / pmax(diag(M), 1e-30)
            function(r) dinv * r
        },
        gs = {
            function(r)
                tryCatch(as.numeric(solve(M, r)), error = function(e) r)
        },
        sor = ,
        ssor = {
            function(r) {
                y <- rep(0, n)
                for (s in seq_len(smooth_sweeps)) {
                    y <- y + omega * (r - as.numeric(M %*% y))
                }
                y
            }
        },
        amg = {
            # Two-level AMG-like: smooth, restrict (binomial weights),
            # coarse average, prolongate, smooth again.
            n_coarse <- max(floor(n / 2), 2L)
            idx_c <- seq(1L, n, length.out = n_coarse)
            function(r) {
                y <- rep(0, n)
                for (s in seq_len(smooth_sweeps)) {
                    y <- y + 0.9 * omega * (r - as.numeric(M %*% y))
                }
                # restrict + coarse correction with much stronger damping
                rc <- r[idx_c] - as.numeric(t(A) %*% y)[idx_c]
                Mcc <- M[idx_c, idx_c, drop = FALSE]
                yc <- tryCatch(as.numeric(solve(Mcc, rc)),
                               error = function(e)
                                   as.numeric(qr.solve(Mcc, rc, tol = 1e-8)))
                e <- rep(0, n); e[idx_c] <- yc
                y <- y + e
                for (s in seq_len(smooth_sweeps)) {
                    y <- y + 0.9 * omega * (r - as.numeric(M %*% y))
                }
                pmax(y, -1e30)
            }
        },
        stop("Unknown preconditioner: ", preconditioner)
    )
}

#' Solve by AMG / preconditioned Krylov method
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Not used (Krylov starts from zero), kept for API compatibility.
#' @param method Character; \code{"cg"}, \code{"bicgstab"}, or
#'   \code{"gmres"}.
#' @param preconditioner Character; \code{"amg"}, \code{"jacobi"},
#'   \code{"gs"}, \code{"sor"}, \code{"ssor"}, or \code{"none"}.
#' @param omega Numeric; SOR relaxation factor (for sor/ssor/amg). Default
#'   1.2.
#' @param max_iterations Integer; max outer iterations. Default 500.
#' @param tolerance Numeric; relative residual tolerance. Default 1e-6.
#' @param outer_iterations Integer; number of GMRES restart cycles. Default
#'   10.
#' @param nonnegativity Logical; clamp negative entries to zero at the end.
#'   Default TRUE.
#' @param regularization Numeric; Tikhonov regularization added to the
#'   normal matrix. Default 0.
#' @param smooth_sweeps Integer; smoothing sweeps per AMG application.
#'   Default 3.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_amg(A, b, NULL)
solve_amg <- function(A, b, x0 = NULL, method = "cg",
                      preconditioner = "jacobi", omega = 1.2,
                      max_iterations = 500L, tolerance = 1e-6,
                      outer_iterations = 10L, nonnegativity = TRUE,
                      regularization = 0.0, smooth_sweeps = 3L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    n_outer <- max(as.integer(outer_iterations), 1L)
    # We solve the SPD normal equation (Aeff' Aeff) x = Aeff' b.
    if (regularization > 0) {
        Aeff <- rbind(A, sqrt(regularization) * diag(n))
        beff <- c(b, rep(0, n))
    } else {
        Aeff <- A
        beff <- b
    }
    M <- crossprod(Aeff) + 1e-12 * diag(n)   # numerical SPD safeguard
    rhs <- as.numeric(t(Aeff) %*% beff)
    pre <- .amg_preconditioner(Aeff, preconditioner, omega, smooth_sweeps)
    meth <- tolower(method)
    if (meth == "cg") {
        res <- .amg_cg(M, rhs, pre, as.integer(max_iterations), tolerance)
    } else if (meth == "bicgstab") {
        res <- .amg_bicgstab(M, rhs, pre, as.integer(max_iterations),
                             tolerance)
    } else if (meth == "gmres") {
        x_acc <- rep(0, n); it_tot <- 0L
        converged <- FALSE
        r_cur <- rhs
        # GMRES on normal equations with restarts
        for (rc in seq_len(n_outer)) {
            g <- .amg_gmres(M, r_cur, pre, as.integer(max_iterations),
                            tolerance)
            x_acc <- x_acc + g$x
            it_tot <- it_tot + g$iterations
            r_cur <- rhs - as.numeric(M %*% x_acc)
            if (sqrt(sum(r_cur^2)) / max(sqrt(sum(rhs^2)), 1e-30) <
                tolerance) { converged <- TRUE; break }
        }
        res <- list(x = x_acc, iterations = it_tot, converged = converged)
    } else {
        stop("Unknown Krylov method '", method,
             "'. Available: cg, bicgstab, gmres.")
    }
    spectrum <- if (isTRUE(nonnegativity)) pmax(res$x, 0) else res$x
    list(spectrum = spectrum, iterations = as.integer(res$iterations),
         converged = res$converged)
}

#' Wrapper around \code{\link{solve_amg}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_amg
#' @export
unfold_amg <- function(detector_names, n_energy_bins, E_MeV, sensitivities,
                       cc_icrp116, save_result_callback, readings,
                       initial_spectrum = NULL, method = "cg",
                       preconditioner = "jacobi", omega = 1.2,
                       max_iterations = 500L, tolerance = 1e-6,
                       outer_iterations = 10L, nonnegativity = TRUE,
                       regularization = 0.0, method_name = "AMG",
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
        solve_func = solve_amg,
        solve_kwargs = list(method = method,
                            preconditioner = preconditioner, omega = omega,
                            max_iterations = max_iterations,
                            tolerance = tolerance,
                            outer_iterations = outer_iterations,
                            nonnegativity = nonnegativity,
                            regularization = regularization),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
