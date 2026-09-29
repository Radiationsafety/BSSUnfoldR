#' AMG / preconditioned Krylov unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_amg.py}.
#' The (Tikhonov-damped) normal equations
#' \deqn{(A^\top A + \lambda I) x = A^\top b}
#' are solved with a preconditioned Krylov method (\code{cg} for the SPD
#' case, \code{bicgstab} and \code{gmres} for the general case) using a
#' preconditioner that approximates \eqn{(A^\top A + \lambda I)^{-1}}:
#' \itemize{
#'   \item \code{"amg"} \eqn{\rightarrow} algebraic multigrid
#'     (smoothed aggregation).  \code{pyamg} is an optional Python
#'     dependency and is not installed in the reference environment, so
#'     \code{bssunfold} emits a \code{RuntimeWarning} and degrades to
#'     \code{"jacobi"}; the port reproduces that fallback exactly.
#'   \item \code{"jacobi"} \eqn{\rightarrow} point-Jacobi
#'     (\eqn{M = D}{M = D}).
#'   \item \code{"gs"} \eqn{\rightarrow} one forward Gauss-Seidel sweep
#'     (\eqn{M = D + L}{M = D + L}).
#'   \item \code{"sor"} \eqn{\rightarrow} one SOR sweep
#'     (\eqn{M = (D + \omega L)/\omega}{M = (D + omega L)/omega}).
#'   \item \code{"ssor"} \eqn{\rightarrow} symmetric SOR sweep
#'     (\eqn{M = (D + \omega L) D^{-1} (D + \omega U) / (\omega (2 - \omega))}).
#'   \item \code{"none"} \eqn{\rightarrow} identity.
#' }
#' Non-negativity is enforced by projected outer restarts: after every
#' Krylov solve the spectrum is clamped to \code{x >= 0} and the
#' iteration is restarted on the residual of the clamped iterate,
#' \code{outer_iterations} times.
#'
#' @name amg-methods
NULL

# Valid argument values, mirroring _VALID_METHODS / _VALID_PRECONDITIONERS
# and _CG_COMPATIBLE of the Python module.
.amg_valid_methods <- c("cg", "bicgstab", "gmres")
.amg_valid_preconditioners <- c("amg", "jacobi", "gs", "sor", "ssor", "none")
## Preconditioners whose application is symmetric positive definite, i.e.
## usable with CG. "gs" and "sor" are non-symmetric.
.amg_cg_compatible <- c("amg", "jacobi", "ssor", "none")
## Auto Tikhonov damping factor (_AUTO_REG_FACTOR): damping is
## 1e-4 * mean(diag(A'A)) when regularization is NULL.
.amg_auto_reg_factor <- 1e-4

# ---- Preconditioner application (Python _stationary_apply) -----------------

.amg_stationary_apply <- function(N, kind, omega, r) {
    D <- diag(N)
    D_safe <- ifelse(abs(D) > 0, D, 1)
    if (identical(kind, "jacobi")) {
        return(r / D_safe)
    }
    lower <- matrix(0, nrow = nrow(N), ncol = ncol(N))
    lower[lower.tri(N)] <- N[lower.tri(N)]
    upper <- matrix(0, nrow = nrow(N), ncol = ncol(N))
    upper[upper.tri(N)] <- N[upper.tri(N)]
    if (identical(kind, "gs")) {
        Mf <- lower + diag(D_safe, nrow = nrow(N))
        return(as.numeric(backsolve(t(Mf), r, transpose = TRUE)))
    }
    if (identical(kind, "sor")) {
        Mf <- omega * lower + diag(D_safe, nrow = nrow(N))
        return(omega * as.numeric(backsolve(t(Mf), r, transpose = TRUE)))
    }
    if (identical(kind, "ssor")) {
        Mf <- omega * lower + diag(D_safe, nrow = nrow(N))
        Mb <- omega * upper + diag(D_safe, nrow = nrow(N))
        t <- as.numeric(backsolve(t(Mf), r, transpose = TRUE))
        t <- D_safe * t
        t <- as.numeric(backsolve(Mb, t))
        return(omega * (2 - omega) * t)
    }
    stop("Unknown stationary kind: ", kind)
}

## build_preconditioner(): a function approximating (A'A + damping I)^{-1} r.
## "amg" degrades to "jacobi", like the Python module without pyamg.
## The returned function counts its applications, mirroring
## _counting_operator(); the counter lives in a shared environment.
.amg_build_preconditioner <- function(Nmat, kind, omega, counter) {
    n <- nrow(Nmat)
    if (identical(kind, "none")) {
        return(function(r) {
            counter$n <- counter$n + 1L
            as.numeric(r)
        })
    }
    if (identical(kind, "amg")) {
        warning("pyamg is not installed -- AMG preconditioner falls back ",
                "to Jacobi.", call. = FALSE)
        kind <- "jacobi"
    }
    function(r) {
        counter$n <- counter$n + 1L
        .amg_stationary_apply(Nmat, kind, omega, as.numeric(r))
    }
}

# ---- Krylov solvers --------------------------------------------------------

## scipy.sparse.linalg.cg mirror.  Stopping test is
## norm(r) < max(atol, rtol * norm(b)) evaluated at the TOP of the loop,
## x0 is always the zero vector here, and info is 0 on success or
## maxiter when the loop is exhausted.
.amg_cg <- function(matvec, psolve, b, rtol, maxiter, atol = 0) {
    x <- numeric(length(b))
    r <- as.numeric(b)                       # x is all zero -> r = b
    bnrm2 <- sqrt(sum(r^2))
    if (bnrm2 == 0) return(list(x = r, info = 0L))
    atol <- max(atol, rtol * bnrm2)
    rho_prev <- NA_real_
    p <- NULL
    for (iteration in seq_len(maxiter)) {
        if (sqrt(sum(r^2)) < atol) return(list(x = x, info = 0L))
        z <- psolve(r)
        rho_cur <- sum(r * z)
        if (iteration > 1L) {
            p <- p * (rho_cur / rho_prev)
            p <- p + z
        } else {
            p <- as.numeric(z)
        }
        q <- matvec(p)
        alpha <- rho_cur / sum(p * q)
        x <- x + alpha * p
        r <- r - alpha * q
        rho_prev <- rho_cur
    }
    list(x = x, info = as.integer(maxiter))
}

## scipy.sparse.linalg.bicgstab mirror (real case, x0 = 0).
.amg_bicgstab <- function(matvec, psolve, b, rtol, maxiter, atol = 0) {
    x <- numeric(length(b))
    r <- as.numeric(b)
    bnrm2 <- sqrt(sum(r^2))
    if (bnrm2 == 0) return(list(x = r, info = 0L))
    atol <- max(atol, rtol * bnrm2)
    if (sqrt(sum(r^2)) < atol) return(list(x = x, info = 0L))
    rhat <- r
    rho <- 1; alpha <- 1; omega_b <- 1
    v <- numeric(length(b))
    for (iteration in seq_len(maxiter)) {
        rho_cur <- sum(rhat * r)
        if (!is.finite(rho_cur) || rho_cur == 0) {
            return(list(x = x, info = as.integer(iteration)))
        }
        if (iteration == 1L) {
            z <- psolve(r)
            v <- as.numeric(matvec(z))
            p <- as.numeric(z)
        } else {
            beta <- (rho_cur / rho) * (alpha / omega_b)
            z_hat <- psolve(r)
            v_hat <- as.numeric(matvec(z_hat))
            p <- z_hat + beta * (p - omega_b * v)
            v <- v_hat + beta * (v - omega_b * v_hat)
        }
        denom <- sum(rhat * v)
        if (!is.finite(denom) || denom == 0) {
            return(list(x = x, info = as.integer(iteration)))
        }
        alpha <- rho_cur / denom
        if (!is.finite(alpha)) return(list(x = x, info = as.integer(iteration)))
        x <- x + alpha * p
        s <- r - alpha * v
        tvec <- as.numeric(matvec(s))
        tt <- sum(tvec * tvec)
        omega_b <- if (is.finite(tt) && tt > 0) sum(tvec * s) / tt else 0
        if (!is.finite(omega_b)) {
            rho <- rho_cur
            next
        }
        x <- x + omega_b * s
        r <- as.numeric(psolve(s - omega_b * tvec))
        if (!all(is.finite(r))) return(list(x = x, info = as.integer(iteration)))
        rho <- rho_cur
        if (sqrt(sum(r^2)) < atol) return(list(x = x, info = 0L))
    }
    list(x = x, info = as.integer(maxiter))
}

## Restarted preconditioned GMRES following scipy's defaults (restart = 20):
## Arnoldi basis, Givens-rotated least squares on the Hessenberg system,
## restarted until maxiter Krylov steps are used.
.amg_gmres <- function(matvec, psolve, b, rtol, maxiter, restart = 20L) {
    n <- length(b)
    x <- numeric(n)
    r <- as.numeric(b)
    bnrm2 <- sqrt(sum(r^2))
    if (bnrm2 == 0) return(list(x = r, info = 0L))
    atol <- rtol * bnrm2
    iteration <- 0L
    while (iteration < maxiter) {
        if (sqrt(sum(r^2)) <= atol) return(list(x = x, info = 0L))
        k <- as.integer(min(restart, maxiter - iteration))
        m <- sqrt(sum(r^2))
        V <- matrix(0, nrow = n, ncol = k + 1L)
        H <- matrix(0, nrow = k + 1L, ncol = k)
        V[, 1L] <- r / m
        jj <- 1L
        while (jj <= k) {
            w <- psolve(as.numeric(matvec(V[, jj])))
            for (i in seq_len(jj)) {
                H[i, jj] <- sum(w * V[, i])
                w <- w - H[i, jj] * V[, i]
            }
            H[jj + 1L, jj] <- sqrt(sum(w^2))
            if (H[jj + 1L, jj] <= 0) break
            V[, jj + 1L] <- w / H[jj + 1L, jj]
            jj <- jj + 1L
        }
        k <- jj - 1L
        if (k < 1L) return(list(x = x, info = as.integer(iteration)))
        Hloc <- H[seq_len(k + 1L), seq_len(k), drop = FALSE]
        g <- numeric(k + 1L)
        g[1L] <- m
        for (i in seq_len(k)) {
            denom <- sqrt(Hloc[i, i]^2 + Hloc[i + 1L, i]^2)
            cs <- if (denom == 0) 1 else Hloc[i, i] / denom
            sn <- if (denom == 0) 0 else Hloc[i + 1L, i] / denom
            for (j in seq_len(k)) {
                temp <- cs * Hloc[i, j] + sn * Hloc[i + 1L, j]
                Hloc[i + 1L, j] <- -sn * Hloc[i, j] + cs * Hloc[i + 1L, j]
                Hloc[i, j] <- temp
            }
            temp <- cs * g[i] + sn * g[i + 1L]
            g[i + 1L] <- -sn * g[i] + cs * g[i + 1L]
            g[i] <- temp
        }
        y <- tryCatch(as.numeric(backsolve(Hloc[seq_len(k), , drop = FALSE],
                                          g[seq_len(k)])),
                      error = function(e) numeric(k))
        x <- x + as.numeric(V[, seq_len(k), drop = FALSE] %*% y)
        r <- as.numeric(b - matvec(x))
        iteration <- iteration + k
    }
    list(x = x, info = as.integer(iteration))
}

# ---- Core solver -----------------------------------------------------------

#' Solve by AMG / preconditioned Krylov method
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial guess (length n) or \code{NULL} for zeros,
#'   matching \code{default_initial = np.zeros(n)} of the Python workflow.
#' @param method Character; \code{"cg"} (default), \code{"bicgstab"}, or
#'   \code{"gmres"}.  With \code{"cg"} the non-symmetric preconditioners
#'   \code{"gs"} and \code{"sor"} are replaced by \code{"ssor"}.
#' @param preconditioner Character; \code{"amg"} (default), \code{"jacobi"},
#'   \code{"gs"}, \code{"sor"}, \code{"ssor"}, or \code{"none"}.
#'   \code{"amg"} falls back to \code{"jacobi"} (warning), as \code{bssunfold}
#'   does when \code{pyamg} is unavailable.
#' @param omega Numeric; SOR/SSOR relaxation factor. Default 1.0,
#'   recommended range \eqn{(0, 2]}{(0, 2]}.
#' @param max_iterations Integer; Krylov iterations per outer restart.
#'   Default 200.
#' @param tolerance Numeric; relative residual tolerance of the normal
#'   equations. Default 1e-10.
#' @param outer_iterations Integer; projected non-negativity restarts.
#'   Default 3.
#' @param nonnegativity Logical; clamp the spectrum to \code{x >= 0}
#'   between restarts. Default TRUE.
#' @param regularization Numeric Tikhonov damping added to the diagonal of
#'   \eqn{A^\top A}{A'A}, or \code{NULL} (default) for the automatic
#'   \code{1e-4 * mean(diag(A'A))}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_amg(A, b, NULL)
solve_amg <- function(A, b, x0 = NULL, method = "cg", preconditioner = "amg",
                      omega = 1.0, max_iterations = 200L, tolerance = 1e-10,
                      outer_iterations = 3L, nonnegativity = TRUE,
                      regularization = NULL) {
    method <- tolower(as.character(method)[1L])
    if (!method %in% .amg_valid_methods) {
        stop("method must be one of ",
             paste(.amg_valid_methods, collapse = ", "),
             ", got '", method, "'")
    }
    preconditioner <- tolower(as.character(preconditioner)[1L])
    if (!preconditioner %in% .amg_valid_preconditioners) {
        stop("preconditioner must be one of ",
             paste(.amg_valid_preconditioners, collapse = ", "),
             ", got '", preconditioner, "'")
    }
    if (identical(method, "cg") && !preconditioner %in% .amg_cg_compatible) {
        warning("preconditioner='", preconditioner,
                "' is nonsymmetric and incompatible with method='cg'; ",
                "switching to 'ssor'", call. = FALSE)
        preconditioner <- "ssor"
    }
    max_iterations <- as.integer(max_iterations)
    if (is.na(max_iterations) || max_iterations < 1L) {
        stop("max_iterations must be a positive integer")
    }
    tolerance <- as.numeric(tolerance)
    if (is.na(tolerance) || tolerance <= 0) {
        stop("tolerance must be a positive number")
    }
    outer_iterations <- as.integer(outer_iterations)
    if (is.na(outer_iterations) || outer_iterations < 1L) {
        stop("outer_iterations must be >= 1, got ", outer_iterations)
    }
    omega <- as.numeric(omega)
    if (!(omega > 0 && omega <= 2)) {
        warning("omega=", omega, " outside the recommended range (0, 2]",
                call. = FALSE)
    }

    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)

    AT_A <- as.matrix(crossprod(A))
    AT_b <- as.numeric(crossprod(A, b))
    if (is.null(regularization)) {
        damping <- .amg_auto_reg_factor * mean(diag(AT_A))
    } else {
        damping <- as.numeric(regularization)
        if (length(damping) != 1L || is.na(damping)) {
            stop("regularization must be a single non-negative number")
        }
        if (damping < 0) {
            stop("regularization must be non-negative, got ", damping)
        }
    }
    Nmat <- AT_A + damping * diag(n)
    matvec <- function(v) as.numeric(Nmat %*% v)

    counter <- new.env(parent = emptyenv())
    counter$n <- 0L
    psolve <- .amg_build_preconditioner(Nmat, preconditioner, omega, counter)

    x <- if (is.null(x0)) numeric(n) else as.numeric(x0)
    if (length(x) != n) x <- numeric(n)
    converged <- FALSE
    inner_converged <- FALSE
    ## cg/gmres apply the preconditioner once per iteration, bicgstab twice.
    apps_per_iteration <- if (identical(method, "bicgstab")) 2L else 1L
    b_norm <- max(sqrt(sum(AT_b^2)), 1e-300)

    for (.outer in seq_len(outer_iterations)) {
        residual <- AT_b - matvec(x)
        if (sqrt(sum(residual^2)) <= tolerance * b_norm) {
            converged <- TRUE
            break
        }
        run <- if (identical(method, "cg")) {
            .amg_cg(matvec, psolve, residual, tolerance, max_iterations)
        } else if (identical(method, "bicgstab")) {
            .amg_bicgstab(matvec, psolve, residual, tolerance, max_iterations)
        } else {
            .amg_gmres(matvec, psolve, residual, tolerance, max_iterations)
        }
        dx <- run$x
        if (any(!is.finite(dx))) dx <- numeric(n)
        x <- x + dx
        if (isTRUE(nonnegativity)) x <- pmax(x, 0)
        if (run$info == 0L) {
            inner_converged <- TRUE
            if (sqrt(sum((AT_b - matvec(x))^2)) <= tolerance * b_norm) {
                converged <- TRUE
                break
            }
        } else if (run$info < 0L) {
            break
        }
    }

    converged <- converged || inner_converged
    total_iterations <- min(counter$n %/% apps_per_iteration,
                            outer_iterations * max_iterations)
    if (isTRUE(nonnegativity)) x <- pmax(x, 0)
    list(spectrum = as.numeric(x), iterations = as.integer(total_iterations),
         converged = converged)
}

#' Wrapper around \code{\link{solve_amg}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_amg
#' @param method_name Character; label stored in the result. Default
#'   \code{"AMG"}.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_amg <- function(detector_names, n_energy_bins, E_MeV, sensitivities,
                       cc_icrp116, save_result_callback, readings,
                       initial_spectrum = NULL, method = "cg",
                       preconditioner = "amg", omega = 1.0,
                       max_iterations = 200L, tolerance = 1e-10,
                       outer_iterations = 3L, nonnegativity = TRUE,
                       regularization = NULL, method_name = "AMG",
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
        default_initial = numeric(n_energy_bins),
        solve_func = make_solve_wrapper(solve_amg,
                                        method = method,
                                        preconditioner = preconditioner,
                                        omega = omega,
                                        max_iterations = max_iterations,
                                        tolerance = tolerance,
                                        outer_iterations = outer_iterations,
                                        nonnegativity = nonnegativity,
                                        regularization = regularization),
        solve_kwargs = list(),
        method_name = method_name,
        extra_output = list(krylov_method = method,
                            preconditioner = preconditioner,
                            omega = omega,
                            outer_iterations = outer_iterations,
                            regularization = regularization),
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
