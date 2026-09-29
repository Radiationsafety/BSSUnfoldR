#' Solve a least-squares problem as NumPy's ``lstsq`` does
#'
#' Private helper used by \code{\link{solve_hybrid_gmres}}.
#'
#' @param Mmat Numeric matrix.
#' @param rhs Numeric right-hand side vector.
#' @return Numeric solution vector.
#' @keywords internal
.hgmres_lstsq <- function(Mmat, rhs) {
    ## numpy.linalg.lstsq(..., rcond = NULL): SVD cut-off eps * max(M, N).
    sv <- svd(Mmat)
    tol <- .Machine$double.eps * max(nrow(Mmat), ncol(Mmat)) * sv$d[1L]
    keep <- sv$d > tol
    dinv <- numeric(length(sv$d))
    dinv[keep] <- 1 / sv$d[keep]
    utb <- as.numeric(crossprod(sv$u, matrix(rhs, ncol = 1L)))
    as.numeric(sv$v %*% (dinv * utb))
}

.hgmres_pinv <- function(Mmat) {
    ## numpy.linalg.pinv default rcond = 1e-15 (relative to the largest s).
    sv <- svd(Mmat)
    cutoff <- 1e-15 * sv$d[1L]
    keep <- sv$d > cutoff
    dinv <- numeric(length(sv$d))
    dinv[keep] <- 1 / sv$d[keep]
    tcrossprod(sv$v * rep(dinv, each = nrow(sv$v)), sv$u)
}

## GCV value of the projected regularized problem (Python _gcv_function).
.hgmres_gcv <- function(lambda_val, B_k, beta_vec) {
    k <- ncol(B_k)
    out <- tryCatch({
        if (lambda_val < 1e-14) {
            bvec <- beta_vec[seq_len(k)]
            x_lambda <- .hgmres_lstsq(B_k, bvec)
            residual <- bvec - as.numeric(B_k %*% x_lambda)
        } else {
            B_reg <- rbind(B_k, lambda_val * diag(k))
            rhs <- c(beta_vec, numeric(k))
            x_lambda <- .hgmres_lstsq(B_reg, rhs)
            residual <- beta_vec - as.numeric(B_k %*% x_lambda)
        }
        numerator <- sum(residual^2)
        if (k < 50L) {
            H <- B_k %*% .hgmres_pinv(B_k)
            denominator <- (length(beta_vec) - sum(diag(H)))^2
        } else {
            denominator <- max(length(beta_vec) * 0.1, 1)
        }
        if (denominator < 1e-10) 1e10 else numerator / denominator
    }, error = function(e) 1e10)
    out
}

#' Hybrid GMRES unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_hybrid_gmres.py}.
#' The method is a \emph{hybrid} GMRES: a Golub-Kahan bidiagonalization
#' builds two orthonormal bases, \code{V} in detector space (length
#' \code{m}) and \code{U} in solution space (length \code{n}), together with
#' the lower bidiagonal projection \code{B_k}.  At every Krylov dimension
#' the projected problem
#' \code{min || B_k y - beta_0 e_1 ||^2 + lambda^2 ||y||^2} is solved and the
#' regularization parameter is chosen by the generalized cross validation
#' (GCV) curve.  The iterate with the smallest GCV value over all Krylov
#' dimensions is returned (this is the classic IRtools
#' \code{IRhybrid_gmres} strategy), and the recovered spectrum is
#' \code{x0 + U_k y}.
#'
#' @section Limitations:
#' The GCV denominator uses \code{trace(B_k B_k^+)} exactly as the Python
#' original does, and the pseudo-inverse cut-off mirrors NumPy's default
#' (\code{rcond = 1e-15}).  The least-squares solves mirror
#' \code{numpy.linalg.lstsq} (SVD with the machine-precision cut-off).
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional numeric initial guess (length n). Default \code{NULL}
#'   = zero spectrum, as in the Python original.
#' @param max_iterations Positive integer; max Krylov dimension. Default 100.
#' @param regularization_method Character; \code{"gcv"} (default) or
#'   \code{"modgcv"} (identical search grid, as in Python), \code{"discrep"}
#'   (discrepancy principle, requires \code{noise_level}), or any other
#'   string to use the fixed \code{regularization}.
#' @param regularization Numeric; fixed regularization parameter, also used
#'   as the starting value of the discrepancy search. Default 0.0.
#' @param noise_level Optional numeric; relative noise for the discrepancy
#'   principle. Default \code{NULL}.
#' @param eta Numeric; safety factor for the discrepancy principle.
#'   Default 1.01.
#' @param reorthogonalization Logical; re-orthogonalize each new basis vector
#'   against the previous ones. Default \code{TRUE}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_hybrid_gmres(A, b, NULL, max_iterations = 3)
solve_hybrid_gmres <- function(A, b, x0 = NULL, max_iterations = 100L,
                                  regularization_method = "gcv",
                                  regularization = 0.0,
                                  noise_level = NULL, eta = 1.01,
                                  reorthogonalization = TRUE) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A); n <- ncol(A)
    max_krylov <- as.integer(min(max_iterations, m, n))
    if (max_krylov < 1L) max_krylov <- 1L
    x0v <- if (is.null(x0)) numeric(n) else as.numeric(x0)
    if (length(x0v) != n) x0v <- numeric(n)

    # ---- Golub-Kahan bidiagonalization (Python lines 207-313) -----------
    r0 <- b - as.numeric(A %*% x0v)
    beta0 <- sqrt(sum(r0^2))
    if (beta0 < 1e-14) {
        return(list(spectrum = pmax(x0v, 0), iterations = 0L,
                    converged = TRUE))
    }

    U <- matrix(0, nrow = n, ncol = max_krylov + 1L)   # solution-space basis
    Vm <- matrix(0, nrow = m, ncol = max_krylov + 1L)  # detector-space basis
    alpha <- numeric(max_krylov)                       # diagonal of B
    beta <- numeric(max_krylov + 1L)                   # sub-diagonal of B
    beta[1L] <- beta0
    Vm[, 1L] <- r0 / beta0
    u1 <- as.numeric(crossprod(A, Vm[, 1L]))
    alpha[1L] <- sqrt(sum(u1^2))
    if (alpha[1L] < 1e-14) {
        return(list(spectrum = pmax(x0v, 0), iterations = 0L,
                    converged = TRUE))
    }
    U[, 1L] <- u1 / alpha[1L]

    gcv_values <- numeric(0)
    best_gcv <- Inf
    best_gcv_val <- Inf
    best_solution <- NULL
    y_lambda <- NULL
    stop_iteration <- max_krylov
    actual_k <- max_krylov

    for (idx in seq_len(max_krylov)) {
        # v_{k+1} = A u_k - alpha_k v_k
        v_new <- as.numeric(A %*% U[, idx]) - alpha[idx] * Vm[, idx]
        if (isTRUE(reorthogonalization) && idx > 1L) {
            for (j in seq_len(idx - 1L)) {
                v_new <- v_new - sum(Vm[, j] * v_new) * Vm[, j]
            }
        }
        beta[idx + 1L] <- sqrt(sum(v_new^2))
        if (beta[idx + 1L] < 1e-14) {   # breakdown
            actual_k <- idx
            break
        }
        Vm[, idx + 1L] <- v_new / beta[idx + 1L]

        # u_{k+1} = A' v_{k+1} - beta_{k+1} u_k
        u_new <- as.numeric(crossprod(A, Vm[, idx + 1L])) -
            beta[idx + 1L] * U[, idx]
        if (isTRUE(reorthogonalization) && idx > 1L) {
            for (j in seq_len(idx - 1L)) {
                u_new <- u_new - sum(U[, j] * u_new) * U[, j]
            }
        }
        if (idx < max_krylov) {
            alpha[idx + 1L] <- sqrt(sum(u_new^2))
            if (alpha[idx + 1L] < 1e-14) {   # alpha breakdown
                actual_k <- idx
                break
            }
            U[, idx + 1L] <- u_new / alpha[idx + 1L]
        }

        current_k <- idx
        # Bidiagonal projection B_k: (current_k + 1) x current_k
        B_k <- matrix(0, nrow = current_k + 1L, ncol = current_k)
        for (i in seq_len(current_k)) {
            B_k[i, i] <- alpha[i]
            B_k[i + 1L, i] <- beta[i + 1L]
        }
        rhs_proj <- numeric(current_k + 1L)
        rhs_proj[1L] <- beta[1L]

        # ---- regularization parameter ------------------------------------
        if (regularization_method %in% c("gcv", "modgcv")) {
            lambdas <- 10^seq(-10, 2, length.out = 50L)
            best_lambda <- regularization
            best_gcv_val <- Inf
            for (lam in lambdas) {
                gcv_val <- .hgmres_gcv(lam, B_k, rhs_proj)
                if (gcv_val < best_gcv_val) {
                    best_gcv_val <- gcv_val
                    best_lambda <- lam
                }
            }
            lambda_k <- best_lambda
            gcv_values <- c(gcv_values, best_gcv_val)
            if (length(gcv_values) >= 3L) {
                recent_min <- min(gcv_values[(length(gcv_values) - 2L):
                                             length(gcv_values)])
                if (best_gcv_val > recent_min * 1.01) {
                    stop_iteration <- current_k
                }
            }
        } else if (identical(regularization_method, "discrep") &&
                   !is.null(noise_level)) {
            threshold <- eta * as.numeric(noise_level) * sqrt(sum(b^2))
            lambda_k <- regularization
            for (.step in seq_len(20L)) {
                L_reg <- rbind(B_k, lambda_k * diag(current_k))
                rhs_reg <- c(rhs_proj, numeric(current_k))
                y_try <- tryCatch(.hgmres_lstsq(L_reg, rhs_reg),
                                  error = function(e) NULL)
                if (is.null(y_try)) break
                x_try <- x0v + as.numeric(U[, seq_len(current_k),
                                            drop = FALSE] %*% y_try)
                res_norm <- sqrt(sum((as.numeric(A %*% x_try) - b)^2))
                lambda_k <- if (res_norm > threshold) {
                    lambda_k * 2
                } else {
                    lambda_k / 2
                }
            }
        } else {
            lambda_k <- regularization
        }

        # ---- solve the regularized projected problem ---------------------
        L_reg <- rbind(B_k, lambda_k * diag(current_k))
        rhs_reg <- c(rhs_proj, numeric(current_k))
        solved <- tryCatch({
            yl <- .hgmres_lstsq(L_reg, rhs_reg)
            xl <- x0v + as.numeric(U[, seq_len(current_k), drop = FALSE] %*% yl)
            list(y = yl, x = xl)
        }, error = function(e) NULL)
        if (is.null(solved)) next
        y_lambda <- solved$y
        x_lambda <- solved$x

        if (best_gcv_val < best_gcv) {
            best_gcv <- best_gcv_val
            best_solution <- x_lambda
        }
    }

    spectrum <- if (is.null(best_solution)) {
        if (is.null(y_lambda)) {
            pmax(x0v, 0)
        } else {
            pmax(x0v + as.numeric(U[, seq_len(actual_k), drop = FALSE] %*%
                                  y_lambda), 0)
        }
    } else {
        pmax(best_solution, 0)
    }

    list(spectrum = as.numeric(spectrum),
         iterations = as.integer(stop_iteration),
         converged = stop_iteration < max_krylov)
}

#' Wrapper around \code{\link{solve_hybrid_gmres}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_hybrid_gmres
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_hybrid_gmres <- function(detector_names, n_energy_bins, E_MeV,
                                   sensitivities, cc_icrp116,
                                   save_result_callback, readings,
                                   initial_spectrum = NULL,
                                   max_iterations = 100L,
                                   regularization_method = "gcv",
                                   regularization = 0.0,
                                   noise_level = NULL, eta = 1.01,
                                   reorthogonalization = TRUE,
                                   calculate_errors = FALSE,
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
        solve_func = make_solve_wrapper(solve_hybrid_gmres,
                                         max_iterations = max_iterations,
                                         regularization_method = regularization_method,
                                         regularization = regularization,
                                         noise_level = noise_level,
                                         eta = eta,
                                         reorthogonalization = reorthogonalization),
        solve_kwargs = list(),
        method_name = "Hybrid_GMRES",
        extra_output = list(regularization_method = regularization_method,
                            regularization = regularization),
        calculate_errors = calculate_errors,
        noise_level = if (is.null(noise_level)) 0.01 else noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
