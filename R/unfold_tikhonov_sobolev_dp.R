#' Tikhonov regularization with the generalized discrepancy principle
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_tikhonov_sobolev_dp.py}.
#' Port of the \code{alfaFinder()} method for the energy-non-invariant
#' apparatus function described in I. N. Ogorodnikov, "Inverse problems of
#' spectroscopy and spectrometry in applied research", Traektoriya
#' Issledovaniy no. 2 (10), pp. 42-83 (2024), sections 3 and 5.
#'
#' The method solves the discretized Fredholm integral equation of the first
#' kind \eqn{A z = b} by minimizing the Tikhonov smoothing functional
#' (article eq. 3.9)
#' \deqn{M[z] = || A z - b ||^2 + \alpha || z ||_W^2}{M[z] = ||Az - b||^2 +
#'   alpha ||z||_W^2}
#' where the Sobolev-space norm \eqn{W_2^1} of the article corresponds, in
#' its discrete form, to the squared norm of the first difference operator
#' \code{L = D1} (the Euler equation 3.10 is the continuous analogue of the
#' normal system \eqn{(A^T A + \alpha L^T L) z = A^T b}).
#'
#' The regularization parameter \code{alpha} is selected by the generalized
#' discrepancy principle (article eqs. 3.8-3.10): \code{alpha*} is the root
#' of the generalized discrepancy
#' \deqn{\rho(\alpha) = || A z_\alpha - b ||^2 - \delta^2}{rho(alpha) =
#'   ||A z_alpha - b||^2 - delta^2}
#' where \code{delta} is the RMS level of the measurement noise.  In the
#' article the root is found by the Newton and chord methods; here a
#' bracketing plus Brent-type scheme is used for robustness on
#' \code{log10(alpha)}.
#'
#' Like the article's \code{alfaFinder()}, the solver reports the achieved
#' discrepancy and diagnostic status codes when no admissible \code{alpha*}
#' exists (data inconsistent with the supplied \code{delta}).
#'
#' @keywords internal
#' @name tikhonov_sobolev_dp
#' @param n Integer; number of energy bins.
#' @param penalty Character; penalty operator \code{"sobolev"} (first
#'   difference, the discrete \eqn{W_2^1}{W_2^1} norm), \code{"curvature"}
#'   (second difference, \eqn{W_2^2}{W_2^2}) or \code{"identity"}
#'   (zeroth-order Tikhonov, standard ridge).
#' @return The discrete penalty operator \code{L} (shape \code{(k, n)}).
#' @examples
#' BSSUnfoldR:::.tikhonov_penalty_matrix(4L, "sobolev")
.tikhonov_penalty_matrix <- function(n, penalty = "sobolev") {
    n <- as.integer(n)
    if (penalty == "sobolev") {
        # First difference operator: discrete Sobolev W_2^1 norm.
        return(as.matrix(create_derivative_matrix(n, 1L)))
    }
    if (penalty == "curvature") {
        # Second difference operator: discrete W_2^2 norm.
        return(as.matrix(create_derivative_matrix(n, 2L)))
    }
    if (penalty == "identity") {
        # Zeroth-order Tikhonov (standard ridge).
        return(diag(n))
    }
    stop("Unsupported penalty: ", deparse(penalty),
         ". Choose from 'sobolev', 'curvature', 'identity'.")
}

# Solve the normal system (N + alpha K) z = rhs, falling back to the
# pseudo-inverse when the system is singular (np.linalg.solve / pinv).
.tikhonov_solve_regularized <- function(N, K, rhs, alpha) {
    matrix_ <- N + alpha * K
    out <- tryCatch(solve(matrix_, rhs), error = function(e) NULL)
    if (is.null(out)) {
        out <- .tikhonov_pinv(matrix_) %*% rhs
    }
    as.numeric(out)
}

# Moore-Penrose pseudo-inverse via SVD, using numpy's default rcond = 1e-15
# singular-value cutoff.
.tikhonov_pinv <- function(M) {
    s <- svd(M)
    tol <- 1e-15 * max(s$d) * max(nrow(M), ncol(M))
    keep <- s$d > tol
    ds <- 1 / s$d[keep]
    n <- length(ds)
    as.matrix(s$v[, seq_len(n), drop = FALSE] %*%
                  (t(s$u[, seq_len(n), drop = FALSE]) * ds))
}

# Generalized discrepancy rho(alpha) (article eq. 3.8):
#   rho(alpha) = ||A z_alpha - b||^2 - delta^2
# where z_alpha solves (A'A + alpha L'L) z = A'b.
.tikhonov_generalized_discrepancy <- function(alpha, N, K, rhs, A, b, delta_sq) {
    z_alpha <- .tikhonov_solve_regularized(N, K, rhs, alpha)
    residual <- as.numeric(A %*% z_alpha) - b
    as.numeric(residual %*% residual) - delta_sq
}

#' @rdname tikhonov_sobolev_dp
#' @param alpha Positive numeric; regularization parameter.
#' @param N,K Numeric matrices \eqn{A^T A}{A'A} and
#'   \eqn{L^T L}{L'L}.
#' @param rhs Numeric right-hand side \eqn{A^T b}{A'b}.
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param delta_sq Numeric; squared noise level \code{delta^2}.
#' @return \code{|| A z_alpha - b ||^2 - delta^2}.
#' @export
generalized_discrepancy <- function(alpha, N, K, rhs, A, b, delta_sq) {
    .tikhonov_generalized_discrepancy(as.numeric(alpha), N, K, rhs, A, b,
                                      as.numeric(delta_sq))
}

#' Find alpha* as the root of the generalized discrepancy
#'
#' Implements the article's step 2 of the regularizing algorithm: the
#' optimal regularization parameter \code{alpha*}, consistent with the data
#' error level, is computed as the root of \code{rho(alpha*) = 0}.  The root
#' is bracketed on the \code{log10} scale and refined with a Brent-type
#' bracketing solver (the combination of bisection, chord/secant and
#' inverse-quadratic methods used in the article).
#'
#' @rdname tikhonov_sobolev_dp
#' @param delta Positive numeric; RMS noise level.  \code{alpha*} satisfies
#'   \eqn{|| A z_{alpha*} - b ||^2 = \delta^2}{||A z - b||^2 = delta^2}.
#' @param L Optional penalty operator (k x n).  Default \code{NULL} = the
#'   first-difference (Sobolev \eqn{W_2^1}) operator.
#' @param alpha_range Numeric vector of length 2; search interval for
#'   \code{alpha}.  Default \code{c(1e-10, 1e10)}.
#' @param max_iter Positive integer; maximum number of root-finder
#'   iterations.  Default 100.
#' @param xtol Positive numeric; absolute tolerance on \code{log10(alpha)}.
#'   Default 1e-8.
#' @return A list with components \code{alpha}, \code{residual_sq},
#'   \code{rho}, \code{status}, \code{converged}, \code{n_iter},
#'   \code{min_residual_sq} and \code{max_residual_sq}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' info <- alpha_finder_generalized_discrepancy(A, b, delta = 0.02 * sqrt(sum(b^2)))
alpha_finder_generalized_discrepancy <- function(A, b, delta, L = NULL,
                                                 alpha_range = c(1e-10, 1e10),
                                                 max_iter = 100L,
                                                 xtol = 1e-8) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)

    if (is.null(L)) {
        L <- .tikhonov_penalty_matrix(n, "sobolev")
    }
    L <- as.matrix(L); storage.mode(L) <- "double"

    delta <- as.numeric(delta)
    if (delta <= 0) stop("delta must be positive, got ", delta)
    delta_sq <- delta * delta

    N <- t(A) %*% A
    K <- t(L) %*% L
    rhs <- as.numeric(t(A) %*% b)

    alpha_min <- as.numeric(alpha_range[1])
    alpha_max <- as.numeric(alpha_range[2])
    if (!(0 < alpha_min && alpha_min < alpha_max)) {
        stop("alpha_range must satisfy 0 < alpha_min < alpha_max, got (",
             alpha_range[1], ", ", alpha_range[2], ")")
    }

    n_calls <- 0L
    rho_fun <- function(alpha) {
        n_calls <<- n_calls + 1L
        .tikhonov_generalized_discrepancy(alpha, N, K, rhs, A, b, delta_sq)
    }
    rho_log10 <- function(t) rho_fun(10^t)

    n_iter <- 0L
    rho_min <- rho_fun(alpha_min)
    rho_max <- rho_fun(alpha_max)
    n_iter <- n_iter + 2L

    # Residual of the plain least-squares end (alpha -> 0): the minimal
    # achievable data misfit, used for diagnostics.
    z_ls <- as.numeric(qr.coef(qr(A), b))
    r_ls <- as.numeric(A %*% z_ls) - b
    min_residual_sq <- as.numeric(r_ls %*% r_ls)

    if (rho_min > 0) {
        # Even the least-regularized solution misfits more than delta:
        # no admissible alpha in the range (article: rho(alpha*) > 0).
        z_alpha <- .tikhonov_solve_regularized(N, K, rhs, alpha_min)
        r_alpha <- as.numeric(A %*% z_alpha) - b
        return(list(alpha = alpha_min,
                    residual_sq = as.numeric(r_alpha %*% r_alpha),
                    rho = rho_min,
                    status = 1L,
                    converged = FALSE,
                    n_iter = n_iter,
                    min_residual_sq = min_residual_sq,
                    max_residual_sq = rho_max + delta_sq))
    }

    if (rho_max < 0) {
        # Even the most-regularized solution fits better than delta:
        # delta is larger than any achievable misfit.
        z_alpha <- .tikhonov_solve_regularized(N, K, rhs, alpha_max)
        r_alpha <- as.numeric(A %*% z_alpha) - b
        return(list(alpha = alpha_max,
                    residual_sq = as.numeric(r_alpha %*% r_alpha),
                    rho = rho_max,
                    status = 2L,
                    converged = FALSE,
                    n_iter = n_iter,
                    min_residual_sq = min_residual_sq,
                    max_residual_sq = rho_max + delta_sq))
    }

    # rho is non-decreasing in alpha: bracket on the log10 grid and refine
    # with the Brent-type root finder.
    t_min <- log10(alpha_min)
    t_max <- log10(alpha_max)
    if (rho_min == 0) {
        root <- t_min
    } else if (rho_max == 0) {
        root <- t_max
    } else {
        fit <- stats::uniroot(rho_log10, c(t_min, t_max), tol = xtol,
                              maxiter = as.integer(max_iter),
                              extendInt = "no", check.conv = TRUE)
        root <- fit$root
    }
    n_iter <- n_iter + n_calls
    alpha_star <- as.numeric(10^root)
    z_alpha <- .tikhonov_solve_regularized(N, K, rhs, alpha_star)
    r_alpha <- as.numeric(A %*% z_alpha) - b
    residual_sq <- as.numeric(r_alpha %*% r_alpha)

    list(alpha = alpha_star,
         residual_sq = residual_sq,
         rho = residual_sq - delta_sq,
         status = 0L,
         converged = TRUE,
         n_iter = n_iter,
         min_residual_sq = min_residual_sq,
         max_residual_sq = rho_max + delta_sq)
}

#' Solve unfolding with Tikhonov + generalized discrepancy principle
#'
#' Core solver mirroring \code{solve_tikhonov_sobolev_dp} in
#' \code{bssunfold/src/bssunfold/core/unfold_tikhonov_sobolev_dp.py}.
#' Minimizes \eqn{||A z - b||^2 + \alpha ||L z||^2} with \code{alpha*} chosen
#' so that the data misfit matches the noise level
#' (\eqn{||A z - b||^2 = \delta^2}), following the article's
#' \code{alfaFinder()}.  By default \code{L} is the first-difference operator
#' (discrete Sobolev \eqn{W_2^1} penalty of the article, Euler equation 3.10)
#' and \eqn{\delta = noise\_level \cdot ||b||}{delta = noise_level * ||b||}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial spectrum guess (accepted for API
#'   compatibility; not used by the linear regularized solver).
#'   Default \code{NULL}.
#' @param noise_level Numeric; relative noise level used to derive
#'   \code{delta} when it is not given explicitly.  Default 0.02 (i.e. 2 %
#'   as in the article).
#' @param delta Optional explicit RMS noise level; overrides
#'   \code{noise_level}.  Default \code{NULL}.
#' @param penalty Character; \code{"sobolev"} (first difference, default),
#'   \code{"curvature"} (second difference) or \code{"identity"}.
#' @param alpha_range Numeric vector of length 2; search interval for the
#'   regularization parameter.  Default \code{c(1e-10, 1e10)}.
#' @param max_iter Positive integer; maximum number of root-finder
#'   iterations.  Default 100.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_tikhonov_sobolev_dp(A, b, rep(0, 3))
solve_tikhonov_sobolev_dp <- function(A, b, x0 = NULL, noise_level = 0.02,
                                      delta = NULL, penalty = "sobolev",
                                      alpha_range = c(1e-10, 1e10),
                                      max_iter = 100L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)

    if (is.null(delta)) {
        delta <- as.numeric(noise_level) * sqrt(sum(b^2))
    }
    delta <- as.numeric(delta)
    if (delta <= 0) {
        stop("delta must be positive, got ", delta,
             "; provide delta or a positive noise_level")
    }

    n <- ncol(A)
    L <- .tikhonov_penalty_matrix(n, penalty)

    info <- alpha_finder_generalized_discrepancy(
        A, b, delta = delta, L = L, alpha_range = alpha_range,
        max_iter = max_iter)

    N <- t(A) %*% A
    K <- t(L) %*% L
    rhs <- as.numeric(t(A) %*% b)
    spectrum <- .tikhonov_solve_regularized(N, K, rhs, info$alpha)

    list(spectrum = spectrum, iterations = as.integer(info$n_iter),
         converged = as.logical(info$converged))
}

#' Tikhonov + discrepancy-principle unfolding (unified workflow wrapper)
#'
#' Thin wrapper around \code{\link{solve_tikhonov_sobolev_dp}} for the
#' unified workflow, mirroring \code{unfold_tikhonov_sobolev_dp} in
#' \code{bssunfold/src/bssunfold/core/unfold_tikhonov_sobolev_dp.py}.
#'
#' Port of the article's \code{alfaFinder()} for multisphere data: Tikhonov
#' regularization with the discrete Sobolev \eqn{W_2^1} penalty and the
#' regularization parameter selected by the generalized discrepancy
#' principle \eqn{||A z - b||^2 = \delta^2}.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_tikhonov_sobolev_dp
#' @param calculate_errors Logical; if \code{TRUE}, run Monte-Carlo
#'   uncertainty estimation.  Default \code{FALSE}.
#' @param n_montecarlo Integer; number of Monte-Carlo samples.
#'   Default 100.
#' @param random_state Optional integer seed for Monte-Carlo.
#' @param max_neutron_energy Optional numeric energy cutoff in MeV.
#'   Default \code{NULL} = no cutoff.
#' @return A result list as produced by \code{\link{run_unfolding}}, with
#'   additional entries \code{alpha}, \code{delta},
#'   \code{discrepancy_status} and \code{dp_converged}.
#' @export
unfold_tikhonov_sobolev_dp <- function(detector_names, n_energy_bins, E_MeV,
                                       sensitivities, cc_icrp116,
                                       save_result_callback, readings,
                                       initial_spectrum = NULL,
                                       noise_level = 0.02, delta = NULL,
                                       penalty = "sobolev",
                                       alpha_range = c(1e-10, 1e10),
                                       max_iter = 100L,
                                       calculate_errors = FALSE,
                                       n_montecarlo = 100L,
                                       save_result = FALSE,
                                       random_state = NULL,
                                       max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)

    holder <- new.env(parent = emptyenv())

    b_selected <- as.numeric(readings[detector_names[detector_names %in%
                                                         names(readings)]])

    solve_func <- function(A, b, x0 = NULL, ...) {
        res <- solve_tikhonov_sobolev_dp(A, b, x0 = NULL,
                                         noise_level = noise_level,
                                         delta = delta, penalty = penalty,
                                         alpha_range = alpha_range,
                                         max_iter = max_iter)
        spectrum <- res$spectrum
        n_iter <- res$iterations
        converged <- res$converged
        # Record only the clean-fit metadata (Monte-Carlo replicates use
        # perturbed readings and would overwrite them).
        if (identical(as.numeric(b), b_selected)) {
            used_delta <- delta
            if (is.null(used_delta)) {
                used_delta <- as.numeric(noise_level) * sqrt(sum(b^2))
            }
            info <- alpha_finder_generalized_discrepancy(
                A, b, delta = used_delta,
                L = .tikhonov_penalty_matrix(ncol(A), penalty),
                alpha_range = alpha_range, max_iter = max_iter)
            holder$alpha <- info$alpha
            holder$discrepancy_status <- info$status
            holder$dp_converged <- info$converged
            holder$residual_sq <- info$residual_sq
        }
        list(spectrum = spectrum, iterations = n_iter, converged = converged)
    }

    result <- run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solve_func,
        solve_kwargs = list(),
        method_name = "TikhonovSobolevDP",
        extra_output = list(penalty = penalty,
                            delta = if (is.null(delta)) NULL else as.numeric(delta),
                            noise_level = as.numeric(noise_level)),
        calculate_errors = calculate_errors,
        noise_level = if (is.null(noise_level) || noise_level == 0)
                          0.01 else noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)

    for (key in c("alpha", "discrepancy_status", "dp_converged", "residual_sq")) {
        if (!is.null(holder[[key]])) result[[key]] <- holder[[key]]
    }

    result
}
