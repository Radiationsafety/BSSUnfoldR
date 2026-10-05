#' NNQP-based unfolding method for neutron spectrum reconstruction
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_nnqp.py}.
#' Method built on the NNQP (Non-Negative Quadratic Programming) solver of
#' Giovannucci & Pehlevan (https://github.com/simonsfoundation/NNQP).
#'
#' NNQP solves, by coordinate descent, the convex program
#' \deqn{\minimize \frac{1}{2} x' Q x + f' x \quad \text{s.t. } x \ge 0}
#' For BSS unfolding the regularised least-squares problem is recast as
#' \deqn{\minimize \frac{1}{2} ||A x - b||^2 + \frac{\alpha}{2} ||L x||^2
#'   + \frac{\alpha_0}{2} ||x||^2 \quad \text{s.t. } x \ge 0}{min 1/2 ||Ax -
#'   b||^2 + alpha/2 ||Lx||^2 + alpha0/2 ||x||^2, x >= 0}
#' which has the QP form above with
#' \deqn{Q = A' A + \alpha L' L + \alpha_0 I, \qquad f = -A' b.}{Q = A'A +
#'   alpha L'L + alpha0 I, f = -A'b.}
#' The coordinate-descent update for coordinate \eqn{i} is
#' \deqn{x_i \leftarrow \max(0, -(M_i \cdot x + f_i) / Q_{ii}), \qquad
#'   M = Q - \operatorname{diag}(Q)}{x_i = max(0, -(M_i.x + f_i)/Q_ii)}
#' which is the classic NNQP update.  The diagonal entries \eqn{Q_{ii}} are
#' always strictly positive because \eqn{Q} is the sum of two
#' positive-semidefinite matrices plus \eqn{\alpha_0 I}{alpha0 * I}, so the
#' regularisation floor \eqn{\alpha_0}{alpha0} (default 1e-6) guarantees
#' strict positive-definiteness even when \eqn{A'A} is rank-deficient —
#' which is the norm for BSS problems (few detectors, many energy bins).
#'
#' The implementation is a self-contained base-R port of the original
#' \code{nnqp.py} (which used \code{numba}): the same algorithm, with the
#' Gauss-Seidel sweep executed as an explicit loop because BSS problems
#' typically have many energy bins.
#'
#' @param Q Numeric symmetric positive-definite Hessian (n x n).
#' @param f Numeric linear term (length n).
#' @param x0 Optional warm-start spectrum (length n).  If \code{NULL} a
#'   uniform \code{[0, 1)} random vector is used (matching the original
#'   implementation).  Default \code{NULL}.
#' @param tol Positive numeric; convergence tolerance on the relative change
#'   of the active components.  Default 1e-6.
#' @param max_iterations Positive integer; iteration cap.  Default 10000.
#' @param random_state Optional integer seed for the initial guess when
#'   \code{x0} is \code{NULL}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @keywords internal
#' @examples
#' Q <- matrix(c(4, 1, 1, 3), nrow = 2)
#' f <- c(-2, -1)
#' BSSUnfoldR:::.nnqp(Q, f, x0 = c(0.5, 0.5))
.nnqp <- function(Q, f, x0 = NULL, tol = 1e-6, max_iterations = 10000L,
                  random_state = NULL) {
    Q <- as.matrix(Q); storage.mode(Q) <- "double"
    f <- as.numeric(f)
    if (nrow(Q) != ncol(Q)) {
        stop("Q must be a square matrix, got shape ", nrow(Q), "x", ncol(Q))
    }
    n <- nrow(Q)
    if (length(f) != n) {
        stop("f length (", length(f), ") must match Q size (", n, ")")
    }

    # Symmetrise defensively (the algorithm relies on the symmetry of Q).
    Q <- 0.5 * (Q + t(Q))

    qdg <- diag(Q)
    # The diagonal must be strictly positive for the coordinate update to be
    # well-defined.  If the user passes a rank-deficient Q, we add a tiny
    # floor -- but only here at the lowest level so callers can override.
    bad <- qdg <= 0
    if (any(bad)) {
        qdg[bad] <- pmax(qdg[bad], 1e-12)
    }
    Dinv <- 1 / qdg
    M <- Q - diag(qdg)  # off-diagonal part

    if (!is.null(x0)) {
        x <- pmax(as.numeric(x0), 0)
        if (length(x) != n) {
            stop("x0 length (", length(x), ") must match Q size (", n, ")")
        }
    } else {
        if (!is.null(random_state)) set.seed(as.integer(random_state))
        x <- stats::runif(n)
    }

    converged <- FALSE
    iterations <- 0L
    zero_count <- 0L
    for (it in seq_len(as.integer(max_iterations))) {
        iterations <- it
        x_prev <- x
        # Coordinate-descent sweep -- update in place, picking up the latest
        # values of already-updated coordinates (Gauss-Seidel flavour).
        # Since M has zero diagonal, M[i, ] %*% x depends only on the other
        # coordinates, so a Gauss-Seidel sweep is correct.
        for (i in seq_len(n)) {
            dum <- Dinv[i] * (-(as.numeric(M[i, , drop = TRUE] %*% x)) - f[i])
            x[i] <- if (dum > 0) dum else 0
        }

        # Convergence test: relative change on the strictly positive
        # components of x_prev (as in the original NNQP).
        active <- x_prev > 0
        if (!any(active)) {
            # First iteration can land on all-zero x; give it one more
            # sweep before declaring convergence.
            if (zero_count == 0L) {
                er <- 1.0
                zero_count <- 1L
            } else {
                er <- 0
            }
        } else {
            denom <- pmax(abs(x_prev[active]), 1e-300)
            er <- max(abs(x[active] - x_prev[active]) / denom)
        }

        if (er <= tol) {
            converged <- TRUE
            break
        }
    }

    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Solve the BSS unfolding problem with NNQP (coordinate descent)
#'
#' Core solver mirroring \code{solve_nnqp} in
#' \code{bssunfold/src/bssunfold/core/unfold_nnqp.py}.  Recasts the
#' regularised non-negative least-squares problem
#' \deqn{\minimize \frac{1}{2} ||A x - b||^2 + \frac{\alpha}{2} ||L x||^2
#'   + \frac{\alpha_0}{2} ||x||^2 \quad \text{s.t. } x \ge 0}{min 1/2 ||Ax -
#'   b||^2 + alpha/2 ||Lx||^2 + alpha0/2 ||x||^2, x >= 0}
#' as the NNQP program with \eqn{Q = A' A + \alpha L' L + \alpha_0 I}{Q =
#' A'A + alpha L'L + alpha0 I} and \eqn{f = -A' b}{f = -A'b}, where \eqn{L}
#' is the finite-difference derivative matrix of order
#' \code{smoothness_order} (0 disables the smoothness term — \eqn{L} is
#' empty) and \eqn{\alpha_0}{alpha0} is the diagonal regularisation floor.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional warm-start spectrum (length n).  If \code{NULL} the
#'   NNQP solver uses a uniform random initial guess (matching the original
#'   implementation).  Default \code{NULL}.
#' @param regularization Numeric; Tikhonov / smoothness regularisation
#'   weight.  Default 1e-4.
#' @param smoothness_order Integer; smoothness penalty order (0, 1 or 2).
#'   Default 0.
#' @param smoothness_weight Numeric; weight for the smoothness term.
#'   Default 1.0.
#' @param tol Positive numeric; convergence tolerance.  Default 1e-6.
#' @param max_iterations Positive integer; iteration cap.  Default 10000.
#' @param floor Numeric; diagonal regularisation floor added to \code{Q} to
#'   guarantee strict positive-definiteness even when \eqn{A'A} is
#'   rank-deficient.  Default 1e-6.
#' @param random_state Optional integer seed for the initial guess when
#'   \code{x0} is \code{NULL}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_nnqp(A, b, rep(0, 3), max_iterations = 200L)
solve_nnqp <- function(A, b, x0 = NULL, regularization = 1e-4,
                       smoothness_order = 0L, smoothness_weight = 1.0,
                       tol = 1e-6, max_iterations = 10000L, floor = 1e-6,
                       random_state = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)

    if (!smoothness_order %in% c(0L, 1L, 2L)) {
        stop("Unsupported smoothness order: ", smoothness_order, ". Use 0, 1 or 2.")
    }

    # Build Q = A'A + alpha * L'L + alpha0 * I
    AtA <- t(A) %*% A
    Q <- AtA + floor * diag(n)
    smoothness_order <- as.integer(smoothness_order)
    if (smoothness_order %in% c(1L, 2L) && regularization > 0) {
        L <- as.matrix(create_derivative_matrix(n, smoothness_order))
        Q <- Q + as.numeric(regularization) * as.numeric(smoothness_weight) *
            (t(L) %*% L)
    }

    f <- -as.numeric(t(A) %*% b)

    .nnqp(Q, f, x0 = x0, tol = tol, max_iterations = max_iterations,
          random_state = random_state)
}

#' NNQP unfolding (unified workflow wrapper)
#'
#' Thin wrapper around \code{\link{solve_nnqp}} for the unified workflow,
#' mirroring \code{unfold_nnqp} in
#' \code{bssunfold/src/bssunfold/core/unfold_nnqp.py}.
#'
#' Solves \eqn{\minimize \frac{1}{2} ||A x - b||^2 + \frac{\alpha}{2} ||L x||^2
#' + \frac{\alpha_0}{2} ||x||^2}{min 1/2 ||Ax-b||^2 + alpha/2 ||Lx||^2 +
#' alpha0/2 ||x||^2} subject to \eqn{x \ge 0}{x >= 0}, where \eqn{L} is the
#' finite-difference derivative matrix of order \code{smoothness_order},
#' using the coordinate-descent NNQP solver of Giovannucci & Pehlevan.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_nnqp
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
unfold_nnqp <- function(detector_names, n_energy_bins, E_MeV,
                        sensitivities, cc_icrp116, save_result_callback,
                        readings, initial_spectrum = NULL,
                        regularization = 1e-4, smoothness_order = 0L,
                        smoothness_weight = 1.0, tol = 1e-6,
                        max_iterations = 10000L, floor = 1e-6,
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
        solve_func = make_solve_wrapper(solve_nnqp,
                                        regularization = regularization,
                                        smoothness_order = smoothness_order,
                                        smoothness_weight = smoothness_weight,
                                        tol = tol,
                                        max_iterations = max_iterations,
                                        floor = floor,
                                        random_state = random_state),
        solve_kwargs = list(),
        method_name = "NNQP",
        extra_output = list(regularization = regularization,
                            smoothness_order = as.integer(smoothness_order),
                            smoothness_weight = smoothness_weight,
                            tol = tol,
                            max_iterations = as.integer(max_iterations),
                            floor = floor),
        calculate_errors = calculate_errors,
        noise_level = noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
