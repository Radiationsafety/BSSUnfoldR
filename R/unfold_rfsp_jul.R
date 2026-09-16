#' RFSP-JUL unfolding (damped least squares)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_rfsp_jul.py}.
#' RFSP-JUL is an iterative, damped least-squares method. At each iteration
#' it minimises
#' \deqn{S^{(k)} = \sum_i W_i [(b_i - (A \phi)_i) / b_i]^2 + \sum_j [(\phi_j - \phi_{prev,j}) / \phi_{prev,j}]^2}
#' Both terms are quadratic in \eqn{\phi}, so the minimiser is found in
#' closed form from the normal equations.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial guess (length n).
#' @param max_iterations Positive integer; default 200.
#' @param tolerance Positive numeric; default 1e-4.
#' @param weights Optional per-detector weights.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_rfsp_jul(A, b, rep(1, 3), max_iterations = 50)
solve_rfsp_jul <- function(A, b, x0, max_iterations = 200L, tolerance = 1e-4,
                             weights = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    if (nrow(A) == 0L || length(b) == 0L) {
        stop("Response matrix and measurements must be non-empty")
    }
    if (all(b <= 0)) stop("All measurements are zero or negative")
    m <- nrow(A); n <- ncol(A)
    W <- if (is.null(weights)) rep(1.0, m) else as.numeric(weights)
    W <- pmax(W, 0.0)
    pos <- b > 0
    if (!any(pos)) stop("All measurements are zero or negative")
    A_pos <- A[pos, , drop = FALSE]
    b_pos <- b[pos]
    W_pos <- W[pos]
    wb_inv <- W_pos / b_pos^2
    rhs_data <- as.numeric(t(A_pos) %*% (wb_inv * b_pos))
    M <- t(A_pos) %*% (wb_inv * A_pos)
    x <- pmax(x0, 1e-12)
    converged <- FALSE; iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        phi_prev <- pmax(x, 1e-12)
        inv_prev2 <- 1.0 / phi_prev^2
        lhs <- M + diag(inv_prev2, nrow = n, ncol = n)
        rhs <- rhs_data + inv_prev2 * phi_prev
        x_new <- tryCatch(as.numeric(qr.solve(lhs, rhs)),
                          error = function(e) {
            sv <- svd(lhs, nu = n, nv = n)
            s_inv <- ifelse(sv$d > 1e-10 * max(sv$d), 1 / sv$d, 0)
            as.numeric(sv$v %*% (s_inv * (t(sv$u) %*% rhs)))
        })
        x_new <- pmax(x_new, 0.0)
        rel_change <- max(abs(x_new - x) / pmax(x, 1e-12))
        x <- x_new
        if (rel_change < tolerance) { converged <- TRUE; break }
    }
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_rfsp_jul}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_rfsp_jul
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_rfsp_jul <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              max_iterations = 200L, tolerance = 1e-4,
                              weights = NULL,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_rfsp_jul,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         weights = weights),
        solve_kwargs = list(),
        method_name = "RFSP-JUL",
        extra_output = list(tolerance = tolerance),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
