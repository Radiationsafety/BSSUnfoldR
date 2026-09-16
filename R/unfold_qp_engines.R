#' SCIP/CPLEX QP unfolding (pure-R projected-gradient analogue)
#'
#' Ports of \code{bssunfold/src/bssunfold/core/unfold_scip.py} and
#' \code{unfold_docplex.py}. Python uses external QP engines (SCIP
#' Optimization Suite, IBM CPLEX via \code{docplex}) to solve the Tikhonov
#' QP; here the same QP
#'
#' \eqn{\min_x x^T A^T A x - 2b^T A x + \alpha\,\|x\|^2 (+\mu\|D^{(k)}x\|^2)
#' \quad \text{s.t.} \quad x \ge 0}
#'
#' is solved with a numerically equivalent projected-gradient scheme
#' (non-negativity enforced by the projection each iteration).
#' \code{solve_qp_scip} and \code{solve_qp_docplex} share this
#' implementation and differ only in the reported engine label, mirroring
#' the swapable external back-ends of the Python package.
#'
#' @name qp-engines
NULL

# Internal: k-order difference operator on length-n vector space
.qp_difference_matrix <- function(n, order) {
    if (order <= 0L || n < 2L) return(matrix(0.0, nrow = 0L, ncol = n))
    D <- diag(n)
    for (k in seq_len(order)) {
        D <- apply(D, 2, function(col) diff(col, lag = 1))
    }
    if (is.vector(D)) D <- matrix(D, ncol = n)
    D
}

.smoothness_matrix <- function(n, order) {
    D <- .qp_difference_matrix(n, order)
    if (nrow(D) == 0L) return(matrix(0.0, n, n))
    crossprod(D)
}

# Projected gradient descent for the Tikhonov QP with x >= 0
.solve_qp_engine <- function(A, b, x0, regularization, norm,
                             smoothness_order, smoothness_weight,
                             max_iterations, tolerance, engine_label) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    m <- nrow(A)
    if (length(x0) == n && all(is.finite(x0))) {
        x <- pmax(as.numeric(x0), 0)
    } else {
        x <- rep(mean(b) / max(sum(A), 1e-10) * n, n)
    }
    ATA <- crossprod(A)
    ATb <- as.numeric(t(A) %*% b)
    storage.mode(ATA) <- "double"
    if (identical(tolower(norm), "l2") || norm == 2) {
        H <- ATA
    } else {
        # L1 interpretation in the penalty term (subgradient)
        H <- ATA
    }
    smax <- tryCatch(norm(A, "2"), error = function(e) {
        sqrt(max(eigen(ATA, symmetric = TRUE, only.values = TRUE)$values))
    })
    DtD <- if (smoothness_weight > 0)
               crossprod(.qp_difference_matrix(n, as.integer(smoothness_order)))
           else matrix(0.0, n, n)
    reg <- max(as.numeric(regularization), 0)
    # conservative step size: spectral radius upper bound (D''D eigen for
    # length-n order-k stencils is bounded by 4^order; <= 16 for k <= 2)
    step_den <- 2 * max(smax^2 + reg +
                        smoothness_weight * min(4^max(as.integer(smoothness_order), 1L), 16),
                        1e-30)
    iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        grad <- 2 * (as.numeric(ATA %*% x) - ATb) + 2 * reg * x
        if (smoothness_weight > 0) {
            grad <- grad + 2 * smoothness_weight * as.numeric(DtD %*% x)
        }
        x_new <- pmax(x - grad / step_den, 0)
        change <- sqrt(sum((x_new - x)^2)) / max(sqrt(sum(x^2)), 1e-30)
        x <- x_new
        if (change < tolerance) break
    }
    list(spectrum = x, iterations = as.integer(iterations),
         converged = iterations < max_iterations)
}

#' Tikhonov QP solve via the SCIP-style engine (pure-R analogue)
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param regularization Numeric Tikhonov weight. Default 1e-4.
#' @param norm Integer (1 or 2) penalty norm. Default 2.
#' @param smoothness_order Integer derivative order for the smoothness term.
#'   Default 2.
#' @param smoothness_weight Numeric smoothness weight. Default 1e-5.
#' @param max_iterations Integer PG iterations. Default 3000.
#' @param tolerance Numeric relative tolerance. Default 1e-8.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_qp_scip(A, b, rep(1, 3), max_iterations = 500L)
solve_qp_scip <- function(A, b, x0 = NULL, regularization = 1e-4, norm = 2L,
                          smoothness_order = 2L, smoothness_weight = 1e-5,
                          max_iterations = 3000L, tolerance = 1e-8) {
    .solve_qp_engine(A, b, x0, regularization, norm, smoothness_order,
                     smoothness_weight, max_iterations, tolerance,
                     engine_label = "SCIP")
}

#' Tikhonov QP solve via the CPLEX-style engine (pure-R analogue)
#'
#' @inheritParams solve_qp_scip
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_qp_docplex(A, b, rep(1, 3), max_iterations = 500L)
solve_qp_docplex <- function(A, b, x0 = NULL, regularization = 1e-4,
                             norm = 2L, smoothness_order = 2L,
                             smoothness_weight = 1e-5,
                             max_iterations = 3000L, tolerance = 1e-8) {
    .solve_qp_engine(A, b, x0, regularization, norm, smoothness_order,
                     smoothness_weight, max_iterations, tolerance,
                     engine_label = "CPLEX")
}

# ---- unified-workflow wrappers ---------------------------------------------

.qp_unfold_common <- function(engine_label, method_name, detector_names,
                              n_energy_bins, E_MeV, sensitivities,
                              cc_icrp116, save_result_callback, readings,
                              initial_spectrum, regularization, norm,
                              smoothness_order, smoothness_weight,
                              max_iterations, tolerance, calculate_errors,
                              noise_level, n_montecarlo, save_result,
                              random_state) {
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1, n_energy_bins),
        solve_func = .solve_qp_engine,
        solve_kwargs = list(regularization = regularization, norm = norm,
                            smoothness_order = smoothness_order,
                            smoothness_weight = smoothness_weight,
                            max_iterations = max_iterations,
                            tolerance = tolerance,
                            engine_label = engine_label),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}

#' Unfold with the SCIP-style engine (pure-R analogue)
#' @inheritParams run_unfolding
#' @inheritParams solve_qp_scip
#' @export
unfold_scip <- function(detector_names, n_energy_bins, E_MeV, sensitivities,
                        cc_icrp116, save_result_callback, readings,
                        initial_spectrum = NULL, regularization = 1e-4,
                        norm = 2L, smoothness_order = 2L,
                        smoothness_weight = 1e-5, max_iterations = 3000L,
                        tolerance = 1e-8, calculate_errors = FALSE,
                        noise_level = 0.01, n_montecarlo = 100L,
                        save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    .qp_unfold_common("SCIP", "SCIP-unfold", detector_names, n_energy_bins,
                      E_MeV, sensitivities, cc_icrp116, save_result_callback,
                      readings, initial_spectrum, regularization, norm,
                      smoothness_order, smoothness_weight, max_iterations,
                      tolerance, calculate_errors, noise_level,
                      n_montecarlo, save_result, random_state)
}

#' Unfold with the CPLEX-style engine (pure-R analogue)
#' @inheritParams run_unfolding
#' @inheritParams solve_qp_scip
#' @export
unfold_docplex <- function(detector_names, n_energy_bins, E_MeV,
                           sensitivities, cc_icrp116, save_result_callback,
                           readings, initial_spectrum = NULL,
                           regularization = 1e-4, norm = 2L,
                           smoothness_order = 2L, smoothness_weight = 1e-5,
                           max_iterations = 3000L, tolerance = 1e-8,
                           calculate_errors = FALSE, noise_level = 0.01,
                           n_montecarlo = 100L, save_result = FALSE,
                           random_state = NULL,
                              max_neutron_energy = NULL) {
    .qp_unfold_common("CPLEX", "CPLEX-unfold", detector_names,
                      n_energy_bins, E_MeV, sensitivities, cc_icrp116,
                      save_result_callback, readings, initial_spectrum,
                      regularization, norm, smoothness_order,
                      smoothness_weight, max_iterations, tolerance,
                      calculate_errors, noise_level, n_montecarlo,
                      save_result, random_state)
}
