#' Hybrid parametric unfolding (parametric model + iterative refinement)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_hybrid_parametric.py}.
#' Hybrid approach: first fit a parametric model (Maxwellian + 1/E + evaporation)
#' to the readings, then refine the residual with an iterative solver.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid (length n).
#' @param log_steps Numeric log-energy bin widths (length n).
#' @param initial_params Optional named list with \code{A_th}, \code{T_th},
#'   \code{A_epi}, \code{A_f}, \code{T_ev}.
#' @param refinement_iterations Positive integer; number of MLEM iterations for
#'   residual refinement. Default 100.
#' @param tolerance Numeric; convergence tolerance. Default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_hybrid_parametric(A, b, E, compute_log_steps(E) * log(10),
#'                              refinement_iterations = 30)
solve_hybrid_parametric <- function(A, b, E, log_steps,
                                       initial_params = NULL,
                                       refinement_iterations = 100L,
                                       tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E); log_steps <- as.numeric(log_steps)
    n <- ncol(A)
    # Stage 1: parametric fit
    res_para <- solve_fruit_like(A, b, E, log_steps,
                                    initial_params = initial_params)
    spectrum <- res_para$spectrum
    # Stage 2: MLEM refinement of residual
    computed <- as.numeric(A %*% spectrum)
    residual_b <- pmax(b - computed, 0)
    # MLEM on the residual: solve A x_res ~= residual_b, add to parametric
    x_res <- rep(0.0, n)
    AT <- t(A)
    eps <- 1e-10
    converged <- FALSE; iterations <- as.integer(res_para$nfev)
    for (i in seq_len(refinement_iterations)) {
        iterations <- iterations + 1L
        Ax <- as.numeric(A %*% x_res)
        Ax_safe <- pmax(Ax, eps)
        ratio <- residual_b / Ax_safe
        correction <- as.numeric(AT %*% ratio)
        x_new <- x_res * correction
        diff <- sqrt(sum((x_new - x_res)^2)) / (sqrt(sum(x_res^2)) + eps)
        x_res <- pmax(x_new, 0)
        if (diff < tolerance) { converged <- TRUE; break }
    }
    spectrum_final <- pmax(spectrum + x_res, 0)
    list(spectrum = spectrum_final,
         iterations = iterations,
         converged = converged,
         parametric_spectrum = spectrum,
         residual_spectrum = x_res,
         params = res_para$params)
}

#' Wrapper around \code{\link{solve_hybrid_parametric}} for the unified
#' workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_hybrid_parametric
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_hybrid_parametric <- function(detector_names, n_energy_bins, E_MeV,
                                         sensitivities, cc_icrp116,
                                         save_result_callback, readings,
                                         initial_spectrum = NULL,
                                         initial_params = NULL,
                                         refinement_iterations = 100L,
                                         tolerance = 1e-6,
                                         calculate_errors = FALSE,
                                         noise_level = 0.01,
                                         n_montecarlo = 100L,
                                         save_result = FALSE,
                                         random_state = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    log_steps <- compute_log_steps(E_MeV) * log(10.0)
    res <- solve_hybrid_parametric(A, b, E_MeV, log_steps,
                                      initial_params = initial_params,
                                      refinement_iterations = refinement_iterations,
                                      tolerance = tolerance)
    spectrum <- res$spectrum
    computed_readings <- as.numeric(A %*% spectrum)
    residual <- b - computed_readings
    result <- list(
        energy = E_MeV,
        spectrum = spectrum,
        spectrum_absolute = spectrum,
        effective_readings = stats::setNames(computed_readings, selected),
        residual = residual,
        residual_norm = sqrt(sum(residual^2)),
        method = "HybridParametric",
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = res$iterations,
        converged = res$converged,
        parametric_spectrum = res$parametric_spectrum,
        residual_spectrum = res$residual_spectrum,
        params = res$params
    )
    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }
    result
}
