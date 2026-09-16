#' Parametric unfolding engine variants (FRUIT-style)
#'
#' Thin variants of \code{\link{unfold_parametric}} matching the Python
#' wrappers \code{unfold_parametric_cvxpy} (SQP with a cvxpy QP backend),
#' \code{unfold_parametric_qpsolvers} (SQP through qpsolvers) and
#' \code{unfold_parametric_combined} (lmfit first pass + QP refinement).
#' The pure-R package routes all three through the same directional
#' divergence/kernel-fitting engine of \code{\link{solve_parametric}} and
#' differ only in the engine label reported in \code{method}, exactly like
#' the Python backends swap. Use the \code{engine} argument to choose the
#' label: \code{"cvxpy"} (SQP), \code{"qpsolvers"} (SQP) or
#' \code{"combined"} (lmfit + QP refinement).
#'
#' @name parametric-variants
NULL

#' Named-engine parametric unfolding
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_parametric
#' @param engine Character; \code{"cvxpy"} (default), \code{"qpsolvers"} or
#'   \code{"combined"}.
#' @export
unfold_parametric_engine <- function(detector_names, n_energy_bins, E_MeV,
                                     sensitivities, cc_icrp116,
                                     save_result_callback, readings,
                                     initial_spectrum = NULL,
                                     initial_params = NULL, engine = "cvxpy",
                                     max_iterations = 100L,
                                     tolerance = 1e-6,
                                     save_result = FALSE,
                                     random_state = NULL,
                              max_neutron_energy = NULL) {
    eng <- tolower(engine)
    method_name <- switch(eng,
        cvxpy = "Parametric (cvxpy SQP)",
        qpsolvers = "Parametric (qpsolvers SQP)",
        combined = "Parametric (combined lmfit + SQP)",
        stop("Unknown engine '", engine,
             "'. Available: cvxpy, qpsolvers, combined."))
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    log_steps <- compute_log_steps(E_MeV) * log(10.0)
    res <- solve_parametric(A, b, E_MeV, log_steps,
                            initial_params = initial_params,
                            max_iterations = max_iterations,
                            tolerance = tolerance)
    spectrum <- res$spectrum
    computed_readings <- as.numeric(A %*% spectrum)
    residual <- b - computed_readings
    result <- list(
        energy = E_MeV, spectrum = spectrum, spectrum_absolute = spectrum,
        effective_readings = stats::setNames(computed_readings, selected),
        residual = residual, residual_norm = sqrt(sum(residual^2)),
        method = method_name,
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = res$iterations, converged = res$converged,
        params = res$params
    )
    if (isTRUE(save_result) && is.function(save_result_callback))
        save_result_callback(result)
    result
}

#' Parametric unfolding via the cvxpy-SQP label (engine = "cvxpy")
#' @inheritParams run_unfolding
#' @inheritParams solve_parametric
#' @export
unfold_parametric_cvxpy <- function(detector_names, n_energy_bins, E_MeV,
                                    sensitivities, cc_icrp116,
                                    save_result_callback, readings,
                                    initial_spectrum = NULL,
                                    initial_params = NULL,
                                    max_iterations = 100L, tolerance = 1e-6,
                                    save_result = FALSE,
                                    random_state = NULL,
                              max_neutron_energy = NULL) {
    unfold_parametric_engine(
        detector_names, n_energy_bins, E_MeV, sensitivities, cc_icrp116,
        save_result_callback, readings,
        initial_spectrum = initial_spectrum,
        initial_params = initial_params, engine = "cvxpy",
        max_iterations = max_iterations, tolerance = tolerance,
        save_result = save_result,
        random_state = random_state)
}

#' Parametric unfolding via the qpsolvers-SQP label (engine =
#' "qpsolvers")
#' @inheritParams run_unfolding
#' @inheritParams solve_parametric
#' @export
unfold_parametric_qpsolvers <- function(detector_names, n_energy_bins, E_MeV,
                                        sensitivities, cc_icrp116,
                                        save_result_callback, readings,
                                        initial_spectrum = NULL,
                                        initial_params = NULL,
                                        max_iterations = 100L,
                                        tolerance = 1e-6,
                                        save_result = FALSE,
                                        random_state = NULL,
                              max_neutron_energy = NULL) {
    unfold_parametric_engine(
        detector_names, n_energy_bins, E_MeV, sensitivities, cc_icrp116,
        save_result_callback, readings,
        initial_spectrum = initial_spectrum,
        initial_params = initial_params, engine = "qpsolvers",
        max_iterations = max_iterations, tolerance = tolerance,
        save_result = save_result,
        random_state = random_state)
}

#' Parametric unfolding via the combined lmfit+SQP label (engine =
#' "combined")
#' @inheritParams run_unfolding
#' @inheritParams solve_parametric
#' @export
unfold_parametric_combined <- function(detector_names, n_energy_bins, E_MeV,
                                       sensitivities, cc_icrp116,
                                       save_result_callback, readings,
                                       initial_spectrum = NULL,
                                       initial_params = NULL,
                                       max_iterations = 100L,
                                       tolerance = 1e-6,
                                       save_result = FALSE,
                                       random_state = NULL,
                              max_neutron_energy = NULL) {
    unfold_parametric_engine(
        detector_names, n_energy_bins, E_MeV, sensitivities, cc_icrp116,
        save_result_callback, readings,
        initial_spectrum = initial_spectrum,
        initial_params = initial_params, engine = "combined",
        max_iterations = max_iterations, tolerance = tolerance,
        save_result = save_result,
        random_state = random_state)
}
