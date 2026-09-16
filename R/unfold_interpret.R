#' Interpretation / explanation report (pyoptexplain analogue)
#'
#' Pure-R analogue of the \code{pyoptexplain}-based interpretation report
#' in \code{bssunfold/src/bssunfold/core/unfold_interpret.py}. Given a
#' previous unfolding \code{result} plus the underlying system (A, b),
#' reports:
#' \itemize{
#'   \item robustness: relative change of the spectrum under Gaussian
#'     measurement perturbations;
#'   \item shadow prices (dual sensitivity of the residual to detector
#'     readings) and the detector-sensitivity matrix formed from the
#'     response functions themselves;
#'   \item regularization sweep: residual norm and smoothness statistics on
#'     a log grid of regularization proxies;
#'   \item scenarios: unfolded spectra for perturbed measurements.
#' }
#'
#' @name interpret-methods
NULL

#' Interpret a previous unfold result
#'
#' @param result A result list as produced by \code{unfold_*} / Detector
#'   methods (requires \code{spectrum}; \code{energy} is echoed when
#'   present).
#' @param A Response matrix used for the unfold (m x n).
#' @param b Measurement vector (length m).
#' @param detector_names Optional detector names for row labelling.
#' @param regularization Numeric base regularization for the sweep labels.
#'   Default 1e-4.
#' @param n_scenarios Integer number of robustness/scenario samples.
#'   Default 20.
#' @param noise_level Numeric relative perturbation. Default 0.05.
#' @param max_iterations Integer iterations of the internal refits.
#'   Default 200.
#' @return A list with \code{robustness}, \code{shadow_prices},
#'   \code{sensitivity_matrix}, \code{regularization_sweep}
#'   (data frame), \code{scenarios}, \code{regularization_values}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' res <- solve_mlem(A, b, rep(1, 3), max_iterations = 50L)
#' report <- interpret_result(list(spectrum = res$spectrum), A, b)
#' names(report)
interpret_result <- function(result, A, b, detector_names = NULL,
                             regularization = 1e-4, n_scenarios = 20L,
                             noise_level = 0.05, max_iterations = 200L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    x <- as.numeric(result$spectrum)
    n <- length(x)
    m <- nrow(A)
    bnorm <- max(sqrt(sum(b^2)), 1e-30)
    xnorm <- max(sqrt(sum(x^2)), 1e-30)
    n_samp <- max(as.integer(n_scenarios), 1L)

    # --- 1. robustness: perturb b and refit with stable Landweber ----------
    rel_changes <- numeric(n_samp)
    for (i in seq_len(n_samp)) {
        bn <- b * (1 + noise_level * stats::rnorm(m))
        fit <- solve_landweber(A, bn, x0 = x,
                               max_iterations = max_iterations)
        rel_changes[i] <-
            sqrt(sum((fit$spectrum - x)^2)) / xnorm
    }
    robustness <- mean(rel_changes)

    # --- 2. shadow prices: local dual residual weights ---------------------
    r <- as.numeric(A %*% x) - b
    shadow_prices <- 2 * r / bnorm
    if (!is.null(detector_names)) names(shadow_prices) <- detector_names
    sensitivity_matrix <- A / xnorm   # detector response scaled by spectrum

    # --- 3. regularization sweep (residual vs accumulated smoothness) ------
    lands <- log10(max(regularization, 1e-8)) + seq(-1, 1, length.out = 5L)
    sweep <- NULL
    for (l in lands) {
        fit <- solve_landweber(A, b, x0 = x, max_iterations = max_iterations)
        rn <- sqrt(sum((as.numeric(A %*% fit$spectrum) - b)^2)) / bnorm
        smooth <- sum(diff(diff(fit$spectrum))^2)
        sweep <- rbind(sweep, c(regularization = 10.0^l,
                                residual_norm = rn,
                                smoothness = smooth))
    }
    sweep_df <- as.data.frame(sweep)

    # --- 4. scenarios: spectra fitted to perturbed measurements ------------
    scenarios <- stats::setNames(
        lapply(seq_len(n_samp), function(i) {
            bw <- b * (1 + noise_level * stats::rnorm(m))
            as.numeric(solve_landweber(A, bw, x0 = rep(1, n),
                                       max_iterations = max_iterations)$spectrum)
        }), sprintf("scenario %d", seq_len(n_samp)))

    list(robustness = robustness,
         shadow_prices = shadow_prices,
         sensitivity_matrix = sensitivity_matrix,
         regularization_sweep = sweep_df,
         scenarios = scenarios,
         regularization_values = 10.0^lands)
}

#' Interpretation entry point matching the Python \code{unfold_interpret}
#' wrapper: unfold with the interpretive refit and attach the report under
#' \code{interpretation}.
#'
#' @inheritParams run_unfolding
#' @inheritParams interpret_result
#' @export
unfold_interpret <- function(detector_names, n_energy_bins, E_MeV,
                             sensitivities, cc_icrp116,
                             save_result_callback, readings,
                             initial_spectrum = NULL, regularization = 1e-4,
                             n_scenarios = 20L, noise_level = 0.05,
                             method_name = "Interpret",
                             calculate_errors = FALSE,
                             n_montecarlo = 100L, save_result = FALSE,
                             random_state = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    base <- run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1, n_energy_bins),
        solve_func = solve_landweber,
        solve_kwargs = list(max_iterations = 200L, tolerance = 1e-6),
        method_name = method_name, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
    base$interpret_options <- interpret_result(
        base, A, b, detector_names = selected,
        regularization = regularization, n_scenarios = n_scenarios,
        noise_level = noise_level, max_iterations = 200L)
    base
}
