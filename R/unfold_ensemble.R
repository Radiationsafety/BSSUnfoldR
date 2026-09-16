#' Ensemble / Cascade / Composite unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_ensemble.py},
#' \code{unfold_composite.py}, and \code{unfold_cascade.py} (combined into
#' one file for brevity).
#'
#' \describe{
#'   \item{Ensemble}{Runs multiple underlying solvers and averages their
#'     spectra (with optional weights).}
#'   \item{Cascade}{Runs a sequence of solvers, each starting from the
#'     result of the previous one.}
#'   \item{Composite}{Runs multiple solvers and picks the result with the
#'     smallest residual norm.}
#' }
#'
#' @name ensemble-methods
NULL

#' Solve by ensemble (weighted average of multiple solvers)
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param solvers List of solver functions, each with signature
#'   \code{solver(A, b, x0, ...)}.
#' @param weights Optional numeric vector of solver weights. Default
#'   \code{NULL} = uniform.
#' @param ... Extra arguments forwarded to each solver.
#' @return A list \code{list(spectrum, iterations = length(solvers),
#'   converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' solvers <- list(solve_mlem, solve_gravel)
#' r <- solve_ensemble(A, b, rep(1, 3), solvers,
#'                     max_iterations = 50L)
solve_ensemble <- function(A, b, x0, solvers, weights = NULL, ...) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    n <- ncol(A)
    if (length(solvers) == 0L) stop("solvers must be a non-empty list")
    if (is.null(weights)) {
        weights <- rep(1.0 / length(solvers), length(solvers))
    } else {
        weights <- weights / sum(weights)
    }
    spectra <- matrix(0.0, nrow = length(solvers), ncol = n)
    for (i in seq_along(solvers)) {
        res <- do.call(solvers[[i]], list(A = A, b = b, x0 = x0, ...))
        spec <- if (is.list(res) && !is.null(res$spectrum)) res$spectrum
                else as.numeric(res)
        spectra[i, ] <- pmax(spec, 0)
    }
    spectrum <- as.numeric(weights %*% spectra)
    list(spectrum = pmax(spectrum, 0.0),
         iterations = as.integer(length(solvers)),
         converged = TRUE,
         ensemble_spectra = spectra,
         ensemble_weights = weights)
}

#' Solve by cascade (sequential solvers, each starting from previous result)
#'
#' @inheritParams solve_ensemble
#' @param solver_kwargs_list Optional list of per-solver kwargs. Each entry
#'   is a list forwarded to the corresponding solver.
#' @param ... Extra arguments forwarded to each solver.
#' @return A list \code{list(spectrum, iterations = length(solvers),
#'   converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' solvers <- list(solve_mlem, solve_landweber)
#' r <- solve_cascade(A, b, rep(1, 3), solvers,
#'                     max_iterations = 50L)
solve_cascade <- function(A, b, x0, solvers, solver_kwargs_list = NULL, ...) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    n <- ncol(A)
    if (length(solvers) == 0L) stop("solvers must be a non-empty list")
    if (!is.null(solver_kwargs_list) &&
        length(solver_kwargs_list) != length(solvers)) {
        stop("solver_kwargs_list length must match solvers length")
    }
    x_cur <- x0
    cascade_spectra <- vector("list", length(solvers))
    for (i in seq_along(solvers)) {
        kwargs <- if (!is.null(solver_kwargs_list)) solver_kwargs_list[[i]]
                   else list()
        args <- c(list(A = A, b = b, x0 = x_cur), kwargs, list(...))
        res <- do.call(solvers[[i]], args)
        spec <- if (is.list(res) && !is.null(res$spectrum)) res$spectrum
                else as.numeric(res)
        x_cur <- pmax(spec, 0)
        cascade_spectra[[i]] <- x_cur
    }
    list(spectrum = as.numeric(x_cur),
         iterations = as.integer(length(solvers)),
         converged = TRUE,
         cascade_spectra = cascade_spectra)
}

#' Solve by composite (pick best residual norm from multiple solvers)
#'
#' @inheritParams solve_ensemble
#' @param ... Extra arguments forwarded to each solver.
#' @return A list \code{list(spectrum, iterations = which_best,
#'   converged = TRUE, best_solver_index = which_best)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' solvers <- list(solve_mlem, solve_gravel, solve_sandii)
#' r <- solve_composite(A, b, rep(1, 3), solvers, max_iterations = 50L)
solve_composite <- function(A, b, x0, solvers, ...) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    n <- ncol(A)
    if (length(solvers) == 0L) stop("solvers must be a non-empty list")
    spectra <- matrix(0.0, nrow = length(solvers), ncol = n)
    residuals <- numeric(length(solvers))
    for (i in seq_along(solvers)) {
        res <- do.call(solvers[[i]], list(A = A, b = b, x0 = x0, ...))
        spec <- if (is.list(res) && !is.null(res$spectrum)) res$spectrum
                else as.numeric(res)
        spec <- pmax(spec, 0)
        spectra[i, ] <- spec
        residuals[i] <- sqrt(sum((as.numeric(A %*% spec) - b)^2))
    }
    best <- which.min(residuals)
    list(spectrum = spectra[best, ],
         iterations = as.integer(best),
         converged = TRUE,
         best_solver_index = as.integer(best),
         composite_spectra = spectra,
         composite_residuals = residuals)
}

# --- Detector-facing wrappers ---

#' Wrapper around \code{\link{solve_ensemble}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_ensemble
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_ensemble <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              solvers = list(solve_mlem, solve_gravel),
                              weights = NULL, ...,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    x0 <- if (is.null(initial_spectrum)) rep(0.5, n_energy_bins)
           else as.numeric(initial_spectrum)
    res <- solve_ensemble(A, b, x0, solvers, weights = weights, ...)
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
        method = "Ensemble",
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = res$iterations,
        converged = res$converged,
        ensemble_spectra = res$ensemble_spectra,
        ensemble_weights = res$ensemble_weights
    )
    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }
    result
}

#' Wrapper around \code{\link{solve_cascade}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_cascade
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_cascade <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              solvers = list(solve_mlem, solve_gravel),
                              solver_kwargs_list = NULL, ...,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    x0 <- if (is.null(initial_spectrum)) rep(0.5, n_energy_bins)
           else as.numeric(initial_spectrum)
    res <- solve_cascade(A, b, x0, solvers,
                          solver_kwargs_list = solver_kwargs_list, ...)
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
        method = "Cascade",
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = res$iterations,
        converged = res$converged,
        cascade_spectra = res$cascade_spectra
    )
    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }
    result
}

#' Wrapper around \code{\link{solve_composite}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_composite
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_composite <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              solvers = list(solve_mlem, solve_gravel),
                              ...,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    x0 <- if (is.null(initial_spectrum)) rep(0.5, n_energy_bins)
           else as.numeric(initial_spectrum)
    res <- solve_composite(A, b, x0, solvers, ...)
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
        method = "Composite",
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = res$iterations,
        converged = res$converged,
        best_solver_index = res$best_solver_index,
        composite_spectra = res$composite_spectra,
        composite_residuals = res$composite_residuals
    )
    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }
    result
}
