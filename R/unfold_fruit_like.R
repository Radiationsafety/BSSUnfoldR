#' FRUIT-like parametric unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_fruit_like.py}.
#' Implements a parametric unfolding method inspired by FRUIT (Fast Real-time
#' Unfolding of neutron spectra with Iterative parameterNization Technique,
#' Bedogni et al., NIM A 580 (2007)).
#'
#' The parametric model consists of:
#' \itemize{
#'   \item Maxwellian thermal component: \eqn{A_{th} \sqrt{E} \exp(-E/T_{th})}
#'   \item 1/E epithermal component: \eqn{A_{epi} / E}
#'   \item Evaporation fast component: \eqn{A_f \exp(-E/T_{ev})}
#' }
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid (length n).
#' @param log_steps Numeric log-energy bin widths (length n).
#' @param initial_params Optional named list with \code{A_th}, \code{T_th},
#'   \code{A_epi}, \code{A_f}, \code{T_ev}.
#' @param method Character; optimizer method passed to \code{\link[stats]{optim}}.
#'   Default \code{"L-BFGS-B"}.
#' @return A list \code{list(spectrum, success, nfev, converged)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_fruit_like(A, b, E, compute_log_steps(E))
solve_fruit_like <- function(A, b, E, log_steps,
                                initial_params = NULL, method = "L-BFGS-B") {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E)
    if (is.null(initial_params)) {
        p0 <- c(A_th = 1e-6, T_th = 0.025e-6, A_epi = 1e-6,
                 A_f = 1e-6, T_ev = 2.0)
    } else {
        p0 <- c(A_th = initial_params$A_th %||% 1e-6,
                 T_th = initial_params$T_th %||% 0.025e-6,
                 A_epi = initial_params$A_epi %||% 1e-6,
                 A_f = initial_params$A_f %||% 1e-6,
                 T_ev = initial_params$T_ev %||% 2.0)
    }
    lower <- c(0, 1e-9, 0, 0, 0.1)
    upper <- c(1e-3, 1e-3, 1e-3, 1e-3, 20.0)

    .residuals <- function(p) {
        spectrum <- .parametric_model(E, p[1], p[2], p[3], p[4], p[5])
        spectrum_steps <- spectrum * log_steps
        as.numeric(A %*% spectrum_steps) - b
    }
    .objective <- function(p) sum(.residuals(p)^2)
    if (method == "L-BFGS-B") {
        result <- stats::optim(p0, .objective, method = "L-BFGS-B",
                                lower = lower, upper = upper,
                                control = list(maxit = 1000))
    } else {
        result <- stats::optim(p0, .objective, method = method,
                                control = list(maxit = 1000))
    }
    p_opt <- result$par
    spectrum <- .parametric_model(E, p_opt[1], p_opt[2], p_opt[3], p_opt[4], p_opt[5])
    spectrum <- spectrum * log_steps
    list(spectrum = pmax(as.numeric(spectrum), 0),
         success = (result$convergence == 0),
         nfev = as.integer(result$counts[1L]),
         converged = (result$convergence == 0),
         params = setNames(as.list(p_opt), names(p0)))
}

.parametric_model <- function(E, A_th, T_th, A_epi, A_f, T_ev,
                                epi_max = 0.1) {
    spectrum <- numeric(length(E))
    thermal <- E < 0.4e-6
    epithermal <- (E >= 0.4e-6) & (E < epi_max)
    fast <- E >= epi_max
    spectrum[thermal] <- spectrum[thermal] +
        (A_th * sqrt(E[thermal]) * exp(-E[thermal] / T_th))
    spectrum[epithermal] <- spectrum[epithermal] +
        (A_epi / (E[epithermal] + 1e-15))
    spectrum[fast] <- spectrum[fast] +
        (A_f * exp(-E[fast] / T_ev))
    spectrum
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Wrapper around \code{\link{solve_fruit_like}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_fruit_like
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_fruit_like <- function(detector_names, n_energy_bins, E_MeV,
                                 sensitivities, cc_icrp116, save_result_callback,
                                 readings, initial_spectrum = NULL,
                                 initial_params = NULL, method = "L-BFGS-B",
                                 calculate_errors = FALSE,
                                 noise_level = 0.01, n_montecarlo = 100L,
                                 save_result = FALSE, random_state = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    log_steps <- compute_log_steps(E_MeV) * log(10.0)
    x0_default <- rep(mean(b) / max(mean(rowSums(A)), 1e-30), n_energy_bins)
    solver <- function(A, b, x0 = NULL, ...) {
        res <- solve_fruit_like(A, b, E_MeV, log_steps,
                                initial_params = initial_params, method = method)
        list(spectrum = res$spectrum, iterations = res$nfev,
             converged = res$success)
    }
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver,
        solve_kwargs = list(),
        method_name = "FRUIT-like",
        extra_output = list(initial_params = initial_params, method = method),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
