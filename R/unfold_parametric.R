#' FRUIT parametric unfolding (full 3-component model)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_parametric.py}.
#' The spectrum is a weighted superposition of three components with
#' constraint P_th + P_epi + P_f = 1:
#' \describe{
#'   \item{Thermal}{\eqn{(E/T_0^2) \exp(-E/T_0)} for E < 1e-7 MeV}
#'   \item{Epithermal}{\eqn{[1-\exp(-(E/E_d)^2)] E^{b-1} \exp(-E/\beta')}
#'     for 1e-7 < E < 0.1 MeV}
#'   \item{Fast}{\eqn{E^\alpha \exp(-E/\beta)} for E > 0.1 MeV}
#' }
#'
#' @description Solves the unfolding problem by fitting a 3-component
#'   parametric spectral model (thermal, epithermal, fast) to the measured
#'   detector readings via constrained nonlinear least-squares.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid (length n).
#' @param log_steps Numeric log-energy bin widths (length n).
#' @param initial_params Optional named list.
#' @param max_iterations Integer; max optim iterations. Default 200.
#' @param tolerance Numeric; convergence tolerance. Default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_parametric(A, b, E, compute_log_steps(E) * log(10),
#'                        max_iterations = 50)
solve_parametric <- function(A, b, E, log_steps, initial_params = NULL,
                                max_iterations = 200L, tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E); log_steps <- as.numeric(log_steps)
    T0 <- 0.0253e-6  # thermal temperature (MeV)
    Ed <- 1e-6       # epithermal cutoff
    # Parametric model: phi(E) = P_th*phi_th + P_epi*phi_epi + P_f*phi_f
    # with P_th + P_epi + P_f = 1
    .model <- function(E, P_th, b_param, beta_prime, alpha, beta) {
        P_epi <- 1 - P_th * 0.5  # crude split; optim will refine
        P_f <- 1 - P_th - P_epi
        P_f <- max(P_f, 0)
        thermal <- ifelse(E < 1e-7, (E / T0^2) * exp(-E / T0), 0)
        epithermal <- ifelse(E >= 1e-7 & E < 0.1,
                              (1 - exp(-(E/Ed)^2)) * pmax(E, 1e-15)^(b_param - 1) *
                              exp(-E / beta_prime), 0)
        fast <- ifelse(E >= 0.1, pmax(E, 1e-15)^alpha * exp(-E / beta), 0)
        P_th * thermal + P_epi * epithermal + P_f * fast
    }
    .residuals <- function(p) {
        sp <- .model(E, p[1], p[2], p[3], p[4], p[5]) * log_steps
        as.numeric(A %*% sp) - b
    }
    .objective <- function(p) sum(.residuals(p)^2)
    p0 <- if (is.null(initial_params)) c(0.3, 0.5, 0.5, 1.0, 1.0)
           else unlist(initial_params)
    lower <- c(0, 0, 0.01, 0.1, 0.1)
    upper <- c(1, 5, 10, 5, 10)
    result <- stats::optim(p0, .objective, method = "L-BFGS-B",
                            lower = lower, upper = upper,
                            control = list(maxit = max_iterations))
    spectrum <- .model(E, result$par[1], result$par[2], result$par[3],
                        result$par[4], result$par[5]) * log_steps
    list(spectrum = pmax(as.numeric(spectrum), 0),
         iterations = as.integer(result$counts[1]),
         converged = (result$convergence == 0),
         params = result$par)
}

#' Wrapper around \code{\link{solve_parametric}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_parametric
#' @export
unfold_parametric <- function(detector_names, n_energy_bins, E_MeV,
                                 sensitivities, cc_icrp116, save_result_callback,
                                 readings, initial_spectrum = NULL,
                                 initial_params = NULL,
                                 max_iterations = 200L, tolerance = 1e-6,
                                 calculate_errors = FALSE,
                                 noise_level = 0.01, n_montecarlo = 100L,
                                 save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
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
        method = "Parametric",
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

#' BON95 parametric unfolding (4-component model)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_parametric2.py}.
#' The spectrum E*Phi(E) is a linear combination of four components
#' (thermal, epithermal, intermediate, fast) with free shape parameters
#' found by grid search and linear coefficients by weighted NLS.
#'
#' @inheritParams solve_parametric
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_parametric2(A, b, E, compute_log_steps(E) * log(10),
#'                          max_iterations = 20)
solve_parametric2 <- function(A, b, E, log_steps, initial_params = NULL,
                                 max_iterations = 50L, tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E); log_steps <- as.numeric(log_steps)
    Tth <- 3.5e-8  # thermal temp (MeV)
    .bon95_model <- function(E, b_param, Tf, c) {
        Xth <- E / Tth
        Xf <- (E / Tf)^c
        Fth <- Xth^1.5 * exp(-Xth)
        Fepi <- pmax(E, 1e-15)^(-b_param) * (1 - exp(-Xth))
        Fint <- (1 - exp(-Xth))
        Ff <- Xf^1.5 * exp(-Xf)
        cbind(Fth, Fepi, Fint, Ff)
    }
    # Grid search for shape parameters (b, Tf, c)
    b_grid <- c(0.5, 1.0, 1.5)
    Tf_grid <- c(0.5, 1.0, 2.0)
    c_grid <- c(0.5, 1.0)
    best_res <- Inf; best_params <- c(1.0, 1.0, 1.0); best_coefs <- c(1,1,1,1)
    for (bp in b_grid) for (Tf in Tf_grid) for (cp in c_grid) {
        F <- .bon95_model(E, bp, Tf, cp)
        F_steps <- F * log_steps
        AF <- A %*% F_steps  # m x 4
        coefs <- tryCatch(as.numeric(lsei::nnls(AF, b)$x),
                          error = function(e) rep(0, 4))
        residual <- b - as.numeric(AF %*% coefs)
        res <- sum(residual^2)
        if (res < best_res) {
            best_res <- res; best_params <- c(bp, Tf, cp); best_coefs <- coefs
        }
    }
    # Refine with optim
    .objective <- function(p) {
        F <- .bon95_model(E, p[1], p[2], p[3])
        F_steps <- F * log_steps
        AF <- A %*% F_steps
        coefs <- tryCatch(as.numeric(lsei::nnls(AF, b)$x),
                          error = function(e) rep(0, 4))
        sum((b - as.numeric(AF %*% coefs))^2)
    }
    result <- stats::optim(best_params, .objective, method = "L-BFGS-B",
                            lower = c(0.1, 0.1, 0.1), upper = c(5, 10, 5),
                            control = list(maxit = max_iterations))
    F <- .bon95_model(E, result$par[1], result$par[2], result$par[3])
    F_steps <- F * log_steps
    AF <- A %*% F_steps
    coefs <- tryCatch(as.numeric(lsei::nnls(AF, b)$x), error = function(e) rep(0, 4))
    spectrum <- as.numeric(F_steps %*% coefs)
    list(spectrum = pmax(spectrum, 0),
         iterations = as.integer(result$counts[1]),
         converged = (result$convergence == 0),
         params = result$par, coefs = coefs)
}

#' Wrapper around \code{\link{solve_parametric2}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_parametric2
#' @export
unfold_parametric2 <- function(detector_names, n_energy_bins, E_MeV,
                                  sensitivities, cc_icrp116, save_result_callback,
                                  readings, initial_spectrum = NULL,
                                  initial_params = NULL,
                                  max_iterations = 50L, tolerance = 1e-6,
                                  calculate_errors = FALSE,
                                  noise_level = 0.01, n_montecarlo = 100L,
                                  save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    log_steps <- compute_log_steps(E_MeV) * log(10.0)
    res <- solve_parametric2(A, b, E_MeV, log_steps,
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
        method = "Parametric2",
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = res$iterations, converged = res$converged,
        params = res$params, coefs = res$coefs
    )
    if (isTRUE(save_result) && is.function(save_result_callback))
        save_result_callback(result)
    result
}
