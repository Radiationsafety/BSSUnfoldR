#' EPIC Tikhonov regularization unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_epic.py}.
#' Implements the Equal Posterior Information Condition (EPIC) Tikhonov
#' regularization (Ortega-Culaciati et al., JGR Solid Earth, 2021).
#' Regularization weights are chosen so that posterior variances of model
#' parameters match user-supplied target variances.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Unused.
#' @param sigma_frac Numeric; fraction of max LS solution as target sigma. Default 0.1.
#' @param max_iterations Integer; max EPIC iterations. Default 30.
#' @param tolerance Numeric; convergence on beta weights. Default 1e-4.
#' @param nonneg Logical; enforce non-negativity via NNLS. Default TRUE.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_epic(A, b, NULL, max_iterations = 10)
solve_epic <- function(A, b, x0 = NULL, sigma_frac = 0.1,
                        max_iterations = 30L, tolerance = 1e-4,
                        nonneg = TRUE) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    AT_A <- crossprod(A)
    AT_b <- as.numeric(t(A) %*% b)
    # Initial LS solution
    x_ls <- tryCatch(as.numeric(qr.solve(A, b, tol = 1e-10)),
                     error = function(e) as.numeric(ginv(A) %*% b))
    scale <- max(abs(x_ls)); if (!is.finite(scale) || scale <= 0) scale <- max(abs(b))
    target_sigmas <- rep(sigma_frac * scale, n)
    # Initial beta (log-weights) = 0 (uniform regularization)
    betas <- rep(0.0, n)
    converged <- FALSE; iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        betas_old <- betas
        weights <- exp(betas)  # regularization weights
        # Solve (A'A + diag(weights)) x = A'b
        lhs <- AT_A + diag(weights, nrow = n)
        if (nonneg) {
            # Use NNLS on augmented system
            Aw <- rbind(A, diag(sqrt(weights), nrow = n))
            bw <- c(b, rep(0, n))
            x <- tryCatch(as.numeric(lsei::nnls(Aw, bw)$x),
                          error = function(e) as.numeric(qr.solve(lhs, AT_b)))
        } else {
            x <- tryCatch(as.numeric(qr.solve(lhs, AT_b)),
                          error = function(e) as.numeric(solve(lhs, AT_b)))
        }
        # Compute posterior variances: diag of (A'A + diag(weights))^{-1}
        cov_inv <- lhs
        posterior_vars <- tryCatch(diag(solve(cov_inv)),
                                    error = function(e) rep(1, n))
        # Update betas to match target sigmas
        ratios <- sqrt(pmax(target_sigmas^2 / pmax(posterior_vars, 1e-30), 1e-30))
        betas <- betas + log(pmax(ratios, 1e-10))
        # Damp the update
        betas <- 0.5 * betas + 0.5 * betas_old
        if (max(abs(betas - betas_old)) < tolerance) { converged <- TRUE; break }
    }
    list(spectrum = pmax(as.numeric(x), 0),
         iterations = iterations, converged = converged,
         betas = betas, posterior_vars = posterior_vars)
}

#' Wrapper around \code{\link{solve_epic}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_epic
#' @export
unfold_epic <- function(detector_names, n_energy_bins, E_MeV,
                           sensitivities, cc_icrp116, save_result_callback,
                           readings, initial_spectrum = NULL,
                           sigma_frac = 0.1, max_iterations = 30L,
                           tolerance = 1e-4, nonneg = TRUE,
                           calculate_errors = FALSE,
                           noise_level = 0.01, n_montecarlo = 100L,
                           save_result = FALSE, random_state = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_epic,
                                         sigma_frac = sigma_frac,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         nonneg = nonneg),
        solve_kwargs = list(),
        method_name = "EPIC",
        extra_output = list(sigma_frac = sigma_frac),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
