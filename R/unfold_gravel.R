#' GRAVEL unfolding
#'
#' Port of the GRAVEL iterative multiplicative algorithm used by the UMG
#' package. Each iteration updates the spectrum by
#' \eqn{x_j <- x_j * exp( sum_i W_{ij} * ln(b_i / E_i) / sum_i W_{ij} )}
#' where \eqn{W_{ij} = b_i * A_{ij} * x_j / E_i}{W_ij = b_i * A_ij * x_j / E_i}
#' and \eqn{E_i = sum_j A_{ij} x_j}{E_i = sum_j A_ij * x_j}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum guess (length n).
#' @param tolerance Positive numeric; default 1e-8.
#' @param max_iterations Positive integer; default 1000.
#' @param regularization Non-negative numeric; default 0.0.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_gravel(A, b, rep(1, 3), max_iterations = 50)
#' r$converged
solve_gravel <- function(A, b, x0, tolerance = 1e-8,
                         max_iterations = 1000L, regularization = 0.0) {
    v <- validate_system(A, b, x0 = x0, max_iterations = max_iterations,
                        tolerance = tolerance)
    A <- v$A; b <- v$b; x0 <- v$x0
    x <- x0
    valid <- b > 0
    if (!any(valid)) stop("All measurements are zero or negative")
    A_v  <- A[valid, , drop = FALSE]
    b_v  <- b[valid]
    eps  <- 1e-10
    J_prev <- 0.0; dJ_prev <- 1.0

    for (it in seq_len(max_iterations)) {
        computed <- as.numeric(A_v %*% x)
        csafe <- pmax(computed, eps)
        xsafe <- pmax(x, 0.0)
        # W_ij = b_i * A_ij * x_j / csafe_i
        W <- sweep(A_v, 1L, b_v / csafe, "*") *
            matrix(xsafe, nrow = nrow(A_v), ncol = ncol(A_v), byrow = TRUE)
        # Zero-out where csafe == 0 or xsafe == 0 (kept naturally because W == 0)
        log_ratio <- (b_v / csafe)
        log_ratio[!is.finite(log_ratio)] <- 0
        log_ratio <- log(b_v / csafe)
        log_ratio[b_v <= 0 | csafe <= eps | computed <= 0] <- 0
        numerator <- as.numeric(t(W) %*% log_ratio)
        denom     <- as.numeric(colSums(W))
        update_mask <- denom > 0
        if (any(update_mask)) {
            reg_term <- regularization * log(x[update_mask] + eps)
            upd <- exp((numerator[update_mask] - reg_term) / denom[update_mask])
            x[update_mask] <- x[update_mask] * upd
        }
        computed_final <- as.numeric(A_v %*% x)
        chi_sq <- sum((computed_final - b_v)^2 / pmax(b_v, eps))
        J <- chi_sq / sum(computed_final)
        dJ <- J_prev - J
        ddJ <- abs(dJ - dJ_prev)
        if (it > 1L && ddJ <= tolerance) {
            return(list(spectrum = as.numeric(x), iterations = it,
                        converged = TRUE))
        }
        J_prev <- J; dJ_prev <- dJ
    }
    list(spectrum = as.numeric(x), iterations = max_iterations,
         converged = FALSE)
}

#' Wrapper around \code{\link{solve_gravel}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_gravel
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_gravel <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, initial_spectrum = NULL,
                          tolerance = 1e-8, max_iterations = 1000L,
                          regularization = 0.0,
                          calculate_errors = FALSE,
                          noise_level = 0.01, n_montecarlo = 100L,
                          save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.5, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_gravel,
                                         tolerance = tolerance,
                                         max_iterations = max_iterations,
                                         regularization = regularization),
        solve_kwargs = list(),
        method_name = "GRAVEL",
        extra_output = list(tolerance = tolerance,
                            regularization = regularization),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
