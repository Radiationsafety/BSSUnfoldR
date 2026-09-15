#' SAND-II unfolding
#'
#' Port of the SAND-II multiplicative ratio algorithm. Each iteration corrects
#' every spectrum bin by the weighted geometric mean of the
#' measured-to-calculated count-rate ratios.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum guess (length n).
#' @param max_iterations Positive integer; default 50.
#' @param tolerance Positive numeric; default 1e-3.
#' @param chi_fac Integer; 1 = stop on chi-square, 0 = stop on relative change.
#'   Default 1.
#' @param relative_uncertainty Numeric; default 0.1.
#' @param sigma Optional explicit per-detector uncertainties. Default \code{NULL}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_sandii(A, b, rep(1, 3), max_iterations = 20)
solve_sandii <- function(A, b, x0, max_iterations = 50L, tolerance = 1e-3,
                         chi_fac = 1L, relative_uncertainty = 0.1,
                         sigma = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    x0 <- as.numeric(x0)
    x <- pmax(x0, 0)
    if (!is.null(sigma)) {
        sigma <- pmax(as.numeric(sigma), 1e-12)
    } else {
        sigma <- relative_uncertainty * pmax(b, 1e-12)
    }
    valid <- b > 0
    if (!any(valid)) stop("All measurements are zero or negative")
    A_v     <- A[valid, , drop = FALSE]
    b_v     <- b[valid]
    sigma_v <- sigma[valid]
    m_v     <- length(b_v)
    eps <- 1e-12

    converged <- FALSE; iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        E <- as.numeric(A_v %*% x)
        Esafe <- pmax(E, 1e-300)
        R <- b_v / Esafe
        W <- A_v *
            sweep(matrix(x, nrow = nrow(A_v), ncol = ncol(A_v), byrow = TRUE),
                  1L, 1 / Esafe, "*")
        denom <- colSums(W)
        denom_safe <- ifelse(denom <= eps, eps, denom)
        numer <- as.numeric(t(W) %*% log(R))
        x_new <- x * exp(numer / denom_safe)
        if (chi_fac == 1L) {
            chi2 <- sum(((b_v - as.numeric(A_v %*% x_new)) / sigma_v)^2)
            if (chi2 <= m_v) { x <- x_new; converged <- TRUE; break }
        } else {
            rel <- abs(x_new - x) / pmax(x, eps)
            if (max(rel) < tolerance) { x <- x_new; converged <- TRUE; break }
        }
        x <- x_new
    }
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_sandii}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_sandii
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_sandii <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, initial_spectrum = NULL,
                          max_iterations = 50L, tolerance = 1e-3,
                          chi_fac = 1L, relative_uncertainty = 0.1,
                          calculate_errors = FALSE,
                          noise_level = 0.01, n_montecarlo = 100L,
                          save_result = FALSE, random_state = NULL) {
    x0_default <- rep(1, n_energy_bins)
    x0_default[1L] <- 0.0
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_sandii,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         chi_fac = chi_fac,
                                         relative_uncertainty = relative_uncertainty),
        solve_kwargs = list(),
        method_name = "SAND-II",
        extra_output = list(chi_fac = as.integer(chi_fac),
                            relative_uncertainty = relative_uncertainty),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
