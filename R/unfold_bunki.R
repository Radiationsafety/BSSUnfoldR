#' BUNKI (SPUNIT) unfolding
#'
#' Port of the BUNKI/SPUNIT algorithm from the BUMS2 package. Each iteration
#' applies the SPUNIT multiplicative update on the lethargy-weighted response
#' matrix followed by 3-point smoothing of the spectrum.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum guess (length n).
#' @param smoothing Numeric; 3-point smoothing weight. Default 0.1.
#' @param max_iterations Positive integer; default 1000.
#' @param tolerance Positive numeric; relative change tolerance. Default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_bunki(A, b, rep(1, 3), max_iterations = 50)
solve_bunki <- function(A, b, x0, smoothing = 0.1,
                        max_iterations = 1000L, tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    n <- ncol(A)
    if (any(b < 0)) stop("BUNKI requires strictly positive measurements")
    if (any(b == 0)) {
        keep <- b > 0
        A <- A[keep, , drop = FALSE]
        b <- b[keep]
        if (length(b) == 0L) {
            stop("BUNKI requires strictly positive measurements")
        }
    }
    x0_safe <- pmax(x0, 0.0)
    aleth <- A * matrix(x0_safe, nrow = nrow(A), ncol = n, byrow = TRUE)
    spl <- rep(1, n)
    bcc <- as.numeric(aleth %*% spl)
    inv_b <- ifelse(b > 0, 1.0 / b, 0.0)
    ss <- as.numeric(t(aleth) %*% inv_b)
    inv_ss <- ifelse(ss > 0, 1.0 / pmax(ss, 1e-37), 0.0)
    denom_s <- 1.0 + 2.0 * smoothing
    converged <- FALSE; iterations <- 0L

    for (k in seq_len(max_iterations)) {
        iterations <- k
        inv_bcc <- ifelse(bcc > 0, 1.0 / pmax(bcc, 1e-37), 0.0)
        spll <- spl * as.numeric(t(aleth) %*% inv_bcc) * inv_ss
        spll <- ifelse(spll < 1e-37, 0.0, spll)
        spll <- ifelse(spl <= 0.0, 0.0, spll)
        new_spl <- spll
        if (n > 2) {
            new_spl[2:(n - 1L)] <- (smoothing * spll[1:(n - 2L)] +
                                    spll[2:(n - 1L)] +
                                    smoothing * spll[3:n]) / denom_s
        }
        bcc <- as.numeric(aleth %*% new_spl)
        rel <- sqrt(sum((new_spl - spl)^2)) / (sqrt(sum(spl^2)) + 1e-12)
        spl <- new_spl
        if (rel < tolerance) { converged <- TRUE; break }
    }
    spectrum <- spl * x0_safe
    list(spectrum = as.numeric(spectrum), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_bunki}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_bunki
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_bunki <- function(detector_names, n_energy_bins, E_MeV,
                         sensitivities, cc_icrp116, save_result_callback,
                         readings, initial_spectrum = NULL,
                         smoothing = 0.1, max_iterations = 1000L,
                         tolerance = 1e-6,
                         calculate_errors = FALSE,
                         noise_level = 0.01, n_montecarlo = 100L,
                         save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(1, n_energy_bins)
    x0_default[1L] <- 0.0
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_bunki,
                                         smoothing = smoothing,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "BUNKI",
        extra_output = list(smoothing = smoothing),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
