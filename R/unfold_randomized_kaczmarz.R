#' Randomized Kaczmarz unfolding
#'
#' R port of \code{bssunfold/core/unfold_randomized_kaczmarz.py}. The
#' randomized variant selects rows probabilistically with probability
#' proportional to their squared norms (Strohmer & Vershynin, 2009), achieving
#' faster convergence than the deterministic cyclic variant on ill-conditioned
#' systems.
#'
#' Row selection is driven by the numpy-legacy \code{RandomState} emulation
#' shared with \code{\link{solve_eki}}, not by R's own RNG: Python draws one
#' \code{random_sample()} per iteration and maps it through
#' \code{np.searchsorted(cum_dist, u)} (the exact behaviour of
#' \code{RandomState.choice(m, p = probabilities)}), so
#' \code{random_state = k} reproduces the Python row sequence for the same
#' \code{k}. The cumulative distribution is built the same way as numpy builds
#' it — a left-to-right running sum of \code{row_norms_sq / sum(row_norms_sq)},
#' the \code{< 1} renormalisation loop and the \code{1e-10} bump on the final
#' element — and the generator is private, so R's global RNG state is neither
#' read nor reseeded.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial guess (length n).
#' @param max_iterations Positive integer; default 1000.
#' @param omega Numeric relaxation in (0, 2]; default 1.0.
#' @param tolerance Positive numeric; default 1e-6.
#' @param random_state Optional integer seed for reproducibility.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_randomized_kaczmarz(A, b, rep(0, 3), max_iterations = 30,
#'                                random_state = 7)
#' length(r$spectrum)

solve_randomized_kaczmarz <- function(A, b, x0, max_iterations = 1000L,
                                      omega = 1.0, tolerance = 1e-6,
                                      random_state = NULL) {
    v <- validate_system(A, b, x0 = x0,
                        max_iterations = max_iterations,
                        tolerance = tolerance)
    A <- v$A; b <- v$b; x0 <- v$x0
    m <- nrow(A); n <- ncol(A)
    rng <- .np_random_state(random_state)
    x <- as.numeric(x0)
    row_norms_sq <- as.numeric(rowSums(A * A))
    total_norm_sq <- sum(row_norms_sq)
    if (total_norm_sq == 0) {
        return(list(spectrum = x, iterations = 0L, converged = TRUE))
    }
    probabilities <- row_norms_sq / total_norm_sq

    # ---- RandomState.choice(m, p = probabilities) --------------------------
    cum_dist <- numeric(m)
    run <- 0
    for (i in seq_len(m)) {
        run <- run + probabilities[i]
        cum_dist[i] <- run
    }
    if (cum_dist[m] < 1) {
        # numpy pads the upper half of the strata by the (tiny) deficit
        half <- m %/% 2L
        if (half > 0L) {
            missing_value <- 1 - cum_dist[m]
            to_add <- missing_value / half
            for (i in seq_len(half)) {
                idx <- i + half
                cum_dist[idx] <- min(cum_dist[idx] + to_add, 1)
            }
        }
    }
    cum_dist[m] <- cum_dist[m] + 1e-10

    converged <- FALSE
    iterations <- 0L
    x_old <- x
    for (k in seq_len(as.integer(max_iterations))) {
        # searchsorted(cum_dist, u, side = "left") -> number of entries < u
        i <- sum(cum_dist < rng$uniform()) + 1L
        if (row_norms_sq[i] > 0) {
            Ai <- A[i, ]
            update <- (b[i] - sum(Ai * x)) / row_norms_sq[i]
            x <- x + omega * update * Ai
            x[x < 0] <- 0
        }
        if (k %% m == 0L) {
            if (sqrt(sum((x - x_old)^2)) < tolerance) {
                converged <- TRUE
                iterations <- k
                break
            }
            x_old <- x
        }
    }
    if (!converged) iterations <- as.integer(max_iterations)
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_randomized_kaczmarz}} for the unified
#' workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_randomized_kaczmarz
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_randomized_kaczmarz <- function(detector_names, n_energy_bins, E_MeV,
                                         sensitivities, cc_icrp116,
                                         save_result_callback, readings,
                                         initial_spectrum = NULL,
                                         max_iterations = 1000L,
                                         omega = 1.0, tolerance = 1e-6,
                                         calculate_errors = FALSE,
                                         noise_level = 0.01,
                                         n_montecarlo = 100L,
                                         save_result = FALSE,
                                         random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_randomized_kaczmarz,
                                         max_iterations = max_iterations,
                                         omega = omega,
                                         tolerance = tolerance,
                                         random_state = random_state),
        solve_kwargs = list(),
        method_name = "Randomized Kaczmarz",
        extra_output = list(tolerance = tolerance, omega = omega,
                            random_state = random_state),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
