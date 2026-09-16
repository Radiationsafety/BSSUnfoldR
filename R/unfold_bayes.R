#' Bayesian iterative unfolding (D'Agostini)
#'
#' R port of \code{bssunfold/core/unfold_bayes.py}. Pure-R implementation of
#' the D'Agostini algorithm. The response matrix is column-normalised so each
#' column sums to 1 (conditional probability \eqn{P(D_j | E_i)}), then the
#' result is rescaled to physical units via division by the column sums.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional prior spectrum (length n). Default \code{NULL} = uniform.
#' @param max_iterations Positive integer; default 4000.
#' @param tolerance Positive numeric; relative L2 convergence tolerance.
#'   Default 1e-3.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_bayes(A, b, rep(1, 3), max_iterations = 50)
solve_bayes <- function(A, b, x0 = NULL, max_iterations = 4000L,
                         tolerance = 1e-3) {
    v <- validate_system(A, b)
    A <- v$A; b <- v$b
    n <- ncol(A)
    column_sums <- colSums(A)
    column_sums_safe <- ifelse(column_sums > 0, column_sums, 1.0)
    P <- sweep(A, 2L, column_sums_safe, "/")
    zero_sens <- column_sums <= 0

    if (!is.null(x0) && sum(x0) > 0) {
        prior <- as.numeric(x0) / sum(x0)
    } else {
        prior <- rep(1.0, n) / n
    }
    total_counts <- sum(b)
    y <- total_counts * prior
    converged <- FALSE
    iterations <- 0L

    for (i in seq_len(max_iterations)) {
        iterations <- i
        f_norm <- as.numeric(P %*% y)
        f_norm_safe <- ifelse(f_norm > 0, f_norm, 1e-300)
        weight <- sweep(P, 1L, b / f_norm_safe, "*")
        y_new <- y * colSums(weight)
        if (any(zero_sens)) {
            y_new[zero_sens] <- prior[zero_sens] * total_counts
        }
        denom <- max(1.0, sqrt(sum(y^2)))
        if (sqrt(sum((y_new - y)^2)) / denom < tolerance) {
            y <- y_new
            converged <- TRUE
            break
        }
        y <- y_new
    }
    x <- y / column_sums_safe
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_bayes}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_bayes
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_bayes <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, initial_spectrum = NULL,
                          max_iterations = 4000L, tolerance = 1e-3,
                          calculate_errors = FALSE,
                          noise_level = 0.01, n_montecarlo = 100L,
                          save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_bayes,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "Bayes_D'Agostini",
        extra_output = list(),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
