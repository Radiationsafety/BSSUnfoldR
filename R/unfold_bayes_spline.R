#' Bayesian iterative unfolding with spline regularization
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_bayes_spline_regularization.py}.
#' D'Agostini iterative Bayesian unfolding with a smoothing spline applied
#' to the physical spectrum between iterations. The smoother is applied in
#' log10-space to handle the large dynamic range of physical spectra and to
#' prevent edge blow-up in low-sensitivity bins.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional prior spectrum. Default \code{NULL} = uniform.
#' @param max_iterations Positive integer; default 4000.
#' @param tolerance Positive numeric; default 1e-3.
#' @param spline_degree Integer; spline degree (1-5). Default 3.
#' @param spline_smooth Numeric; smoothing parameter passed to
#'   \code{\link[stats]{smooth.spline}}. Default 1e-2.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_bayes_spline(A, b, rep(1, 3), max_iterations = 50)
solve_bayes_spline <- function(A, b, x0 = NULL, max_iterations = 4000L,
                                  tolerance = 1e-3, spline_degree = 3L,
                                  spline_smooth = 1e-2) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n_energy <- ncol(A)
    column_sums <- colSums(A)
    column_sums_safe <- ifelse(column_sums > 0, column_sums, 1.0)
    P <- sweep(A, 2L, column_sums_safe, "/")
    zero_sens <- column_sums <= 0

    if (!is.null(x0) && sum(x0) > 0) {
        prior <- as.numeric(x0) / sum(x0)
    } else {
        prior <- rep(1.0, n_energy) / n_energy
    }
    total_counts <- sum(b)
    y <- total_counts * prior
    x_indices <- seq_len(n_energy)
    x_smooth <- rep(0.0, n_energy)
    converged <- FALSE; iterations <- 0L

    for (i in seq_len(max_iterations)) {
        iterations <- i
        y_old <- y
        f_norm <- as.numeric(P %*% y)
        f_norm_safe <- ifelse(f_norm > 0, f_norm, 1e-300)
        weight <- sweep(P, 1L, b / f_norm_safe, "*")
        y <- y_old * colSums(weight)
        if (any(zero_sens)) {
            y[zero_sens] <- prior[zero_sens] * total_counts
        }
        x <- y / column_sums_safe
        # Spline in log10-space
        if (n_energy > spline_degree + 1L) {
            log_x <- log10(pmax(x, 1e-300))
            sp <- tryCatch(
                stats::smooth.spline(x_indices, log_x,
                                     df = spline_degree + 2L,
                                     spar = spline_smooth),
                error = function(e) NULL
            )
            log_x_smooth <- if (is.null(sp)) log_x else as.numeric(predict(sp, x_indices)$y)
            x_smooth <- 10.0^log_x_smooth
        } else {
            x_smooth <- x
        }
        x_smooth <- pmax(x_smooth, 0)
        y <- x_smooth * column_sums_safe
        denom <- max(1.0, sqrt(sum(y_old^2)))
        if (sqrt(sum((y - y_old)^2)) / denom < tolerance) {
            converged <- TRUE; break
        }
    }
    list(spectrum = as.numeric(x_smooth), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_bayes_spline}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_bayes_spline
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_bayes_spline_regularization <- function(detector_names, n_energy_bins,
                                                 E_MeV, sensitivities,
                                                 cc_icrp116, save_result_callback,
                                                 readings,
                                                 initial_spectrum = NULL,
                                                 max_iterations = 4000L,
                                                 tolerance = 1e-3,
                                                 spline_degree = 3L,
                                                 spline_smooth = 1e-2,
                                                 calculate_errors = FALSE,
                                                 noise_level = 0.01,
                                                 n_montecarlo = 100L,
                                                 save_result = FALSE,
                                                 random_state = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_bayes_spline,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         spline_degree = spline_degree,
                                         spline_smooth = spline_smooth),
        solve_kwargs = list(),
        method_name = "Bayes_Spline",
        extra_output = list(spline_degree = spline_degree,
                            spline_smooth = spline_smooth),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
