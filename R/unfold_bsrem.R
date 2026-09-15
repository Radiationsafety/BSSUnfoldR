#' BSREM (block sequential regularised EM) unfolding
#'
#' R port of \code{bssunfold/core/unfold_bsrem.py}. Penalised EM with a
#' user-supplied relaxation sequence \code{alpha(n)}, ordered subsets, and a
#' floor clamp that prevents spectrum bins from being locked at zero.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param prior Character prior type (see \code{\link{prior_gradient}}).
#'   Default \code{"none"}.
#' @param beta Numeric prior weight. Default 1e-3.
#' @param prior_delta Numeric width/floor parameter. Default 1.0.
#' @param gamma Numeric edge-preservation parameter. Default 1.0.
#' @param max_iterations Positive integer; default 50.
#' @param n_subsets Positive integer; number of ordered subsets. Default 1.
#' @param tolerance Positive numeric; default 1e-6.
#' @param relaxation Numeric constant relaxation parameter (default 1.0) or a
#'   function \code{function(n)} returning the relaxation at iteration n.
#'   Default \code{NULL} = constant 1.
#' @param addition_after_iteration Numeric floor value clamped after every
#'   sub-iteration. Default 1e-4.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_bsrem(A, b, rep(1, 3), n_subsets = 3)
solve_bsrem <- function(A, b, x0, prior = "none", beta = 1e-3,
                        prior_delta = 1.0, gamma = 1.0,
                        max_iterations = 50L, n_subsets = 1L,
                        tolerance = 1e-6,
                        relaxation = NULL,
                        addition_after_iteration = 1e-4) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    m <- nrow(A); n <- ncol(A)
    if (n_subsets < 1L) stop("n_subsets must be >= 1")
    if (n_subsets > m) stop("n_subsets (", n_subsets,
                            ") must not exceed number of detectors (", m, ")")

    prior <- tolower(as.character(prior))
    if (!(prior %in% c("none", "quadratic", "logcosh", "relative_difference"))) {
        stop("Unknown prior '", prior, "'. Choose from 'none', 'quadratic', ",
             "'logcosh', 'relative_difference'.")
    }

    if (is.null(relaxation)) {
        relax_seq <- function(n) 1.0
    } else if (is.function(relaxation)) {
        relax_seq <- relaxation
    } else {
        relax_val <- as.numeric(relaxation)
        relax_seq <- function(n) relax_val
    }

    eps <- 1e-11
    groups <- (seq_len(m) - 1L) %% n_subsets
    idx_list <- split(seq_len(m), groups)
    idx_list <- idx_list[lengths(idx_list) > 0L]
    norm_all <- colSums(A)
    x <- pmax(x0, 0)
    converged <- FALSE
    iterations <- 0L

    for (it in seq_len(max_iterations)) {
        iterations <- it
        x_old <- x
        alpha <- relax_seq(it)
        for (idx in idx_list) {
            omega <- length(idx) / m
            A_sub <- A[idx, , drop = FALSE]
            b_sub <- b[idx]
            norm_sub <- colSums(A_sub)
            Ax <- as.numeric(A_sub %*% x) + eps
            ratio <- b_sub / Ax
            correction <- as.numeric(t(A_sub) %*% ratio)
            if (prior == "none") {
                grad <- rep(0.0, n)
            } else {
                grad <- omega * prior_gradient(x, prior, beta,
                                              prior_delta, gamma)
            }
            update <- correction - norm_sub - grad
            step <- ifelse(omega * norm_all > eps,
                           alpha / (omega * norm_all + eps), 0.0)
            x <- x + x * step * update
            x[x <= addition_after_iteration] <- addition_after_iteration
        }
        rel <- sqrt(sum((x - x_old)^2)) / (sqrt(sum(x_old^2)) + eps)
        if (rel < tolerance) { converged <- TRUE; break }
    }
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_bsrem}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_bsrem
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_bsrem <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, initial_spectrum = NULL,
                          prior = "none", beta = 1e-3,
                          prior_delta = 1.0, gamma = 1.0,
                          max_iterations = 50L, n_subsets = 1L,
                          tolerance = 1e-6, relaxation = NULL,
                          addition_after_iteration = 1e-4,
                          calculate_errors = FALSE,
                          noise_level = 0.01, n_montecarlo = 100L,
                          save_result = FALSE, random_state = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    x0_default[1L] <- 0.0
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_bsrem,
                                         prior = prior,
                                         beta = beta,
                                         prior_delta = prior_delta,
                                         gamma = gamma,
                                         max_iterations = max_iterations,
                                         n_subsets = n_subsets,
                                         tolerance = tolerance,
                                         relaxation = relaxation,
                                         addition_after_iteration = addition_after_iteration),
        solve_kwargs = list(),
        method_name = "BSREM",
        extra_output = list(
            prior = prior, beta = beta, n_subsets = as.integer(n_subsets)
        ),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
