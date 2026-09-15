#' MAP-EM (penalised expectation maximisation) unfolding
#'
#' R port of \code{bssunfold/core/unfold_mapem.py} (one-step-late OSMAPOSL).
#' The update reads
#' \deqn{x^{n+1} = x^n * A^T ( b / (A x^n + eps) ) / ( A^T 1 + beta * grad V(x^n) )}
#' with \eqn{V} a nearest-neighbour prior over the energy axis. Available
#' priors: \code{quadratic}, \code{logcosh}, \code{relative_difference}. Set
#' \code{prior = "none"} to recover plain MLEM.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param prior Character: \code{"none"}, \code{"quadratic"},
#'   \code{"logcosh"} or \code{"relative_difference"}. Default
#'   \code{"quadratic"}.
#' @param beta Numeric prior weight. Default 1e-3.
#' @param prior_delta Numeric width/floor parameter. Default 1.0.
#' @param gamma Numeric edge-preservation parameter. Default 1.0.
#' @param max_iterations Positive integer; default 50.
#' @param tolerance Positive numeric; default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_mapem(A, b, rep(1, 3), prior = "quadratic")
solve_mapem <- function(A, b, x0, prior = "quadratic", beta = 1e-3,
                        prior_delta = 1.0, gamma = 1.0,
                        max_iterations = 50L, tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    prior <- tolower(as.character(prior))
    if (!(prior %in% c("none", "quadratic", "logcosh", "relative_difference"))) {
        stop("Unknown prior '", prior, "'. Choose from 'none', 'quadratic', ",
             "'logcosh', 'relative_difference'.")
    }
    eps <- 1e-11
    x <- pmax(x0, 0)
    norm <- colSums(A)
    converged <- FALSE
    iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        x_old <- x
        Ax <- as.numeric(A %*% x) + eps
        ratio <- b / Ax
        correction <- as.numeric(t(A) %*% ratio)
        if (prior == "none") {
            x <- pmax(x * correction / (norm + eps), 0.0)
        } else {
            grad <- prior_gradient(x, prior, beta, prior_delta, gamma)
            x <- pmax(x * correction / (norm + grad + eps), 0.0)
        }
        rel <- sqrt(sum((x - x_old)^2)) / (sqrt(sum(x_old^2)) + eps)
        if (rel < tolerance) { converged <- TRUE; break }
    }
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_mapem}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_mapem
#' @return A result list as produced by \code{\link{run_unfolding}}. If a
#'   non-\code{"none"} prior is selected, the list also carries a
#'   \code{prior_value} entry giving \eqn{V(x)} at the final spectrum.
#' @export
unfold_mapem <- function(detector_names, n_energy_bins, E_MeV,
                         sensitivities, cc_icrp116, save_result_callback,
                         readings, initial_spectrum = NULL,
                         prior = "quadratic", beta = 1e-3,
                         prior_delta = 1.0, gamma = 1.0,
                         max_iterations = 50L, tolerance = 1e-6,
                         calculate_errors = FALSE,
                         noise_level = 0.01, n_montecarlo = 100L,
                         save_result = FALSE, random_state = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    x0_default[1L] <- 0.0
    result <- run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_mapem,
                                         prior = prior,
                                         beta = beta,
                                         prior_delta = prior_delta,
                                         gamma = gamma,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "MAP-EM",
        extra_output = list(
            prior = prior, beta = beta,
            prior_delta = prior_delta, gamma = gamma
        ),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = FALSE
    )
    if (tolower(prior) != "none") {
        result$prior_value <- prior_value(result$spectrum, prior, beta,
                                          prior_delta, gamma)
    }
    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }
    result
}
