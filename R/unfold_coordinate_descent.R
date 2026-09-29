#' Coordinate-descent unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_coordinate_descent.py}.
#' Block-free coordinate descent for the non-negative least-squares problem
#' \deqn{\min_{x\ge0} \frac12\|Ax-b\|^2 + \lambda_1\|x\|_1 + \frac{\lambda_2}{2}\|x\|^2.}
#' Each coordinate is updated in closed form
#' \deqn{x_j \leftarrow \max(0, (a_j^T r + \|a_j\|^2 x_j - \lambda_1)/(\|a_j\|^2 + \lambda_2))}
#' where \eqn{r = b - Ax} is the running residual, so a full sweep costs
#' \eqn{O(mn)}{O(m n)} with no matrix products.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Initial spectrum (length n); projected onto the nonnegative orthant.
#' @param max_iterations Positive integer; maximum full sweeps; default 2000L.
#' @param tolerance Positive numeric relative change tolerance per sweep; default 1e-8.
#' @param l1_penalty Numeric L1 penalty weight; default 0.
#' @param l2_penalty Numeric ridge penalty weight; default 0.
#' @param selection Character; \code{"cyclic"} or \code{"random"} coordinate
#'   order; default \code{"cyclic"}.
#' @param random_state Optional integer seed for the random coordinate order.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
solve_coordinate_descent <- function(A, b, x0, max_iterations = 2000L,
                                     tolerance = 1e-8, l1_penalty = 0.0,
                                     l2_penalty = 0.0, selection = "cyclic",
                                     random_state = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)

    selection <- tolower(as.character(selection))
    if (!selection %in% c("cyclic", "random")) {
        stop("selection must be 'cyclic' or 'random'")
    }

    x <- pmax(as.numeric(x0), 0.0)
    l1_penalty <- max(as.numeric(l1_penalty), 0.0)
    l2_penalty <- max(as.numeric(l2_penalty), 0.0)

    col_sq <- colSums(A * A)                    # ||a_j||^2
    residual <- b - as.numeric(A %*% x)         # running residual r = b - A x

    if (!is.null(random_state)) set.seed(as.integer(random_state))

    converged <- FALSE
    iterations <- 0L
    for (k in seq_len(max_iterations)) {
        x_prev <- x
        order <- if (selection == "cyclic") seq_len(n) else sample.int(n)

        for (j in order) {
            c <- col_sq[j]
            if (c <= 0) next
            # partial correlation: a_j^T (r + a_j x_j) = a_j^T r + c x_j
            rho_j <- sum(A[, j] * residual) + c * x[j]
            x_new <- (rho_j - l1_penalty) / (c + l2_penalty)
            x_new <- max(x_new, 0.0)
            delta <- x_new - x[j]
            if (delta != 0.0) {
                residual <- residual - delta * A[, j]
                x[j] <- x_new
            }
        }

        iterations <- k
        rel_change <- sqrt(sum((x - x_prev)^2)) /
            max(sqrt(sum(x^2)), 1e-30)
        if (rel_change < tolerance) {
            converged <- TRUE
            break
        }
    }

    list(spectrum = as.numeric(x), iterations = iterations, converged = converged)
}

#' @rdname solve_coordinate_descent
#' @inheritParams run_unfolding
#' @param initial_spectrum Optional initial spectrum guess.
#' @param selection Character; \code{"cyclic"} or \code{"random"}.
#' @param calculate_errors Logical; run Monte-Carlo uncertainty.
#' @param noise_level Numeric; MC noise level.
#' @param n_montecarlo Integer; MC samples.
#' @param save_result Logical; call the save callback.
#' @param random_state Optional integer seed (also used by 'random' selection).
#' @param max_neutron_energy Optional energy cutoff in MeV.
#' @return A result list produced by \code{\link{run_unfolding}}.
#' @export
unfold_coordinate_descent <- function(detector_names, n_energy_bins, E_MeV,
                                      sensitivities, cc_icrp116,
                                      save_result_callback, readings,
                                      initial_spectrum = NULL,
                                      max_iterations = 2000L, tolerance = 1e-8,
                                      l1_penalty = 0.0, l2_penalty = 0.0,
                                      selection = "cyclic",
                                      calculate_errors = FALSE,
                                      noise_level = 0.01, n_montecarlo = 100L,
                                      save_result = FALSE, random_state = NULL,
                                      max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_coordinate_descent,
                                        max_iterations = max_iterations,
                                        tolerance = tolerance,
                                        l1_penalty = l1_penalty,
                                        l2_penalty = l2_penalty,
                                        selection = selection,
                                        random_state = random_state),
        solve_kwargs = list(),
        method_name = "Coordinate Descent",
        extra_output = list(l1_penalty = l1_penalty,
                            l2_penalty = l2_penalty,
                            coordinate_selection = selection),
        calculate_errors = calculate_errors,
        noise_level = noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
