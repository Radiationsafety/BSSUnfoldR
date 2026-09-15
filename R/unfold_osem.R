#' OSEM (ordered-subset expectation maximisation) unfolding
#'
#' R port of \code{bssunfold/core/unfold_osem.py}. OSEM generalises MLEM by
#' updating the spectrum with one subset of detectors at a time, accelerating
#' convergence:
#' \deqn{x^{n+1} = x^n * A_m^T ( b_m / (A_m x^n + eps) ) / ( A_m^T 1 + eps)}
#' With \code{n_subsets=1} the update reduces to standard MLEM.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param max_iterations Positive integer; default 50.
#' @param n_subsets Positive integer; number of detector subsets. Default 1.
#' @param tolerance Positive numeric; default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_osem(A, b, rep(1, 3), n_subsets = 3)
solve_osem <- function(A, b, x0, max_iterations = 50L, n_subsets = 1L,
                       tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    m <- nrow(A); n <- ncol(A)
    if (n_subsets < 1L) stop("n_subsets must be >= 1")
    if (n_subsets > m) stop("n_subsets (", n_subsets,
                            ") must not exceed number of detectors (", m, ")")

    eps <- 1e-11
    # split seq_len(m) into n_subsets chunks
    # Split m rows into n_subsets groups using even modulo distribution
    # (base::cut() fails when n_subsets == m and is awkward for small inputs.)
    groups <- (seq_len(m) - 1L) %% n_subsets
    idx_list <- split(seq_len(m), groups)
    # remove empty chunks (happens when n_subsets > m, but we already checked)
    idx_list <- idx_list[lengths(idx_list) > 0L]

    x <- pmax(x0, 0)
    converged <- FALSE
    iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        x_old <- x
        for (idx in idx_list) {
            A_sub <- A[idx, , drop = FALSE]
            b_sub <- b[idx]
            norm <- colSums(A_sub)
            Ax <- as.numeric(A_sub %*% x) + eps
            ratio <- b_sub / Ax
            correction <- as.numeric(t(A_sub) %*% ratio)
            x <- pmax(x * correction / (norm + eps), 0.0)
        }
        rel <- sqrt(sum((x - x_old)^2)) / (sqrt(sum(x_old^2)) + eps)
        if (rel < tolerance) { converged <- TRUE; break }
    }
    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' Wrapper around \code{\link{solve_osem}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_osem
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_osem <- function(detector_names, n_energy_bins, E_MeV,
                         sensitivities, cc_icrp116, save_result_callback,
                         readings, initial_spectrum = NULL,
                         max_iterations = 50L, n_subsets = 1L,
                         tolerance = 1e-6,
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
        solve_func = make_solve_wrapper(solve_osem,
                                         max_iterations = max_iterations,
                                         n_subsets = n_subsets,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "OSEM",
        extra_output = list(n_subsets = as.integer(n_subsets)),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
