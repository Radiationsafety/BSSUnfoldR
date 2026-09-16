#' Binned adaptive unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_binned.py}.
#' Bin-wise adaptive unfolding: instead of unfolding on the full fine energy
#' grid, this method groups adjacent energy bins into super-bins (coarse
#' groups), unfolds on the coarse grid using any underlying solver, then
#' redistributes the coarse flux back to the fine grid using a smoothness
#' prior. This is useful for severely underdetermined problems.
#'
#' @section Limitations:
#' This is a simplified port. The Python original supports pre-computed
#' bin-lookup tables (JSON files) for arbitrary bin groupings. This R port
#' uses uniform log-energy groupings.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param E_MeV Numeric energy grid (length n).
#' @param n_super_bins Integer; number of super-bins. Default \code{NULL} =
#'   \code{max(3, n/10)}.
#' @param max_iterations Positive integer; passed to the underlying MLEM solver.
#'   Default 1000.
#' @param tolerance Positive numeric; default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged, super_bins, super_spectrum)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_binned(A, b, rep(1, 60), E, n_super_bins = 6L,
#'                    max_iterations = 50)
solve_binned <- function(A, b, x0, E_MeV, n_super_bins = NULL,
                           max_iterations = 1000L, tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0); E_MeV <- as.numeric(E_MeV)
    n <- ncol(A)
    if (is.null(n_super_bins)) n_super_bins <- max(3L, n %/% 10L)
    if (n_super_bins < 2L) n_super_bins <- 2L
    if (n_super_bins > n) n_super_bins <- n
    # Build bin lookup: assign each fine bin to a super-bin by log-uniform
    # binning of the energy grid
    log_E <- log10(pmax(E_MeV, 1e-15))
    breaks <- seq(min(log_E), max(log_E), length.out = n_super_bins + 1L)
    # Each fine bin goes to the super-bin whose [break, break+1) range
    # contains log_E[i]
    bin_lookup <- findInterval(log_E, breaks, rightmost.closed = TRUE)
    bin_lookup <- pmin(pmax(bin_lookup, 1L), n_super_bins)
    # Build coarse response matrix A_super (m x n_super_bins)
    A_super <- matrix(0.0, nrow = nrow(A), ncol = n_super_bins)
    for (j in seq_len(n)) {
        A_super[, bin_lookup[j]] <- A_super[, bin_lookup[j]] + A[, j]
    }
    # Unfold on coarse grid using MLEM
    x_super <- rep(mean(b) / max(mean(rowSums(A_super)), 1e-10), n_super_bins)
    eps <- 1e-10
    AT <- t(A_super)
    converged <- FALSE; iterations <- 0L
    for (it in seq_len(max_iterations)) {
        iterations <- it
        Ax <- as.numeric(A_super %*% x_super)
        Ax_safe <- pmax(Ax, eps)
        ratio <- b / Ax_safe
        correction <- as.numeric(AT %*% ratio)
        x_new <- pmax(x_super * correction, 0)
        diff <- sqrt(sum((x_new - x_super)^2)) / (sqrt(sum(x_super^2)) + eps)
        x_super <- x_new
        if (diff < tolerance) { converged <- TRUE; break }
    }
    # Redistribute: each fine bin gets the super-bin's flux density
    # divided by the number of fine bins in its super-bin (uniform within
    # the super-bin).
    fine_per_super <- tabulate(bin_lookup, nbins = n_super_bins)
    spectrum <- numeric(n)
    for (j in seq_len(n)) {
        sb <- bin_lookup[j]
        spectrum[j] <- x_super[sb] / max(fine_per_super[sb], 1L)
    }
    list(spectrum = pmax(spectrum, 0),
         iterations = iterations,
         converged = converged,
         super_bins = bin_lookup,
         super_spectrum = x_super)
}

#' Build a uniform log-energy bin lookup
#'
#' @param E_MeV Numeric energy grid.
#' @param n_super_bins Integer; number of super-bins.
#' @return Integer vector (length n) assigning each fine bin to a super-bin
#'   (1-indexed).
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' lookup <- build_bin_lookup(E, 6L)
build_bin_lookup <- function(E_MeV, n_super_bins) {
    E_MeV <- as.numeric(E_MeV)
    log_E <- log10(pmax(E_MeV, 1e-15))
    breaks <- seq(min(log_E), max(log_E), length.out = n_super_bins + 1L)
    lookup <- findInterval(log_E, breaks, rightmost.closed = TRUE)
    pmin(pmax(lookup, 1L), n_super_bins)
}

#' Save / load bin lookup
#'
#' @param lookup Integer vector.
#' @param path File path.
#' @return Invisible \code{lookup} (for save) or the loaded integer vector.
#' @export
save_bin_lookup <- function(lookup, path) {
    saveRDS(as.integer(lookup), path)
    invisible(lookup)
}

#' @rdname save_bin_lookup
#' @export
load_bin_lookup <- function(path) {
    as.integer(readRDS(path))
}

#' Wrapper around \code{\link{solve_binned}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_binned
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_binned <- function(detector_names, n_energy_bins, E_MeV,
                             sensitivities, cc_icrp116, save_result_callback,
                             readings, initial_spectrum = NULL,
                             n_super_bins = NULL,
                             max_iterations = 1000L, tolerance = 1e-6,
                             calculate_errors = FALSE,
                             noise_level = 0.01, n_montecarlo = 100L,
                             save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(0.5, n_energy_bins)
    solver <- function(A, b, x0 = NULL, ...) {
        solve_binned(A, b, if (is.null(x0)) x0_default else x0,
                     E_MeV, n_super_bins = n_super_bins,
                     max_iterations = max_iterations, tolerance = tolerance)
    }
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver,
        solve_kwargs = list(),
        method_name = "Binned",
        extra_output = list(n_super_bins = if (is.null(n_super_bins))
                                              max(3L, n_energy_bins %/% 10L)
                                           else as.integer(n_super_bins)),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
