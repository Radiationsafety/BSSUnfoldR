#' LOUHI unfolding method (Routti & Sandberg 1980)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_louhi.py}.
#' Port of the LOUHI78 general purpose unfolding program (J. T. Routti and
#' V. Sandberg, "General purpose unfolding program LOUHI78 with linear and
#' nonlinear regressions", Computer Physics Communications 21 (1980)
#' 119-135, doi:10.1016/0010-4655(80)90021-4).
#'
#' LOUHI formulates Bonner-sphere unfolding as a constrained weighted
#' least-squares problem with generalized smoothing:
#' \deqn{\chi^2(\phi) = \sum_i [(b_i - (A \phi)_i) / \sigma_i]^2
#'   + \lambda^2 \sum_k [L (\phi - \phi_0)]_k^2}{
#'   chi2(phi) = sum_i [(b_i - (A phi)_i)/sigma_i]^2 + lambda^2 sum_k
#'   [L (phi - phi0)]_k^2}
#' where \code{phi0} is the default (a-priori) spectrum, \code{sigma_i} the
#' per-detector measurement uncertainties and \code{L} the smoothing
#' operator.  The quadratic form is minimized by quadratic programming with
#' non-negativity constraints (LOUHI's LSI step) using Hildreth's iterative
#' coordinate algorithm, and the smoothing weight \code{lambda} can either
#' be fixed (linear mode) or adjusted automatically by a nonlinear
#' regression on the total chi-square so that the data misfit reaches its
#' expected value (nonlinear mode, LOUHI's smoothing-parameter search).
#'
#' @param n Integer; number of energy bins.
#' @param smooth_order Integer; order of the smoothing functional: \code{0}
#'   shrinks the solution toward the default spectrum (identity operator),
#'   \code{1} penalizes first differences of the deviation from the default
#'   spectrum and \code{2} penalizes second differences.  Default 1.
#' @return The \code{(n, n)} (order 0), \code{(n-1, n)} (order 1) or
#'   \code{(n-2, n)} (order 2) smoothing matrix \code{L}.
#' @keywords internal
#' @examples
#' \dontrun{
#' louhi_smoothing_matrix(4L, 1L)
#' }
louhi_smoothing_matrix <- function(n, smooth_order = 1L) {
    n <- as.integer(n)
    smooth_order <- as.integer(smooth_order)
    if (!smooth_order %in% c(0L, 1L, 2L)) {
        stop("smooth_order must be one of 0, 1, 2, got ", smooth_order)
    }
    if (n < 1L) stop("n must be positive, got ", n)
    identity <- diag(n)
    if (smooth_order == 0L || n == 1L) return(identity)
    first <- matrix(0, nrow = n - 1L, ncol = n)
    idx <- seq_len(n - 1L)
    first[cbind(idx, idx)] <- -1.0
    first[cbind(idx, idx + 1L)] <- 1.0
    if (smooth_order == 1L || n == 2L) return(first)
    second <- matrix(0, nrow = n - 2L, ncol = n)
    idx2 <- seq_len(n - 2L)
    second[cbind(idx2, idx2)] <- 1.0
    second[cbind(idx2, idx2 + 1L)] <- -2.0
    second[cbind(idx2, idx2 + 2L)] <- 1.0
    second
}

# Minimize 1/2 x' H x - g' x subject to x >= 0 with Hildreth's iterative
# quadratic-programming algorithm as used by the LSI step of LOUHI78: cyclic
# coordinate minimization with projection of every coordinate onto its
# non-negativity constraint.  The sweep stops when the relative change of the
# quadratic objective between two consecutive sweeps drops below `tolerance`
# (the ill-conditioning of Bonner-sphere response matrices makes the
# objective a far more robust convergence indicator than the raw iterates).
#
#   H          symmetric positive-definite Hessian (n x n)
#   g          linear term (n,)
#   x0         starting point (clipped to x >= 0)
#   max_iterations  maximum number of full coordinate sweeps (default 500)
#   tolerance       relative objective change for convergence (default 1e-6)
#
# Returns list(spectrum, sweeps, converged).
.louhi_hildreth_qp <- function(H, g, x0, max_iterations = 100L, tolerance = 1e-6) {
    x <- pmax(as.numeric(x0), 0)
    n <- length(x)
    diag_ <- pmax(diag(H), 1e-12)
    converged <- FALSE
    sweeps <- 0L
    prev_obj <- Inf
    for (sweep in seq_len(as.integer(max_iterations))) {
        sweeps <- sweep
        for (j in seq_len(n)) {
            grad_j <- as.numeric(H[j, , drop = TRUE] %*% x) - g[j]
            x[j] <- max(0, x[j] - grad_j / diag_[j])
        }
        Hx <- as.numeric(H %*% x)
        obj <- 0.5 * as.numeric(x %*% Hx) - as.numeric(g %*% x)
        if (abs(prev_obj - obj) <= tolerance * max(1, abs(obj))) {
            converged <- TRUE
            break
        }
        prev_obj <- obj
    }
    list(spectrum = x, sweeps = sweeps, converged = converged)
}

# Assemble the normal equations of the LOUHI quadratic program.  Returns
# list(H, g) for `min 1/2 x' H x - g' x` of the objective
# chi2_stat + lambda^2 * ||L (x - x0)||^2.
.louhi_weighted_normal_equations <- function(A, b, sigma, x0, L, smoothness) {
    w <- 1 / sigma^2
    ata <- t(A) %*% (A * w)
    atb <- as.numeric(t(A) %*% (b * w))
    ltl <- t(L) %*% L
    H <- 2 * (ata + smoothness^2 * ltl)
    g <- 2 * (atb + smoothness^2 * as.numeric(ltl %*% x0))
    list(H = H, g = g)
}

# Weighted chi-square of the data term for spectrum x.
.data_chi2 <- function(A, b, sigma, x) {
    resid <- (b - as.numeric(A %*% x)) / sigma
    as.numeric(resid %*% resid)
}

# Bracket the crossing of a monotone objective on [lo, hi]: scan a uniform
# grid on log10(lambda) in [lo, hi] and return the first adjacent pair of
# grid points where the objective changes sign.  Falls back to the full
# interval when no crossing exists.
.smoothness_bracket <- function(objective, lo = -6.0, hi = 6.0,
                                n_scan = 25L) {
    grid <- seq(lo, hi, length.out = as.integer(n_scan))
    values <- vapply(grid, function(t) as.numeric(objective(t)), numeric(1))
    for (k in seq_len(length(grid) - 1L)) {
        if (values[k] == 0 || values[k] * values[k + 1L] <= 0) {
            return(c(grid[k], grid[k + 1L]))
        }
    }
    c(lo, hi)
}

#' Solve the unfolding problem with the LOUHI78 algorithm
#'
#' Core solver mirroring \code{solve_louhi} in
#' \code{bssunfold/src/bssunfold/core/unfold_louhi.py}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional default (a-priori) spectrum (length n) used as the
#'   starting point and as the reference of the generalized smoothing term.
#'   Default \code{NULL} = flat unit spectrum.
#' @param smoothness Numeric; smoothing weight \code{lambda} of the
#'   generalized smoothing term.  Default 1.0.  Ignored when
#'   \code{auto_smooth = TRUE}.
#' @param smooth_order Integer; order of the smoothing operator: \code{0}
#'   (identity), \code{1} (first differences) or \code{2} (second
#'   differences).  Default 1.
#' @param max_iterations Positive integer; maximum number of Hildreth
#'   coordinate sweeps.  Default 500.
#' @param tolerance Positive numeric; maximum relative objective change
#'   between sweeps for convergence.  Default 1e-6.
#' @param relative_uncertainty Numeric; relative measurement uncertainty
#'   used to derive detector sigma values when \code{sigma} is not supplied.
#'   Default 0.1.  (The Python implementation derives the sigmas inside
#'   \code{_validate_inputs} with a hard-coded 0.1 factor and never reads
#'   this argument; the port reproduces that behaviour exactly.)
#' @param sigma Optional explicit per-detector measurement uncertainties
#'   (length m).  When given, overrides \code{relative_uncertainty}.
#'   Default \code{NULL}.
#' @param auto_smooth Logical; nonlinear regression mode of LOUHI78: adjust
#'   the smoothing weight by a golden-section search on \code{log10(lambda)}
#'   so that the data chi-square reaches \code{chi2_target}.  Default
#'   \code{FALSE}.
#' @param chi2_target Optional target data chi-square for
#'   \code{auto_smooth}.  Default \code{NULL} = the number of detectors (the
#'   expected value of the chi-square).
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_louhi(A, b, rep(1, 3), smoothness = 1.0)
solve_louhi <- function(A, b, x0 = NULL, smoothness = 1.0,
                        smooth_order = 1L, max_iterations = 500L,
                        tolerance = 1e-6, relative_uncertainty = 0.1,
                        sigma = NULL, auto_smooth = FALSE,
                        chi2_target = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    if (length(b) == 0L || nrow(A) == 0L) {
        stop("Response matrix must be a non-empty 2-D array, got shape (",
             nrow(A), ", ", ncol(A), ")")
    }
    if (length(b) != nrow(A)) {
        stop("Measurement vector length ", length(b),
             " does not match response matrix with ", nrow(A), " rows")
    }
    if (!any(b > 0)) {
        stop("At least one positive measurement is required; all readings are zero or negative")
    }
    if (!is.null(sigma)) {
        sigma <- as.numeric(sigma)
        if (length(sigma) != length(b)) {
            stop("sigma must have length ", length(b), ", got ", length(sigma))
        }
        sigma <- pmax(sigma, 1e-12)
    } else {
        # Python _validate_inputs() hard-codes the 0.1 relative factor and
        # never receives `relative_uncertainty`; mirrored here verbatim.
        sigma <- 0.1 * pmax(b, 1e-12)
    }

    n <- ncol(A)
    if (is.null(x0)) {
        x0 <- rep(1, n)
    } else {
        x0 <- as.numeric(x0)
    }
    x0 <- pmax(x0, 0)
    if (length(x0) != n) {
        stop("Default spectrum length ", length(x0),
             " does not match the number of energy bins ", n)
    }
    if (smoothness < 0) stop("smoothness must be non-negative, got ", smoothness)
    L <- louhi_smoothing_matrix(n, smooth_order)

    solve_with <- function(lam) {
        ng <- .louhi_weighted_normal_equations(A, b, sigma, x0, L, lam)
        .louhi_hildreth_qp(ng$H, ng$g, x0, max_iterations, tolerance)$spectrum
    }

    lambda_used <- as.numeric(smoothness)
    target <- if (!is.null(chi2_target)) as.numeric(chi2_target) else as.numeric(nrow(A))

    misfit <- function(log_lam) {
        lam <- 10^log_lam
        chi2 <- .data_chi2(A, b, sigma, solve_with(lam))
        abs(chi2 - target)
    }

    if (isTRUE(auto_smooth) && max_iterations > 0) {
        bracket <- .smoothness_bracket(function(t) {
            .data_chi2(A, b, sigma, solve_with(10^t)) - target
        })
        best <- .bss_golden_section(misfit, bracket[1], bracket[2])
        lambda_used <- 10^best$t_opt
    }

    ng <- .louhi_weighted_normal_equations(A, b, sigma, x0, L, lambda_used)
    res <- .louhi_hildreth_qp(ng$H, ng$g, x0,
                        max_iterations = max_iterations, tolerance = tolerance)
    list(spectrum = as.numeric(res$spectrum), iterations = res$sweeps,
         converged = res$converged)
}

#' LOUHI unfolding (unified workflow wrapper)
#'
#' Thin wrapper around \code{\link{solve_louhi}} for the unified workflow,
#' mirroring \code{unfold_louhi} in
#' \code{bssunfold/src/bssunfold/core/unfold_louhi.py}.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_louhi
#' @param calculate_errors Logical; if \code{TRUE}, run Monte-Carlo
#'   uncertainty estimation.  Default \code{FALSE}.
#' @param noise_level Numeric; relative Gaussian noise level for Monte-Carlo.
#'   Default 0.01.
#' @param n_montecarlo Integer; number of Monte-Carlo samples.  Default 100.
#' @param random_state Optional integer seed for Monte-Carlo.
#' @param max_neutron_energy Optional numeric energy cutoff in MeV.
#'   Default \code{NULL} = no cutoff.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_louhi <- function(detector_names, n_energy_bins, E_MeV,
                         sensitivities, cc_icrp116, save_result_callback,
                         readings, initial_spectrum = NULL,
                         smoothness = 1.0, smooth_order = 1L,
                         auto_smooth = FALSE, chi2_target = NULL,
                         max_iterations = 500L, tolerance = 1e-6,
                         relative_uncertainty = 0.1,
                         calculate_errors = FALSE,
                         noise_level = 0.01,
                         n_montecarlo = 100L,
                         save_result = FALSE, random_state = NULL,
                         max_neutron_energy = NULL) {
    default_initial <- rep(1.0, n_energy_bins)

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = default_initial,
        solve_func = make_solve_wrapper(solve_louhi,
                                        smoothness = smoothness,
                                        smooth_order = smooth_order,
                                        auto_smooth = auto_smooth,
                                        chi2_target = chi2_target,
                                        max_iterations = max_iterations,
                                        tolerance = tolerance,
                                        relative_uncertainty = relative_uncertainty),
        solve_kwargs = list(),
        method_name = "LOUHI",
        extra_output = list(smoothness = as.numeric(smoothness),
                            smooth_order = as.integer(smooth_order),
                            auto_smooth = as.logical(auto_smooth)),
        calculate_errors = calculate_errors,
        noise_level = noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
