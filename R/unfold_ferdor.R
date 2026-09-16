#' FERDOR (Ferret) unfolding
#'
#' Classic weighted least-squares with second-difference smoothing. The
#' smoothing weight \eqn{alpha} is adjusted iteratively (bisection) so the
#' reduced chi-square approaches \code{chi_squared_target} (discrepancy
#' principle). Non-negativity is enforced via the \pkg{lsei} NNLS solver.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum guess (length n) -- accepted for API
#'   compatibility; the solution does not depend on it.
#' @param max_iterations Positive integer; default 100.
#' @param tolerance Positive numeric; default 1e-3.
#' @param smoothing Numeric; initial smoothing weight alpha. Default 1e-3.
#' @param chi_squared_target Numeric; default 1.0.
#' @param relative_uncertainty Numeric; default 0.1.
#' @param sigma Optional per-detector uncertainties. Default \code{NULL}.
#' @param min_alpha Numeric; default 1e-12.
#' @param max_alpha Numeric; default 1e12.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_ferdor(A, b, rep(1, 3))
solve_ferdor <- function(A, b, x0, max_iterations = 100L, tolerance = 1e-3,
                          smoothing = 1e-3, chi_squared_target = 1.0,
                          relative_uncertainty = 0.1, sigma = NULL,
                          min_alpha = 1e-12, max_alpha = 1e12) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    m <- nrow(A); n <- ncol(A)
    if (m == 0L) stop("Measurement vector b is empty")
    if (all(b <= 0)) {
        stop("FERDOR requires at least one strictly positive measurement")
    }
    if (!is.null(sigma)) {
        sigma <- pmax(as.numeric(sigma), 1e-12)
    } else {
        sigma <- relative_uncertainty * pmax(abs(b), 1e-12)
    }
    if (length(sigma) != m) {
        stop("sigma must have length ", m)
    }
    Wsqrt <- 1.0 / sigma
    Aw <- A * matrix(Wsqrt, nrow = m, ncol = n, byrow = FALSE)
    bw <- b * Wsqrt
    ATA <- crossprod(Aw)
    ATb <- as.numeric(crossprod(Aw, bw))

    if (n > 2L) {
        L <- as.matrix(create_derivative_matrix(n, 2L))
        LTL <- crossprod(L)
    } else {
        L <- NULL
        LTL <- matrix(0.0, n, n)
    }
    dof <- max(m - 1L, 1L)

    # Try unregularised NNLS first
    x_init <- .ferdor_solve(ATA, ATb, LTL, 0.0, Aw, bw, L)
    if (!is.null(x_init)) {
        r0 <- as.numeric(A %*% x_init) - b
        chi2_0 <- sum((r0 / sigma)^2)
        ratio_0 <- chi2_0 / dof
        if (ratio_0 <= chi_squared_target) {
            return(list(spectrum = as.numeric(x_init), iterations = 1L,
                        converged = TRUE))
        }
    } else {
        x_init <- pmax(x0, 0.0)
    }

    lo <- min_alpha; hi <- max_alpha
    alpha <- min(max(smoothing, lo), hi)
    spectrum <- pmax(x0, 0.0)
    converged <- FALSE; iterations <- 0L

    for (it in seq_len(max_iterations)) {
        iterations <- it
        x <- .ferdor_solve(ATA, ATb, LTL, alpha, Aw, bw, L)
        if (is.null(x)) break
        spectrum <- x
        r <- as.numeric(A %*% x) - b
        chi2 <- sum((r / sigma)^2)
        ratio <- chi2 / dof
        if (abs(ratio - chi_squared_target) <=
            tolerance * max(abs(chi_squared_target), 1.0)) {
            converged <- TRUE; break
        }
        if (ratio > chi_squared_target) {
            hi <- alpha
        } else {
            lo <- alpha
        }
        if (hi <= lo) { converged <- TRUE; break }
        new_alpha <- sqrt(lo * hi)
        if (abs(new_alpha - alpha) <= tolerance * max(abs(alpha), 1e-30)) {
            converged <- TRUE; break
        }
        alpha <- new_alpha
    }
    list(spectrum = as.numeric(spectrum), iterations = iterations,
         converged = converged)
}

.ferdor_solve <- function(ATA, ATb, LTL, alpha, Aw = NULL, bw = NULL,
                          L = NULL) {
    nnls_extract <- function(A, b) {
        r <- tryCatch(lsei::nnls(A, b), error = function(e) NULL)
        if (is.null(r)) return(NULL)
        as.numeric(r$x)
    }
    if (alpha < 1e-20) {
        if (!is.null(Aw) && !is.null(bw)) {
            x <- nnls_extract(Aw, bw)
            if (!is.null(x)) return(x)
        }
        x <- tryCatch(as.numeric(qr.solve(ATA, ATb)),
                      error = function(e) NULL)
        if (is.null(x)) return(NULL)
        return(pmax(x, 0.0))
    }
    P <- ATA + alpha * LTL
    x <- tryCatch(as.numeric(qr.solve(P, ATb)), error = function(e) NULL)
    if (!is.null(x)) return(pmax(x, 0.0))
    # Augmented NNLS fallback
    if (!is.null(Aw) && !is.null(bw) && !is.null(L) && alpha > 0 &&
        any(LTL != 0)) {
        L_aug <- tryCatch({
            chol_L <- tryCatch(chol(LTL), error = function(e) NULL)
            if (is.null(chol_L)) NULL else t(chol_L) * sqrt(alpha)
        }, error = function(e) NULL)
        if (!is.null(L_aug)) {
            Aw_aug <- rbind(Aw, L_aug)
            bw_aug <- c(bw, rep(0, nrow(L_aug)))
            x <- nnls_extract(Aw_aug, bw_aug)
            if (!is.null(x)) return(x)
        }
    }
    if (!is.null(Aw) && !is.null(bw)) {
        x <- nnls_extract(Aw, bw)
        if (!is.null(x)) return(x)
    }
    NULL
}

#' Wrapper around \code{\link{solve_ferdor}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_ferdor
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_ferdor <- function(detector_names, n_energy_bins, E_MeV,
                           sensitivities, cc_icrp116, save_result_callback,
                           readings, initial_spectrum = NULL,
                           max_iterations = 100L, tolerance = 1e-3,
                           smoothing = 1e-3, chi_squared_target = 1.0,
                           relative_uncertainty = 0.1,
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
        solve_func = make_solve_wrapper(solve_ferdor,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance,
                                         smoothing = smoothing,
                                         chi_squared_target = chi_squared_target,
                                         relative_uncertainty = relative_uncertainty),
        solve_kwargs = list(),
        method_name = "FERDOR",
        extra_output = list(chi_squared_target = chi_squared_target,
                            relative_uncertainty = relative_uncertainty,
                            smoothing = smoothing),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
