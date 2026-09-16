#' Iterative refinement unfolding (two-pass)
#'
#' R port of \code{bssunfold/core/unfold_iterative_refinement.py}.
#' Two-pass method combining a fast first-pass solver (e.g. MLEM) with a
#' residual-correction second-pass solver (e.g. Landweber or CGLS).
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial guess (length n). Default \code{NULL} = flat 0.5.
#' @param first_pass_solver Optional function \code{(A, b, x0, ...)} for the
#'   first pass. Default \code{NULL} = \code{\link{solve_mlem}}.
#' @param second_pass_solver Optional function for the residual-correction
#'   second pass. Default \code{NULL} = \code{\link{solve_landweber}}.
#' @param first_pass_kwargs Optional list of kwargs for \code{first_pass_solver}.
#' @param second_pass_kwargs Optional list of kwargs for
#'   \code{second_pass_solver}.
#' @param alpha Optional numeric blending factor in
#'   \code{x_final = x1 + alpha * x2}. Default \code{NULL} = auto-line search.
#' @param max_alpha_search Integer; line-search grid size. Default 20.
#' @return A list \code{list(spectrum = ..., info = list(...))}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_iterative_refinement(A, b, rep(1, 3))
solve_iterative_refinement <- function(A, b, x0 = NULL,
                                         first_pass_solver = NULL,
                                         second_pass_solver = NULL,
                                         first_pass_kwargs = NULL,
                                         second_pass_kwargs = NULL,
                                         alpha = NULL,
                                         max_alpha_search = 20L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (is.null(x0)) x0 <- rep(0.5, n)
    if (is.null(first_pass_solver)) first_pass_solver <- solve_mlem
    if (is.null(second_pass_solver)) second_pass_solver <- solve_landweber
    if (is.null(first_pass_kwargs)) {
        first_pass_kwargs <- list(max_iterations = 150L, tolerance = 1e-4)
    }
    if (is.null(second_pass_kwargs)) {
        second_pass_kwargs <- list(max_iterations = 100L, tolerance = 1e-5)
    }

    res1 <- do.call(first_pass_solver,
                    c(list(A = A, b = b, x0 = x0), first_pass_kwargs))
    x1 <- if (is.list(res1) && !is.null(res1$spectrum)) res1$spectrum
          else as.numeric(res1)
    x1 <- pmax(as.numeric(x1), 0)
    r <- b - as.numeric(A %*% x1)
    res2 <- do.call(second_pass_solver,
                    c(list(A = A, b = r, x0 = rep(0.0, n)),
                      second_pass_kwargs))
    x2 <- if (is.list(res2) && !is.null(res2$spectrum)) res2$spectrum
          else as.numeric(res2)
    x2 <- as.numeric(x2)
    if (!is.null(alpha)) {
        best_alpha <- alpha
    } else {
        candidates <- seq(0.0, 2.0, length.out = max_alpha_search)
        best_alpha <- 0.0
        best_res <- sqrt(sum((as.numeric(A %*% x1) - b)^2))
        for (a in candidates) {
            x_cand <- x1 + a * x2
            res_cand <- sqrt(sum((as.numeric(A %*% x_cand) - b)^2))
            if (res_cand < best_res) {
                best_res <- res_cand
                best_alpha <- a
            }
        }
    }
    spectrum <- pmax(x1 + best_alpha * x2, 0)
    info <- list(
        first_pass_residual = sqrt(sum(r^2)),
        second_pass_correction_norm = sqrt(sum(x2^2)),
        alpha = best_alpha,
        final_residual = sqrt(sum((as.numeric(A %*% spectrum) - b)^2))
    )
    list(spectrum = spectrum, info = info)
}

#' Wrapper around \code{\link{solve_iterative_refinement}} for the unified
#' workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_iterative_refinement
#' @return A result list as produced by \code{\link{run_unfolding}} plus a
#'   \code{parameters} entry with the two-pass diagnostics.
#' @export
unfold_iterative_refinement <- function(detector_names, n_energy_bins, E_MeV,
                                          sensitivities, cc_icrp116,
                                          save_result_callback, readings,
                                          initial_spectrum = NULL,
                                          first_pass_kwargs = NULL,
                                          second_pass_kwargs = NULL,
                                          alpha = NULL,
                                          max_alpha_search = 20L,
                                          calculate_errors = FALSE,
                                          noise_level = 0.01,
                                          n_montecarlo = 100L,
                                          save_result = FALSE,
                                          random_state = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    x0_default <- rep(0.5, n_energy_bins)
    x0 <- if (is.null(initial_spectrum)) x0_default else as.numeric(initial_spectrum)
    res <- solve_iterative_refinement(
        A, b, x0,
        first_pass_kwargs = first_pass_kwargs,
        second_pass_kwargs = second_pass_kwargs,
        alpha = alpha, max_alpha_search = max_alpha_search
    )
    spectrum <- pmax(res$spectrum, 0)
    computed_readings <- as.numeric(A %*% spectrum)
    residual <- b - computed_readings
    result <- list(
        energy = E_MeV,
        spectrum = spectrum,
        spectrum_absolute = spectrum,
        effective_readings = stats::setNames(computed_readings, selected),
        residual = residual,
        residual_norm = sqrt(sum(residual^2)),
        method = "IterativeRefinement",
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = 0L,
        parameters = c(
            list(first_pass_kwargs = first_pass_kwargs,
                 second_pass_kwargs = second_pass_kwargs),
            res$info
        )
    )
    if (isTRUE(calculate_errors) && n_montecarlo > 0L) {
        if (!is.null(random_state)) set.seed(as.integer(random_state))
        spectra_mc <- vector("list", n_montecarlo)
        for (i in seq_len(n_montecarlo)) {
            b_pert <- b * (1.0 + noise_level * rnorm(length(b)))
            r_mc <- tryCatch(
                solve_iterative_refinement(
                    A, pmax(b_pert, 0), x0,
                    first_pass_kwargs = first_pass_kwargs,
                    second_pass_kwargs = second_pass_kwargs,
                    alpha = alpha, max_alpha_search = max_alpha_search
                ),
                error = function(e) NULL
            )
            if (!is.null(r_mc)) spectra_mc[[i]] <- r_mc$spectrum
        }
        ok <- !vapply(spectra_mc, is.null, logical(1L))
        if (any(ok)) {
            mat <- do.call(rbind, spectra_mc[ok])
            result$spectrum_uncertainty <- matrixStats_col_sd(mat)
            result$calculate_errors <- TRUE
            result$n_montecarlo <- sum(ok)
        }
    }
    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }
    result
}
