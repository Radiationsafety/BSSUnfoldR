#' Hybrid parametric unfolding (parametric model + iterative refinement)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_hybrid_parametric.py}.
#' Hybrid approach: the FRUIT parametric model (see
#' \code{\link{solve_parametric}}) provides a physically motivated initial
#' guess, which is then refined with a nonparametric iterative solver
#' (Landweber or MLEM) so that the readings are fitted better.
#'
#' @description The parametric starting point is generated exactly like the
#' Python original: a brute-force scan over \eqn{P_{th}}/\eqn{P_{epi}}
#' (\code{.fruit_find_initial_params}) with the bin widths \emph{recomputed
#' internally} as \code{compute_log_steps(E) * log(10)}, i.e. the passed
#' \code{log_steps} is not used for the guess.  The parametric spectrum is
#' weighted by those steps once more (lethargy density on top of the already
#' lethargy-weighted response matrix), floor-clamped at \code{1e-30} and then
#' refined.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param E Numeric energy grid (length n).
#' @param log_steps Numeric log-energy bin widths (length n).  Kept for API
#'   compatibility with Python; the guess is recomputed internally from
#'   \code{E}, exactly like \code{unfold_hybrid_parametric.py}.
#' @param refinement_method Character; \code{"landweber"} (default) or
#'   \code{"mlem"}.
#' @param max_iterations Integer; maximum refinement iterations. Default 100.
#' @param tolerance Numeric; convergence tolerance. Default 1e-6.
#' @param step_size Numeric; Landweber relaxation parameter. Default 0.01.
#' @param initial_params Optional named list overriding the FRUIT parameters
#'   \code{b}, \code{beta_prime}, \code{alpha}, \code{beta}, \code{P_th},
#'   \code{P_epi}.  When \code{NULL} the Python grid scan is used.
#' @return A list \code{list(spectrum, iterations, converged, message,
#'   n_iter)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_hybrid_parametric(A, b, E, compute_log_steps(E) * log(10),
#'                              max_iterations = 30)
solve_hybrid_parametric <- function(A, b, E, log_steps,
                                    refinement_method = "landweber",
                                    max_iterations = 100L,
                                    tolerance = 1e-6,
                                    step_size = 0.01,
                                    initial_params = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); E <- as.numeric(E)
    n_energy <- ncol(A)
    refinement_method <- match.arg(tolower(refinement_method),
                                   c("landweber", "mlem"))

    # ---- Stage 1: parametric initial guess (FRUIT grid scan) ----
    parametric_guess <- tryCatch({
        # Python recomputes the steps internally and IGNORES the passed value.
        # A grid/matrix mismatch makes the broadcast fail there, which the
        # except clause below mirrors.
        if (length(E) != n_energy) {
            stop("energy grid length does not match the response matrix")
        }
        steps <- .compute_log_steps_n(E, n_energy) * log(10)
        best_params <- if (is.null(initial_params)) {
            .fruit_find_initial_params(A, b, E, steps, n_grid = 5L,
                                       return_top = 1L)
        } else {
            .fruit_get_initial_params(initial_params)
        }
        guess <- .fruit_model_vec(E, best_params) * steps
        pmax(as.numeric(guess), 1e-30)
    }, error = function(e) {
        # Fallback to a flat spectrum, as in Python
        rep(mean(b) / max(mean(rowSums(A)), .Machine$double.xmin), n_energy)
    })

    # ---- Stage 2: iterative refinement ----
    if (refinement_method == "landweber") {
        ref <- .landweber_iteration(parametric_guess, A, b, step_size,
                                    max_iterations, tolerance)
    } else {
        ref <- .mlem_iteration(parametric_guess, A, b, max_iterations,
                               tolerance)
    }
    n_iter <- ref$n_iter
    success <- n_iter < as.integer(max_iterations)
    message <- if (success) {
        sprintf("Converged in %d iterations", n_iter)
    } else {
        "Max iterations reached"
    }
    list(spectrum = ref$x, iterations = as.integer(n_iter),
         converged = success, message = message, n_iter = n_iter,
         refinement_method = refinement_method,
         parametric_guess = parametric_guess)
}

#' Log10 bin widths for the first \code{n} energy grid points
#'
#' Exact port of \code{_matrix_utils.compute_log_steps(E_MeV, n_energy_bins)}:
#' the array is sized by \code{n_energy_bins}, so a truncated response matrix
#' changes the last (edge) step while the package-level
#' \code{\link{compute_log_steps}} always uses the full grid.
#' @keywords internal
#' @noRd
.compute_log_steps_n <- function(E_MeV, n_energy_bins) {
    log_e <- log10(as.numeric(E_MeV) + 1e-15)
    out <- rep(0, as.integer(n_energy_bins))
    if (n_energy_bins > 1L) {
        out[1L] <- log_e[2L] - log_e[1L]
        out[n_energy_bins] <- log_e[n_energy_bins] - log_e[n_energy_bins - 1L]
        if (n_energy_bins > 2L) {
            idx <- seq.int(2L, n_energy_bins - 1L)
            out[idx] <- (log_e[idx + 1L] - log_e[idx - 1L]) / 2
        }
    } else {
        out[1L] <- 1
    }
    out
}

#' Landweber refinement step (port of _landweber_iteration)
#' @keywords internal
#' @noRd
.landweber_iteration <- function(spectrum, A, b, step_size, max_iter,
                                 tolerance) {
    x <- as.numeric(spectrum)
    max_iter <- as.integer(max_iter)
    for (i in seq_len(max_iter)) {
        residual <- as.numeric(b) - as.numeric(A %*% x)
        gradient <- as.numeric(crossprod(A, residual))
        x_new <- pmax(x + step_size * gradient, 0)
        if (sqrt(sum((x_new - x)^2)) < tolerance) {
            return(list(x = x_new, n_iter = i))
        }
        x <- x_new
    }
    list(x = x, n_iter = max_iter)
}

#' MLEM refinement step (port of _mlem_iteration)
#' @keywords internal
#' @noRd
.mlem_iteration <- function(spectrum, A, b, max_iter, tolerance) {
    x <- pmax(as.numeric(spectrum), 1e-15)
    b <- as.numeric(b)
    max_iter <- as.integer(max_iter)
    for (i in seq_len(max_iter)) {
        computed <- pmax(as.numeric(A %*% x), 1e-15)
        ratio <- b / computed
        correction <- as.numeric(crossprod(A, ratio))
        x_new <- pmax(x * correction, 0)
        if (sqrt(sum((x_new - x)^2)) / (sqrt(sum(x^2)) + 1e-15) < tolerance) {
            return(list(x = x_new, n_iter = i))
        }
        x <- x_new
    }
    list(x = x, n_iter = max_iter)
}

#' Wrapper around \code{\link{solve_hybrid_parametric}} for the unified
#' workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_hybrid_parametric
#' @param ln_steps Optional numeric \eqn{d(\ln E)} vector; accepted for API
#'   parity with Python (only used for the dose calculation of the standard
#'   output).
#' @param calculate_errors Logical; run Monte-Carlo uncertainty estimation.
#' @param noise_level Numeric; relative noise for the Monte-Carlo part.
#' @param n_montecarlo Integer; Monte-Carlo sample count.
#' @param save_result Logical; hand the result to \code{save_result_callback}.
#' @param random_state Optional integer seed for the Monte-Carlo part.
#' @param max_neutron_energy Optional energy cutoff in MeV.
#' @param refinement_iterations Integer; deprecated alias of
#'   \code{max_iterations}. Default \code{NULL}.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_hybrid_parametric <- function(detector_names, n_energy_bins, E_MeV,
                                         sensitivities, cc_icrp116,
                                         save_result_callback, readings,
                                         ln_steps = NULL,
                                         initial_spectrum = NULL,
                                         refinement_method = "landweber",
                                         max_iterations = 100L,
                                         tolerance = 1e-6,
                                         step_size = 0.01,
                                         initial_params = NULL,
                                         refinement_iterations = NULL,
                                         calculate_errors = FALSE,
                                         noise_level = 0.01,
                                         n_montecarlo = 100L,
                                         save_result = FALSE,
                                         random_state = NULL,
                              max_neutron_energy = NULL) {
    if (!is.null(refinement_iterations) && missing(max_iterations)) {
        max_iterations <- refinement_iterations
    }
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b

    solver <- function(A, b, x0 = NULL, ...) {
        res <- solve_hybrid_parametric(
            A, b, E_MeV, compute_log_steps(E_MeV),
            refinement_method = refinement_method,
            max_iterations = max_iterations, tolerance = tolerance,
            step_size = step_size, initial_params = initial_params)
        list(spectrum = res$spectrum, iterations = res$iterations,
             converged = res$converged)
    }

    x0_default <- rep(mean(b) / max(mean(rowSums(A)), .Machine$double.xmin),
                      n_energy_bins)

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver, solve_kwargs = list(),
        method_name = sprintf("hybrid_parametric (%s)", refinement_method),
        extra_output = list(refinement_method = refinement_method,
                            max_iterations = max_iterations,
                            tolerance = tolerance, step_size = step_size),
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
