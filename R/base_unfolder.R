#' Base unfolding workflow
#'
#' \code{run_unfolding} is the unified pipeline used by every
#' \code{unfold_*} method in this package. It centralises:
#'
#' \itemize{
#'   \item building the response matrix \eqn{A} and measurement vector \eqn{b}
#'         from the user-supplied readings and detector sensitivities;
#'   \item normalising the initial spectrum guess (or falling back to the
#'         method default);
#'   \item calling the core solver function (which may return either a single
#'         spectrum or a \code{list(spectrum, iterations, converged)});
#'   \item standardising the output dictionary (spectrum, residual,
#'         effective readings, dose rates, ...);
#'   \item optionally running Monte-Carlo uncertainty estimation;
#'   \item optionally saving the result to a history list.
#' }
#'
#' @param detector_names Character vector of detector names.
#' @param n_energy_bins Integer; number of energy bins in the output spectrum.
#' @param E_MeV Numeric energy grid (length \code{n_energy_bins}).
#' @param sensitivities Named list of numeric sensitivity vectors, one per
#'   detector (each of length \code{n_energy_bins}).
#' @param cc_icrp116 Named list of ICRP-116 conversion coefficients (or
#'   \code{NULL} to skip dose-rate computation).
#' @param save_result_callback Function of one argument (the result list)
#'   invoked when \code{save_result = TRUE}. Default \code{NULL} = no saving.
#' @param readings Named numeric vector of detector readings.
#' @param initial_spectrum Optional numeric spectrum guess (length
#'   \code{n_energy_bins}) or \code{NULL} to use \code{default_initial}.
#' @param default_initial Numeric default spectrum (length
#'   \code{n_energy_bins}).
#' @param solve_func Solver function with signature
#'   \code{solve_func(A, b, x0 = ..., ...)} returning a numeric spectrum
#'   or a list \code{list(spectrum, iterations, converged)}.
#' @param solve_kwargs Named list of extra keyword arguments forwarded to
#'   \code{solve_func}.
#' @param method_name Character string used as the \code{method} field of the
#'   result.
#' @param extra_output Optional named list of extra entries to merge into the
#'   result.
#' @param calculate_errors Logical; if \code{TRUE}, run Monte-Carlo
#'   uncertainty estimation. Default \code{FALSE}.
#' @param noise_level Numeric; relative Gaussian noise level for Monte-Carlo.
#'   Default 0.01.
#' @param n_montecarlo Integer; number of Monte-Carlo samples. Default 100.
#' @param random_state Optional integer seed for Monte-Carlo.
#' @param save_result Logical; if \code{TRUE}, call
#'   \code{save_result_callback}. Default \code{FALSE}.
#' @param max_neutron_energy Optional numeric energy cutoff in MeV: bins
#'   above this energy are removed from the response matrix during the
#'   solve and the returned spectrum is expanded back to the full grid
#'   with exact zeros above the cutoff (Python \code{max_neutron_energy}
#'   parameter). Default \code{NULL} = no cutoff.
#' @return A list with at minimum the following components:
#' \describe{
#'   \item{energy}{copy of E_MeV.}
#'   \item{spectrum}{unfolded, non-negative spectrum.}
#'   \item{spectrum_absolute}{same as \code{spectrum}.}
#'   \item{effective_readings}{named numeric vector of \code{A * spectrum}.}
#'   \item{residual}{\code{b - A * spectrum}.}
#'   \item{residual_norm}{sqrt(sum(residual^2)).}
#'   \item{method}{the \code{method_name} string.}
#'   \item{doserates}{output of \code{\link{calculate_dose_rates}}.}
#'   \item{iterations}{if the solver returned a list with \code{iterations}.}
#'   \item{converged}{if the solver returned a list with \code{converged}.}
#'   \item{montecarlo_samples}{present iff \code{calculate_errors = TRUE}.}
#'   \item{noise_level}{present iff \code{calculate_errors = TRUE}.}
#' }
#' @export
#' @examples
#' # Tiny synthetic MLEM-like example.
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' det_names <- c("d1", "d2", "d3")
#' sens <- lapply(setNames(det_names, det_names), function(n) A[which(det_names == n), ])
#' E <- c(1e-9, 1e-6, 1e-3)
#' solver <- function(A, b, x0 = NULL, ...) {
#'     x <- if (is.null(x0)) rep(1, ncol(A)) else x0
#'     for (i in 1:5) x <- x * (t(A) %*% (b / (A %*% x) + 1e-10))
#'     list(spectrum = as.numeric(x), iterations = 5L, converged = TRUE)
#' }
#' res <- run_unfolding(
#'     detector_names = det_names, n_energy_bins = 3L, E_MeV = E,
#'     sensitivities = sens, cc_icrp116 = NULL, save_result_callback = NULL,
#'     readings = c(d1 = 1, d2 = 0.6, d3 = 0.4),
#'     initial_spectrum = NULL, default_initial = rep(1, 3),
#'     solve_func = solver, solve_kwargs = list(),
#'     method_name = "demo"
#' )
#' str(res)
run_unfolding <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, initial_spectrum, default_initial,
                          solve_func, solve_kwargs = list(),
                          method_name = "unfold",
                          extra_output = NULL,
                          calculate_errors = FALSE,
                          noise_level = 0.01, n_montecarlo = 100L,
                          random_state = NULL,
                          save_result = FALSE,
                          max_neutron_energy = NULL) {
    # ---- 0. Validate inputs ----
    if (!is.numeric(readings) || length(readings) == 0L) {
        stop("'readings' must be a non-empty named numeric vector.")
    }
    if (!is.character(detector_names) || length(detector_names) == 0L) {
        stop("'detector_names' must be a non-empty character vector.")
    }
    n_energy_bins <- as.integer(n_energy_bins)
    if (n_energy_bins <= 0L) {
        stop("'n_energy_bins' must be a positive integer, got ", n_energy_bins)
    }
    E_MeV <- as.numeric(E_MeV)
    if (length(E_MeV) != n_energy_bins) {
        stop("Length of E_MeV (", length(E_MeV),
             ") must match n_energy_bins (", n_energy_bins, ")")
    }
    if (!is.numeric(noise_level) || length(noise_level) != 1L ||
        noise_level <= 0 || noise_level > 1) {
        stop("'noise_level' must be a number in (0, 1], got ", noise_level)
    }
    if (!is.numeric(n_montecarlo) || length(n_montecarlo) != 1L ||
        n_montecarlo < 0) {
        stop("'n_montecarlo' must be a non-negative integer, got ",
             n_montecarlo)
    }
    n_montecarlo <- as.integer(n_montecarlo)

    # ---- 1. Build the (A, b) system ----
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected

    # ---- 1b. Optional energy cutoff (Python max_neutron_energy) ----
    E_full <- as.numeric(E_MeV)
    n_full <- n_energy_bins
    A_kept <- A
    if (!is.null(max_neutron_energy)) {
        cutoff <- as.numeric(max_neutron_energy)
        if (is.finite(cutoff) && cutoff > 0) {
            keep <- E_full <= cutoff
            if (sum(keep) == 0L) {
                stop("max_neutron_energy (", cutoff,
                     ") is below the lowest energy bin.")
            }
            A_kept <- A[, keep, drop = FALSE]
        }
    }

    # ---- 2. Normalize the initial spectrum (validated on the full grid) --
    x0 <- .normalize_initial(initial_spectrum, default_initial, n_full)

    # ---- 3. Solve (on the trimmed matrix when a cutoff is set) -----------
    kwargs <- c(list(x0 = x0), solve_kwargs)
    use_cutoff <- (!is.null(max_neutron_energy) &&
                       is.finite(as.numeric(max_neutron_energy)) &&
                       as.numeric(max_neutron_energy) > 0 &&
                       ncol(A_kept) < ncol(A))
    if (!isTRUE(use_cutoff)) {
        solve_result <- do.call(solve_func, c(list(A = A, b = b), kwargs))
    } else {
        kwargs_kept <- kwargs
        x0_kept <- as.numeric(x0)[keep_idx <- which(E_full <= as.numeric(max_neutron_energy))]
        kwargs_kept$x0 <- x0_kept
        solve_result <- do.call(solve_func,
                                c(list(A = A_kept, b = b), kwargs_kept))
    }

    extra_meta <- list()
    if (is.list(solve_result) &&
        !is.null(solve_result$spectrum) &&
        (is.null(dim(solve_result$spectrum)) || is.numeric(solve_result$spectrum))) {
        spectrum <- as.numeric(solve_result$spectrum)
        if (!is.null(solve_result$iterations)) {
            extra_meta$iterations <- as.integer(solve_result$iterations)
        }
        if (!is.null(solve_result$converged)) {
            extra_meta$converged <- as.logical(solve_result$converged)
        }
    } else {
        spectrum <- as.numeric(solve_result)
    }
    # Expand a cut-off solve back to the full energy grid with exact zeros
    # above the cutoff, and standardize on the full grid.
    if (isTRUE(use_cutoff)) {
        spectrum_full <- numeric(length(E_full))
        spectrum_full[keep_idx] <- spectrum
        spectrum <- spectrum_full
    }
    if (length(spectrum) != n_full) {
        stop("Solver returned a spectrum of length ", length(spectrum),
             " but n_energy_bins = ", n_full)
    }
    E_MeV <- E_full
    if (!is.null(extra_output)) {
        extra_output <- c(extra_output, extra_meta)
    } else {
        extra_output <- extra_meta
    }
    # The computed model must be evaluated on the FULL response matrix:
    output <- .standardize_output(
        spectrum = spectrum, A = A, b = b, E_MeV = E_MeV,
        selected = selected, cc_icrp116 = cc_icrp116,
        method = method_name, extra = extra_output
    )

    # ---- 5. Monte-Carlo uncertainty ----
    if (isTRUE(calculate_errors) && n_montecarlo > 0L) {
        mc <- .montecarlo_for_run(
            solve_func = solve_func, readings = readings,
            noise_level = noise_level, n_montecarlo = n_montecarlo,
            n_energy_bins = n_energy_bins, random_state = random_state,
            solve_kwargs = solve_kwargs, detector_names = detector_names,
            sensitivities = sensitivities, x0 = x0
        )
        for (n in names(mc)) output[[n]] <- mc[[n]]
        output$montecarlo_samples <- n_montecarlo
        output$noise_level <- noise_level
    }

    # ---- 6. Save result ----
    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(output)
    }

    output
}

# --- internal helpers ---

.build_system <- function(readings, detector_names, sensitivities) {
    selected <- detector_names[detector_names %in% names(readings)]
    if (length(selected) == 0L) {
        stop("None of the detector_names appear in 'readings'.")
    }
    b <- as.numeric(readings[selected])
    A <- do.call(rbind, lapply(selected, function(n) {
        if (!n %in% names(sensitivities)) {
            stop("Detector '", n, "' missing from 'sensitivities'.")
        }
        as.numeric(sensitivities[[n]])
    }))
    list(A = A, b = b, selected = selected)
}

.normalize_initial <- function(initial_spectrum, default_initial, n_energy_bins) {
    if (!is.null(initial_spectrum)) {
        if (is.list(initial_spectrum)) {
            if (!is.null(initial_spectrum$spectrum)) {
                initial_spectrum <- initial_spectrum$spectrum
            } else {
                return(as.numeric(default_initial))
            }
        }
        spectrum <- as.numeric(initial_spectrum)
        if (length(spectrum) != n_energy_bins) {
            stop("Initial spectrum length (", length(spectrum),
                 ") must match number of energy bins (", n_energy_bins, ")")
        }
        return(pmax(spectrum, 0))
    }
    as.numeric(default_initial)
}

.standardize_output <- function(spectrum, A, b, E_MeV, selected, cc_icrp116,
                                method, extra = NULL) {
    spectrum_nonneg <- pmax(as.numeric(spectrum), 0)
    computed_readings <- as.numeric(A %*% spectrum_nonneg)
    residual <- b - computed_readings
    out <- list(
        energy = E_MeV,
        spectrum = spectrum_nonneg,
        spectrum_absolute = spectrum_nonneg,
        effective_readings = stats::setNames(computed_readings, selected),
        residual = residual,
        residual_norm = sqrt(sum(residual^2)),
        method = method,
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum_nonneg, cc_icrp116)
                    else numeric(0L)
    )
    if (!is.null(extra)) {
        for (n in names(extra)) out[[n]] <- extra[[n]]
    }
    out
}

.montecarlo_for_run <- function(solve_func, readings, noise_level,
                                n_montecarlo, n_energy_bins, random_state,
                                solve_kwargs, detector_names, sensitivities,
                                x0) {
    mc_solver <- function(noisy_readings, ...) {
        sys <- .build_system(noisy_readings, detector_names, sensitivities)
        A_noisy <- sys$A; b_noisy <- sys$b
        kw <- list(...)
        # filter kwargs used by montecarlo dispatcher but not the solver
        kw <- kw[!names(kw) %in% c("detector_names", "sensitivities")]
        kw$x0 <- x0
        result <- do.call(solve_func,
                          c(list(A = A_noisy, b = b_noisy), kw))
        if (is.list(result) && !is.null(result$spectrum)) {
            return(as.numeric(result$spectrum))
        }
        as.numeric(result)
    }
    mc_kwargs <- c(solve_kwargs,
                   list(detector_names = detector_names,
                        sensitivities = sensitivities))
    do.call(monte_carlo_uncertainty,
            c(list(func = mc_solver, readings = readings,
                   noise_level = noise_level, n_samples = n_montecarlo,
                   n_energy_bins = n_energy_bins,
                   random_state = random_state), mc_kwargs))
}
