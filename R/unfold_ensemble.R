#' Ensemble / Cascade / Composite unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_ensemble.py},
#' \code{unfold_composite.py}, and \code{unfold_cascade.py} (combined into
#' one file for brevity).
#'
#' \describe{
#'   \item{Ensemble}{Runs multiple underlying solvers and combines their
#'     spectra (weighted average with inverse-residual weights, or median,
#'     trimmed mean, best-residual selection).}
#'   \item{Cascade}{Runs a sequence of Detector-level unfolding stages
#'     (TSVD -> MLEM -> Bayes-Spline by default), each stage potentially
#'     reusing the previous stage's spectrum as initial guess or prior.}
#'   \item{Composite}{Runs a pool of Detector-level methods and combines
#'     their results with confidence-weighted averaging (mean cosine
#'     similarity to the other members).}
#' }
#'
#' @name ensemble-methods
NULL

# --- internal helpers mirroring the Python module machinery ---

#' Cosine similarity, mirroring bssunfold.utils.comparison.cosine_similarity.
#' Returns 0 when either vector has zero norm.
#' @param p,q Numeric vectors of equal length.
#' @return Numeric cosine similarity of \code{p} and \code{q}.
#' @keywords internal
.bss_cosine_similarity <- function(p, q) {
    p <- as.numeric(p); q <- as.numeric(q)
    nrm_p <- sqrt(sum(p^2))
    nrm_q <- sqrt(sum(q^2))
    if (nrm_p == 0 || nrm_q == 0) return(0)
    sum(p * q) / (nrm_p * nrm_q)
}

#' Confidence weight of one solution against the rest of the ensemble,
#' mirroring _confidence_weight() in unfold_composite.py.
#' @param spectrum Numeric candidate spectrum.
#' @param others List of the other ensemble spectra.
#' @return Numeric weight between 0 and 1: the mean cosine similarity to
#'   \code{others}, clipped to that range (1 for a lone member).
#' @keywords internal
.bss_confidence_weight <- function(spectrum, others) {
    if (length(others) == 0L) return(1)
    sims <- vapply(others, function(o) .bss_cosine_similarity(spectrum, o),
                   numeric(1L))
    sims <- sims[is.finite(sims)]
    if (length(sims) == 0L) return(1)
    mean_s <- mean(sims)
    min(max(mean_s, 0), 1)
}

#' Quality metrics, mirroring compute_quality_metrics() in unfold_cascade.py.
#' @param spectrum Numeric unfolded spectrum.
#' @param reconstructed_readings Numeric readings folded from \code{spectrum}.
#' @param measured_readings Numeric measured readings of the same detectors.
#' @param energy Numeric energy grid in MeV, used for the hardness ratio.
#' @return List with \code{chi_square}, \code{smoothness}, \code{flux_error},
#'   \code{negativity_count}, \code{hardness_ratio}, \code{peak_count} and
#'   \code{overall_quality}.
#' @keywords internal
.bss_compute_quality_metrics <- function(spectrum, reconstructed_readings,
                                          measured_readings, energy) {
    eps <- 1e-10
    spectrum <- as.numeric(spectrum)
    reconstructed_readings <- as.numeric(reconstructed_readings)
    measured_readings <- as.numeric(measured_readings)
    energy <- as.numeric(energy)

    residuals <- (measured_readings - reconstructed_readings) /
        (reconstructed_readings + eps)
    chi_square <- sum(residuals^2)

    log_spectrum <- log(spectrum + eps)
    n_spec <- length(log_spectrum)
    smoothness <- if (n_spec > 2L) {
        second_deriv <- diff(log_spectrum, differences = 2L)
        # np.std: population standard deviation
        sd_p <- sqrt(mean((second_deriv - mean(second_deriv))^2))
        1 / (1 + sd_p)
    } else 1
    total_flux_spec <- sum(spectrum)
    total_flux_readings <- sum(measured_readings)
    flux_error <- abs(total_flux_spec - total_flux_readings) /
        (total_flux_readings + eps)
    negativity_count <- sum(spectrum < 0)

    hardness_ratio <- if (length(spectrum) > 10L) {
        thermal_region <- energy < 0.5
        fast_region <- energy > 5.0
        thermal_fraction <- sum(spectrum[thermal_region]) / (total_flux_spec + eps)
        fast_fraction <- sum(spectrum[fast_region]) / (total_flux_spec + eps)
        fast_fraction / (thermal_fraction + eps)
    } else 0

    peak_count <- if (length(spectrum) > 3L) {
        x <- spectrum
        sum(x[-c(1L, length(x))] > x[-c(length(x) - 1L, length(x))] &
            x[-c(1L, length(x))] > x[-c(1L, 2L)])
    } else 0

    list(
        chi_square = chi_square,
        smoothness = smoothness,
        flux_error = flux_error,
        negativity_count = as.integer(negativity_count),
        hardness_ratio = hardness_ratio,
        peak_count = as.integer(peak_count),
        overall_quality = smoothness / (1 + chi_square + flux_error * 10)
    )
}

.bss_numeric <- function(x) suppressWarnings(as.numeric(x))

#' Build a default-ensemble member: solver wrapper that injects Python's
#' conservative default kwargs (max_iterations = 200, tolerance = 1e-4),
#' letting caller-supplied dots override them (keep-last deduplication).
#' @param solver Function with signature \code{(A, b, x0, ...)}.
#' @param kwargs Named list of defaults injected ahead of the caller's
#'   \code{...} arguments.
#' @return A solver function suitable for \code{run_unfolding}.
#' @keywords internal
.bss_default_member <- function(solver, kwargs) {
    force(solver)
    function(A, b, x0, ...) {
        dots <- list(...)
        call_args <- c(list(A = A, b = b, x0 = x0), kwargs, dots)
        keep <- !duplicated(names(call_args), fromLast = TRUE)
        do.call(solver, call_args[keep])
    }
}

#' Default ensemble members, mirroring _ensure_default_methods() in
#' unfold_ensemble.py: MLEM, Bayes, Landweber, CGLS, GRAVEL with
#' max_iterations = 200 and tolerance = 1e-4 each.
#' @keywords internal
.bss_default_ensemble_solvers <- function() {
    kwargs <- list(max_iterations = 200L, tolerance = 1e-4)
    list(
        .bss_default_member(solve_mlem, kwargs),
        .bss_default_member(solve_bayes, kwargs),
        .bss_default_member(solve_landweber, kwargs),
        .bss_default_member(solve_cgls, kwargs),
        .bss_default_member(solve_gravel, kwargs)
    )
}

#' METHOD_DISPATCH from unfold_cascade.py (short name -> Detector.unfold_*).
#' @keywords internal
.bss_cascade_dispatch <- c(
    tsvd = "unfold_tsvd",
    bayes = "unfold_bayes",
    cvxpy = "unfold_cvxpy",
    qpsolvers = "unfold_qpsolvers",
    statreg = "unfold_statreg",
    landweber = "unfold_landweber",
    mlem = "unfold_mlem",
    bayes_spline = "unfold_bayes_spline_regularization",
    cgls = "unfold_cgls",
    hybrid_gmres = "unfold_hybrid_gmres",
    parametric2 = "unfold_parametric2",
    hybrid_parametric = "unfold_hybrid_parametric",
    tikhonov_tv = "unfold_tikhonov_tv",
    gravel = "unfold_gravel",
    kaczmarz = "unfold_kaczmarz",
    genetic = "unfold_genetic",
    mystic = "unfold_mystic",
    mystic_hybrid = "unfold_mystic_hybrid",
    scip = "unfold_scip",
    docplex = "unfold_docplex",
    epic = "unfold_epic",
    cs = "unfold_cs",
    lanczos = "unfold_lanczos"
)

#' METHOD_DISPATCH from unfold_composite.py.
#' @keywords internal
.bss_composite_dispatch <- c(
    tsvd = "unfold_tsvd",
    bayes = "unfold_bayes",
    cvxpy = "unfold_cvxpy",
    statreg = "unfold_statreg",
    lanczos = "unfold_lanczos",
    mlem = "unfold_mlem",
    landweber = "unfold_landweber",
    bayes_spline = "unfold_bayes_spline_regularization",
    gravel = "unfold_gravel",
    qpsolvers = "unfold_qpsolvers",
    hybrid_parametric = "unfold_hybrid_parametric",
    parametric2 = "unfold_parametric2",
    genetic = "unfold_genetic",
    interpret = "unfold_interpret",
    maeo_ensemble = "unfold_maeo",
    mystic = "unfold_mystic",
    mystic_hybrid = "unfold_mystic_hybrid",
    cs = "unfold_cs",
    scip = "unfold_scip",
    docplex = "unfold_docplex",
    epic = "unfold_epic",
    kaczmarz = "unfold_kaczmarz"
)

#' GENERAL_METHODS from unfold_composite.py (curated fallback pool).
#' @keywords internal
.bss_general_methods <- c("tsvd", "mlem", "cvxpy", "qpsolvers", "bayes_spline")

#' DEFAULT_ENSEMBLE_WEIGHTS from unfold_composite.py.
#' @keywords internal
.bss_composite_base_weights <- c(
    tsvd = 1.0, bayes = 1.0, cvxpy = 1.0, statreg = 1.0, lanczos = 1.0,
    mlem = 1.0, landweber = 1.0, bayes_spline = 1.0, gravel = 1.0,
    qpsolvers = 1.0, hybrid_parametric = 1.0, parametric2 = 1.0,
    genetic = 0.8, interpret = 0.8, maeo_ensemble = 0.8, mystic = 0.8,
    mystic_hybrid = 0.85, cs = 0.8, scip = 0.8, docplex = 0.8,
    epic = 0.8, kaczmarz = 1.0
)

#' The QP backends that Python's unfold_qpsolvers delegates to are not
#' available in this environment (Python raises ImportError; the bssunfold
#' qpsolvers method is then skipped by the composite/cascade dispatchers).
#' @keywords internal
.bss_qp_backend_available <- function() {
    requireNamespace("quadprog", quietly = TRUE) ||
        requireNamespace("osqp", quietly = TRUE)
}

.bss_resolve_method <- function(dispatch, name, fallback = TRUE) {
    known <- name %in% names(dispatch)
    attr <- if (known) unname(dispatch[name])
            else if (fallback) paste0("unfold_", name) else ""
    if (!nzchar(attr)) return(NULL)
    get0(attr, inherits = TRUE, mode = "function")
}

#' Default "general" cascade stages, mirroring
#' create_default_cascade("general") in unfold_cascade.py.
#' @keywords internal
.bss_default_cascade_stages <- function() {
    list(
        list(method = "tsvd",
             params = list(k = 15, method = "discrepancy"),
             use_as_initial = FALSE, use_as_prior = FALSE,
             store_intermediate = TRUE, quality_threshold = 0.3,
             max_iterations = NULL),
        list(method = "mlem",
             params = list(max_iterations = 150L),
             use_as_initial = TRUE, use_as_prior = FALSE,
             store_intermediate = FALSE, quality_threshold = NULL,
             max_iterations = NULL),
        list(method = "bayes_spline",
             params = list(spline_smooth = 0.3),
             use_as_initial = TRUE, use_as_prior = TRUE,
             store_intermediate = FALSE, quality_threshold = NULL,
             max_iterations = NULL)
    )
}

#' Stage-level cascade driver mirroring the module-level unfold_cascade()
#' in bssunfold/core/unfold_cascade.py (default "general" stage sequence,
#' quality-threshold early stop, Detector-level member calls).
#' @inheritParams run_unfolding
#' @return List with the final \code{spectrum}, \code{stages_run},
#'   \code{method_sequence}, \code{cascade_spectra},
#'   \code{intermediate_results}, \code{quality_metrics}, \code{status} and
#'   \code{message}.
#' @keywords internal
.bss_run_cascade_stages <- function(detector_names, n_energy_bins, E_MeV,
                                     sensitivities, cc_icrp116,
                                     save_result_callback, readings,
                                     calculate_errors = FALSE,
                                     save_result = FALSE) {
    stages <- .bss_default_cascade_stages()

    # _build_response_matrix(): stack every detector sensitivity.
    A_stack <- do.call(rbind, lapply(detector_names,
                                      function(d) as.numeric(sensitivities[[d]])))
    measured <- as.numeric(readings[detector_names])

    current_spectrum <- NULL
    intermediate_results <- list()
    convergence_history <- list()
    stage_spectra <- list()
    stages_run <- 0L
    method_sequence <- character(0L)

    for (stage_idx in seq_along(stages)) {
        stage <- stages[[stage_idx]]
        method_name <- stage$method
        unfold_func <- .bss_resolve_method(.bss_cascade_dispatch, method_name)
        if (is.null(unfold_func)) next   # method not found, skipping

        params <- stage$params
        params$save_result <- save_result

        if (!is.null(current_spectrum) && isTRUE(stage$use_as_initial)) {
            if (length(current_spectrum) == n_energy_bins) {
                params$initial_spectrum <- current_spectrum
            }
        }

        if (!is.null(current_spectrum) && isTRUE(stage$use_as_prior)) {
            if (length(current_spectrum) == n_energy_bins) {
                if (method_name %in% c("bayes", "bayes_spline")) {
                    params$initial_spectrum <- current_spectrum
                } else if ("reference_spectrum" %in%
                           names(formals(unfold_func))) {
                    params$reference_spectrum <- current_spectrum
                }
            }
        }

        if (!is.null(stage$max_iterations)) {
            cur <- params$max_iterations
            params$max_iterations <- if (!is.null(cur)) {
                min(as.numeric(cur), as.numeric(stage$max_iterations))
            } else stage$max_iterations
        }

        params$calculate_errors <- (stage_idx == length(stages)) &&
            isTRUE(calculate_errors)

        args <- c(list(detector_names = detector_names,
                        n_energy_bins = n_energy_bins, E_MeV = E_MeV,
                        sensitivities = sensitivities,
                        cc_icrp116 = cc_icrp116,
                        save_result_callback = save_result_callback,
                        readings = readings), params)
        result <- tryCatch(do.call(unfold_func, args),
                            error = function(e) NULL)
        if (is.null(result)) next        # stage failed, continue cascade
        stages_run <- stages_run + 1L
        method_sequence <- c(method_sequence, method_name)

        spec <- result$spectrum
        if (!is.null(spec)) {
            current_spectrum <- as.numeric(spec)
            reconstructed <- as.numeric(A_stack %*% current_spectrum)
            metrics <- .bss_compute_quality_metrics(
                current_spectrum, reconstructed, measured, E_MeV)
            convergence_history[[length(convergence_history) + 1L]] <-
                c(list(stage = stage_idx - 1L, method = method_name), metrics)
            stage_spectra[[length(stage_spectra) + 1L]] <- current_spectrum

            if (isTRUE(stage$store_intermediate)) {
                key <- paste0("stage_", stage_idx - 1L, "_", method_name)
                intermediate_results[[key]] <- list(
                    spectrum = current_spectrum, metrics = metrics)
            }

            if (!is.null(stage$quality_threshold) &&
                metrics$overall_quality >= stage$quality_threshold) {
                break   # quality threshold met, stop the cascade
            }
        }
    }

    if (is.null(current_spectrum)) {
        return(list(spectrum = NULL, stages_run = stages_run,
                    intermediate_results = intermediate_results,
                    quality_metrics = list(),
                    method_sequence = method_sequence,
                    convergence_history = convergence_history,
                    cascade_spectra = list(),
                    status = "ERROR", message = "No successful stages"))
    }

    reconstructed <- as.numeric(A_stack %*% current_spectrum)
    final_metrics <- .bss_compute_quality_metrics(
        current_spectrum, reconstructed, measured, E_MeV)

    list(spectrum = current_spectrum, stages_run = stages_run,
         intermediate_results = intermediate_results,
         quality_metrics = final_metrics,
         method_sequence = method_sequence,
         convergence_history = convergence_history,
         cascade_spectra = stage_spectra,
         status = "OK",
         message = sprintf("Successfully completed %d cascade stages",
                            stages_run))
}

#' Pool-level composite driver mirroring the module-level unfold_composite()
#' in bssunfold/core/unfold_composite.py: run the general method pool with
#' detector defaults, drop invalid outputs, combine with confidence-weighted
#' averaging (base weight * mean cosine similarity to the other members).
#' @inheritParams run_unfolding
#' @param n_methods Integer; number of leading general methods to run.
#'   Default 5.
#' @return List with the combined \code{spectrum}, \code{successful_methods},
#'   \code{consistency}, \code{weights}, \code{individual_spectra},
#'   \code{composite_spectra}, \code{status} and \code{message}.
#' @keywords internal
.bss_run_composite_pool <- function(detector_names, n_energy_bins, E_MeV,
                                     sensitivities, cc_icrp116,
                                     save_result_callback, readings,
                                     save_result = FALSE, n_methods = 5L) {
    candidates <- .bss_general_methods[seq_len(min(n_methods,
                                                    length(.bss_general_methods)))]

    names_ok <- character(0L)
    spectra_list <- list()
    messages <- list()

    for (name in candidates) {
        # Python: qpsolvers raises ImportError when the backend package is
        # missing, so the dispatcher records the failure and skips it.
        if (name == "qpsolvers" && !.bss_qp_backend_available()) {
            messages[[name]] <- "ImportError: qpsolvers backend unavailable"
            next
        }
        unfold_func <- .bss_resolve_method(.bss_composite_dispatch, name,
                                           fallback = FALSE)
        if (is.null(unfold_func)) {
            messages[[name]] <- "unknown method"
            next
        }
        args <- list(detector_names = detector_names,
                      n_energy_bins = n_energy_bins, E_MeV = E_MeV,
                      sensitivities = sensitivities, cc_icrp116 = cc_icrp116,
                      save_result_callback = save_result_callback,
                      readings = readings, save_result = save_result)
        result <- tryCatch(do.call(unfold_func, args),
                            error = function(e) {
                                conditionMessage(e)
                            })
        if (is.character(result)) {
            messages[[name]] <- result
            next
        }
        spec <- if (is.list(result)) result$spectrum else NULL
        if (is.null(spec) || any(is.nan(as.numeric(spec))) ||
            sum(as.numeric(spec)) <= 0) {
            messages[[name]] <- "invalid output"
            next
        }
        names_ok <- c(names_ok, name)
        spectra_list[[length(spectra_list) + 1L]] <- as.numeric(spec)
    }

    if (length(spectra_list) == 0L) {
        return(list(spectrum = NULL, successful_methods = character(0L),
                    consistency = 0, weights = list(),
                    individual_spectra = list(),
                    status = "ERROR",
                    message = "No method succeeded"))
    }

    nm <- length(names_ok)
    stacked <- do.call(rbind, spectra_list)

    combined <- numeric(ncol(stacked))
    total_weight <- 0
    used_weights <- numeric(nm)
    names(used_weights) <- names_ok
    for (i in seq_len(nm)) {
        others <- spectra_list[-i]
        conf <- .bss_confidence_weight(stacked[i, ], others)
        base <- unname(.bss_composite_base_weights[names_ok[i]])
        if (is.na(base)) base <- 1.0
        w <- base * conf
        combined <- combined + w * stacked[i, ]
        total_weight <- total_weight + w
        used_weights[i] <- w
    }
    if (total_weight > 0) combined <- combined / total_weight

    consistency <- 0
    n_pairs <- 0L
    for (i in seq_len(nm - 1L)) {
        for (j in seq.int(i + 1L, nm)) {
            s <- .bss_cosine_similarity(stacked[i, ], stacked[j, ])
            if (is.finite(s)) {
                consistency <- consistency + s
                n_pairs <- n_pairs + 1L
            }
        }
    }
    consistency <- if (n_pairs > 0L) consistency / n_pairs else 0

    list(spectrum = combined, successful_methods = names_ok,
         consistency = consistency, weights = as.list(used_weights),
         individual_spectra = spectra_list,
         composite_spectra = stacked,
         status = "OK",
         message = sprintf("Combined %d/%d methods with consistency %.3f",
                            nm, length(candidates), consistency))
}

# --- exported solver-level helpers -------------------------------------

#' Solve by ensemble (combination of multiple solvers)
#'
#' R port of \code{solve_ensemble} in
#' \code{bssunfold/src/bssunfold/core/unfold_ensemble.py}. Runs each solver
#' on \code{(A, b, x0)}, clips every member spectrum at zero, then combines
#' them according to \code{combination}. For \code{"weighted_average"}, a
#' \code{NULL} \code{weights} argument derives inverse-residual weights
#' \code{1 / ||A x_i - b||} (normalised to sum 1), exactly like Python.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param solvers List of solver functions, each with signature
#'   \code{solver(A, b, x0, ...)}.
#' @param weights Optional numeric vector of solver weights. Default
#'   \code{NULL} = inverse-residual weights (as in Python).
#' @param combination Combination strategy: \code{"weighted_average"},
#'   \code{"median"}, \code{"trimmed_mean"} or \code{"best_residual"}.
#' @param trim_fraction Fraction of extreme values to discard per bin for
#'   \code{"trimmed_mean"}.
#' @param ... Extra arguments forwarded to each solver.
#' @return A list \code{list(spectrum, iterations, converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' solvers <- list(solve_mlem, solve_gravel)
#' r <- solve_ensemble(A, b, rep(1, 3), solvers,
#'                     max_iterations = 50L)
solve_ensemble <- function(A, b, x0, solvers, weights = NULL,
                            combination = "weighted_average",
                            trim_fraction = 0.2, ...) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    n <- ncol(A)
    if (length(solvers) == 0L) stop("solvers must be a non-empty list")
    valid_combinations <- c("weighted_average", "median", "trimmed_mean",
                             "best_residual")
    combination <- match.arg(combination, valid_combinations)

    dots <- list(...)
    spectra <- list()
    residuals <- numeric(0L)
    names_ok <- character(0L)

    for (idx in seq_along(solvers)) {
        member_name <- paste0("method_", idx - 1L)
        res <- tryCatch({
            call_args <- c(list(A = A, b = b, x0 = x0), dots)
            keep <- !duplicated(names(call_args), fromLast = TRUE)
            do.call(solvers[[idx]], call_args[keep])
        }, error = function(e) {
            warning(sprintf("Ensemble method %s failed: %s", member_name,
                             conditionMessage(e)))
            NULL
        })
        if (is.null(res)) next
        spec <- if (is.list(res) && !is.null(res$spectrum)) res$spectrum
                else as.numeric(res)
        spec <- pmax(as.numeric(spec), 0)
        spectra[[length(spectra) + 1L]] <- spec
        residuals[length(residuals) + 1L] <-
            sqrt(sum((as.numeric(A %*% spec) - b)^2))
        names_ok <- c(names_ok, member_name)
    }

    if (length(spectra) == 0L) stop("All ensemble methods failed")

    spectra_arr <- do.call(rbind, spectra)
    w_used <- NULL
    info_str <- ""

    if (combination == "best_residual") {
        best_idx <- which.min(residuals)[1L]
        spectrum <- spectra_arr[best_idx, ]
        info_str <- sprintf("best=%s (res=%.4e)", names_ok[best_idx],
                             residuals[best_idx])
    } else if (combination == "median") {
        spectrum <- apply(spectra_arr, 2L, stats::median)
        info_str <- sprintf("median of %d methods", length(spectra))
    } else if (combination == "trimmed_mean") {
        nm <- nrow(spectra_arr)
        k <- max(1L, as.integer(floor(trim_fraction * nm)))
        spectrum <- if (k < nm) {
            sorted <- apply(spectra_arr, 2L, sort)
            rows <- seq.int(k + 1L, nm - k)
            if (length(rows) > 0L) {
                colMeans(sorted[rows, , drop = FALSE])
            } else {
                rep(NaN, n)
            }
        } else {
            colMeans(spectra_arr)
        }
        info_str <- sprintf("trimmed_mean (trim=%s) of %d methods",
                             format(trim_fraction), length(spectra))
    } else {  # weighted_average
        if (is.null(weights)) {
            w <- 1 / (residuals + 1e-30)
            w <- w / sum(w)
        } else {
            w <- as.numeric(weights)
            w <- w / sum(w)
        }
        w_used <- w
        spectrum <- as.numeric(w %*% spectra_arr)
        info_str <- sprintf("weighted_average of %d methods",
                             length(spectra))
    }

    list(spectrum = as.numeric(spectrum),
         iterations = as.integer(length(spectra)),
         converged = TRUE,
         ensemble_spectra = spectra_arr,
         ensemble_weights = if (!is.null(w_used)) w_used else NULL,
         combination = combination,
         method_names = names_ok,
         residuals = residuals,
         info_str = info_str)
}

#' Solve by cascade (sequential solvers, each starting from previous result)
#'
#' @inheritParams solve_ensemble
#' @param solver_kwargs_list Optional list of per-solver kwargs. Each entry
#'   is a list forwarded to the corresponding solver.
#' @param ... Extra arguments forwarded to each solver.
#' @return A list \code{list(spectrum, iterations = length(solvers),
#'   converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' solvers <- list(solve_mlem, solve_landweber)
#' r <- solve_cascade(A, b, rep(1, 3), solvers,
#'                     max_iterations = 50L)
solve_cascade <- function(A, b, x0, solvers, solver_kwargs_list = NULL, ...) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    n <- ncol(A)
    if (length(solvers) == 0L) stop("solvers must be a non-empty list")
    if (!is.null(solver_kwargs_list) &&
        length(solver_kwargs_list) != length(solvers)) {
        stop("solver_kwargs_list length must match solvers length")
    }
    x_cur <- x0
    cascade_spectra <- vector("list", length(solvers))
    for (i in seq_along(solvers)) {
        kwargs <- if (!is.null(solver_kwargs_list)) solver_kwargs_list[[i]]
                   else list()
        args <- c(list(A = A, b = b, x0 = x_cur), kwargs, list(...))
        res <- do.call(solvers[[i]], args)
        spec <- if (is.list(res) && !is.null(res$spectrum)) res$spectrum
                else as.numeric(res)
        x_cur <- pmax(spec, 0)
        cascade_spectra[[i]] <- x_cur
    }
    list(spectrum = as.numeric(x_cur),
         iterations = as.integer(length(solvers)),
         converged = TRUE,
         cascade_spectra = cascade_spectra)
}

#' Solve by composite (pick best residual norm from multiple solvers)
#'
#' @inheritParams solve_ensemble
#' @param ... Extra arguments forwarded to each solver.
#' @return A list \code{list(spectrum, iterations = which_best,
#'   converged = TRUE, best_solver_index = which_best)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' solvers <- list(solve_mlem, solve_gravel, solve_sandii)
#' r <- solve_composite(A, b, rep(1, 3), solvers, max_iterations = 50L)
solve_composite <- function(A, b, x0, solvers, ...) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b); x0 <- as.numeric(x0)
    n <- ncol(A)
    if (length(solvers) == 0L) stop("solvers must be a non-empty list")
    spectra <- matrix(0.0, nrow = length(solvers), ncol = n)
    residuals <- numeric(length(solvers))
    for (i in seq_along(solvers)) {
        res <- do.call(solvers[[i]], list(A = A, b = b, x0 = x0, ...))
        spec <- if (is.list(res) && !is.null(res$spectrum)) res$spectrum
                else as.numeric(res)
        spec <- pmax(spec, 0)
        spectra[i, ] <- spec
        residuals[i] <- sqrt(sum((as.numeric(A %*% spec) - b)^2))
    }
    best <- which.min(residuals)
    list(spectrum = spectra[best, ],
         iterations = as.integer(best),
         converged = TRUE,
         best_solver_index = as.integer(best),
         composite_spectra = spectra,
         composite_residuals = residuals)
}

# --- Detector-facing wrappers ---

#' Wrapper around \code{\link{solve_ensemble}} for the unified workflow.
#'
#' Mirrors \code{unfold_ensemble} from
#' \code{bssunfold/core/unfold_ensemble.py}: when \code{solvers} is not
#' supplied, the Python default ensemble (MLEM, Bayes, Landweber, CGLS,
#' GRAVEL; \code{max_iterations = 200}, \code{tolerance = 1e-4}) is used and
#' the members are combined with inverse-residual weights.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_ensemble
#' @param combination Combination strategy (default \code{"weighted_average"}).
#' @param trim_fraction Trim fraction for \code{trimmed_mean}.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_ensemble <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              solvers = list(solve_mlem, solve_gravel),
                              weights = NULL, ...,
                              combination = "weighted_average",
                              trim_fraction = 0.2,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    if (!is.null(random_state)) set.seed(random_state)

    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    x0 <- if (is.null(initial_spectrum)) rep(0.5, n_energy_bins)
           else as.numeric(initial_spectrum)

    if (missing(solvers)) solvers <- .bss_default_ensemble_solvers()

    res <- solve_ensemble(A, b, x0, solvers, weights = weights,
                           combination = combination,
                           trim_fraction = trim_fraction, ...)
    spectrum <- res$spectrum
    computed_readings <- as.numeric(A %*% spectrum)
    residual <- b - computed_readings
    result <- list(
        energy = E_MeV,
        spectrum = spectrum,
        spectrum_absolute = spectrum,
        effective_readings = stats::setNames(computed_readings, selected),
        residual = residual,
        residual_norm = sqrt(sum(residual^2)),
        method = "Ensemble",
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = 0L,
        converged = res$converged,
        ensemble_spectra = res$ensemble_spectra,
        ensemble_weights = res$ensemble_weights,
        combination = res$combination,
        method_names = res$method_names,
        member_residuals = res$residuals,
        info_str = res$info_str
    )

    # Monte-Carlo uncertainty (mirrors the calculate_errors branch).
    if (isTRUE(calculate_errors)) {
        if (!is.null(random_state)) set.seed(random_state)
        spectra_mc <- list()
        for (i in seq_len(as.integer(n_montecarlo))) {
            b_pert <- pmax(b * (1 + noise_level * stats::rnorm(length(b))), 0)
            x_mc <- tryCatch(
                solve_ensemble(A, b_pert, x0, solvers, weights = weights,
                                combination = combination,
                                trim_fraction = trim_fraction, ...)$spectrum,
                error = function(e) NULL)
            if (!is.null(x_mc)) spectra_mc[[length(spectra_mc) + 1L]] <- x_mc
        }
        if (length(spectra_mc) > 0L) {
            mc_arr <- do.call(rbind, spectra_mc)
            result$spectrum_uncertainty <- apply(mc_arr, 2L, function(v)
                sqrt(mean((v - mean(v))^2)))
            result$calculate_errors <- TRUE
            result$n_montecarlo <- length(spectra_mc)
        }
    }

    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }
    result
}

#' Wrapper around the stage-based cascade for the unified workflow.
#'
#' Mirrors \code{unfold_cascade} from
#' \code{bssunfold/core/unfold_cascade.py}: by default it runs the
#' "general" cascade (TSVD with k = 15 discrepancy, early stop at quality
#' 0.3 -> MLEM with 150 iterations seeded by TSVD -> Bayes-Spline with
#' spline_smooth = 0.3 seeded by MLEM as initial guess and prior). Passing
#' an explicit \code{solvers} list keeps the legacy solver-level cascade.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_cascade
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_cascade <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              solvers = list(solve_mlem, solve_gravel),
                              solver_kwargs_list = NULL, ...,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    x0 <- if (is.null(initial_spectrum)) rep(0.5, n_energy_bins)
           else as.numeric(initial_spectrum)

    if (missing(solvers)) {
        res <- .bss_run_cascade_stages(detector_names, n_energy_bins, E_MeV,
                                        sensitivities, cc_icrp116,
                                        save_result_callback, readings,
                                        calculate_errors = calculate_errors,
                                        save_result = save_result)
    } else {
        stage_res <- solve_cascade(A, b, x0, solvers,
                                    solver_kwargs_list = solver_kwargs_list,
                                    ...)
        res <- list(spectrum = stage_res$spectrum,
                    stages_run = length(solvers),
                    intermediate_results = list(),
                    quality_metrics = list(),
                    method_sequence = paste0("solver_", seq_len(length(solvers))),
                    convergence_history = list(),
                    cascade_spectra = stage_res$cascade_spectra,
                    status = "OK", message = "legacy solver cascade")
    }

    spectrum <- res$spectrum
    if (is.null(spectrum)) {
        result <- list(
            energy = E_MeV,
            spectrum = NULL,
            spectrum_absolute = NULL,
            method = "Cascade",
            doserates = numeric(0L),
            iterations = as.integer(res$stages_run),
            converged = FALSE,
            cascade_spectra = res$cascade_spectra,
            stages_run = res$stages_run,
            method_sequence = res$method_sequence,
            quality_metrics = res$quality_metrics,
            convergence_history = res$convergence_history,
            status = res$status,
            message = res$message
        )
        if (isTRUE(save_result) && is.function(save_result_callback)) {
            save_result_callback(result)
        }
        return(result)
    }

    computed_readings <- as.numeric(A %*% spectrum)
    residual <- b - computed_readings
    result <- list(
        energy = E_MeV,
        spectrum = spectrum,
        spectrum_absolute = spectrum,
        effective_readings = stats::setNames(computed_readings, selected),
        residual = residual,
        residual_norm = sqrt(sum(residual^2)),
        method = "Cascade",
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = as.integer(res$stages_run),
        converged = TRUE,
        cascade_spectra = res$cascade_spectra,
        stages_run = res$stages_run,
        method_sequence = res$method_sequence,
        quality_metrics = res$quality_metrics,
        convergence_history = res$convergence_history,
        status = res$status,
        message = res$message
    )
    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }
    result
}

#' Wrapper around the confidence-weighted method pool for the unified
#' workflow.
#'
#' Mirrors \code{unfold_composite} from
#' \code{bssunfold/core/unfold_composite.py}: by default it runs the general
#' method pool (TSVD, MLEM, CVXPY, QPsolvers, Bayes-Spline with detector
#' defaults), skips members that fail or return invalid spectra, and
#' combines the rest with confidence-weighted averaging. Passing an
#' explicit \code{solvers} list keeps the legacy best-residual behaviour.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_composite
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_composite <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              solvers = list(solve_mlem, solve_gravel),
                              ...,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    sys <- .build_system(readings, detector_names, sensitivities)
    A <- sys$A; b <- sys$b; selected <- sys$selected
    x0 <- if (is.null(initial_spectrum)) rep(0.5, n_energy_bins)
           else as.numeric(initial_spectrum)

    if (missing(solvers)) {
        res <- .bss_run_composite_pool(detector_names, n_energy_bins, E_MeV,
                                        sensitivities, cc_icrp116,
                                        save_result_callback, readings,
                                        save_result = save_result)
    } else {
        best_res <- solve_composite(A, b, x0, solvers, ...)
        res <- list(spectrum = best_res$spectrum,
                    successful_methods = paste0("solver_",
                        seq_len(ncol(best_res$composite_spectra))),
                    consistency = NA_real_,
                    weights = list(),
                    individual_spectra = list(),
                    composite_spectra = best_res$composite_spectra,
                    composite_residuals = best_res$composite_residuals,
                    best_solver_index = best_res$best_solver_index,
                    status = "OK", message = "legacy solver composite")
    }

    spectrum <- res$spectrum
    if (is.null(spectrum)) {
        result <- list(
            energy = E_MeV,
            spectrum = NULL,
            spectrum_absolute = NULL,
            method = "Composite",
            doserates = numeric(0L),
            iterations = 0L,
            converged = FALSE,
            successful_methods = res$successful_methods,
            consistency = res$consistency,
            composite_weights = res$weights,
            composite_spectra = res$individual_spectra,
            status = res$status,
            message = res$message
        )
        if (isTRUE(save_result) && is.function(save_result_callback)) {
            save_result_callback(result)
        }
        return(result)
    }

    computed_readings <- as.numeric(A %*% spectrum)
    residual <- b - computed_readings
    result <- list(
        energy = E_MeV,
        spectrum = spectrum,
        spectrum_absolute = spectrum,
        effective_readings = stats::setNames(computed_readings, selected),
        residual = residual,
        residual_norm = sqrt(sum(residual^2)),
        method = "Composite",
        doserates = if (!is.null(cc_icrp116))
                        calculate_dose_rates(spectrum, cc_icrp116)
                    else numeric(0L),
        iterations = length(res$successful_methods),
        converged = TRUE,
        successful_methods = res$successful_methods,
        consistency = res$consistency,
        composite_weights = res$weights,
        composite_spectra = res$composite_spectra,
        individual_spectra = res$individual_spectra,
        status = res$status,
        message = res$message
    )
    if (!is.null(res$composite_residuals)) {
        result$composite_residuals <- res$composite_residuals
        result$best_solver_index <- res$best_solver_index
    }
    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }
    result
}
