#' Method-name pipeline unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_combined.py}.
#' Runs a sequential multi-method pipeline: each stage is specified by a
#' \code{list(method = "mlem", params = list(max_iterations = 500L))} entry.
#' The first stage starts from the initial spectrum; every subsequent stage
#' starts from the result of the previous stage. This is the Python package's
#' "combined approach" for chaining multiple unfolding methods.
#'
#' Recognized method names (dispatched to \code{solve_*} functions):
#' "mlem", "mlem_stop", "gravel", "doroshenko", "kaczmarz",
#' "randomized_kaczmarz", "landweber", "bayes", "bayes_spline", "maxed",
#' "imaxed", "amaxed", "sandii", "bunki", "bunkiut", "ferdor", "tsvd",
#' "cgls", "gks", "lanczos", "fista", "osem", "mapem", "bsrem", "sart",
#' "crystal_ball", "rfsp_jul", "staysl", "ensemble", "cascade", "composite",
#' "pspline_reml", "amg", "pdhg", "douglas_rachford".
#'
#' @name combined-methods
NULL

# Lookup table: canonical method name -> internal representative solver
# function name string (all exported as solve_*).
.combined_solver_registry <- function() {
    c(
        mlem = "solve_mlem",
        mlem_stop = "solve_mlem_stop",
        gravel = "solve_gravel",
        doroshenko = "solve_doroshenko",
        kaczmarz = "solve_kaczmarz",
        randomized_kaczmarz = "solve_randomized_kaczmarz",
        landweber = "solve_landweber",
        bayes = "solve_bayes",
        bayes_spline = "solve_bayes_spline",
        maxed = "solve_maxed",
        imaxed = "solve_imaxed",
        amaxed = "solve_amaxed",
        tsvd = "solve_tsvd",
        cgls = "solve_cgls",
        gks = "solve_gks",
        lanczos = "solve_lanczos",
        fista = "solve_fista",
        osem = "solve_osem",
        mapem = "solve_mapem",
        bsrem = "solve_bsrem",
        sart = "solve_sart",
        crystal_ball = "solve_crystal_ball",
        rfsp_jul = "solve_rfsp_jul",
        staysl = "solve_staysl",
        ensemble = "solve_ensemble",
        cascade = "solve_cascade",
        composite = "solve_composite",
        pspline_reml = "solve_pspline_reml",
        amg = "solve_amg",
        pdhg = "solve_pdhg",
        douglas_rachford = "solve_douglas_rachford",
        sandii = "solve_sandii",
        bunki = "solve_bunki",
        bunkiut = "solve_bunkiut",
        ferdor = "solve_ferdor",
        reconst = "solve_reconst",
        statreg = "solve_statreg",
        express = "solve_express",
        directed_divergence = "solve_directed_divergence"
    )
}

#' Solve by named-method pipeline
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param pipeline List of stages; each stage is a list with fields
#'   \code{method} (character) and optional \code{params} (list of extra
#'   arguments for the stage solver).
#' @param ... Extra arguments forwarded to every stage solver.
#' @return A list \code{list(spectrum, iterations, converged,
#'   pipeline_spectra)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_combined(A, b, rep(1, 3),
#'                     pipeline = list(
#'                         list(method = "mlem",
#'                              params = list(max_iterations = 30L)),
#'                         list(method = "landweber",
#'                              params = list(max_iterations = 10L))))
solve_combined <- function(A, b, x0 = NULL, pipeline = NULL, ...) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (is.null(x0)) x0 <- rep(mean(b) / max(sum(A), 1e-10) * n, n)
    x0 <- as.numeric(x0)
    if (length(pipeline) == 0L) {
        stop("'pipeline' must be a non-empty list of stages.")
    }
    registry <- .combined_solver_registry()
    x_cur <- pmax(x0, 1e-30)
    stage_spectra <- vector("list", length(pipeline))
    converged <- TRUE
    iterations <- 0L
    for (i in seq_along(pipeline)) {
        stage <- pipeline[[i]]
        if (!is.list(stage) || is.null(stage$method)) {
            stop("pipeline[[", i, "]] must be a list with a 'method' field.")
        }
        mname <- tolower(stage$method)
        if (!(mname %in% names(registry))) {
            stop("Unknown pipeline method '", stage$method,
                 "'. Known: ", paste(names(registry), collapse = ", "))
        }
        solver <- get(registry[[mname]], mode = "function")
        params <- if (is.null(stage$params)) list() else stage$params
        args <- c(list(A = A, b = b, x0 = x_cur), params, list(...))
        res <- do.call(solver, args)
        if (is.list(res) && !is.null(res$spectrum)) {
            spec <- as.numeric(res$spectrum)
            if (!is.null(res$converged)) converged <- converged &&
                isTRUE(res$converged)
        } else {
            spec <- as.numeric(res)
        }
        x_cur <- pmax(spec, 0)
        stage_spectra[[i]] <- x_cur
    }
    list(spectrum = x_cur,
         iterations = as.integer(length(pipeline)),
         converged = converged,
         stage_spectra = stage_spectra,
         stage_methods = vapply(pipeline, function(s) as.character(s$method),
                                character(1L)))
}

#' Wrapper around \code{\link{solve_combined}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_combined
#' @param ... Extra arguments forwarded to \code{\link{solve_combined}}.
#' @export
unfold_combined <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            pipeline = list(
                                list(method = "mlem",
                                     params = list(max_iterations = 200L)),
                                list(method = "landweber",
                                     params = list(max_iterations = 20L))),
                            ..., method_name = "Combined",
                            calculate_errors = FALSE,
                            noise_level = 0.01, n_montecarlo = 100L,
                            save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1, n_energy_bins),
        solve_func = solve_combined, solve_kwargs = list(pipeline = pipeline),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
