#' NSDUAZ unfolding (catalogue + SPUNIT)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_nsduaz.py}.
#' NSDUAZ (Ortiz-Rodriguez & Vega-Carrillo, 2012) selects the initial
#' guess spectrum from a catalogue of standard neutron spectra by a
#' statistical test on count-rate ratios relative to the reference sphere,
#' then runs SPUNIT (BUNKI) iteration with a ~1% convergence tolerance.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n). Typically selected by
#'   \code{\link{select_catalogue_initial}}.
#' @param smoothing Numeric; 3-point smoothing factor. Default 0.1.
#' @param max_iterations Positive integer; default 1000.
#' @param tolerance Positive numeric; default 0.01.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_nsduaz(A, b, rep(1, 3), max_iterations = 50)
solve_nsduaz <- function(A, b, x0, smoothing = 0.1,
                          max_iterations = 1000L, tolerance = 0.01) {
    v <- validate_system(A, b, x0 = x0, max_iterations = max_iterations,
                         tolerance = tolerance)
    solve_bunki(A = v$A, b = v$b, x0 = v$x0, smoothing = smoothing,
                max_iterations = max_iterations, tolerance = tolerance)
}

#' Built-in catalogue of standard neutron spectra
#'
#' Returns a list of normalised analytic standard spectra on the energy
#' grid \code{E_MeV}: \code{ambe} (241Am/9Be alpha-n), \code{cf252}
#' (Watt fission), and \code{reactor} (thermal + 1/E + fast).
#'
#' @param E_MeV Numeric energy grid in MeV.
#' @return A named list of numeric vectors.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' cat <- builtin_catalogue(E)
#' names(cat)
builtin_catalogue <- function(E_MeV) {
    E <- pmax(as.numeric(E_MeV), 1e-9)
    # Watt fission spectrum (252Cf): exp(-E/a) sinh(sqrt(b*E))
    watt <- exp(-E / 1.025) * sinh(sqrt(2.926 * E))
    if (sum(watt) <= 0) watt <- rep(1, length(E))
    watt <- watt / sum(watt)
    # 241Am/9Be: evaporation continuum + 4.2 MeV peak
    ambe <- exp(-E / 2.0) + 3.5 * exp(-0.5 * ((E - 4.2) / 1.2)^2)
    if (sum(ambe) <= 0) ambe <- rep(1, length(E))
    ambe <- ambe / sum(ambe)
    # Reactor-like: thermal Maxwellian + 1/E + fast fission
    kT <- 0.0253e-6  # 0.0253 eV in MeV
    thermal <- (E / kT) * exp(-E / kT)
    epithermal <- ifelse(E > 1e-6, 1.0 / pmax(E, 1e-9), 0.0)
    reactor <- 1e-3 * thermal + 0.1 * epithermal + watt
    if (sum(reactor) <= 0) reactor <- rep(1, length(E))
    reactor <- reactor / sum(reactor)
    list(ambe = ambe, cf252 = watt, reactor = reactor)
}

#' Select initial spectrum from catalogue
#'
#' The experimental count rates are normalised to the reading of the
#' reference sphere (typically the 20.32 cm / 8-inch sphere) and compared
#' with the relative count-rate pattern predicted by each catalogue
#' spectrum folded through the response matrix. The entry minimising the
#' weighted chi-square of the relative ratios is chosen and rescaled so
#' that its predicted reference reading matches the measured one.
#'
#' @param readings Named numeric vector of detector readings.
#' @param detector_names Character vector of detector names.
#' @param sensitivities Named list of numeric sensitivity vectors.
#' @param catalogue Optional named list of spectra on the energy grid.
#'   Default \code{NULL} = \code{\link{builtin_catalogue}}.
#' @param reference_name Optional character; name of the reference detector.
#'   Default \code{NULL} = auto-detect the 8-inch / 20-cm sphere.
#' @param E_MeV Optional energy grid. Default \code{NULL} = use length of
#'   sensitivities to build a representative log grid.
#' @return A list \code{list(spectrum, label)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' sens <- setNames(lapply(1:5, function(i) A[i, ]),
#'                   c("0in","3in","5in","8in","12in"))
#' readings <- c("0in"=100, "3in"=80, "5in"=60, "8in"=40, "12in"=10)
#' sel <- select_catalogue_initial(readings, names(sens), sens, E_MeV = E)
#' sel$label
select_catalogue_initial <- function(readings, detector_names, sensitivities,
                                        catalogue = NULL,
                                        reference_name = NULL,
                                        E_MeV = NULL) {
    selected <- detector_names[detector_names %in% names(readings)]
    if (length(selected) == 0L) {
        stop("No detector readings available for catalogue selection")
    }
    b <- as.numeric(readings[selected])
    A <- do.call(rbind, lapply(selected, function(n) as.numeric(sensitivities[[n]])))
    if (!is.null(reference_name)) {
        if (!(reference_name %in% selected)) {
            stop("reference_name '", reference_name,
                 "' is not present in readings")
        }
        ref_idx <- which(selected == reference_name)
    } else {
        ref_idx <- .nsduaz_find_reference_index(selected, A)
    }
    n_bins <- ncol(A)
    if (is.null(catalogue)) {
        if (!is.null(E_MeV) && length(E_MeV) == n_bins) {
            catalogue <- builtin_catalogue(E_MeV)
        } else {
            E_rep <- 10^seq(log10(1e-9), log10(1e2), length.out = n_bins)
            catalogue <- builtin_catalogue(E_rep)
        }
    }
    b_ref <- b[ref_idx]
    if (b_ref <= 0) {
        stop("Reference sphere reading must be strictly positive")
    }
    r_ratio <- b / b_ref
    best_label <- NULL; best_chi <- Inf; best_scale <- 1.0; best_spec <- NULL
    for (label in names(catalogue)) {
        spec <- as.numeric(catalogue[[label]])
        if (length(spec) != n_bins) {
            stop("Catalogue spectrum '", label, "' has length ", length(spec),
                 ", expected ", n_bins)
        }
        if (!any(spec > 0)) next
        c <- as.numeric(A %*% pmax(spec, 0))
        c_ref <- c[ref_idx]
        if (c_ref <= 0) next
        s_ratio <- c / c_ref
        denom <- pmax(s_ratio, 1e-12)
        chi <- sum(((r_ratio - s_ratio) / denom)^2)
        if (chi < best_chi) {
            best_chi <- chi; best_label <- label
            best_scale <- b_ref / c_ref; best_spec <- spec
        }
    }
    if (is.null(best_spec)) {
        stop("Catalogue is empty or has no usable spectrum")
    }
    list(spectrum = pmax(best_scale * best_spec, 0.0), label = best_label)
}

.nsduaz_find_reference_index <- function(detector_names, A) {
    for (i in seq_along(detector_names)) {
        nm <- tolower(detector_names[i])
        if (any(grepl(x = nm, pattern = "20.32|20in|8in|8 in"))) return(i)
    }
    which.max(rowSums(abs(A)))
}

#' Wrapper around \code{\link{solve_nsduaz}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_nsduaz
#' @param catalogue Optional named list passed to
#'   \code{\link{select_catalogue_initial}}.
#' @param reference_name Optional character passed to
#'   \code{\link{select_catalogue_initial}}.
#' @return A result list as produced by \code{\link{run_unfolding}} with
#'   an extra \code{catalogue_label} field.
#' @export
unfold_nsduaz <- function(detector_names, n_energy_bins, E_MeV,
                            sensitivities, cc_icrp116, save_result_callback,
                            readings, initial_spectrum = NULL,
                            smoothing = 0.1, max_iterations = 1000L,
                            tolerance = 0.01,
                            catalogue = NULL, reference_name = NULL,
                            calculate_errors = FALSE,
                            noise_level = 0.01, n_montecarlo = 100L,
                            save_result = FALSE, random_state = NULL) {
    if (is.null(initial_spectrum)) {
        sel <- select_catalogue_initial(
            readings = readings, detector_names = detector_names,
            sensitivities = sensitivities, catalogue = catalogue,
            reference_name = reference_name, E_MeV = E_MeV
        )
        x0_init <- sel$spectrum
        catalogue_label <- sel$label
    } else {
        x0_init <- as.numeric(initial_spectrum)
        catalogue_label <- NA_character_
    }
    result <- run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = x0_init,
        default_initial = rep(1.0, n_energy_bins),
        solve_func = make_solve_wrapper(solve_nsduaz,
                                         smoothing = smoothing,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "NSDUAZ",
        extra_output = list(smoothing = smoothing,
                            catalogue_label = catalogue_label),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
    result
}
