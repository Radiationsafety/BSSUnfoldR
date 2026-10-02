#' Binned adaptive unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_binned.py}.
#' Bin-wise adaptive unfolding: for every energy bin the best unfolding method
#' is taken from a pre-computed benchmark lookup table and the final spectrum
#' is assembled bin by bin, picking the winning method's value at that bin.
#' This exploits the empirical observation that different unfolding methods
#' excel in different energy regions.
#'
#' The lookup shipped with the Python package (\code{data/bin_lookup.json},
#' 60 bins, top 5 methods per bin) is embedded verbatim below, so no external
#' data file is required.  When the energy grid has a different number of bins
#' the table is used as-is (bins \code{1..n_bins}), exactly like the Python
#' implementation.
#'
#' @section Limitations:
#' The Python original builds lookup tables from benchmark output directories
#' (\code{tools/build_bin_lookup.py}); that analysis needs pandas and
#' \code{.npz} files and is not ported.  \code{\link{build_bin_lookup}} keeps
#' the simple uniform log-energy grouping of the R port instead.
#'
#' @name unfold_binned_lookup
NULL

# ---- Shipped per-bin method ranking ---------------------------------------- #
# Values copied verbatim from the Python package data file
# (bssunfold/src/bssunfold/data/bin_lookup.json, n_bins = 60) so that the R
# port assembles the same per-bin winner. RANKS are 1-based indices into
# METHODS and SCORES are the reported per-bin mean absolute errors.
.BINNED_LOOKUP_METHODS <- c(
    "amaxed_regularization", "bayes", "bayes_spline_regularization", "bayesian_parametric",
    "bsrem", "bunki", "bunkiut", "combined",
    "composite", "crystal_ball", "cvxpy", "docplex",
    "doroshenko", "ferdor", "fruit_like", "genetic",
    "gks", "hybrid_gmres", "hybrid_parametric", "imaxed",
    "interpret", "kaczmarz", "landweber", "lmfit",
    "lmfit_ic", "mapem", "mcmc", "mlem",
    "mlem_odl", "mlem_stop", "mystic_hybrid", "odl_pdhg",
    "osem", "parametric", "parametric2", "parametric_combined",
    "parametric_qpsolvers", "qubo", "reconst", "rfsp_jul",
    "sandii", "sart", "scip", "scipy_direct_method",
    "smt", "statreg", "staysl", "tikhonov_legendre",
    "tikhonov_tv", "tsvd", "zfit"
)
.BINNED_LOOKUP_RANKS <- list(
    c(7L, 28L, 41L, 42L, 6L),   # bin 0
    c(36L, 16L, 17L, 19L, 50L),   # bin 1
    c(16L, 17L, 36L, 19L, 50L),   # bin 1
    c(16L, 10L, 44L, 18L, 50L),   # bin 3
    c(36L, 9L, 25L, 23L, 42L),   # bin 4
    c(36L, 29L, 46L, 25L, 27L),   # bin 5
    c(29L, 21L, 46L, 25L, 36L),   # bin 6
    c(29L, 21L, 9L, 25L, 18L),   # bin 7
    c(21L, 29L, 9L, 18L, 10L),   # bin 8
    c(18L, 29L, 10L, 44L, 9L),   # bin 9
    c(46L, 9L, 18L, 44L, 10L),   # bin 10
    c(29L, 25L, 16L, 9L, 42L),   # bin 11
    c(18L, 10L, 44L, 29L, 25L),   # bin 12
    c(1L, 29L, 25L, 9L, 6L),   # bin 13
    c(29L, 1L, 25L, 6L, 20L),   # bin 14
    c(29L, 2L, 6L, 25L, 1L),   # bin 15
    c(2L, 13L, 29L, 27L, 3L),   # bin 16
    c(2L, 13L, 29L, 3L, 27L),   # bin 17
    c(2L, 27L, 29L, 3L, 13L),   # bin 18
    c(2L, 27L, 29L, 3L, 6L),   # bin 19
    c(2L, 29L, 27L, 6L, 3L),   # bin 20
    c(2L, 29L, 27L, 6L, 3L),   # bin 21
    c(2L, 27L, 29L, 9L, 25L),   # bin 22
    c(25L, 29L, 9L, 2L, 37L),   # bin 23
    c(25L, 29L, 9L, 37L, 34L),   # bin 24
    c(25L, 29L, 18L, 44L, 10L),   # bin 25
    c(25L, 29L, 16L, 34L, 37L),   # bin 26
    c(16L, 25L, 29L, 34L, 37L),   # bin 27
    c(16L, 34L, 29L, 25L, 37L),   # bin 28
    c(16L, 29L, 37L, 34L, 25L),   # bin 29
    c(29L, 37L, 16L, 34L, 25L),   # bin 30
    c(29L, 37L, 16L, 25L, 10L),   # bin 31
    c(29L, 25L, 37L, 16L, 18L),   # bin 32
    c(29L, 25L, 16L, 27L, 37L),   # bin 33
    c(29L, 25L, 46L, 45L, 16L),   # bin 34
    c(25L, 46L, 29L, 23L, 16L),   # bin 35
    c(45L, 25L, 46L, 29L, 16L),   # bin 36
    c(45L, 29L, 16L, 25L, 23L),   # bin 37
    c(25L, 29L, 10L, 44L, 23L),   # bin 38
    c(45L, 18L, 10L, 44L, 23L),   # bin 39
    c(45L, 29L, 16L, 10L, 25L),   # bin 40
    c(45L, 25L, 46L, 39L, 29L),   # bin 41
    c(45L, 18L, 44L, 10L, 25L),   # bin 42
    c(45L, 29L, 25L, 46L, 6L),   # bin 43
    c(45L, 29L, 25L, 46L, 34L),   # bin 44
    c(45L, 25L, 29L, 21L, 46L),   # bin 45
    c(45L, 29L, 46L, 25L, 39L),   # bin 46
    c(45L, 29L, 46L, 25L, 39L),   # bin 47
    c(45L, 25L, 29L, 1L, 39L),   # bin 48
    c(45L, 25L, 21L, 46L, 11L),   # bin 49
    c(45L, 1L, 36L, 25L, 14L),   # bin 50
    c(45L, 25L, 12L, 11L, 14L),   # bin 51
    c(45L, 4L, 13L, 15L, 31L),   # bin 52
    c(14L, 4L, 45L, 15L, 34L),   # bin 53
    c(13L, 14L, 4L, 15L, 34L),   # bin 54
    c(4L, 34L, 15L, 16L, 45L),   # bin 55
    c(14L, 4L, 34L, 15L, 37L),   # bin 56
    c(14L, 38L, 4L, 34L, 15L),   # bin 57
    c(13L, 32L, 38L, 47L, 4L),   # bin 58
    c(32L, 38L, 47L, 4L, 34L)   # bin 59
)
.BINNED_LOOKUP_SCORES <- list(
    c(0.00478276, 0.00478276, 0.00478276, 0.00478276, 0.00478276),
    c(0.00617382, 0.00692595, 0.00747123, 0.00752685, 0.00839962),
    c(0.00751391, 0.00794948, 0.00847423, 0.00851606, 0.00867272),
    c(0.00644984, 0.00686703, 0.00698476, 0.00703179, 0.00787985),
    c(0.01065601, 0.01327324, 0.01411693, 0.01482919, 0.01587195),
    c(0.02505947, 0.02595429, 0.03062535, 0.03128008, 0.03463535),
    c(0.03048539, 0.03139974, 0.03144028, 0.0318068, 0.03192583),
    c(0.02459476, 0.02466655, 0.02508558, 0.02551272, 0.02566129),
    c(0.02288877, 0.02339166, 0.02344035, 0.02353862, 0.02368889),
    c(0.01822504, 0.01833786, 0.01840098, 0.01843203, 0.01879111),
    c(0.01197523, 0.01231975, 0.01361766, 0.01394233, 0.0144592),
    c(0.01108216, 0.01252723, 0.01302215, 0.01335229, 0.01363607),
    c(0.01162239, 0.01164407, 0.01173345, 0.01174847, 0.01197299),
    c(0.00808386, 0.00826448, 0.00919753, 0.00945073, 0.00984584),
    c(0.00689174, 0.00715466, 0.00844061, 0.00880502, 0.00909549),
    c(0.00650001, 0.00775425, 0.00835871, 0.00854649, 0.0086008),
    c(0.00638005, 0.00712555, 0.00719647, 0.00820132, 0.00848765),
    c(0.00575697, 0.00692963, 0.00761288, 0.00814204, 0.00885104),
    c(0.00508393, 0.00700035, 0.00700535, 0.00790755, 0.00825242),
    c(0.00546988, 0.00729633, 0.00780612, 0.0081861, 0.00837667),
    c(0.00603403, 0.0077536, 0.00796929, 0.00839647, 0.00865718),
    c(0.00634339, 0.00789479, 0.00811802, 0.00833489, 0.00896149),
    c(0.00672339, 0.00716722, 0.00753917, 0.0083602, 0.00847965),
    c(0.00585221, 0.00595516, 0.00682346, 0.00683792, 0.00787536),
    c(0.00425451, 0.00461081, 0.00642317, 0.006647, 0.00689847),
    c(0.00344575, 0.00396678, 0.00427944, 0.0044097, 0.00454855),
    c(0.00373139, 0.0039869, 0.00437344, 0.00483072, 0.00512108),
    c(0.00244942, 0.00407611, 0.00421999, 0.00426101, 0.0046919),
    c(0.0023623, 0.00388104, 0.00425582, 0.00435878, 0.0043799),
    c(0.00383676, 0.00439479, 0.00444188, 0.00505631, 0.00564023),
    c(0.00478619, 0.00502639, 0.00626809, 0.00673161, 0.00761941),
    c(0.00521644, 0.0062544, 0.0071636, 0.00775631, 0.00887135),
    c(0.00566453, 0.00770478, 0.00772011, 0.00784889, 0.00801616),
    c(0.00525241, 0.00700363, 0.00791123, 0.00894447, 0.00911363),
    c(0.0054297, 0.00551146, 0.00740844, 0.00805257, 0.00873965),
    c(0.00522115, 0.00574493, 0.00642512, 0.00878453, 0.00972023),
    c(0.00612319, 0.00766479, 0.00784457, 0.00932353, 0.01081992),
    c(0.00576357, 0.00959757, 0.01032019, 0.01055764, 0.01062555),
    c(0.00934573, 0.01058615, 0.01112554, 0.01167803, 0.0118062),
    c(0.00059277, 0.0049544, 0.00604629, 0.00611298, 0.00943516),
    c(0.00085913, 0.01010861, 0.0136153, 0.01385846, 0.01389752),
    c(0.00103983, 0.00992067, 0.01036939, 0.01069273, 0.01099005),
    c(0.00016419, 0.01084931, 0.01106263, 0.01136863, 0.0155784),
    c(0.00010719, 0.01421451, 0.01575523, 0.0194332, 0.02146954),
    c(5.016e-05, 0.01012536, 0.01116138, 0.01735593, 0.01811945),
    c(2.4e-05, 0.01220935, 0.01255788, 0.01682709, 0.0169632),
    c(7.54e-06, 0.01146082, 0.01345257, 0.01377408, 0.01523824),
    c(0.0, 0.00871272, 0.01058755, 0.01085525, 0.01125678),
    c(0.0, 0.00684669, 0.00710567, 0.00975001, 0.01092524),
    c(0.0, 0.00543307, 0.00724935, 0.0073203, 0.00784015),
    c(0.0, 0.00252939, 0.00276888, 0.00305088, 0.00305649),
    c(0.0, 0.00130285, 0.0016108, 0.00167145, 0.00168009),
    c(0.0, 0.0001127, 0.00011362, 0.00012424, 0.00012651),
    c(0.0, 0.0, 0.0, 3.3e-07, 6.7e-07),
    c(0.0, 0.0, 0.0, 0.0, 0.0),
    c(0.0, 0.0, 0.0, 0.0, 0.0),
    c(0.0, 0.0, 0.0, 0.0, 0.0),
    c(0.0, 0.0, 0.0, 0.0, 0.0),
    c(0.0, 0.0, 0.0, 0.0, 0.0),
    c(0.0, 0.0, 0.0, 0.0, 0.0)
)

# Mapping short name -> Detector unfold_* function, mirroring METHOD_DISPATCH.
.BINNED_METHOD_DISPATCH <- c(
    tsvd = "unfold_tsvd", bayes = "unfold_bayes", cvxpy = "unfold_cvxpy",
    statreg = "unfold_statreg", lanczos = "unfold_lanczos", mlem = "unfold_mlem",
    landweber = "unfold_landweber",
    bayes_spline = "unfold_bayes_spline_regularization", gravel = "unfold_gravel",
    qpsolvers = "unfold_qpsolvers", hybrid_parametric = "unfold_hybrid_parametric",
    parametric2 = "unfold_parametric2",
    genetic = "unfold_genetic", interpret = "unfold_interpret",
    maeo_ensemble = "unfold_maeo_ensemble",
    mystic = "unfold_mystic", mystic_hybrid = "unfold_mystic_hybrid", cs = "unfold_cs",
    scip = "unfold_scip", docplex = "unfold_docplex", epic = "unfold_epic",
    kaczmarz = "unfold_kaczmarz", sart = "unfold_sart", osem = "unfold_osem",
    bsrem = "unfold_bsrem", mapem = "unfold_mapem", ferdor = "unfold_ferdor",
    rebunki = "unfold_rebunki", nsduaz = "unfold_nsduaz",
    doroshenko = "unfold_doroshenko",
    sandii = "unfold_sandii", bunki = "unfold_bunki", bunkiut = "unfold_bunkiut",
    reconst = "unfold_reconst", amaxed = "unfold_amaxed",
    amaxed_regularization = "unfold_amaxed_regularization",
    imaxed = "unfold_imaxed", maxed = "unfold_maxed", mlem_odl = "unfold_mlem_odl",
    mlem_stop = "unfold_mlem_stop", cgls = "unfold_cgls", gks = "unfold_gks",
    hybrid_gmres = "unfold_hybrid_gmres",
    tikhonov_legendre = "unfold_tikhonov_legendre", tikhonov_tv = "unfold_tikhonov_tv",
    fista = "unfold_fista", crystal_ball = "unfold_crystal_ball",
    rfsp_jul = "unfold_rfsp_jul",
    staysl = "unfold_staysl", parametric = "unfold_parametric",
    parametric_cvxpy = "unfold_parametric",
    parametric_qpsolvers = "unfold_parametric",
    parametric_combined = "unfold_parametric", lmfit = "unfold_lmfit",
    lmfit_ic = "unfold_lmfit",
    scipy_direct_method = "unfold_scipy_direct_method", qubo = "unfold_qubo",
    zfit = "unfold_zfit", mcmc = "unfold_mcmc",
    bayesian_parametric = "unfold_bayesian_parametric",
    eki = "unfold_eki", maeo = "unfold_maeo", odl_pdhg = "unfold_odl_pdhg",
    odl_douglas_rachford = "unfold_odl_douglas_rachford",
    combined = "unfold_combined", cascade = "unfold_cascade",
    composite = "unfold_composite"
)

# Aliases that map to a base method with fixed extra params.  The keys are
# Detector attribute names, while the lookup values are resolved through
# \code{.BINNED_METHOD_DISPATCH}, so this table is copied verbatim from Python
# (including the fact that it never fires).
.BINNED_ALIASES <- list(
    `unfold_parametric_cvxpy` = list(base = "unfold_parametric",
                                     fixed = list(optimizer = "cvxpy")),
    `unfold_parametric_qpsolvers` = list(base = "unfold_parametric",
                                         fixed = list(optimizer = "qpsolvers")),
    `unfold_parametric_combined` = list(base = "unfold_parametric",
                                        fixed = list(optimizer = "combined")),
    `unfold_lmfit_ic` = list(base = "unfold_lmfit", fixed = list())
)

# Per-method wall-clock timeout, mirroring _run_with_timeout().  R cannot read
# back the previous limits (setTimeLimit() returns NULL before R 4.4), so the
# limit is cleared with Inf when the guarded call returns.
.binned_run_with_timeout <- function(thunk, timeout) {
    if (is.null(timeout) || !is.finite(timeout) || timeout <= 0) {
        return(thunk())
    }
    setTimeLimit(elapsed = as.numeric(timeout), transient = FALSE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE), add = TRUE)
    tryCatch(thunk(), error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("elapsed time limit", msg)) {
            stop("MethodTimeout", call. = FALSE)
        }
        stop(msg, call. = FALSE)
    })
}

# Adapt an unfold_* free function to the (A, b, x0) interface used by
# solve_binned(), mirroring _make_detector_solver().
.make_binned_solver <- function(fn, selected, detector_names, n_energy_bins,
                                E_MeV, sensitivities, cc_icrp116,
                                save_result_callback) {
    force(fn); force(selected)
    function(A_mat, b_vec, x0 = NULL, ...) {
        readings_vec <- stats::setNames(as.numeric(b_vec), selected)
        call_kwargs <- list(...)
        if (!is.null(x0)) {
            call_kwargs <- c(call_kwargs, list(initial_spectrum = as.numeric(x0)))
        }
        accepted <- names(formals(fn))
        if (!("..." %in% accepted)) {
            call_kwargs <- call_kwargs[names(call_kwargs) %in% accepted]
        }
        res <- do.call(fn, c(list(detector_names = detector_names,
                                  n_energy_bins = n_energy_bins,
                                  E_MeV = E_MeV,
                                  sensitivities = sensitivities,
                                  cc_icrp116 = cc_icrp116,
                                  save_result_callback = save_result_callback,
                                  readings = readings_vec), call_kwargs))
        if (is.list(res) && !is.null(res$spectrum)) {
            return(as.numeric(res$spectrum))
        }
        as.numeric(res)
    }
}

#' Load a bin lookup table
#'
#' Load a pre-computed bin lookup table, mirroring
#' \code{bssunfold.core.unfold_binned.load_bin_lookup}.
#'
#' @param path Path to a JSON lookup file. \code{NULL} (the default) returns
#'   the lookup shipped with the package (embedded in this file).
#' @return A list with \code{bin_to_methods} (a named list, one entry per bin,
#'   each holding \code{methods} and \code{scores}), \code{unique_methods}
#'   (character vector) and \code{n_bins} (integer).
#' @export
#' @examples
#' lk <- load_bin_lookup()
#' lk$n_bins
#' lk$bin_to_methods[["0"]]
load_bin_lookup <- function(path = NULL) {
    if (is.null(path)) {
        b2m <- lapply(seq_along(.BINNED_LOOKUP_RANKS), function(i) {
            list(methods = .BINNED_LOOKUP_METHODS[.BINNED_LOOKUP_RANKS[[i]]],
                 scores = .BINNED_LOOKUP_SCORES[[i]])
        })
        names(b2m) <- as.character(seq_len(length(b2m)) - 1L)
        return(list(bin_to_methods = b2m,
                    unique_methods = .BINNED_LOOKUP_METHODS,
                    n_bins = length(.BINNED_LOOKUP_RANKS)))
    }
    if (!file.exists(path)) {
        stop("Bin lookup not found at ", path,
             ".  Run ``python tools/build_bin_lookup.py`` to generate it.")
    }
    raw <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    keys <- as.integer(names(raw$bin_to_methods))
    o <- order(keys)
    keys <- keys[o]
    b2m <- lapply(raw$bin_to_methods[o], function(v) {
        if (is.data.frame(v)) {
            list(methods = as.character(v[[1]]), scores = as.numeric(v[[2]]))
        } else {
            list(methods = vapply(v, function(p) as.character(p[[1]]), ""),
                 scores = vapply(v, function(p) as.numeric(p[[2]]), 0))
        }
    })
    names(b2m) <- as.character(keys)
    list(bin_to_methods = b2m,
         unique_methods = as.character(raw$unique_methods),
         n_bins = if (is.null(raw$n_bins)) length(b2m) else as.integer(raw$n_bins))
}

#' Save a bin lookup table
#'
#' Persist a bin lookup table to JSON, mirroring \code{save_bin_lookup} in
#' \code{bssunfold.core.unfold_binned}.
#'
#' @param lookup A lookup list as returned by \code{\link{load_bin_lookup}}.
#' @param path File path; parent directories are created if needed.
#' @return Invisible \code{lookup}.
#' @export
#' @examples
#' \dontrun{
#' save_bin_lookup(load_bin_lookup(), tempfile(fileext = ".json"))
#' }
save_bin_lookup <- function(lookup, path) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    b2m <- lookup$bin_to_methods
    out <- lapply(seq_along(b2m), function(i) {
        lapply(seq_along(b2m[[i]]$methods), function(j) {
            c(b2m[[i]]$methods[j], b2m[[i]]$scores[j])
        })
    })
    names(out) <- names(b2m)
    jsonlite::write_json(list(bin_to_methods = out,
                              unique_methods = lookup$unique_methods,
                              n_bins = lookup$n_bins),
                         path, auto_unbox = TRUE, pretty = 2)
    invisible(lookup)
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

#' Bin-wise adaptive solver (low level)
#'
#' Run every candidate method and assemble a spectrum bin by bin, mirroring
#' \code{bssunfold.core.unfold_binned.solve_binned}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n) forwarded to every candidate.
#' @param E_MeV Numeric energy grid (length n). Kept for API compatibility;
#'   the bin-wise assembly does not use it.
#' @param n_super_bins,max_iterations,tolerance Legacy arguments kept for API
#'   compatibility; unused by the bin-wise assembly.
#' @param bin_lookup A lookup table from \code{\link{load_bin_lookup}}, or
#'   \code{NULL} for the shipped table.
#' @param methods Named list; \code{methods[[name]]} must be a function
#'   accepting \code{(A, b, x0 = NULL, ...)} and returning a length-n numeric
#'   spectrum. \code{NULL} means no candidate is available.
#' @param timeout_per_method Wall-clock timeout per candidate, in seconds.
#' @return A list \code{spectrum, method_map, candidate_methods,
#'   successful_methods, individual_spectra, errors, n_bins}.
#' @export
#' @examples
#' \dontrun{
#' r <- solve_binned(A, b, rep(1 / 60, 60), E)
#' }
solve_binned <- function(A, b, x0, E_MeV, n_super_bins = NULL,
                           max_iterations = 1000L, tolerance = 1e-6,
                           bin_lookup = NULL, methods = NULL,
                           timeout_per_method = 30) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n_bins <- ncol(A)
    if (is.null(bin_lookup)) bin_lookup <- load_bin_lookup()
    bin_to_methods <- bin_lookup$bin_to_methods

    candidate_names <- as.character(bin_lookup$unique_methods)
    if (length(candidate_names) == 0L) {
        seen <- unique(unlist(lapply(bin_to_methods, function(r) r$methods)))
        candidate_names <- sort(seen)
    }

    # ---- Run each candidate ----
    spectra <- list()
    successes <- character()
    errors <- list()

    for (name in candidate_names) {
        if (is.null(methods) || !name %in% names(methods)) {
            errors[[name]] <- "not provided"
            next
        }
        solver_fn <- methods[[name]]
        spec <- tryCatch({
            raw <- .binned_run_with_timeout(
                function() solver_fn(A, b, x0 = x0), timeout_per_method)
            raw <- as.numeric(raw)
            if (!is.null(dim(raw))) {
                NULL
            } else if (length(raw) != n_bins || !all(is.finite(raw)) ||
                       sum(raw) <= 0) {
                NULL
            } else {
                pmax(raw, 0)
            }
        }, error = function(e) {
            errors[[name]] <<- conditionMessage(e)
            NULL
        })
        if (!is.null(spec)) {
            spectra[[name]] <- spec
            successes <- c(successes, name)
        } else if (is.null(errors[[name]])) {
            errors[[name]] <- "invalid output"
        }
    }

    # ---- Assemble spectrum bin by bin ----
    assembled <- numeric(n_bins)
    method_map <- rep(-1L, n_bins)

    for (b_idx in seq_len(n_bins)) {
        key <- as.character(b_idx - 1L)
        ranking <- if (!is.null(bin_to_methods[[key]])) {
            bin_to_methods[[key]]$methods
        } else {
            character()
        }
        picked <- FALSE
        for (method_name in ranking) {
            if (!is.null(spectra[[method_name]])) {
                assembled[b_idx] <- spectra[[method_name]][b_idx]
                idx <- match(method_name, candidate_names)
                method_map[b_idx] <- if (is.na(idx)) -1L else idx - 1L
                picked <- TRUE
                break
            }
        }
        if (!picked) {
            vals <- vapply(spectra, function(s) s[b_idx], 0)
            assembled[b_idx] <- if (length(vals) > 0L) stats::median(vals) else 0
        }
    }

    list(spectrum = pmax(assembled, 0),
         method_map = method_map,
         candidate_methods = candidate_names,
         successful_methods = successes,
         individual_spectra = spectra,
         errors = errors,
         n_bins = n_bins)
}

#' Wrapper around \code{\link{solve_binned}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @param bin_lookup Pre-computed lookup table, or \code{NULL} to load
#'   \code{lookup_path} (default: the shipped table).
#' @param lookup_path Path to a JSON lookup file; ignored when \code{bin_lookup}
#'   is supplied.
#' @param timeout_per_method Wall-clock timeout per candidate method, seconds.
#' @param n_super_bins,max_iterations,tolerance Legacy arguments kept for API
#'   parity with the Python original; ignored by the binned solver.
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
                              max_neutron_energy = NULL,
                             bin_lookup = NULL, lookup_path = NULL,
                             timeout_per_method = 30) {
    if (is.null(bin_lookup)) bin_lookup <- load_bin_lookup(lookup_path)

    n <- as.integer(n_energy_bins)
    default_initial <- rep(1, n) / n
    x0 <- default_initial

    selected <- detector_names[detector_names %in% names(readings)]

    # Build the candidate solvers, mirroring METHOD_DISPATCH / _ALIASES.
    candidate_names <- as.character(bin_lookup$unique_methods)
    solver_dict <- list()
    for (name in candidate_names) {
        if (!name %in% names(.BINNED_METHOD_DISPATCH)) next
        dispatch_name <- .BINNED_METHOD_DISPATCH[[name]]
        fn <- NULL
        if (exists(dispatch_name, mode = "function")) {
            fn <- get(dispatch_name, mode = "function")
        }
        if (is.null(fn)) next
        alias <- .BINNED_ALIASES[[dispatch_name]]
        if (is.null(alias)) {
            alias_fixed <- list()
        } else {
            alias_fixed <- alias$fixed
        }
        kw <- c(alias_fixed,
                list(save_result = FALSE, calculate_errors = FALSE,
                     verbose = FALSE))
        solver_dict[[name]] <- list(
            solver = .make_binned_solver(fn, selected, detector_names, n,
                                         E_MeV, sensitivities, cc_icrp116,
                                         save_result_callback),
            kw = kw)
    }

    meta <- new.env(parent = emptyenv())
    solver <- function(A, b, x0 = NULL, ...) {
        methods <- lapply(solver_dict, function(e) {
            function(A_mat, b_vec, x0 = NULL) {
                do.call(e$solver, c(list(A_mat, b_vec, x0 = x0), e$kw))
            }
        })
        # Python always forwards its own ones/n default as x0, ignoring any
        # caller-supplied initial spectrum.
        out <- solve_binned(A, b, default_initial, E_MeV,
                            bin_lookup = bin_lookup,
                            methods = methods,
                            timeout_per_method = timeout_per_method)
        meta$result <- out
        out$spectrum
    }
    out <- run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = default_initial,
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
    if (!is.null(meta$result)) {
        out$method_map <- meta$result$method_map
        out$candidate_methods <- meta$result$candidate_methods
        out$successful_methods <- meta$result$successful_methods
        out$individual_spectra <- meta$result$individual_spectra
        out$bin_lookup <- bin_lookup
    }
    out
}
