test_that("PTB unfolding runs end-to-end for every batch-1 method", {
    det <- make_ptb_detector()
    rds <- make_ptb_readings(det, target_0in = 1000)
    common <- list(max_iterations = 50L)
    method_calls <- list(
        unfold_mlem     = common,
        unfold_gravel   = common,
        unfold_maxed    = common,
        unfold_sandii   = c(common, list(chi_fac = 1L)),
        unfold_bunki    = common,
        unfold_ferdor   = common,
        unfold_staysl   = list(),  # one-step, no max_iterations
        unfold_landweber = common,
        unfold_cgls     = list(max_iterations = 30L),
        unfold_tsvd     = list(method = "gcv")
    )
    for (fn_name in names(method_calls)) {
        fn <- get(fn_name)
        r <- do.call(fn, c(list(
            detector_names = det$detector_names,
            n_energy_bins = det$n_energy_bins, E_MeV = det$E_MeV,
            sensitivities = det$sensitivities, cc_icrp116 = det$cc_icrp116,
            save_result_callback = function(o) det$save_result(o),
            readings = rds
        ), method_calls[[fn_name]]))
        expect_length(r$spectrum, 60L)
        expect_true(all(r$spectrum >= 0))
        expect_named(r$doserates,
                     c("AP", "PA", "LLAT", "RLAT", "ROT", "ISO"))
    }
})

test_that("PTB unfolding runs end-to-end for every batch-2 method", {
    det <- make_ptb_detector()
    rds <- make_ptb_readings(det, target_0in = 1000)
    common_args <- list(max_iterations = 50L)
    method_calls <- list(
        unfold_osem  = c(common_args, list(n_subsets = 5L)),
        unfold_mapem = c(common_args, list(prior = "quadratic", beta = 1e-2)),
        unfold_bsrem = c(common_args, list(n_subsets = 5L)),
        unfold_sart  = common_args,
        unfold_kaczmarz = common_args,
        unfold_randomized_kaczmarz = c(common_args, list(random_state = 42L)),
        unfold_lanczos = list(),
        unfold_tikhonov_legendre = list(delta = 0.05, n_polynomials = 15L),
        unfold_rebunki  = c(common_args, list(tolerance = 0.05)),
        unfold_doroshenko = common_args
    )
    for (fn_name in names(method_calls)) {
        fn <- get(fn_name)
        r <- do.call(fn, c(list(
            detector_names = det$detector_names,
            n_energy_bins = det$n_energy_bins, E_MeV = det$E_MeV,
            sensitivities = det$sensitivities, cc_icrp116 = det$cc_icrp116,
            save_result_callback = function(o) det$save_result(o),
            readings = rds
        ), method_calls[[fn_name]]))
        expect_length(r$spectrum, 60L)
        expect_true(all(r$spectrum >= 0))
    }
})

test_that("Detector R6 convenience wrappers save results to history", {
    det <- make_ptb_detector()
    rds <- make_ptb_readings(det, target_0in = 1000)
    det$history <- list()
    invisible(det$unfold_mlem(rds, max_iterations = 20L, save_result = TRUE))
    invisible(det$unfold_gravel(rds, max_iterations = 20L, save_result = TRUE))
    invisible(det$unfold_rebunki(rds, max_iterations = 20L, save_result = TRUE))
    expect_length(det$history, 3L)
    methods_seen <- sapply(det$history, function(h) h$method)
    expect_equal(methods_seen, c("MLEM", "GRAVEL", "ReBUNKI"))
})

test_that("Monte-Carlo uncertainty estimation works for a batch-1 method", {
    det <- make_ptb_detector()
    rds <- make_ptb_readings(det, target_0in = 1000)
    r <- unfold_mlem(det$detector_names, det$n_energy_bins, det$E_MeV,
                      det$sensitivities, det$cc_icrp116,
                      function(o) det$save_result(o), rds,
                      max_iterations = 30L,
                      calculate_errors = TRUE,
                      n_montecarlo = 5L, noise_level = 0.05,
                      random_state = 42L)
    expect_equal(r$montecarlo_samples, 5L)
    expect_equal(r$noise_level, 0.05)
    expect_equal(dim(r$spectrum_uncert_all), c(5L, 60L))
})
