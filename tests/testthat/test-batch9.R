# Unit tests for the batch-9 ports: adaptive cascade, MAEO-ensemble,
# fission GA (BonnerFinder) and the pure-R CUQI Bayesian samplers.

test_that("unfold_adaptive_cascade selects methods greedily and returns the cascade contract", {
    r <- unfold_adaptive_cascade(
        detector_names = c("D1", "D2", "D3"), n_energy_bins = 3L,
        E_MeV = tiny_E,
        sensitivities = list(D1 = tiny_A[1, ], D2 = tiny_A[2, ],
                             D3 = tiny_A[3, ]),
        cc_icrp116 = NULL,
        save_result_callback = function(o) invisible(o),
        readings = setNames(as.numeric(tiny_A %*% c(2, 1, 0.5)),
                            c("D1", "D2", "D3")),
        max_stages = 2L, initial_method = "tsvd")
    expect_type(r$spectrum, "double")
    expect_length(r$spectrum, 3L)
    expect_equal(r$method, "Adaptive Cascade")
    expect_true(length(r$method_sequence) >= 1L)
    expect_true(r$status %in% c("OK", "PARTIAL_SUCCESS", "ERROR"))
    expect_false(r$status == "ERROR")
    expect_type(r$quality_metrics, "list")
})

test_that("solve_maeo_ensemble passes through labels and is reproducible with a seed", {
    r1 <- solve_maeo_ensemble(tiny_A, tiny_b, rep(1, 3),
                              n_cycles = 2L, n_gen_per_cycle = 3L,
                              pop_size = 8L, seed = 42L)
    r2 <- solve_maeo_ensemble(tiny_A, tiny_b, rep(1, 3),
                              n_cycles = 2L, n_gen_per_cycle = 3L,
                              pop_size = 8L, seed = 42L)
    expect_named(r1, c("spectrum", "iterations", "converged", "knee_index",
                       "maeo_spectra", "maeo_objs", "migration_method",
                       "parallel_enabled", "algorithms"),
                 ignore.order = TRUE)
    expect_equal(r1$spectrum, r2$spectrum)
    expect_true(all(r1$spectrum >= 0))
    expect_equal(r1$algorithms, c("nsga3", "ctea", "agemoea2", "spea2"))
    expect_false(r1$parallel_enabled)
})

test_that("solve_fission_ga returns params, validation and a reproducible spectrum", {
    ln_steps <- compute_log_steps(tiny_E) * log(10)
    r1 <- solve_fission_ga(tiny_A, tiny_b, tiny_E, ln_steps,
                           ga_popsize = 4L, ga_maxiter = 8L,
                           random_state = 5L)
    r2 <- solve_fission_ga(tiny_A, tiny_b, tiny_E, ln_steps,
                           ga_popsize = 4L, ga_maxiter = 8L,
                           random_state = 5L)
    expect_length(r1$spectrum, 3L)
    expect_true(all(r1$spectrum >= 0))
    expect_true(is.finite(r1$params$cost))
    expect_equal(r1$spectrum, r2$spectrum)
    expect_equal(names(r1$params),
                 c("a1", "a2", "a3", "b", "beta", "alpha", "TF",
                   "phi_scale", "weight_fractions", "cost"))
})

test_that("unfold_fission_ga attaches model params and validation to the result", {
    dn <- c("D1", "D2", "D3")
    sens <- list(D1 = tiny_A[1, ], D2 = tiny_A[2, ], D3 = tiny_A[3, ])
    readings <- setNames(as.numeric(tiny_A %*% c(2, 1, 0.5)), dn)
    r <- suppressWarnings(unfold_fission_ga(
        dn, 3L, tiny_E, sens, NULL, function(o) invisible(o), readings,
        ga_popsize = 4L, ga_maxiter = 8L, random_state = 9L))
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
    expect_true("model_params" %in% names(r))
    expect_true("validation" %in% names(r))
    expect_true(all(c("fom_percent", "relative_uncertainties",
                      "residual_sign_changes", "norm_ok", "passed") %in%
                        names(r$validation)))
})

test_that("solve_cuqi_bayesian runs every documented sampler and validates inputs", {
    for (s in c("pcn", "cwmh", "mala", "ula", "gibbs", "gibbs_nuts")) {
        r <- suppressWarnings(solve_cuqi_bayesian(
            tiny_A, tiny_b, sampler = s, n_samples = 60L,
            n_burnin = 30L, chains = 1L, random_state = 1L))
        expect_equal(length(r$spectrum), 3L, info = s)
        expect_true(all(r$spectrum > 0), info = s)
        expect_true(is.finite(r$stats$acc_rate) || s == "ula", info = s)
    }
    # NUTS falls back to CWMH with a warning
    expect_warning(solve_cuqi_bayesian(tiny_A, tiny_b, sampler = "nuts",
                                       n_samples = 20L, n_burnin = 10L,
                                       chains = 1L),
                   "NUTS")
    expect_error(solve_cuqi_bayesian(tiny_A, tiny_b, sampler = "zzz"))
    expect_error(solve_cuqi_bayesian(tiny_A, tiny_b, prior = "zzz"))
    expect_error(solve_cuqi_bayesian(tiny_A, tiny_b, hierarchical = TRUE,
                                     prior = "ou"))
    expect_error(solve_cuqi_bayesian(tiny_A, tiny_b, credible_level = 0))
})

test_that("solve_cuqi_bayesian is seed-reproducible and reports diagnostics", {
    r1 <- solve_cuqi_bayesian(tiny_A, tiny_b, sampler = "pcn",
                              n_samples = 80L, n_burnin = 40L,
                              chains = 2L, random_state = 17L)
    r2 <- solve_cuqi_bayesian(tiny_A, tiny_b, sampler = "pcn",
                              n_samples = 80L, n_burnin = 40L,
                              chains = 2L, random_state = 17L)
    expect_equal(r1$spectrum, r2$spectrum)
    expect_equal(dim(r1$stats$samples), c(160L, 3L))
    expect_equal(dim(r1$stats$theta_samples), c(160L, 3L))
    expect_length(r1$stats$rhat, 3L)
    expect_true(all(r1$stats$rhat > 0.9))
    expect_length(r1$stats$ess, 3L)
    expect_true(all(r1$stats$hpd_upper >= r1$stats$hpd_lower))
    expect_null(r1$stats$delta_samples)
    # hierarchical Gibbs carries hyperparameter draws
    rg <- solve_cuqi_bayesian(tiny_A, tiny_b, sampler = "gibbs",
                              n_samples = 60L, n_burnin = 30L,
                              chains = 2L, random_state = 4L)
    expect_length(rg$stats$delta_samples, 120L)
    expect_true(all(rg$stats$delta_samples > 0))
})

test_that("unfold_cuqi wraps run_unfolding and merges the credible band", {
    dn <- c("D1", "D2", "D3")
    sens <- list(D1 = tiny_A[1, ], D2 = tiny_A[2, ], D3 = tiny_A[3, ])
    readings <- setNames(as.numeric(tiny_A %*% c(2, 1, 0.5)), dn)
    r <- unfold_cuqi(dn, 3L, tiny_E, sens, NULL,
                     function(o) invisible(o), readings,
                     sampler = "pcn", n_samples = 60L, n_burnin = 30L,
                     chains = 1L, random_state = 2L)
    expect_equal(r$method, "CUQI-Bayesian")
    expect_length(r$spectrum, 3L)
    expect_length(r$spectrum_uncertainty, 3L)
    expect_true(all(r$spectrum_lower >= 0))
    expect_true(all(r$spectrum_upper >= r$spectrum_lower))
    expect_true("cuqi_stats" %in% names(r))
    expect_equal(r$cuqi_stats$backend, "R-pure (CUQI-style)")
})
