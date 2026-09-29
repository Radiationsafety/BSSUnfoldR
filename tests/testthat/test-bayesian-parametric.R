test_that("solve_bayesian_parametric follows the Python chain conventions", {
    E <- 10^seq(-9, 2, length.out = 40)
    ls <- compute_log_steps(E)
    set.seed(1)
    A <- matrix(runif(4 * 40, 0, 1), nrow = 4) * 0.1
    truth <- c(A_th = 4e-5, T_th = 0.025e-6, A_epi = 3e-5, A_f = 2e-5,
               T_ev = 2.5)
    spec <- .parametric_model(E, truth[[1]], truth[[2]], truth[[3]],
                              truth[[4]], truth[[5]]) * ls
    b <- as.numeric(A %*% spec)

    expect_named(formals(solve_bayesian_parametric),
                 c("A", "b", "E", "log_steps", "sigma", "initial_params",
                   "n_samples", "burn_in", "proposal_scale", "random_state"))
    expect_equal(formals(solve_bayesian_parametric)$sigma, 0.02)

    r <- solve_bayesian_parametric(A, b, E, ls, sigma = 1e-3,
                                   n_samples = 400L, burn_in = 100L,
                                   random_state = 3L)
    expect_true(all(is.finite(r$spectrum)))
    expect_gte(min(r$spectrum), 0)
    expect_equal(length(r$spectrum), length(E))
    # The chain is seeded from PCG64 in Python, so only the posterior mean has
    # to come back to the generating parameters, not the individual draws.
    expect_lt(norm(r$spectrum - spec, "2") / norm(spec, "2"), 0.2)
    expect_gt(r$acceptance_rate, 0)
})

test_that("unfold_bayesian_parametric is reachable from the Detector", {
    det <- make_ptb_detector()
    readings <- make_ptb_readings(det)
    r <- det$unfold_bayesian_parametric(readings, n_samples = 60L,
                                        burn_in = 20L, random_state = 5L)
    expect_equal(r$method, "bayesian_parametric")
    expect_length(r$spectrum, det$n_energy_bins)
    expect_true(all(is.finite(r$spectrum)))
    expect_equal(r$sigma, 0.02)
})

test_that("compute_log_steps returns base-10 widths, not natural-log ones", {
    # unfold_fruit_like converts to natural log; unfold_bayesian_parametric
    # deliberately does not, so the two conventions must stay distinguishable.
    E <- 10^seq(-9, 2, length.out = 40)
    ls <- compute_log_steps(E)
    step <- diff(log10(E))[1]
    expect_equal(ls, rep(step, 40L), tolerance = 1e-6)
    expect_equal(max(ls) * log(10), step * log(10), tolerance = 1e-6)
})
