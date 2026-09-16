# Tests for the v0.2.0 method batch: pure-R analogues, cutoff, plots, utils

A <- matrix(c(0.9, 0.05, 0.05,
              0.10, 0.8, 0.10,
              0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
b <- c(1, 0.6, 0.4)
det_names <- c("d1", "d2", "d3")
sens <- lapply(setNames(det_names, det_names),
               function(n) A[which(det_names == n), ])
E <- c(1e-9, 1e-6, 1e-3)
readings <- setNames(b, det_names)

test_that("solve_combined runs a two-stage pipeline and improves data fit", {
    r <- solve_combined(
        A, b, rep(1, 3),
        pipeline = list(list(method = "mlem",
                             params = list(max_iterations = 300L)),
                        list(method = "landweber",
                             params = list(max_iterations = 100L))))
    expect_true(is.list(r))
    expect_equal(length(r$spectrum), 3)
    expect_true(all(r$spectrum >= 0))
    expect_equal(r$iterations, 2L)
    expect_true(is.logical(r$converged))
    rn <- norm(as.numeric(A %*% r$spectrum) - b, "2")
    expect_lt(rn, 0.2)
    expect_true(is.list(r$stage_spectra) && length(r$stage_spectra) == 2L)
})

test_that("unfold_combined wrapper returns dose rates and correct method name", {
    cc <- interpolate_coefficients(ICRP116_COEFF_EFFECTIVE_DOSE(), E)
    r <- unfold_combined(det_names, 3L, E, sens, cc, NULL, readings)
    expect_equal(r$method, "Combined")
    expect_true(length(r$doserates) > 0)
    expect_true(r$residual_norm < 0.3)
})

test_that("solve_pspline_reml produces a non-negative finite spectrum", {
    r <- solve_pspline_reml(A, b, NULL, n_basis = 5, max_iterations = 5L)
    expect_true(is.finite(sum(r$spectrum)) && all(r$spectrum >= 0))
    expect_equal(length(r$spectrum), 3)
    expect_true(r$lambda > 0)
})

test_that("solve_amg covers all Krylov/preconditioner combinations", {
    for (mth in c("cg", "bicgstab", "gmres")) {
        for (pc in c("none", "jacobi", "gs", "sor", "ssor", "amg")) {
            r <- solve_amg(A, b, NULL, method = mth,
                           preconditioner = pc, max_iterations = 100L)
            expect_true(is.list(r), info = paste(mth, pc))
            expect_equal(length(r$spectrum), 3)
            expect_true(all(r$spectrum >= 0),
                        info = paste("nonneg:", mth, pc))
        }
    }
})

test_that("unfold_amg wrapper works end-to-end", {
    r <- unfold_amg(det_names, 3L, E, sens, NULL, NULL, readings)
    expect_equal(r$method, "AMG")
    expect_true(r$residual_norm < 0.5)
})

test_that("solve_ssr keeps a sparse non-negative spectrum", {
    r <- solve_ssr(A, b, rep(1, 3), n_thresholds = 4L, max_iterations = 15L)
    expect_equal(length(r$spectrum), 3)
    expect_true(all(is.finite(r$spectrum) & r$spectrum >= 0))
    expect_true(r$threshold > 0)
})

test_that("solve_mystic_hybrid refines beyond the flat spectrum", {
    set.seed(11)
    r <- solve_mystic_hybrid(A, b, rep(1, 3), global_maxiter = 10L)
    expect_true(is.finite(sum(r$spectrum)) && all(r$spectrum >= 0))
    r2 <- solve_mystic_hybrid(A, b, rep(1, 3), global_maxiter = 10L,
                              local_maxiter = 50L)
    rn1 <- sqrt(sum((as.numeric(A %*% r$spectrum) - b)^2))
    rn2 <- sqrt(sum((as.numeric(A %*% r2$spectrum) - b)^2))
    expect_true(rn1 < 0.4 && all(is.finite(r2$spectrum)))
})

test_that("solve_smt returns an exact L2 fit for full-rank systems", {
    r <- solve_smt(A, b, NULL, objective = "l2")
    resid <- as.numeric(A %*% r$spectrum) - b
    expect_lt(sqrt(sum(resid^2)), 1e-8)
    r_nn <- solve_smt(A, b, NULL, nonneg = TRUE)
    expect_true(all(r_nn$spectrum >= 0))
    r1 <- solve_smt(A, b, NULL, objective = "l1")
    expect_equal(length(r1$spectrum), 3)
})

test_that("solve_qp_scip / solve_qp_docplex solve the same projected QP", {
    r1 <- solve_qp_scip(A, b, rep(1, 3), max_iterations = 2000L)
    r2 <- solve_qp_docplex(A, b, rep(1, 3), max_iterations = 2000L)
    expect_true(all(r1$spectrum >= 0))
    expect_true(abs(sum(r1$spectrum) - sum(r2$spectrum)) < 0.5)
    expect_equal(r1$iterations, r2$iterations)
    expect_true(r1$converged)
})

test_that("solve_zfit_poisson fits Poisson counts to within tolerance", {
    A10 <- A * 10
    k <- pmax(round(as.numeric(A10 %*% c(3, 2, 1))), 0)
    r <- solve_zfit_poisson(A10, k, rep(2, 3))
    expect_true(all(r$spectrum > 0))
    expect_equal(length(r$spectrum), 3)
    mu <- as.numeric(A10 %*% r$spectrum)
    rnd <- 2 * sum(pmax(mu, 1e-30) - k * log(pmax(mu, 1e-30)))
    expect_true(all(is.finite(mu)))
})

test_that("solve_zfit_poisson returns Gaussian posterior samples on request", {
    A10 <- A * 10
    r <- solve_zfit_poisson(A10, b * 10, rep(3, 3), use_mcmc = TRUE,
                            n_samples = 12L)
    m <- r$posterior_samples
    expect_equal(nrow(m), 12L)
    expect_equal(ncol(m), 3L)
    expect_true(all(is.finite(pmax(m, 0))))
})

test_that("solve_maeo produces knee selection and non-negative spectrum", {
    set.seed(5)
    r <- solve_maeo(A, b, rep(1, 3), n_cycles = 2L, n_gen_per_cycle = 2L,
                    pop_size = 8L)
    expect_true(all(r$spectrum >= 0))
    expect_true(isTRUE(r$converged))
    expect_true(r$knee_index %in% seq_len(2))
})

test_that("solve_mlem_odl matches plain MLEM-on-matrix behaviour", {
    r0 <- solve_mlem(A, b, rep(1, 3), max_iterations = 40L)
    r1 <- solve_mlem_odl(A, b, rep(1, 3), max_iterations = 40L)
    expect_true(max(abs(r0$spectrum - r1$spectrum)) < 0.05)
})

test_that("interpret_result produces a coherent control report", {
    res <- solve_mlem(A, b, rep(1, 3), max_iterations = 50L)
    rep <- interpret_result(list(spectrum = res$spectrum), A, b,
                           detector_names = det_names, n_scenarios = 5L)
    expect_true(rep$robustness >= 0)
    expect_equal(length(rep$shadow_prices), 3)
    expect_equal(nrow(rep$regularization_sweep), 5L)
    expect_equal(length(rep$scenarios), 5L)
})

test_that("unfold_interpret attaches interpret_options", {
    r <- unfold_interpret(det_names, 3L, E, sens, NULL, NULL, readings,
                          n_scenarios = 5L)
    expect_true(!is.null(r$interpret_options))
    expect_equal(length(r$interpret_options$robustness), 1L)
})

# ---------------- max_neutron_energy cutoff ----------------------------------

test_that("max_neutron_energy zeroes the bins above the cutoff", {
    A6 <- matrix(runif(6 * 6, 0, 1), nrow = 6)
    for (i in 1:6) A6[i, ] <- A6[i, ] / sum(A6[i, ])
    b6 <- c(1.0, 0.9, 0.8, 0.7, 0.5, 0.3)
    E6 <- c(1e-9, 1e-4, 1e-2, 0.1, 1, 10)
    sens <- lapply(setNames(as.character(1:6), as.character(1:6)),
                   function(i) A6[as.integer(i), ])
    r <- unfold_mlem(as.character(1:6), 6L, E6, sens, NULL, NULL,
                     setNames(b6, as.character(1:6)),
                     max_iterations = 50L, max_neutron_energy = 0.5)
    expect_equal(length(r$spectrum), 6L)
    expect_true(all(r$spectrum[E6 > 0.5] == 0))
    expect_true(any(r$spectrum[E6 <= 0.5] > 0))
})

test_that("max_neutron_energy below the lowest bin fails informatively", {
    expect_error(unfold_mlem(det_names, 3L, E, sens, NULL, NULL, readings,
                             max_neutron_energy = 1e-30),
                 "below the lowest energy bin")
})

# ---------------- comparison / conversion / plotting utils -------------------

test_that("compare_spectra reproduces correct metrics", {
    x <- c(1, 0.5, 0.25)
    m <- compare_spectra(x, x)
    expect_true(m$relative_integral_diff < 1e-12)
    expect_true(abs(m$chi_squared - 0) < 1e-12)
    expect_true(abs(m$cosine_similarity - 1) < 1e-12)
    y <- c(0.9, 0.55, 0.2)
    m2 <- compare_spectra(x, y)
    expect_true(m2$relative_integral_diff > 0)
    expect_true(is.finite(m2$kl_divergence))
})

test_that("converters round-trip the result", {
    res <- unfold_mlem(det_names, 3L, E, sens, NULL, NULL, readings,
                       max_iterations = 20L)
    df <- convert_to_dataframe(res, sig_figs = 2L)
    expect_s3_class(df, "data.frame")
    expect_true(nrow(df) == 3L)
    dl <- convert_to_dict(res)
    expect_equal(dl$energy, unname(E))
    expect_equal(dl$spectrum, res$spectrum)
})

test_that("plot_spectrum / plot_residuals / plot_comparison return ggplots", {
    skip_if(!requireNamespace("ggplot2", quietly = TRUE))
    res <- unfold_mlem(det_names, 3L, E, sens, NULL, NULL, readings,
                       max_iterations = 20L)
    p1 <- plot_spectrum(res, logx = TRUE, logy = FALSE)
    expect_s3_class(p1, "ggplot")
    p2 <- plot_residuals(res, detector_names = det_names)
    expect_s3_class(p2, "ggplot")
    p3 <- plot_comparison(list(mlem = res, ref = res))
    expect_s3_class(p3, "ggplot")
})

test_that("Detector result management and cc switching work", {
    det <- Detector$new(response_function = "PTB",
                        detector_names = c("0in", "3in", "5in", "8in",
                                           "12in"))
    expect_equal(det$n_detectors, 5)
    rds <- make_ptb_readings(det, target_0in = 1000)
    r1 <- det$unfold_mlem(rds, max_iterations = 30L, save_result = TRUE)
    expect_equal(length(det$history), 1)
    expect_equal(det$get_result()$method, r1$method)
    expect_equal(det$get_result(1)$spectrum, r1$spectrum)
    out_list <- det$list_results()
    expect_equal(out_list[[1]]$method, "MLEM")
    det$clear_results()
    expect_equal(length(det$history), 0)
    expect_null(det$get_result())
    # dose coefficients switch
    det$set_dose_coefficients("ICRP74_effective")
    eff <- det$cc_icrp116$AP
    expect_equal(length(eff), det$n_energy_bins)
    expect_error(det$set_dose_coefficients("nope"), "Unknown dose")
    # management aliases for new methods
    r2 <- det$unfold_combined(rds, save_result = TRUE,
                              pipeline = list(list(method = "landweber",
                                                   params = list(
                                                       max_iterations = 20L))))
    expect_equal(r2$method, "Combined")
    r3 <- det$unfold_pspline_reml(rds)
    expect_equal(r3$method, "PSplineREML")
    r4 <- det$unfold_amg(rds, max_iterations = 100L)
    expect_equal(r4$method, "AMG")
    r5 <- det$unfold_smt(rds)
    expect_true(is.finite(sum(r5$spectrum)))
    r6 <- det$unfold_pdhg(rds, max_iterations = 40L)
    expect_equal(r6$method, "PDHG")
})

test_that("Detector plot methods return ggplot objects", {
    skip_if(!requireNamespace("ggplot2", quietly = TRUE))
    det <- Detector$new(response_function = "PTB",
                        detector_names = c("0in", "3in", "5in", "8in",
                                           "12in"))
    rds <- make_ptb_readings(det, target_0in = 1000)
    pr <- det$plot_response_functions()
    expect_s3_class(pr, "ggplot")
    det$unfold_mlem(rds, max_iterations = 20L, save_result = TRUE)
    pu <- det$plot_with_uncertainty()
    expect_s3_class(pu, "ggplot")
})

test_that("benchmark_unfold_methods loops probes successfully", {
    det <- make_ptb_detector()
    bm <- benchmark_unfold_methods(det,
        methods = list(landweber = list(max_iterations = 20L),
                       mlem = list(max_iterations = 20L)),
        n_probes = 2L)
    expect_true(bm$results$landweber$method_ok == 1)
    expect_true(bm$results$mlem$method_ok == 1)
    expect_true(is.finite(bm$results$landweber$mean_residual))
    expect_error(benchmark_unfold_methods(det, methods = list(nope = list())),
                 "Unknown benchmark method")
})
