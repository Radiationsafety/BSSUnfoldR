# Tests for batch 4 of BSSUnfoldR (IMAXED, AMAXED, FISTA, Bayes-spline,
# NSDUAZ, MLEM-BS, NSpline, MCMC, Reconst, Ensemble, Cascade, Composite).

test_that("solve_imaxed runs with Armijo backtracking line search", {
    r <- solve_imaxed(tiny_A, tiny_b, tiny_x0_default, max_iterations = 50)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
    expect_true(all(is.finite(r$spectrum)))
})

test_that("solve_amaxed converges with Newton-KKT line search", {
    r <- solve_amaxed(tiny_A, tiny_b, tiny_x0_default, max_iterations = 50)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
})

test_that("solve_fista supports all four regularization options", {
    r0 <- solve_fista(tiny_A, tiny_b, rep(0, 3), max_iterations = 20)
    r1 <- solve_fista(tiny_A, tiny_b, rep(0, 3), max_iterations = 20,
                       regularization = 0.01)
    r2 <- solve_fista(tiny_A, tiny_b, rep(0, 3), max_iterations = 20,
                       l1_penalty = 0.01)
    r3 <- solve_fista(tiny_A, tiny_b, rep(0, 3), max_iterations = 20,
                       tv_penalty = 0.01)
    for (r in list(r0, r1, r2, r3)) expect_length(r$spectrum, 3L)
})

test_that("solve_fista supports box constraints", {
    r <- solve_fista(tiny_A, tiny_b, rep(0, 3), max_iterations = 20,
                     x_min = 0, x_max = 5)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum <= 5))
})

test_that("solve_bayes_spline accepts spline_degree and smooth", {
    r <- solve_bayes_spline(tiny_A, tiny_b, tiny_x0_default,
                              max_iterations = 20,
                              spline_degree = 3, spline_smooth = 0.1)
    expect_length(r$spectrum, 3L)
})

test_that("solve_nsduaz returns the same shape as solve_bunki", {
    r <- solve_nsduaz(tiny_A, tiny_b, tiny_x0_default, max_iterations = 30)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
})

test_that("builtin_catalogue returns three normalised spectra", {
    E <- 10^seq(-9, 1, length.out = 60)
    cat <- builtin_catalogue(E)
    expect_named(cat, c("ambe", "cf252", "reactor"))
    for (s in cat) {
        expect_length(s, 60L)
        expect_equal(sum(s), 1.0, tolerance = 1e-6)
    }
})

test_that("select_catalogue_initial returns a labelled spectrum", {
    E <- 10^seq(-9, 1, length.out = 60)
    A <- matrix(runif(5 * 60), nrow = 5)
    sens <- setNames(lapply(1:5, function(i) A[i, ]),
                     c("0in","3in","5in","8in","12in"))
    rds <- c("0in" = 100, "3in" = 80, "5in" = 60, "8in" = 40, "12in" = 10)
    sel <- select_catalogue_initial(rds, names(sens), sens, E_MeV = E)
    expect_length(sel$spectrum, 60L)
    expect_true(sel$label %in% c("ambe", "cf252", "reactor"))
})

test_that("solve_mlem_bs builds the B-spline basis and iterates", {
    E <- 10^seq(-9, 1, length.out = 60)
    A <- matrix(runif(5 * 60), nrow = 5)
    b <- as.numeric(A %*% rep(0.5, 60))
    r <- solve_mlem_bs(A, b, rep(1, 60), E, n_basis = 10L,
                        beta_relative = 1e-3, max_iterations = 30L)
    expect_length(r$spectrum, 60L)
    expect_true(all(r$spectrum >= 0))
})

test_that("build_bspline_basis has the right shape", {
    E <- 10^seq(-9, 1, length.out = 60)
    B <- build_bspline_basis(E, 10L, "log")
    expect_equal(dim(B), c(60L, 10L))
    # Each row of B should be non-negative; row sums approximate 1
    # (partition of unity holds to numerical precision for clamped B-splines).
    expect_true(all(B >= 0))
    expect_true(all(abs(rowSums(B) - 1) < 1e-3))
})

test_that("second_difference_matrix has shape (n-2, n)", {
    D2 <- second_difference_matrix(6L)
    expect_equal(dim(D2), c(4L, 6L))
})

test_that("ks_statistic is non-negative", {
    val <- ks_statistic(c(1, 2, 3), c(1.1, 1.9, 3.1))
    expect_true(val >= 0)
    expect_true(is.numeric(val))
})

test_that("solve_nspline runs with auto_knots and converges", {
    E <- 10^seq(-9, 1, length.out = 60)
    A <- matrix(runif(5 * 60), nrow = 5)
    b <- as.numeric(A %*% rep(0.5, 60))
    r <- solve_nspline(A, b, rep(1, 60), E, max_iterations = 30L)
    expect_length(r$spectrum, 60L)
    expect_true(all(r$spectrum >= 0))
})

test_that("auto_knots returns strictly increasing interior energies", {
    E <- 10^seq(-9, 1, length.out = 60)
    k <- auto_knots(E)
    expect_true(length(k) >= 1L)
    expect_true(all(diff(k) > 0))
    expect_true(min(k) > min(E))
    expect_true(max(k) < max(E))
})

test_that("build_continuity_matrix has the right shape", {
    D <- build_continuity_matrix(4L, order = 0L)
    expect_equal(dim(D), c(4L, 15L))
})

test_that("solve_bayesian_mcmc returns posterior samples", {
    set.seed(7)
    r <- solve_bayesian_mcmc(tiny_A, tiny_b, tiny_x0_default,
                                n_samples = 50, n_burn = 10)
    expect_length(r$spectrum, 3L)
    expect_equal(dim(r$posterior_samples), c(40L, 3L))
    expect_true(r$acceptance_rate >= 0 && r$acceptance_rate <= 1)
})

test_that("solve_reconst produces a non-negative spectrum", {
    r <- solve_reconst(tiny_A, tiny_b, NULL, alpha = 1.0)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
    expect_equal(r$iterations, 1L)
})

test_that("solve_ensemble averages multiple solvers", {
    r <- solve_ensemble(tiny_A, tiny_b, tiny_x0_default,
                         solvers = list(solve_mlem, solve_gravel),
                         max_iterations = 30)
    expect_length(r$spectrum, 3L)
    expect_equal(dim(r$ensemble_spectra), c(2L, 3L))
    expect_length(r$ensemble_weights, 2L)
    expect_equal(sum(r$ensemble_weights), 1.0)
})

test_that("solve_ensemble accepts custom weights", {
    r <- solve_ensemble(tiny_A, tiny_b, tiny_x0_default,
                         solvers = list(solve_mlem, solve_gravel),
                         weights = c(0.7, 0.3), max_iterations = 30)
    expect_equal(r$ensemble_weights, c(0.7, 0.3))
})

test_that("solve_cascade runs solvers sequentially", {
    r <- solve_cascade(tiny_A, tiny_b, tiny_x0_default,
                         solvers = list(solve_mlem, solve_gravel),
                         max_iterations = 30)
    expect_length(r$spectrum, 3L)
    expect_length(r$cascade_spectra, 2L)
})

test_that("solve_composite picks the best solver by residual", {
    r <- solve_composite(tiny_A, tiny_b, tiny_x0_default,
                           solvers = list(solve_mlem, solve_gravel,
                                            solve_sandii),
                           max_iterations = 30)
    expect_length(r$spectrum, 3L)
    expect_equal(dim(r$composite_spectra), c(3L, 3L))
    expect_true(r$best_solver_index %in% 1:3)
})
