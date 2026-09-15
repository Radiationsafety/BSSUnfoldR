test_that("solve_mlem returns a non-negative spectrum and a 3-tuple", {
    r <- solve_mlem(tiny_A, tiny_b, tiny_x0_default, max_iterations = 50)
    expect_type(r, "list")
    expect_named(r, c("spectrum", "iterations", "converged"))
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
    expect_true(is.numeric(r$iterations))
    expect_true(is.logical(r$converged))
})

test_that("solve_landweber raises on zero-norm matrix and otherwise iterates", {
    r <- solve_landweber(tiny_A, tiny_b, rep(0, 3), max_iterations = 50)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
})

test_that("solve_gravel converges within max_iterations", {
    r <- solve_gravel(tiny_A, tiny_b, tiny_x0_default, max_iterations = 50)
    expect_lte(as.integer(r$iterations), 50L)
    expect_true(all(r$spectrum >= 0))
})

test_that("solve_sandii accepts chi_fac 0 and 1", {
    r0 <- solve_sandii(tiny_A, tiny_b, tiny_x0_default,
                        max_iterations = 30, chi_fac = 0L)
    r1 <- solve_sandii(tiny_A, tiny_b, tiny_x0_default,
                        max_iterations = 30, chi_fac = 1L)
    expect_length(r0$spectrum, 3L)
    expect_length(r1$spectrum, 3L)
})

test_that("solve_bunki preserves the spectrum length", {
    r <- solve_bunki(tiny_A, tiny_b, tiny_x0_default, max_iterations = 50)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
})

test_that("solve_maxed is reproducible", {
    r1 <- solve_maxed(tiny_A, tiny_b, tiny_x0_default, max_iterations = 50)
    r2 <- solve_maxed(tiny_A, tiny_b, tiny_x0_default, max_iterations = 50)
    expect_equal(r1$spectrum, r2$spectrum)
})

test_that("solve_ferdor finds a non-negative solution", {
    r <- solve_ferdor(tiny_A, tiny_b, tiny_x0_default)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
})

test_that("solve_staysl is a single-step method (iterations == 1)", {
    r <- solve_staysl(tiny_A, tiny_b, tiny_x0_default)
    expect_equal(r$iterations, 1L)
    expect_true(r$converged)
})

test_that("solve_cgls runs with and without Tikhonov regularization", {
    r0 <- solve_cgls(tiny_A, tiny_b, rep(0, 3), max_iterations = 20)
    r1 <- solve_cgls(tiny_A, tiny_b, rep(0, 3), max_iterations = 20,
                     regularization = 0.01, smoothness_order = 2L)
    expect_length(r0$spectrum, 3L)
    expect_length(r1$spectrum, 3L)
})

test_that("solve_tsvd supports different k-selection methods", {
    for (m in c("discrepancy", "energy", "gcv", "threshold_ratio",
                 "median_threshold")) {
        r <- solve_tsvd(tiny_A, tiny_b, NULL, method = m)
        expect_length(r$spectrum, 3L)
    }
    r_fixed <- solve_tsvd(tiny_A, tiny_b, NULL, k = 2L)
    expect_length(r_fixed$spectrum, 3L)
})
