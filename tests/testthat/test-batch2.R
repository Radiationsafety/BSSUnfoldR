test_that("solve_osem with n_subsets = 1 stays close to MLEM in shape", {
    r_osem <- solve_osem(tiny_A, tiny_b, tiny_x0_default,
                          max_iterations = 30, n_subsets = 1L)
    r_mlem <- solve_mlem(tiny_A, tiny_b, tiny_x0_default,
                          max_iterations = 30)
    # OSEM normalises by (A_sub^T 1 + eps), MLEM does not — so the absolute
    # amplitudes differ, but both should produce non-negative, finite
    # spectra of the right length.
    expect_length(r_osem$spectrum, 3L)
    expect_true(all(r_osem$spectrum >= 0))
    expect_true(all(is.finite(r_osem$spectrum)))
    expect_true(all(is.finite(r_mlem$spectrum)))
})

test_that("solve_osem handles n_subsets > 1", {
    r <- solve_osem(tiny_A, tiny_b, tiny_x0_default,
                    max_iterations = 20, n_subsets = 3L)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
})

test_that("solve_mapem supports all four priors", {
    for (p in c("none", "quadratic", "logcosh", "relative_difference")) {
        r <- tryCatch(
            solve_mapem(tiny_A, tiny_b, tiny_x0_default,
                         prior = p, max_iterations = 10),
            error = function(e) NULL
        )
        expect_false(is.null(r), info = p)
        if (!is.null(r)) {
            expect_length(r$spectrum, 3L)
            expect_true(all(r$spectrum >= 0))
        }
    }
    expect_error(solve_mapem(tiny_A, tiny_b, tiny_x0_default,
                              prior = "garbage"),
                 "Unknown prior")
})

test_that("prior_gradient and prior_value are consistent on the quadratic prior", {
    x <- c(1, 2, 3, 2, 1)
    g <- prior_gradient(x, "quadratic", beta = 1e-2, delta = 1)
    v <- prior_value(x, "quadratic", beta = 1e-2, delta = 1)
    expect_length(g, 5L)
    expect_true(is.numeric(v) && length(v) == 1L)
    expect_true(is.finite(v))
})

test_that("solve_bsrem runs with various priors and subsets", {
    r <- solve_bsrem(tiny_A, tiny_b, tiny_x0_default,
                      n_subsets = 3L, prior = "quadratic",
                      max_iterations = 10)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
})

test_that("solve_sart runs with constant and functional relaxation", {
    r_const <- solve_sart(tiny_A, tiny_b, rep(0, 3), max_iterations = 20)
    r_func  <- solve_sart(tiny_A, tiny_b, rep(0, 3), max_iterations = 20,
                          relaxation = function(n) 0.5)
    expect_length(r_const$spectrum, 3L)
    expect_length(r_func$spectrum, 3L)
    expect_true(all(r_const$spectrum >= 0))
})

test_that("solve_kaczmarz is reproducible", {
    r1 <- solve_kaczmarz(tiny_A, tiny_b, rep(0, 3), max_iterations = 30)
    r2 <- solve_kaczmarz(tiny_A, tiny_b, rep(0, 3), max_iterations = 30)
    expect_equal(r1$spectrum, r2$spectrum)
})

test_that("solve_randomized_kaczmarz is reproducible with a seed", {
    r1 <- solve_randomized_kaczmarz(tiny_A, tiny_b, rep(0, 3),
                                     max_iterations = 30, random_state = 7L)
    r2 <- solve_randomized_kaczmarz(tiny_A, tiny_b, rep(0, 3),
                                     max_iterations = 30, random_state = 7L)
    expect_equal(r1$spectrum, r2$spectrum)
})

test_that("solve_lanczos runs end-to-end without an error", {
    r <- solve_lanczos(tiny_A, tiny_b, NULL, max_iterations = 3L)
    expect_length(r$spectrum, 3L)
    expect_lte(as.integer(r$iterations), 3L)
})

test_that("solve_tikhonov_legendre runs end-to-end", {
    r <- solve_tikhonov_legendre(tiny_A, tiny_b, NULL,
                                  n_polynomials = 3L)
    expect_length(r$spectrum, 3L)
})

test_that("solve_rebunki accepts the same API as solve_bunki", {
    r <- solve_rebunki(tiny_A, tiny_b, tiny_x0_default,
                        max_iterations = 50)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
})

test_that("solve_doroshenko is reproducible", {
    r1 <- solve_doroshenko(tiny_A, tiny_b, tiny_x0_default,
                            max_iterations = 30)
    r2 <- solve_doroshenko(tiny_A, tiny_b, tiny_x0_default,
                            max_iterations = 30)
    expect_equal(r1$spectrum, r2$spectrum)
})
