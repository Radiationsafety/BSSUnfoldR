# Tests for batch 3 of BSSUnfoldR (Bayes, Directed divergence, Express,
# Iterative refinement, BUNKI-UT, MLEM-STOP, StatReg, Tikhonov-TV, GKS,
# Crystal Ball).

test_that("solve_bayes converges and produces a non-negative spectrum", {
    r <- solve_bayes(tiny_A, tiny_b, tiny_x0_default, max_iterations = 50)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
    expect_true(all(is.finite(r$spectrum)))
})

test_that("solve_directed_divergence accepts smoothing options", {
    r0 <- solve_directed_divergence(tiny_A, tiny_b, tiny_x0_default,
                                     max_iterations = 20)
    r1 <- solve_directed_divergence(tiny_A, tiny_b, tiny_x0_default,
                                     max_iterations = 20,
                                     smoothness_order = 2L,
                                     smoothness_weight = 0.01)
    expect_length(r0$spectrum, 3L)
    expect_length(r1$spectrum, 3L)
    expect_true(all(r0$spectrum >= 0))
})

test_that("solve_directed_divergence rejects negative measurements", {
    expect_error(solve_directed_divergence(tiny_A, c(-1, 0.6, 0.4),
                                            tiny_x0_default),
                 "non-negative")
})

test_that("solve_express runs end-to-end with a small grid", {
    E <- 10^seq(-9, 1, length.out = 10)
    A <- matrix(runif(3 * 10), nrow = 3)
    b <- as.numeric(A %*% rep(0.5, 10))
    r <- solve_express(A, b, E, n_groups = 4, max_iterations = 3)
    expect_length(r$spectrum, 10L)
    expect_true(all(r$spectrum >= 0))
})

test_that("solve_express rejects non-increasing energy grid", {
    expect_error(solve_express(tiny_A, tiny_b, c(1e-3, 1e-6, 1e-9)),
                 "strictly increasing")
})

test_that("solve_iterative_refinement returns a spectrum and info", {
    r <- solve_iterative_refinement(tiny_A, tiny_b, tiny_x0_default,
                                     first_pass_kwargs = list(max_iterations = 20L),
                                     second_pass_kwargs = list(max_iterations = 20L))
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
    expect_true(is.list(r$info))
    expect_true(is.numeric(r$info$alpha))
    expect_true(is.numeric(r$info$final_residual))
})

test_that("solve_bunkiut preserves spectrum length", {
    r <- solve_bunkiut(tiny_A, tiny_b, tiny_x0_default, max_iterations = 30)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
})

test_that("solve_mlem_stop runs and reports convergence status", {
    r <- solve_mlem_stop(tiny_A, tiny_b, tiny_x0_default,
                          max_iterations = 100, j_threshold = 1e-3)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= 0))
    expect_true(is.logical(r$converged))
})

test_that("solve_statreg supports EmpiricalBayes and User methods", {
    r_eb <- solve_statreg(tiny_A, tiny_b, NULL,
                           unfoldermethod = "EmpiricalBayes")
    r_user <- solve_statreg(tiny_A, tiny_b, NULL,
                              unfoldermethod = "User", regularization = 1e-2)
    expect_length(r_eb$spectrum, 3L)
    expect_length(r_user$spectrum, 3L)
    expect_error(solve_statreg(tiny_A, tiny_b, NULL,
                                unfoldermethod = "garbage"),
                 "Unknown method")
})

test_that("solve_tikhonov_tv supports TT, TV and T types", {
    for (t in c("TT", "TV", "T")) {
        r <- solve_tikhonov_tv(tiny_A, tiny_b, NULL,
                                max_iterations = 20, type_ = t)
        expect_length(r$spectrum, 3L)
        expect_true(all(r$spectrum >= 0), info = t)
    }
    expect_error(solve_tikhonov_tv(tiny_A, tiny_b, NULL, type_ = "X"),
                 "Unsupported type_")
})

test_that("solve_gks runs with manual regularization fallback", {
    r <- solve_gks(tiny_A, tiny_b, NULL, max_iterations = 3L,
                    regularization_method = "manual", regularization = 1e-2)
    expect_length(r$spectrum, 3L)
    expect_lte(as.integer(r$iterations), 3L)
})

test_that("solve_gks supports gcv, dp, lcurve and manual", {
    for (m in c("gcv", "dp", "lcurve", "manual")) {
        r <- solve_gks(tiny_A, tiny_b, NULL, max_iterations = 3L,
                        regularization_method = m,
                        regularization = 1e-2, noise_level = 0.05)
        expect_length(r$spectrum, 3L)
    }
    expect_error(solve_gks(tiny_A, tiny_b, NULL,
                            regularization_method = "garbage"),
                 "Unsupported regularization")
})

test_that("solve_crystal_ball is a single-step method (iterations == 1)", {
    r <- solve_crystal_ball(tiny_A, tiny_b, NULL, regularization = 0.01)
    expect_equal(r$iterations, 1L)
    expect_true(r$converged)
    expect_true(all(r$spectrum >= 0))
})

test_that("solve_crystal_ball rejects all-zero measurements", {
    expect_error(solve_crystal_ball(tiny_A, c(0, 0, 0), NULL),
                 "All measurements are zero or negative")
})
