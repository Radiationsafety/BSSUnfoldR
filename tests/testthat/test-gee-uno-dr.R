# Tests for the additional v0.2.1 methods: GEE, Uno, Douglas-Rachford.

A <- matrix(c(0.9, 0.05, 0.05,
              0.10, 0.8, 0.10,
              0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
b <- c(1, 0.6, 0.4)
det_names <- c("d1", "d2", "d3")
sens <- lapply(setNames(det_names, det_names),
               function(n) A[which(det_names == n), ])
E <- c(1e-9, 1e-6, 1e-3)
readings <- setNames(b, det_names)

# ---------------------------------------------------------------------------
# GEE
# ---------------------------------------------------------------------------
test_that("working_correlation builds the three correlation structures", {
    R_ind <- working_correlation(0.0, 5, "independence")
    expect_equal(R_ind, diag(5))
    R_ex <- working_correlation(0.3, 5, "exchangeable")
    expect_equal(dim(R_ex), c(5, 5))
    expect_equal(diag(R_ex), rep(1, 5))
    expect_true(all(R_ex[upper.tri(R_ex)] == 0.3))
    R_ar1 <- working_correlation(0.5, 4, "ar1")
    expect_equal(R_ar1[1, 4], 0.5^3)
})

test_that("working_correlation validates alpha ranges", {
    expect_error(working_correlation(1.0, 5, "exchangeable"),
                 "alpha must be")
    expect_error(working_correlation(2.0, 5, "ar1"), "alpha must be")
})

test_that("estimate_alpha returns moment estimates for exchangeable/ar1", {
    r <- c(-0.3, 0.1, 0.4, -0.2, 0.1)
    est_ex <- estimate_alpha(r, "exchangeable")
    expect_true(est_ex$alpha >= -0.95 / 4)
    expect_true(est_ex$alpha <= 0.95)
    expect_true(est_ex$phi > 0)
    est_ar1 <- estimate_alpha(r, "ar1")
    expect_true(est_ar1$alpha >= -0.95)
    expect_true(est_ar1$alpha <= 0.95)
    est_ind <- estimate_alpha(r, "independence")
    expect_equal(est_ind$alpha, 0)
})

test_that("solve_gee returns a non-negative spectrum for the tiny system", {
    r <- solve_gee(A, b, max_iterations = 50)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= -1e-10))
    expect_true(is.logical(r$converged))
    expect_true(is.integer(r$iterations))
})

test_that("gee_fit returns rich diagnostics for all families and corstrs", {
    for (fam in c("gaussian", "poisson", "gamma")) {
        for (cs in c("independence", "exchangeable", "ar1")) {
            d <- gee_fit(A, b, family = fam, corstr = cs,
                         max_iterations = 30)
            expect_length(d$spectrum, 3L)
            expect_true(all(is.finite(d$spectrum)))
            expect_length(d$robust_se, 3L)
            expect_length(d$naive_se, 3L)
            expect_true(is.numeric(d$alpha))
            expect_true(is.numeric(d$phi))
            expect_equal(d$family, fam)
            expect_equal(d$corstr, cs)
            expect_true(is.numeric(d$pearson_chi2))
        }
    }
})

test_that("gee_fit rejects invalid family/corstr", {
    expect_error(gee_fit(A, b, family = "binomial"), "family must be")
    expect_error(gee_fit(A, b, corstr = "toeplitz"), "corstr must be")
})

test_that("unfold_gee wrapper returns a standardized result with GEE meta", {
    res <- unfold_gee(det_names, 3L, E, sens, NULL, NULL, readings,
                       family = "gaussian", corstr = "exchangeable",
                       max_iterations = 30)
    expect_length(res$spectrum, 3L)
    expect_true(all(res$spectrum >= -1e-10))
    expect_true(!is.null(res$alpha) || is.null(res$alpha))
    expect_equal(res$method, "GEE")
})

# ---------------------------------------------------------------------------
# Uno
# ---------------------------------------------------------------------------
test_that("uno_objective and uno_gradient are consistent numerically", {
    x <- c(0.5, 0.5, 0.5)
    f0 <- uno_objective(A, b, rep(1, 3), 0.01, x)
    g <- uno_gradient(A, b, rep(1, 3), 0.01, x)
    eps <- 1e-6
    for (i in 1:3) {
        xp <- x; xp[i] <- xp[i] + eps
        xm <- x; xm[i] <- xm[i] - eps
        fd <- (uno_objective(A, b, rep(1, 3), 0.01, xp) -
                   uno_objective(A, b, rep(1, 3), 0.01, xm)) / (2 * eps)
        expect_lt(abs(fd - g[i]), 1e-4)
    }
    expect_true(is.finite(f0))
})

test_that("uno_filter acceptability test", {
    filt <- list(c(0.5, 0.1))
    # better in both: acceptable
    expect_true(uno_filter(filt, 0.4, 0.05))
    # worse in both: not acceptable
    expect_false(uno_filter(filt, 0.6, 0.2))
    # empty filter: always acceptable
    expect_true(uno_filter(list(), 0.4, 0.05))
})

test_that("solve_uno_full returns full diagnostics for filter_sqp", {
    d <- solve_uno_full(A, b, preset = "filter_sqp")
    expect_length(d$spectrum, 3L)
    expect_equal(d$preset, "filter_sqp")
    expect_equal(d$hessian, "exact")
    expect_true(is.numeric(d$objective))
    expect_true(is.numeric(d$constraint_violation))
    expect_true(is.numeric(d$dual_infeasibility))
    expect_true(is.integer(d$n_iterations))
    expect_true(is.logical(d$converged))
})

test_that("solve_uno_full supports the ipopt_like preset (exact & bfgs)", {
    for (hess in c("exact", "bfgs")) {
        d <- solve_uno_full(A, b, preset = "ipopt_like", hessian = hess,
                             max_iterations = 100)
        expect_length(d$spectrum, 3L)
        expect_equal(d$preset, "ipopt_like")
        expect_equal(d$hessian, hess)
        expect_true(all(is.finite(d$spectrum)))
    }
})

test_that("solve_uno_full validates preset name and hessian mode", {
    expect_error(solve_uno_full(A, b, preset = "ipopt"),
                 "preset must be one of")
    expect_error(solve_uno_full(A, b, preset = "ipopt_like", hessian = "gauss"),
                 "hessian must be")
})

test_that("solve_uno supports 'uniform' and 'poisson' weights", {
    r1 <- solve_uno(A, b, weights = "uniform")
    r2 <- solve_uno(A, b, weights = "poisson")
    expect_length(r1$spectrum, 3L)
    expect_length(r2$spectrum, 3L)
})

test_that("solve_uno accepts explicit weight arrays", {
    r <- solve_uno(A, b, weights = c(1, 2, 3))
    expect_length(r$spectrum, 3L)
    expect_error(solve_uno(A, b, weights = c(1, -1, 2)),
                 "positive array")
    expect_error(solve_uno(A, b, weights = c(1, 2)),
                 "positive array of length")
})

test_that("unfold_uno wrapper returns a standardized result with Uno meta", {
    res <- unfold_uno(det_names, 3L, E, sens, NULL, NULL, readings,
                       preset = "filter_sqp")
    expect_length(res$spectrum, 3L)
    expect_true(all(res$spectrum >= -1e-10))
    expect_equal(res$method, "Uno (filter_sqp)")
})

# ---------------------------------------------------------------------------
# Douglas-Rachford
# ---------------------------------------------------------------------------
test_that("solve_douglas_rachford returns a non-negative spectrum", {
    r <- solve_douglas_rachford(A, b, rep(1, 3), max_iterations = 100)
    expect_length(r$spectrum, 3L)
    expect_true(all(r$spectrum >= -1e-10))
    expect_true(is.logical(r$converged))
})

test_that("solve_douglas_rachford supports pure L2 (use_tv = FALSE)", {
    r <- solve_douglas_rachford(A, b, rep(1, 3), use_tv = FALSE,
                                 max_iterations = 100)
    expect_length(r$spectrum, 3L)
    expect_true(all(is.finite(r$spectrum)))
})

test_that("unfold_douglas_rachford wrapper runs end-to-end", {
    res <- unfold_douglas_rachford(det_names, 3L, E, sens, NULL, NULL,
                                     readings, max_iterations = 50)
    expect_length(res$spectrum, 3L)
    expect_equal(res$method, "Douglas-Rachford")
})

# ---------------------------------------------------------------------------
# Detector registration of the new methods
# ---------------------------------------------------------------------------
test_that("Detector class exposes unfold_gee, unfold_uno, unfold_douglas_rachford", {
    det <- Detector$new(response_function = "PTB",
                         detector_names = c("0in", "3in", "5in", "8in", "12in"))
    expect_true("unfold_gee" %in% names(det))
    expect_true("unfold_uno" %in% names(det))
    expect_true("unfold_douglas_rachford" %in% names(det))
})
