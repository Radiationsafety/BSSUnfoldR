test_that("Detector construction from PTB works", {
    det <- make_ptb_detector()
    expect_true(inherits(det, "Detector"))
    expect_equal(det$n_energy_bins, 60L)
    expect_true("0in" %in% det$detector_names)
    expect_length(det$sensitivities[["0in"]], 60L)
    expect_named(det$cc_icrp116,
                 c("E_MeV", "AP", "PA", "LLAT", "RLAT", "ROT", "ISO"))
})

test_that("Detector$save_result appends to history", {
    det <- make_ptb_detector()
    det$history <- list()
    fake_result <- list(method = "TEST", spectrum = rep(1, 60),
                         residual_norm = 0)
    det$save_result(fake_result)
    expect_length(det$history, 1L)
    expect_equal(det$history[[1]]$method, "TEST")
})

test_that("validate_readings filters to detector_names and rejects negatives", {
    rds <- c("0in" = 100, "3in" = 200, "12in" = 50)
    v <- validate_readings(rds, c("0in", "3in", "12in"))
    expect_equal(unname(v), c(100, 200, 50))
    expect_error(validate_readings(c("0in" = -1), "0in"), "negative")
    expect_error(validate_readings(c("0in" = Inf), "0in"), "infinite")
    expect_error(validate_readings(c("0in" = NaN), "0in"), "NA")
    expect_error(validate_readings(c("0in" = 100), c("12in")),
                 "No valid detector readings")
})

test_that("validate_energy_grid checks positivity and strict monotonicity", {
    expect_error(validate_energy_grid(c(-1, 1, 2)), "positive")
    expect_error(validate_energy_grid(c(1, 1, 2)), "strictly increasing")
    expect_silent(validate_energy_grid(c(1, 2, 3), min_points = 3L))
})

test_that("validate_spectrum rejects mismatched length and NaN", {
    E <- c(1e-9, 1e-6, 1e-3)
    expect_error(validate_spectrum(c(1, 2), E), "length")
    expect_error(validate_spectrum(c(1, NaN, 3), E), "NA")
    expect_error(validate_spectrum(c(1, -2, 3), E), "negative")
})

test_that("validate_system matches x0 length to A column count", {
    A <- matrix(c(1, 0.5, 0.2, 0.4), nrow = 2)
    expect_error(validate_system(A, c(1, 0.6), x0 = c(1, 2, 3)),
                 "Length of x0")
    v <- validate_system(A, c(1, 0.6), x0 = c(1, 2),
                         max_iterations = 10L, tolerance = 1e-6)
    expect_equal(v$x0, c(1, 2))
})

test_that("calculate_dose_rates returns one entry per geometry", {
    cc <- ICRP116_COEFF_EFFECTIVE_DOSE()
    spec <- rep(1, length(cc$E_MeV))
    d <- calculate_dose_rates(spec, cc)
    expected_geoms <- setdiff(names(cc), "E_MeV")
    expect_named(d, expected_geoms)
    expect_true(all(d >= 0))
})

test_that("get_coefficients returns the right datasets", {
    cc116 <- get_coefficients("ICRP116")
    expect_true("AP" %in% names(cc116))
    cc74  <- get_coefficients("ICRP74_effective")
    expect_true("ROT" %in% names(cc74))
    expect_error(get_coefficients("garbage"), "Unknown dose coefficient name")
})

test_that("interpolate_coefficients returns vectors on the target grid", {
    cc <- get_coefficients("NRB99_2009_effective")
    E_target <- 10^seq(-9, 2, length.out = 25)
    ic <- interpolate_coefficients(cc, E_target)
    expect_equal(length(ic$E_MeV), 25L)
    expect_equal(length(ic$AP), 25L)
})

test_that("create_derivative_matrix has the right shape for orders 1 and 2", {
    L1 <- create_derivative_matrix(6L, 1L)
    L2 <- create_derivative_matrix(6L, 2L)
    L3 <- create_derivative_matrix(3L, 2L)  # tiny case (was a regression bug)
    expect_equal(dim(L1), c(5L, 6L))
    expect_equal(dim(L2), c(4L, 6L))
    expect_equal(dim(L3), c(1L, 3L))
    # second-difference stencil
    expect_equal(as.numeric(L3[1, ]),
                 c(1, -2, 1))
})

test_that("make_regularization_operator returns identity/dense correctly", {
    expect_equal(make_regularization_operator(3L, 0L, identity_for_zero = TRUE),
                 diag(3))
    expect_null(make_regularization_operator(3L, 0L, identity_for_zero = FALSE))
    expect_equal(dim(make_regularization_operator(5L, 1L)), c(4L, 5L))
    expect_equal(dim(make_regularization_operator(5L, 2L)), c(3L, 5L))
})

test_that("monte_carlo_uncertainty returns expected statistics", {
    # Trivial solver: just return b (length 2) repeated
    toy_solver <- function(rds, ...) as.numeric(rds)
    rds <- c(d1 = 1, d2 = 0.6)
    mc <- monte_carlo_uncertainty(toy_solver, rds, 0.05, 10L, 2L, 42L)
    expected <- c("spectrum_uncert_mean", "spectrum_uncert_std",
                 "spectrum_uncert_min", "spectrum_uncert_max",
                 "spectrum_uncert_median",
                 "spectrum_uncert_percentile_5",
                 "spectrum_uncert_percentile_95",
                 "spectrum_uncert_all")
    expect_true(all(expected %in% names(mc)))
    expect_equal(dim(mc$spectrum_uncert_all), c(10L, 2L))
})
