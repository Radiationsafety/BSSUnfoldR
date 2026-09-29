## ---- numpy legacy RandomState emulation ---------------------------------
## bssunfold's Python EKI drives its ensemble from
## ``numpy.random.RandomState(random_state)``, i.e. the legacy MT19937 with
## the legacy polar-method gaussian and a private cached draw.  R's own
## ``rnorm`` stream is a different generator *and* a different
## element-to-ensemble-member assignment, so an R port built on it can only
## ever agree statistically, never pointwise.  These helpers reproduce the
## legacy stream bit-for-bit so a given ``random_state`` gives the same
## ensemble in both languages.
##
## The MT state is held as doubles in [0, 2^32) because R integers cannot
## represent -2^31 (INT_MIN), which the 0x80000000 mask needs.  Bitwise
## and/xor go through R's vectorised integer operators on two 16-bit halves,
## which is exact for values below 2^32.
.NP_2P32 <- 4294967296          # 2^32
.NP_2P53 <- 9007199254740992    # 2^53
.NP_2P16 <- 65536

.np_bw_and <- function(a, b) {
    lo <- bitwAnd(as.integer(a %% .NP_2P16), as.integer(b %% .NP_2P16))
    hi <- bitwAnd(as.integer(a %/% .NP_2P16), as.integer(b %/% .NP_2P16))
    hi * .NP_2P16 + lo
}
.np_bw_xor <- function(a, b) {
    lo <- bitwXor(as.integer(a %% .NP_2P16), as.integer(b %% .NP_2P16))
    hi <- bitwXor(as.integer(a %/% .NP_2P16), as.integer(b %/% .NP_2P16))
    hi * .NP_2P16 + lo
}
## logical right shift (floor, so it is a true >> for non-negative operands)
.np_sr <- function(x, k) floor(x / 2^k)
## left shift truncated to 32 bits
.np_sl <- function(x, k) (x * 2^k) %% .NP_2P32
## exact (a * b) mod 2^32; a plain double product would exceed 2^53
.np_mm32 <- function(a, b) {
    r0 <- (b * (a %% .NP_2P16)) %% .NP_2P32
    r1 <- (b * (a %/% .NP_2P16)) %% .NP_2P32
    (r0 + r1 * .NP_2P16) %% .NP_2P32
}

.np_mt_seed <- function(seed) {
    key <- numeric(624L)
    key[1L] <- seed %% .NP_2P32
    for (i in 2:624) {
        prev <- key[i - 1L]
        key[i] <- (.np_mm32(.np_bw_xor(prev, .np_sr(prev, 30)), 1812433253) +
                     (i - 1L)) %% .NP_2P32
    }
    key
}
.np_mt_twist <- function(key) {
    ## The MT19937 twist is an in-place sequential recurrence, not a
    ## simultaneous one: past the N - M boundary (0-based index 227) each word
    ## reads the *already updated* word M - N positions earlier. A purely
    ## vectorised key[(i + 397) mod 624] therefore matches numpy only for the
    ## first 227 outputs. See the reference implementation in mt19937ar.c.
    N <- 624L
    M <- 397L
    for (i in seq_len(N - M)) {                    # 0-based 0..226
        j <- i + 1L
        y <- .np_bw_and(key[i], 2147483648) + .np_bw_and(key[j], 2147483647)
        ma <- if (.np_bw_and(y, 1) == 1) 2567483615 else 0
        key[i] <- .np_bw_xor(.np_bw_xor(key[i + M], .np_sr(y, 1)), ma)
    }
    for (i in (N - M + 1):(N - 1)) {               # 0-based 227..622
        j <- i + 1L
        y <- .np_bw_and(key[i], 2147483648) + .np_bw_and(key[j], 2147483647)
        ma <- if (.np_bw_and(y, 1) == 1) 2567483615 else 0
        key[i] <- .np_bw_xor(.np_bw_xor(key[i - (N - M)], .np_sr(y, 1)), ma)
    }
    i <- N                                            # 0-based 623
    y <- .np_bw_and(key[i], 2147483648) + .np_bw_and(key[1L], 2147483647)
    ma <- if (.np_bw_and(y, 1) == 1) 2567483615 else 0
    key[i] <- .np_bw_xor(.np_bw_xor(key[M], .np_sr(y, 1)), ma)
    key
}
.np_mt_temper <- function(y) {
    y <- .np_bw_xor(y, .np_sr(y, 11))
    y <- .np_bw_xor(y, .np_bw_and(.np_sl(y, 7), 2636928640))    # 0x9d2c5680
    y <- .np_bw_xor(y, .np_bw_and(.np_sl(y, 15), 4022730752))   # 0xefc60000
    .np_bw_xor(y, .np_sr(y, 18))
}

## A private numpy-legacy generator.  Never touches R's global RNG, exactly as
## RandomState does not touch numpy's global one.
.np_random_state <- function(random_state) {
    if (is.null(random_state)) {
        seed <- as.numeric(sample.int(.Machine$integer.max, 1L)) - 1
    } else {
        seed <- as.numeric(random_state)[1L]
        if (is.na(seed) || !is.finite(seed) || seed < 0) seed <- 0
    }
    key <- .np_mt_seed(seed)
    out <- numeric(624L)
    pos <- 624L
    has_gauss <- 0L
    cached <- 0

    u32 <- function() {
        if (pos >= 624L) {
            key <<- .np_mt_twist(key)
            out <<- .np_mt_temper(key)
            pos <<- 0L
        }
        v <- out[pos + 1L]
        pos <<- pos + 1L
        v
    }
    ## legacy_double(): 53-bit uniform from two uint32 draws
    dbl <- function() (.np_sr(u32(), 5) * 67108864 + .np_sr(u32(), 6)) / .NP_2P53
    ## legacy_gauss(): polar method, second draw cached across calls
    gauss <- function() {
        if (has_gauss == 1L) {
            has_gauss <<- 0L
            v <- cached
            cached <<- 0
            return(v)
        }
        repeat {
            x1 <- 2 * dbl() - 1
            x2 <- 2 * dbl() - 1
            r2 <- x1 * x1 + x2 * x2
            if (r2 < 1 && r2 != 0) break
        }
        f <- sqrt(-2 * log(r2) / r2)
        cached <<- f * x1
        has_gauss <<- 1L
        f * x2
    }
    gaussn <- function(n) {
        v <- numeric(n)
        for (i in seq_len(n)) v[i] <- gauss()
        v
    }
    ## ``uniform`` is ``random_sample()``: one legacy_double(), i.e. two
    ## uint32 draws.  Additive export so other ports (randomized Kaczmarz)
    ## can consume the same legacy stream without re-implementing MT19937.
    list(gauss = gauss, gaussn = gaussn, uniform = dbl, u32 = u32)
}

## numpy fills a (nrows x ncols) array in C order, so flat index
## (i - 1) * ncols + j carries the gaussian for element (i, j).  Building the
## matrix with nrow = ncols and transposing reproduces that mapping; a plain
## R matrix() would fill column-major and pair the draws with the wrong
## ensemble members.
.np_normal_mat <- function(rng, loc, scale, nrows, ncols) {
    g <- rng$gaussn(nrows * ncols)
    m <- t(matrix(g, nrow = ncols, ncol = nrows))
    as.numeric(loc) + as.numeric(scale) * m
}

#' Ensemble Kalman Inversion (EKI) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_eki.py}.
#' Implements the EKI method of Iglesias et al. (2013) for Bayesian posterior
#' approximation without MCMC. The ensemble of particles is propagated through
#' the forward model and updated via the Kalman gain equation with optional
#' regularization for stability.
#'
#' The ensemble and the perturbed observations are drawn from a bit-exact
#' emulation of numpy's legacy \code{RandomState}, so \code{random_state = k}
#' reproduces the Python result for the same \code{k}. As in Python the
#' generator is private: R's global RNG state is neither read beyond one seed
#' draw when \code{random_state} is \code{NULL}, nor reseeded.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial guess (length n), used as the centre of the
#'   initial ensemble.
#' @param n_ensemble Integer; number of ensemble members. Default 50.
#' @param n_iterations Integer; number of EKI iterations. Default 50.
#' @param regularization Numeric; Tikhonov-style regularization added to the
#'   covariance diagonal for numerical stability. Default 1e-4.
#' @param inflation Numeric; covariance inflation factor to prevent ensemble
#'   collapse. Default 1.02.
#' @param noise_std Optional numeric; standard deviation of measurement noise.
#'   If \code{NULL}, estimated as 5\% of \code{||b|| / sqrt(m)}.
#' @param random_state Optional integer seed.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' set.seed(7)
#' r <- solve_eki(A, b, rep(1, 3), n_ensemble = 20, n_iterations = 30)
solve_eki <- function(A, b, x0, n_ensemble = 50L, n_iterations = 50L,
                        regularization = 1e-4, inflation = 1.02,
                        noise_std = NULL, random_state = NULL) {
    v <- validate_system(A, b, x0 = x0)
    A <- v$A; b <- v$b; x0 <- v$x0
    m <- nrow(A); n <- ncol(A)
    rng <- .np_random_state(random_state)
    if (is.null(noise_std)) {
        noise_std <- if (m > 0) 0.05 * sqrt(sum(b^2)) / sqrt(m) else 1e-6
    }
    noise_var <- noise_std^2
    sigma_prior <- abs(x0) + 1e-6
    ensemble <- .np_normal_mat(rng, x0, sigma_prior, n, n_ensemble)
    denom <- max(as.integer(n_ensemble) - 1L, 1L)
    eye_m <- diag(m)
    for (it in seq_len(as.integer(n_iterations))) {
        predictions <- A %*% ensemble
        pred_mean <- rowMeans(predictions)
        state_mean <- rowMeans(ensemble)
        pred_pert <- predictions - pred_mean
        state_pert <- ensemble - state_mean
        C_dd <- (pred_pert %*% t(pred_pert)) / denom
        C_dd <- C_dd + (noise_var + regularization) * eye_m
        C_md <- (state_pert %*% t(pred_pert)) / denom
        C_d_inv <- tryCatch(solve(C_dd, eye_m), error = function(e) {
            # Pseudo-inverse fallback
            sv <- svd(C_dd, nu = m, nv = 0)
            s_inv <- ifelse(sv$d > 1e-10 * max(sv$d), 1 / sv$d, 0)
            sv$u %*% (s_inv * (t(sv$u)))
        })
        # Innovation: b + noise*perturbation - predictions
        noise <- .np_normal_mat(rng, b, rep(noise_std, m), m, n_ensemble)
        innovation <- noise - predictions
        ensemble <- ensemble + C_md %*% C_d_inv %*% innovation
        ensemble <- ensemble * inflation
        ensemble <- pmax(ensemble, 0)
    }
    mean_spectrum <- rowMeans(ensemble)
    list(spectrum = pmax(as.numeric(mean_spectrum), 0),
         iterations = as.integer(n_iterations),
         converged = TRUE,
         ensemble = ensemble)
}

#' Wrapper around \code{\link{solve_eki}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_eki
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_eki <- function(detector_names, n_energy_bins, E_MeV,
                        sensitivities, cc_icrp116, save_result_callback,
                        readings, initial_spectrum = NULL,
                        n_ensemble = 50L, n_iterations = 50L,
                        regularization = 1e-4, inflation = 1.02,
                        noise_std = NULL,
                        calculate_errors = FALSE,
                        noise_level = 0.01, n_montecarlo = 100L,
                        save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    x0_default <- rep(1.0 / max(n_energy_bins, 1L), n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_eki,
                                         n_ensemble = n_ensemble,
                                         n_iterations = n_iterations,
                                         regularization = regularization,
                                         inflation = inflation,
                                         noise_std = noise_std,
                                         random_state = random_state),
        solve_kwargs = list(),
        method_name = "EKI",
        extra_output = list(
            n_ensemble = as.integer(n_ensemble),
            n_iterations = as.integer(n_iterations),
            regularization = regularization,
            inflation = inflation
        ),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
