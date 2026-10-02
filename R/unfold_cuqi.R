#' CUQI-style Bayesian MCMC unfolding (pure-R port)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_cuqi.py}.  The Python
#' original drives the Bayesian inverse-problem library CUQIpy; this port
#' implements the same hierarchical log-scale posterior with pure-R samplers:
#' \itemize{
#'   \item the spectrum is modelled as \code{x = exp(theta)} with a Gaussian
#'         (GMRF or Ornstein-Uhlenbeck) smoothness prior centered on a
#'         data-driven NNLS guess and a Gaussian likelihood with per-detector
#'         sigma \code{noise_level * |b|};
#'   \item \code{sampler = "pcn"} runs a preconditioned Crude-Lanczos
#'         (pCN) chain on the centered log-spectrum \code{theta - mu};
#'   \item \code{"cwmh"} runs a centered adaptive random-walk Metropolis;
#'   \item \code{"mala"} / \code{"ula"} run Laplace-whitened (Gauss-Newton
#'         MAP) Metropolis-adjusted / unadjusted Langevin chains;
#'   \item \code{"gibbs"} / \code{"gibbs_nuts"} run the hierarchical scheme
#'         where the GMRF precision has a Gamma hyperprior
#'         (\code{delta ~ Gamma(delta_alpha, delta_beta)}) drawn with its
#'         conjugate update and the spectrum with pCN / random-walk MH;
#'   \item \code{"nuts"} is mapped onto \code{"cwmh"} with a warning
#'         (no pure-R NUTS implementation).
#' }
#'
#' @section Limitations vs. Python bssunfold:
#' The CUQIpy backend (and NUTS itself) is replaced by statistically
#' equivalent pure-R samplers, so draws differ chain-by-chain; acceptance
#' rates, credible intervals and the Gelman-Rubin diagnostics are computed
#' in R (CUQIpy's ESS is replaced by a Geyer initial-positive-sequence
#' approximation).
#'
#' @name cuqi-methods
NULL

# Supported sampler identifiers (public API), mirroring _VALID_SAMPLERS.
.cuqi_valid_samplers <- c("pcn", "cwmh", "ula", "mala", "nuts",
                          "gibbs", "gibbs_nuts")
.cuqi_langevin_samplers <- c("ula", "mala")
.cuqi_gibbs_samplers <- c("gibbs", "gibbs_nuts")
# Default proposal scales per sampler (_DEFAULT_SCALES); for the whitened
# Langevin samplers the scale is in whitened units.
.cuqi_default_scales <- c(pcn = 0.05, cwmh = 0.05, nuts = 0.05, mala = 0.3,
                          ula = 0.01, gibbs = 0.05, gibbs_nuts = 0.05)

# ---- priors ---------------------------------------------------------------

# Dense Ornstein-Uhlenbeck correlation C[i, j] = exp(-|i - j| / lengthscale).
.cuqi_ou_correlation <- function(n_bins, lengthscale) {
    idx <- seq_len(n_bins) - 1L
    corr <- exp(-abs(outer(idx, idx, "-")) / max(as.numeric(lengthscale), 1e-9))
    corr + 1e-9 * diag(n_bins)
}

.cuqi_ou_precision <- function(n_bins, lengthscale) {
    solve(.cuqi_ou_correlation(n_bins, lengthscale)) +
        1e-9 * diag(n_bins)
}

# First/second-order difference precision with zero boundary anchors,
# mirroring the CUQIpy GMRF finite-difference operator (_gmrf_precision).
.cuqi_gmrf_precision <- function(n_bins, order) {
    order <- as.integer(order)
    if (!order %in% c(1L, 2L)) {
        stop("gmrf_order must be 1 or 2, got ", order)
    }
    rows <- max(n_bins - order, 1L)
    diff <- matrix(0, nrow = rows, ncol = n_bins)
    for (i in seq_len(rows)) {
        diff[i, i] <- 1
        if (order == 1L) {
            diff[i, i + 1L] <- -1
        } else {
            diff[i, i + 1L] <- -2
            diff[i, i + 2L] <- 1
        }
    }
    precision <- crossprod(diff)
    precision[1, 1] <- precision[1, 1] + 1
    if (order == 2L) precision[2, 2] <- precision[2, 2] + 1
    precision + 1e-9 * diag(n_bins)
}

# Data-driven log-space prior center (_prior_center_nnls): the user
# initial_spectrum when given, otherwise the NNLS solution of A @ x = b.
.cuqi_prior_center <- function(A_matrix, b_readings, initial_spectrum, n_energy) {
    if (!is.null(initial_spectrum)) {
        center <- pmax(as.numeric(initial_spectrum), 0)
        if (!is.vector(initial_spectrum) || length(center) != n_energy) {
            center <- numeric(n_energy)
        }
    } else {
        center <- tryCatch(.bss_nnls(A_matrix, b_readings),
                           error = function(e) NULL)
        if (is.null(center)) {
            center <- pmax(as.numeric(qr.solve(A_matrix, b_readings)), 0)
        }
        center <- pmax(as.numeric(center), 0)
        if (!any(center > 0)) center <- rep(1, n_energy)
    }
    log(pmax(center, 1e-6))
}

# ---- log-posterior building blocks ----------------------------------------

# Log-likelihood of the linear forward model b = A exp(theta) with diagonal
# Gaussian noise sigma2, and the gradient of the full log-posterior
# grad = exp(theta) * (A' (resid / sigma2)) - Qp (theta - mu).
.cuqi_loglik <- function(theta, A_matrix, b_readings, sigma2) {
    resid <- b_readings - as.numeric(A_matrix %*% exp(theta))
    -0.5 * sum(resid^2 / sigma2)
}

.cuqi_grad_logpost <- function(theta, A_matrix, b_readings, sigma2, mu, Qp) {
    e <- exp(theta)
    resid <- b_readings - as.numeric(A_matrix %*% e)
    as.numeric(e * (crossprod(A_matrix, resid / sigma2))) -
        as.numeric(Qp %*% (theta - mu))
}

.cuqi_logpost <- function(theta, A_matrix, b_readings, sigma2, mu, Qp) {
    dev <- theta - mu
    .cuqi_loglik(theta, A_matrix, b_readings, sigma2) -
        0.5 * as.numeric(crossprod(dev, Qp %*% dev))
}

# ---- Gauss-Newton MAP + Laplace whitening (_gauss_newton_map) -------------

.cuqi_gauss_newton_map <- function(A_matrix, b_readings, sigma2, mu, Qp,
                                   n_iter = 25L) {
    n_energy <- ncol(A_matrix)
    theta <- mu
    logp_now <- .cuqi_logpost(theta, A_matrix, b_readings, sigma2, mu, Qp)
    for (it in seq_len(as.integer(n_iter))) {
        e <- exp(theta)
        resid <- b_readings - as.numeric(A_matrix %*% e)
        jac <- A_matrix * matrix(e, nrow = nrow(A_matrix), ncol = n_energy,
                                 byrow = FALSE)
        curvature <- Qp + crossprod(jac, jac / sigma2)
        grad <- as.numeric(e * (crossprod(A_matrix, resid / sigma2))) -
            as.numeric(Qp %*% (theta - mu))
        step <- tryCatch(solve(curvature + 1e-10 * diag(n_energy), grad),
                         error = function(e) NULL)
        if (is.null(step)) break
        alpha <- 1
        improved <- FALSE
        for (ls in seq_len(25L)) {
            candidate <- theta + alpha * step
            cand_logp <- .cuqi_logpost(candidate, A_matrix, b_readings,
                                       sigma2, mu, Qp)
            if (is.finite(cand_logp) &&
                cand_logp >= logp_now - 1e-4 * alpha * sum(grad * step)) {
                improved <- TRUE
                break
            }
            alpha <- alpha * 0.5
        }
        if (!improved || alpha < 1e-9) break
        theta <- candidate
        logp_now <- cand_logp
        if (sqrt(sum((alpha * step)^2)) < 1e-12) break
    }
    e <- exp(theta)
    resid <- b_readings - as.numeric(A_matrix %*% e)
    jac <- A_matrix * matrix(e, nrow = nrow(A_matrix), ncol = n_energy,
                             byrow = FALSE)
    curvature <- Qp + crossprod(jac, jac / sigma2)
    laplace_cov <- solve(curvature + 1e-10 * diag(n_energy))
    whitening <- t(chol(laplace_cov + 1e-10 * diag(n_energy)))
    list(theta = theta, whitening = whitening)
}

# ---- sampler chains (one chain each) ---------------------------------------

# pCN on the centered log-spectrum t = theta - mu: proposal
# t' = sqrt(1 - scale^2) t + scale * xi,  xi ~ N(0, prior_cov) (prior-
# invariant), accepted on the likelihood ratio (CUQIpy PCN).
.cuqi_run_pcn_chain <- function(A_matrix, b_readings, sigma2, mu, prior_cov,
                                scale, n_samples, n_burnin, thin) {
    n_energy <- ncol(A_matrix)
    Lc <- t(chol(as.matrix(prior_cov) + 1e-12 * diag(n_energy)))
    t_cur <- numeric(n_energy)
    ll_cur <- .cuqi_loglik(mu + t_cur, A_matrix, b_readings, sigma2)
    beta <- sqrt(max(1 - scale^2, 0))
    n_iter <- as.integer(n_burnin) + as.integer(n_samples) * max(as.integer(thin), 1L)
    keep_every <- max(as.integer(thin), 1L)
    n_stored <- max(as.integer(n_samples) %/% keep_every, 1L)
    samples <- matrix(0, nrow = n_stored, ncol = n_energy)
    stored <- 0L
    accepted <- 0L
    steps_kept <- 0L
    for (it in seq_len(n_iter)) {
        prop <- beta * t_cur + as.numeric(Lc %*% rnorm(n_energy)) * scale
        ll_prop <- .cuqi_loglik(mu + prop, A_matrix, b_readings, sigma2)
        take <- is.finite(ll_prop) &&
            (ll_prop >= ll_cur || log(runif(1)) < ll_prop - ll_cur)
        if (take) {
            t_cur <- prop
            ll_cur <- ll_prop
            if (it > n_burnin) accepted <- accepted + 1L
        }
        if (it > n_burnin && (it - n_burnin) %% keep_every == 0L) {
            steps_kept <- steps_kept + 1L
            if (steps_kept > n_stored) steps_kept <- n_stored  # tail trim
            samples[steps_kept, ] <- t_cur
        }
    }
    list(samples = samples,
         acc = if (steps_kept > 0L) accepted / steps_kept else NaN)
}

# Centered adaptive random-walk Metropolis (CUQIpy CWMH, simplified:
# step-size feedback to the 0.234 acceptance sweet spot during warmup).
.cuqi_run_cwmh_chain <- function(lp_fn, t0, scale, n_samples, n_burnin, thin,
                                 adapt = TRUE) {
    n_energy <- length(t0)
    t_cur <- as.numeric(t0)
    lp_cur <- lp_fn(t_cur)
    keep_every <- max(as.integer(thin), 1L)
    n_stored <- max(as.integer(n_samples) %/% keep_every, 1L)
    samples <- matrix(0, nrow = n_stored, ncol = n_energy)
    n_iter <- as.integer(n_burnin) + as.integer(n_samples) * keep_every
    accepted <- 0L
    steps_kept <- 0L
    win_acc <- 0L
    win_n <- 0L
    for (it in seq_len(n_iter)) {
        prop <- t_cur + rnorm(n_energy, sd = scale)
        lp_prop <- lp_fn(prop)
        take <- is.finite(lp_prop) &&
            (lp_prop >= lp_cur || log(runif(1)) < lp_prop - lp_cur)
        if (take) {
            t_cur <- prop
            lp_cur <- lp_prop
            if (it > n_burnin) accepted <- accepted + 1L
        }
        if (it > n_burnin && (it - n_burnin) %% keep_every == 0L) {
            steps_kept <- steps_kept + 1L
            if (steps_kept > n_stored) steps_kept <- n_stored
            samples[steps_kept, ] <- t_cur
        }
        if (adapt && it <= n_burnin) {
            win_acc <- win_acc + as.integer(take)
            win_n <- win_n + 1L
            if (win_n >= 50L) {
                rate <- win_acc / win_n
                if (rate > 0.441) scale <- scale * 1.1
                if (rate < 0.441 && rate < 0.114) scale <- scale / 1.1
                else if (rate < 0.234) scale <- scale / 1.05
                win_acc <- 0L
                win_n <- 0L
            }
        }
    }
    list(samples = samples,
         acc = if (steps_kept > 0L) accepted / steps_kept else NaN)
}

# Langevin chain in the whitened coordinate z (theta = center + W %*% z):
# ULA is the plain Euler-Maruyama update, MALA adds the Metropolis
# correction with the asymmetric Gaussian proposal density.
.cuqi_run_langevin_chain <- function(A_matrix, b_readings, sigma2, mu, Qp,
                                     center, whitening, scale, sampler,
                                     n_samples, n_burnin, thin) {
    n_energy <- ncol(A_matrix)
    lp_z <- function(z) {
        .cuqi_logpost(center + as.numeric(whitening %*% z),
                      A_matrix, b_readings, sigma2, mu, Qp)
    }
    grad_z <- function(z) {
        t <- center + as.numeric(whitening %*% z)
        as.numeric(crossprod(whitening,
            .cuqi_grad_logpost(t, A_matrix, b_readings, sigma2, mu, Qp)))
    }
    z_cur <- numeric(n_energy)
    lp_cur <- lp_z(z_cur)
    keep_every <- max(as.integer(thin), 1L)
    n_stored <- max(as.integer(n_samples) %/% keep_every, 1L)
    samples <- matrix(0, nrow = n_stored, ncol = n_energy)
    n_iter <- as.integer(n_burnin) + as.integer(n_samples) * keep_every
    accepted <- 0L
    steps_kept <- 0L
    noise_sd <- sqrt(2 * scale)
    for (it in seq_len(n_iter)) {
        g <- grad_z(z_cur)
        prop <- z_cur + scale * g + rnorm(n_energy, sd = noise_sd)
        lp_prop <- lp_z(prop)
        if (sampler == "mala") {
            gp <- grad_z(prop)
            logq_fwd <- -sum((prop - z_cur - scale * g)^2) / (4 * scale)
            logq_back <- -sum((z_cur - prop - scale * gp)^2) / (4 * scale)
            delta <- lp_prop - lp_cur + logq_back - logq_fwd
        } else {
            delta <- Inf  # ULA: unadjusted, always accept
        }
        take <- is.finite(lp_prop) &&
            (delta >= 0 || log(runif(1)) < delta)
        if (take) {
            z_cur <- prop
            lp_cur <- lp_prop
            if (it > n_burnin) accepted <- accepted + 1L
        }
        if (it > n_burnin && (it - n_burnin) %% keep_every == 0L) {
            steps_kept <- steps_kept + 1L
            if (steps_kept > n_stored) steps_kept <- n_stored
            samples[steps_kept, ] <- as.numeric(whitening %*% z_cur)
        }
    }
    acc <- if (sampler == "ula") NaN else
        if (steps_kept > 0L) accepted / steps_kept else NaN
    list(samples = samples, acc = acc)
}

# Hierarchical Gibbs chain: delta ~ Gamma(delta_alpha + n/2,
# delta_beta + 0.5 * t' Q1 t) conjugate update, theta | delta sampled with
# pCN ('gibbs') or random-walk MH ('gibbs_nuts', the NUTS stand-in).
.cuqi_run_gibbs_chain <- function(A_matrix, b_readings, sigma2, mu, Q1,
                                  Q1_cov_chol, sampler, delta_alpha,
                                  delta_beta, scale, n_samples, n_burnin,
                                  thin) {
    n_energy <- ncol(A_matrix)
    centered <- (sampler == "gibbs")
    beta <- sqrt(max(1 - scale^2, 0))
    t_cur <- numeric(n_energy)
    quad <- function(t) as.numeric(crossprod(t, Q1 %*% t))
    delta <- stats::rgamma(1, shape = delta_alpha + n_energy / 2,
                           rate = delta_beta + quad(t_cur) / 2)
    Lc <- Q1_cov_chol  # already t(chol(Q1^-1)) = lower factor of the cov
    lp_cur <- function(delta_now) {
        .cuqi_loglik(mu + t_cur, A_matrix, b_readings, sigma2) -
            0.5 * delta_now * quad(t_cur)
    }
    keep_every <- max(as.integer(thin), 1L)
    n_stored <- max(as.integer(n_samples) %/% keep_every, 1L)
    theta_samples <- matrix(0, nrow = n_stored, ncol = n_energy)
    delta_samples <- numeric(n_stored)
    n_iter <- as.integer(n_burnin) + as.integer(n_samples) * keep_every
    accepted <- 0L
    steps_kept <- 0L
    for (it in seq_len(n_iter)) {
        # 1. conjugate Gamma update of the precision hyperparameter
        delta <- stats::rgamma(1, shape = delta_alpha + n_energy / 2,
                               rate = delta_beta + quad(t_cur) / 2)
        # 2. one spectral sweep at fixed delta
        if (centered) {
            prop <- beta * t_cur +
                (as.numeric(Lc %*% rnorm(n_energy)) / sqrt(delta)) * scale
        } else {
            # RW step sized to the current prior width (1/sqrt(delta)):
            # symmetric given delta, so MH remains valid.
            prop <- t_cur + rnorm(n_energy, sd = scale / sqrt(delta))
        }
        take <- FALSE
        if (centered) {
            ll_cur <- .cuqi_loglik(mu + t_cur, A_matrix, b_readings, sigma2)
            ll_prop <- .cuqi_loglik(mu + prop, A_matrix, b_readings, sigma2)
            take <- is.finite(ll_prop) &&
                (ll_prop >= ll_cur || log(runif(1)) < ll_prop - ll_cur)
        } else {
            lp_c <- lp_cur(delta)
            prop_lp <- .cuqi_loglik(mu + prop, A_matrix, b_readings, sigma2) -
                0.5 * delta * quad(prop)
            take <- is.finite(prop_lp) &&
                (prop_lp >= lp_c || log(runif(1)) < prop_lp - lp_c)
        }
        if (take) {
            t_cur <- prop
            if (it > n_burnin) accepted <- accepted + 1L
        }
        if (it > n_burnin && (it - n_burnin) %% keep_every == 0L) {
            steps_kept <- steps_kept + 1L
            if (steps_kept > n_stored) steps_kept <- n_stored  # tail trim
            theta_samples[steps_kept, ] <- t_cur
            delta_samples[steps_kept] <- delta
        }
    }
    list(samples = theta_samples + matrix(mu, nrow = n_stored,
                                          ncol = n_energy, byrow = TRUE),
         delta = delta_samples,
         acc = if (steps_kept > 0L) accepted / steps_kept else NaN)
}

# ---- diagnostics -----------------------------------------------------------

# Classic Gelman-Rubin R-hat (_gelman_rubin), input array m x n x p.
.cuqi_gelman_rubin <- function(chains_stack) {
    m <- dim(chains_stack)[1]
    n <- dim(chains_stack)[2]
    p <- dim(chains_stack)[3]
    if (m < 2L || n < 2L) return(rep(1, p))
    means <- apply(chains_stack, c(1, 3), mean)   # m x p
    vars <- apply(chains_stack, c(1, 3), stats::var)
    grand <- colMeans(means)
    between <- n / (m - 1) * colSums(sweep(means, 2, grand)^2)
    within <- colMeans(vars)
    var_hat <- (n - 1) / n * within + between / n
    rhat <- suppressWarnings(sqrt(var_hat / within))
    if (is.null(dim(rhat))) rhat <- as.numeric(rhat)
    ifelse(is.finite(rhat), rhat, 1)
}

# Shortest credible (HPD) interval per energy bin (_hpd_interval).
.cuqi_hpd_interval <- function(samples, prob) {
    xs <- as.matrix(samples)
    n_total <- nrow(xs)
    n_keep <- max(ceiling(prob * n_total), 1L)
    bounds <- apply(xs, 2, function(col) {
        s <- sort(as.numeric(col))
        if (n_keep >= n_total) return(c(s[1], s[n_total]))
        widths <- s[n_keep:n_total] - s[seq_len(n_total - n_keep + 1L)]
        i <- which.min(widths)
        c(s[i], s[i + n_keep - 1L])
    })
    list(lower = bounds[1, ], upper = bounds[2, ])
}

# Effective sample size per column, Geyer initial-positive-sequence
# approximation (CUQIpy's compute_ess is not available in the R port).
.cuqi_ess <- function(theta_all) {
    x <- as.matrix(theta_all)
    n <- nrow(x)
    p <- ncol(x)
    max_pair <- min(100L, max(n %/% 4L, 1L))
    out <- numeric(p)
    for (j in seq_len(p)) {
        v <- x[, j]
        v <- v - mean(v)
        var0 <- sum(v^2) / n
        if (var0 <= 0 || !is.finite(var0)) {
            out[j] <- n
            next
        }
        s <- 0
        for (m_idx in seq_len(max_pair)) {
            k1 <- 2 * m_idx - 1L
            k2 <- 2 * m_idx
            if (k2 > n %/% 2L) break
            pair <- (sum(v[seq_len(n - k1)] * v[(k1 + 1L):n]) +
                     sum(v[seq_len(n - k2)] * v[(k2 + 1L):n])) / n / var0
            if (pair < 0) break
            s <- s + pair
        }
        out[j] <- min(max(n / (1 + 2 * s), 1), n)
    }
    out
}

# ---- public low-level solver -----------------------------------------------

#' Solve the unfolding problem with the pure-R CUQI-style Bayesian samplers
#'
#' R port of \code{bssunfold.core.unfold_cuqi.solve_cuqi_bayesian}.  See
#' \code{\link{cuqi-methods}} for the model and sampler description.
#'
#' @param A_matrix Numeric response matrix (n_detectors x n_energy).
#' @param b_readings Numeric measured readings (length n_detectors).
#' @param E Optional energy grid (unused, kept for API consistency).
#' @param log_steps Optional log energy steps (unused, kept for API
#'   consistency; the forward model follows \code{b = A @ spectrum}).
#' @param sampler One of \code{"pcn"} (default), \code{"cwmh"},
#'   \code{"mala"}, \code{"ula"}, \code{"nuts"} (mapped onto
#'   \code{"cwmh"} with a warning), \code{"gibbs"}, \code{"gibbs_nuts"}.
#' @param noise_level Relative likelihood noise scale; sigma = noise_level *
#'   |b| (default 0.05).
#' @param prior \code{"gmrf"} (default) or \code{"ou"}.
#' @param gmrf_order Order of the GMRF difference operator, 1 or 2.
#' @param lengthscale OU correlation length in bins for \code{prior = "ou"}.
#' @param prec Fixed prior precision scale (inferred from the data by the
#'   hierarchical samplers).
#' @param hierarchical Force the hierarchical Gibbs scheme; NULL derives it
#'   from \code{sampler}.
#' @param delta_alpha,delta_beta Shape / rate of the Gamma hyperprior on the
#'   GMRF precision.
#' @param n_samples,n_burnin Number of posterior samples / warmup iterations
#'   per chain.
#' @param thin Thinning interval.
#' @param chains Number of independent chains.
#' @param scale Proposal step size; NULL uses the per-sampler defaults.
#' @param max_depth,step_size Accepted for API parity with the CUQIpy NUTS
#'   options (unused by the pure-R samplers).
#' @param credible_level Credible mass (percent) of the HPD interval.
#' @param initial_spectrum Optional prior-center guess; NULL uses NNLS.
#' @param random_state Optional integer seed (per chain: seed + chain).
#' @param progressbar Unused, kept for API consistency.
#' @return A list \code{list(spectrum, stats, iterations, converged)} where
#'   \code{stats} carries the posterior samples and diagnostics (mean,
#'   median, std, hpd_lower, hpd_upper, ess, rhat, acc_rate,
#'   delta_samples, sampling metadata).
#' @rdname cuqi-methods
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_cuqi_bayesian(A, b, sampler = "pcn", n_samples = 60L,
#'                          n_burnin = 30L, chains = 1L, random_state = 42)
#' length(r$spectrum)
solve_cuqi_bayesian <- function(A_matrix, b_readings, E = NULL,
                                log_steps = NULL, sampler = "pcn",
                                noise_level = 0.05, prior = "gmrf",
                                gmrf_order = 1L, lengthscale = 3.0,
                                prec = 1.0, hierarchical = NULL,
                                delta_alpha = 1.0, delta_beta = 1e-4,
                                n_samples = 2000L, n_burnin = 1000L,
                                thin = 1L, chains = 2L, scale = NULL,
                                max_depth = 8L, step_size = NULL,
                                credible_level = 95.0,
                                initial_spectrum = NULL,
                                random_state = NULL,
                                progressbar = FALSE) {
    sampler_l <- tolower(as.character(sampler))
    if (!sampler_l %in% .cuqi_valid_samplers) {
        stop("Unknown sampler '", sampler, "'. Valid options: ",
             paste(.cuqi_valid_samplers, collapse = ", "))
    }
    prior_l <- tolower(as.character(prior))
    if (!prior_l %in% c("gmrf", "ou")) {
        stop("Unknown prior '", prior, "'. Valid options: 'gmrf', 'ou'")
    }
    if (is.null(hierarchical)) {
        hierarchical <- sampler_l %in% .cuqi_gibbs_samplers
    }
    hierarchical <- as.logical(hierarchical)
    if (hierarchical && !sampler_l %in% .cuqi_gibbs_samplers) {
        sampler_l <- "gibbs"
    }
    if (hierarchical && prior_l != "gmrf") {
        stop("Hierarchical sampling (Gamma hyperprior on the precision) is ",
             "only implemented for the 'gmrf' prior; set prior='gmrf' or ",
             "use a non-hierarchical sampler.")
    }
    if (sampler_l == "nuts") {
        warning("NUTS is not available in the pure-R CUQI port; falling ",
                "back to the adaptive random-walk sampler (cwmh).")
        sampler_l <- "cwmh"
    }
    if (is.null(scale)) scale <- as.numeric(.cuqi_default_scales[[sampler_l]])
    prob <- as.numeric(credible_level) / 100
    if (!(prob > 0 && prob < 1)) {
        stop("credible_level must be in (0, 100) percent, got ",
             credible_level)
    }

    A_matrix <- as.matrix(A_matrix); storage.mode(A_matrix) <- "double"
    b_readings <- as.numeric(b_readings)
    n_detectors <- nrow(A_matrix)
    n_energy <- ncol(A_matrix)

    mu <- .cuqi_prior_center(A_matrix, b_readings, initial_spectrum, n_energy)

    # Likelihood noise: relative scale sigma = noise_level * |b| per detector
    sigma2 <- (as.numeric(noise_level) * (abs(b_readings) + 1e-6))^2

    n_samples <- as.integer(n_samples)
    n_burnin <- as.integer(n_burnin)
    thin <- max(as.integer(thin), 1L)
    chains <- as.integer(chains)

    theta_per_chain <- vector("list", chains)
    delta_parts <- vector("list", chains)
    acc_rates <- numeric(chains)

    Q1 <- .cuqi_gmrf_precision(n_energy, as.integer(gmrf_order))
    if (prior_l == "gmrf") {
        Qp <- as.numeric(prec) * Q1
        prior_cov <- solve(Q1) / as.numeric(prec)
    } else {
        corr <- .cuqi_ou_correlation(n_energy, lengthscale)
        Qp <- as.numeric(prec) * (solve(corr) + 1e-9 * diag(n_energy))
        prior_cov <- corr / as.numeric(prec)
    }

    for (c in seq_len(chains)) {
        if (!is.null(random_state)) {
            set.seed(as.integer(random_state) + (c - 1L))
        }
        if (hierarchical) {
            sigma2_scalar <- rep(mean(sigma2), n_detectors)
            out <- .cuqi_run_gibbs_chain(
                A_matrix, b_readings, sigma2_scalar, mu, Q1,
                Q1_cov_chol = t(chol(solve(Q1))), sampler = sampler_l,
                delta_alpha = delta_alpha, delta_beta = delta_beta,
                scale = scale, n_samples = n_samples, n_burnin = n_burnin,
                thin = thin)
            theta_per_chain[[c]] <- out$samples
            delta_parts[[c]] <- out$delta
            acc_rates[c] <- out$acc
        } else if (sampler_l %in% .cuqi_langevin_samplers) {
            gn <- .cuqi_gauss_newton_map(A_matrix, b_readings, sigma2, mu, Qp)
            out <- .cuqi_run_langevin_chain(
                A_matrix, b_readings, sigma2, mu, Qp, gn$theta, gn$whitening,
                scale = scale, sampler = sampler_l, n_samples = n_samples,
                n_burnin = n_burnin, thin = thin)
            theta_per_chain[[c]] <- out$samples +
                matrix(gn$theta, nrow = nrow(out$samples), ncol = n_energy,
                       byrow = TRUE)
            acc_rates[c] <- out$acc
        } else if (sampler_l == "pcn") {
            out <- .cuqi_run_pcn_chain(A_matrix, b_readings, sigma2, mu,
                                       prior_cov, scale = scale,
                                       n_samples = n_samples,
                                       n_burnin = n_burnin, thin = thin)
            theta_per_chain[[c]] <- out$samples +
                matrix(mu, nrow = nrow(out$samples), ncol = n_energy,
                       byrow = TRUE)
            acc_rates[c] <- out$acc
        } else {  # cwmh on the centered log-spectrum
            lp_fn <- function(t) {
                .cuqi_logpost(mu + t, A_matrix, b_readings, sigma2, mu, Qp)
            }
            out <- .cuqi_run_cwmh_chain(lp_fn, rep(0, n_energy),
                                        scale = scale, n_samples = n_samples,
                                        n_burnin = n_burnin, thin = thin)
            theta_per_chain[[c]] <- out$samples +
                matrix(mu, nrow = nrow(out$samples), ncol = n_energy,
                       byrow = TRUE)
            acc_rates[c] <- out$acc
        }
    }

    theta_all <- do.call(rbind, theta_per_chain)
    n_stored <- nrow(theta_per_chain[[1]])
    delta_all <- if (hierarchical) unlist(delta_parts) else NULL
    acc_rate <- mean(acc_rates[is.finite(acc_rates)])
    if (length(acc_rate) == 0L) acc_rate <- NaN

    x_samples <- exp(theta_all)
    mean_spectrum <- colMeans(x_samples)
    median_spectrum <- apply(x_samples, 2, stats::median)
    # np.std: population standard deviation
    std_spectrum <- apply(x_samples, 2,
                          function(col) sqrt(mean((col - mean(col))^2)))
    hpd <- .cuqi_hpd_interval(x_samples, prob)

    ess <- .cuqi_ess(theta_all)
    rhat <- NULL
    if (chains > 1L && n_stored > 1L) {
        stack <- array(as.numeric(theta_all), dim = c(chains, n_stored,
                                                      n_energy))
        rhat <- .cuqi_gelman_rubin(stack)
    }

    stats <- list(
        samples = x_samples,
        theta_samples = theta_all,
        mean = mean_spectrum,
        median = median_spectrum,
        std = std_spectrum,
        hpd_lower = hpd$lower,
        hpd_upper = hpd$upper,
        ess = ess,
        rhat = rhat,
        acc_rate = acc_rate,
        delta_samples = delta_all,
        sampler = sampler_l,
        prior = prior_l,
        hierarchical = hierarchical,
        n_samples_total = nrow(theta_all),
        n_chains = chains,
        n_samples = n_samples,
        n_burnin = n_burnin,
        thin = thin,
        noise_level = as.numeric(noise_level),
        credible_level = as.numeric(credible_level),
        gmrf_order = as.integer(gmrf_order),
        lengthscale = as.numeric(lengthscale),
        prec = as.numeric(prec),
        delta_alpha = as.numeric(delta_alpha),
        delta_beta = as.numeric(delta_beta),
        scale = as.numeric(scale),
        prior_center = exp(mu),
        backend = "R-pure (CUQI-style)"
    )

    list(spectrum = as.numeric(mean_spectrum), stats = stats,
         iterations = nrow(theta_all), converged = TRUE)
}

# ---- public workflow wrapper ------------------------------------------------

#' Wrapper around \code{\link{solve_cuqi_bayesian}} for the unified workflow
#'
#' Mirrors \code{bssunfold.core.unfold_cuqi.unfold_cuqi}: builds the system
#' matrix from the detector readings, runs the Bayesian MCMC model and
#' returns the standardized unfolding result enriched with
#' \code{spectrum_uncertainty}, \code{spectrum_lower}/\code{spectrum_upper}
#' (HPD credible band) and the full \code{cuqi_stats} diagnostics.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_cuqi_bayesian
#' @param ln_steps Numeric log energy steps. Unused by the Bayesian model
#'   (the forward model follows \code{b = A @ spectrum}); accepted for API
#'   consistency with the other workflows.
#' @param reading_uncertainties,reading_covariance,noise_model,measurement_time
#'   Accepted for API parity with the Monte-Carlo error model; unused.
#' @param mc_noise_level Relative Gaussian noise level for the Monte-Carlo
#'   uncertainty pass when \code{calculate_errors = TRUE}.
#' @return A result list as produced by \code{\link{run_unfolding}} plus
#'   \code{spectrum_uncertainty}, \code{spectrum_lower},
#'   \code{spectrum_upper} and \code{cuqi_stats}.
#' @rdname cuqi-methods
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' detector_names <- c("D1", "D2", "D3")
#' sens <- list(D1 = A[1, ], D2 = A[2, ], D3 = A[3, ])
#' readings <- setNames(as.numeric(A %*% rep(1, 3)), detector_names)
#' r <- unfold_cuqi(detector_names, 3L, c(0.1, 1, 5), sens, NULL,
#'                  function(out) invisible(out), readings,
#'                  sampler = "gibbs", n_samples = 60L, n_burnin = 30L,
#'                  chains = 1L, random_state = 7)
#' length(r$spectrum)
unfold_cuqi <- function(detector_names, n_energy_bins, E_MeV,
                        sensitivities, cc_icrp116, save_result_callback,
                        readings, ln_steps = NULL,
                        initial_spectrum = NULL,
                        sampler = "pcn",
                        noise_level = 0.05,
                        prior = "gmrf",
                        gmrf_order = 1L,
                        lengthscale = 3.0,
                        prec = 1.0,
                        hierarchical = NULL,
                        delta_alpha = 1.0,
                        delta_beta = 1e-4,
                        n_samples = 2000L,
                        n_burnin = 1000L,
                        thin = 1L,
                        chains = 2L,
                        scale = NULL,
                        max_depth = 8L,
                        step_size = NULL,
                        credible_level = 95.0,
                        calculate_errors = FALSE,
                        mc_noise_level = 0.01,
                        n_montecarlo = 100L,
                        save_result = FALSE,
                        random_state = NULL,
                        progressbar = FALSE,
                        reading_uncertainties = NULL,
                        reading_covariance = NULL,
                        noise_model = "gaussian",
                        measurement_time = NULL,
                        max_neutron_energy = NULL) {
    # The main solve is captured in a holder so its posterior statistics
    # (samples, diagnostics) can be merged into the standardized output.
    holder <- new.env(parent = emptyenv())

    # NOTE: unlike iterative methods, the Bayesian prior center must NOT
    # fall back to run_unfolding's default x0: when the user provides no
    # initial_spectrum the data-driven NNLS center is used instead.
    solve_wrapper <- function(A, b, x0 = NULL, ...) {
        out <- solve_cuqi_bayesian(
            A_matrix = A, b_readings = b, E = E_MeV,
            log_steps = rep(1, n_energy_bins),
            sampler = sampler, noise_level = noise_level, prior = prior,
            gmrf_order = gmrf_order, lengthscale = lengthscale, prec = prec,
            hierarchical = hierarchical, delta_alpha = delta_alpha,
            delta_beta = delta_beta, n_samples = n_samples,
            n_burnin = n_burnin, thin = thin, chains = chains,
            scale = scale, max_depth = max_depth, step_size = step_size,
            credible_level = credible_level,
            initial_spectrum = initial_spectrum,
            random_state = random_state, progressbar = progressbar)
        if (is.null(holder$stats)) holder$stats <- out$stats
        out
    }

    result <- run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1, n_energy_bins),
        solve_func = solve_wrapper, solve_kwargs = list(),
        method_name = "CUQI-Bayesian",
        extra_output = list(sampler = sampler, prior = prior,
                            n_samples = n_samples, n_burnin = n_burnin,
                            chains = chains),
        calculate_errors = calculate_errors,
        noise_level = mc_noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = FALSE,
        max_neutron_energy = max_neutron_energy)

    # Merge the posterior statistics into the standardized result *before*
    # it is handed to the history callback so saved results carry the full
    # output (pattern shared with unfold_mcmc).
    if (!is.null(holder$stats)) {
        stats <- holder$stats
        result$spectrum_uncertainty <- as.numeric(stats$std)
        result$spectrum_lower <- pmax(stats$hpd_lower, 0)
        result$spectrum_upper <- as.numeric(stats$hpd_upper)
        result$cuqi_stats <- stats
    }

    if (isTRUE(save_result) && is.function(save_result_callback)) {
        save_result_callback(result)
    }

    result
}
