#' GEE (Generalized Estimation Equation) unfolding with robust inference
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_gee.py}.
#'
#' The unfolded spectrum is obtained from the generalized estimating
#' equations
#' \deqn{U(x) = A^T R(\alpha)^{-1} (b - A x) - \lambda G x = 0,}
#' where the m detector spheres are treated as a correlated cluster of
#' repeated measurements with an m x m working correlation matrix
#' \code{R(alpha)}. Three correlation structures are supported:
#' \describe{
#'   \item{independence}{\eqn{R = I} (Liang & Zeger's first working model)}
#'   \item{exchangeable}{\eqn{R_{ij} = \alpha} for \eqn{i \neq j},
#'     with the classical Liang-Zeger moment estimator of \eqn{\alpha}
#'     from the Pearson residuals}
#'   \item{ar1}{AR-1 structure \eqn{R_{ij} = \alpha^{|i-j|}} estimated
#'     from the lag-1 products of the Pearson residuals}
#' }
#'
#' Family / link pairings follow the \code{gee()} semantics
#' (\code{family=gaussian} by default, \code{poisson} and \code{gamma}
#' variants), implemented as quasi-likelihood variance functions:
#' \describe{
#'   \item{gaussian}{\eqn{v(\mu) = 1}}
#'   \item{poisson}{\eqn{v(\mu) = \mu}}
#'   \item{gamma}{\eqn{v(\mu) = \mu^2}}
#' }
#'
#' Because the number of energy bins exceeds the number of spheres, the
#' score equations alone do not define a unique solution; a ridge on the
#' second-difference roughness penalty \eqn{G = D_2^T D_2} (regularisation
#' \code{lam} relative to the mean diagonal of \eqn{A^T R^{-1} A}) is used.
#'
#' Statistical inference (the selling point of \code{gee} over plain GLS)
#' is delivered by the two canonical sandwich covariance estimators of
#' Liang & Zeger: robust ("sandwich") and naive (model-based) variance.
#'
#' @name gee-methods
#' @rdname gee-methods
NULL

.TINY <- 1e-300

.gee_validate_family_corstr <- function(family, corstr) {
    family <- tolower(as.character(family))
    corstr <- tolower(as.character(corstr))
    families <- c("gaussian", "poisson", "gamma")
    corstrings <- c("independence", "exchangeable", "ar1")
    if (!(family %in% families))
        stop("family must be one of ", paste(families, collapse = ", "),
             ", got '", family, "'")
    if (!(corstr %in% corstrings))
        stop("corstr must be one of ", paste(corstrings, collapse = ", "),
             ", got '", corstr, "'")
    list(family = family, corstr = corstr)
}

#' Build the m x m working correlation matrix
#'
#' @param alpha Numeric; working correlation parameter.
#' @param m Integer; cluster size (number of detector spheres).
#' @param kind Character: \code{"exchangeable"}, \code{"ar1"} or
#'   \code{"independence"}.
#' @return Symmetric positive-definite working correlation matrix.
#' @export
#' @examples
#' R1 <- working_correlation(0.3, 5, "exchangeable")
#' R2 <- working_correlation(0.5, 5, "ar1")
#' R3 <- working_correlation(0.0, 5, "independence")
working_correlation <- function(alpha, m, kind = "exchangeable") {
    kind <- tolower(as.character(kind))
    if (!is.numeric(m) || length(m) != 1L || m < 1)
        stop("m must be a positive integer, got ", m)
    m <- as.integer(m)
    if (kind == "independence") return(diag(m))
    if (kind == "exchangeable") {
        lo <- -1.0 / (m - 1)
        if (!(lo < alpha && alpha < 1.0))
            stop("exchangeable alpha must be in (", lo, ", 1) for m=",
                 m, ", got ", alpha)
        R <- matrix(alpha, m, m); diag(R) <- 1.0
        return(R)
    }
    if (kind == "ar1") {
        if (!(-1.0 < alpha && alpha < 1.0))
            stop("ar1 alpha must be in (-1, 1), got ", alpha)
        idx <- abs(outer(seq_len(m) - 1L, seq_len(m) - 1L, "-"))
        return(alpha^idx)
    }
    stop("unknown working correlation kind '", kind, "'")
}

#' Moment estimates of the working correlation parameter
#'
#' @param r_pearson Numeric vector of Pearson residuals.
#' @param kind Character: \code{"exchangeable"} (default), \code{"ar1"}
#'   or \code{"independence"}.
#' @return A list \code{list(alpha, phi)}.
#' @export
#' @examples
#' r <- c(-0.3, 0.1, 0.4, -0.2, 0.1)
#' est <- estimate_alpha(r, "exchangeable")
estimate_alpha <- function(r_pearson, kind = "exchangeable") {
    r <- as.numeric(r_pearson)
    m <- length(r)
    kind <- tolower(as.character(kind))
    if (m < 2L) return(list(alpha = 0.0, phi = 1.0))
    phi <- sum(r * r) / m
    if (kind == "independence") return(list(alpha = 0.0, phi = phi))
    if (kind == "exchangeable") {
        R0 <- outer(r, r)
        mask <- upper.tri(R0) | lower.tri(R0)
        num <- sum(R0[mask])
        den <- m * (m - 1)
        if (den <= .Machine$double.eps) return(list(alpha = 0.0, phi = phi))
        a <- num / den
    } else if (kind == "ar1") {
        num <- sum(r[-1] * r[-m])
        den <- sum(r[-m] * r[-m])
        if (den <= .Machine$double.eps) return(list(alpha = 0.0, phi = phi))
        a <- num / den
    } else {
        stop("unknown working correlation kind '", kind, "'")
    }
    lo <- -0.95 / max(m - 1L, 1L)
    a <- min(max(a, lo), 0.95)
    list(alpha = a, phi = phi)
}

.gee_variance_mu <- function(mu, family) {
    if (family == "gaussian") return(rep(1.0, length(mu)))
    if (family == "poisson") return(pmax(mu, .TINY))
    if (family == "gamma") return(pmax(mu * mu, .TINY))
    stop("unknown family '", family, "'")
}

#' Fit the unfolding model by generalized estimating equations
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial spectrum guess.
#' @param family Character: \code{"gaussian"} (default), \code{"poisson"}
#'   or \code{"gamma"}.
#' @param corstr Character: \code{"exchangeable"} (default), \code{"ar1"}
#'   or \code{"independence"}.
#' @param regularization Numeric; relative ridge on the difference
#'   penalty. Default 1e-4. \code{0} disables it.
#' @param max_iterations Integer; max GEE iterations. Default 100.
#' @param tolerance Numeric; relative convergence tolerance. Default 1e-6.
#' @param diff_order Integer; order of the roughness penalty. Default 2.
#' @return A diagnostics list with \code{spectrum}, \code{cov_robust},
#'   \code{cov_naive}, \code{robust_se}, \code{naive_se}, \code{alpha},
#'   \code{phi}, \code{residuals}, \code{pearson_residuals},
#'   \code{pearson_chi2}, \code{df}, \code{iterations}, \code{converged},
#'   \code{family}, \code{corstr}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' d <- gee_fit(A, b, family = "gaussian", corstr = "exchangeable",
#'              max_iterations = 50)
gee_fit <- function(A, b, x0 = NULL, family = "gaussian",
                     corstr = "exchangeable", regularization = 1e-4,
                     max_iterations = 100L, tolerance = 1e-6,
                     diff_order = 2L) {
    fc <- .gee_validate_family_corstr(family, corstr)
    family <- fc$family; corstr <- fc$corstr
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A); n <- ncol(A)
    if (regularization < 0)
        stop("regularization must be non-negative, got ", regularization)

    if (is.null(x0)) {
        x <- if (family %in% c("poisson", "gamma")) rep(1, n) else rep(0, n)
    } else {
        x <- pmax(as.numeric(x0), 0.0)
    }

    # Roughness penalty, scale-normalised.
    if (regularization > 0 && diff_order > 0 && n > diff_order) {
        D <- as.matrix(create_derivative_matrix(n, diff_order))
        G <- t(D) %*% D
        g_norm <- max(mean(diag(G)), 1.0)
        G <- G / g_norm
    } else {
        G <- matrix(0.0, n, n)
    }

    alpha <- 0.0; phi <- 1.0
    converged <- FALSE
    it <- 0L
    delta_prev <- Inf
    repeat {
        it <- it + 1L
        if (it > as.integer(max_iterations)) { it <- it - 1L; break }
        Ax <- A %*% x
        mu <- pmax(as.numeric(Ax), .TINY)
        v <- .gee_variance_mu(mu, family)
        r_raw <- b - mu
        r_pear <- r_raw / sqrt(v)
        alpha_phi <- estimate_alpha(r_pear, corstr)
        alpha_new <- alpha_phi$alpha
        phi_new <- alpha_phi$phi

        R <- working_correlation(alpha_new, m, corstr)
        Rinv <- tryCatch(solve(R), error = function(e) qr.solve(R))
        AW <- t(A) %*% Rinv
        H <- AW %*% A
        h_norm <- max(mean(diag(H)), 1.0)
        damp <- regularization * h_norm
        lhs <- H + damp * G
        rhs <- AW %*% b
        x_new <- tryCatch(as.numeric(solve(lhs, rhs)),
                          error = function(e)
                              as.numeric(qr.solve(lhs, rhs, tol = 1e-10)))
        if (!all(is.finite(x_new))) {
            x_new <- x; converged <- FALSE; break
        }
        delta <- sqrt(sum((x_new - x)^2)) /
                  max(sqrt(sum(x_new^2)), sqrt(sum(x^2)), 1.0)
        x <- pmax(x_new, 0.0)
        alpha <- alpha_new; phi <- phi_new
        stalled <- delta >= delta_prev * (1.0 - 1e-12) && delta < 1e10
        if (delta <= tolerance || (stalled && it > 1L)) {
            converged <- TRUE; break
        }
        delta_prev <- delta
    }

    # Final diagnostics on the converged estimate.
    Ax <- A %*% x
    mu <- pmax(as.numeric(Ax), .TINY)
    v <- .gee_variance_mu(mu, family)
    r_raw <- b - mu
    r_pear <- r_raw / sqrt(v)
    R <- working_correlation(alpha, m, corstr)
    Rinv <- tryCatch(solve(R), error = function(e) qr.solve(R))
    AW <- t(A) %*% Rinv
    H <- AW %*% A
    h_norm <- max(mean(diag(H)), 1.0)
    damp <- regularization * h_norm
    bread <- H + damp * G
    Ninv <- tryCatch(solve(bread), error = function(e) qr.solve(bread))

    # Robust ("sandwich") covariance.
    grad_sq <- pmax(r_pear^2, 0.0)
    Au <- A * grad_sq
    meat <- t(A) %*% Au
    cov_robust <- Ninv %*% meat %*% Ninv
    # Naive (model-based) covariance.
    cov_naive <- Ninv %*% H %*% Ninv
    cov_robust <- 0.5 * (cov_robust + t(cov_robust))
    cov_naive <- 0.5 * (cov_naive + t(cov_naive))

    robust_se <- sqrt(pmax(diag(cov_robust), 0.0))
    naive_se <- sqrt(pmax(diag(cov_naive), 0.0))
    pearson_chi2 <- sum(r_pear^2)
    df <- max(m - n, 1L)

    list(spectrum = x,
         cov_robust = cov_robust, cov_naive = cov_naive,
         robust_se = robust_se, naive_se = naive_se,
         alpha = alpha, phi = phi,
         residuals = r_raw, pearson_residuals = r_pear,
         pearson_chi2 = pearson_chi2, df = df,
         iterations = as.integer(it), converged = converged,
         family = family, corstr = corstr)
}

#' Solve unfolding problem by generalized estimating equations
#'
#' @inheritParams gee_fit
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_gee(A, b, max_iterations = 50)
solve_gee <- function(A, b, x0 = NULL, family = "gaussian",
                       corstr = "exchangeable", regularization = 1e-4,
                       max_iterations = 100L, tolerance = 1e-6) {
    diag <- gee_fit(A, b, x0 = x0, family = family, corstr = corstr,
                     regularization = regularization,
                     max_iterations = max_iterations,
                     tolerance = tolerance)
    list(spectrum = diag$spectrum,
         iterations = diag$iterations,
         converged = diag$converged)
}

#' Wrapper around \code{\link{solve_gee}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams gee_fit
#' @export
unfold_gee <- function(detector_names, n_energy_bins, E_MeV, sensitivities,
                        cc_icrp116, save_result_callback, readings,
                        initial_spectrum = NULL, family = "gaussian",
                        corstr = "exchangeable", regularization = 1e-4,
                        max_iterations = 100L, tolerance = 1e-6,
                        calculate_errors = FALSE, noise_level = 0.01,
                        n_montecarlo = 100L, save_result = FALSE,
                        random_state = NULL,
                        max_neutron_energy = NULL) {
    x0_default <- rep(1, n_energy_bins)

    extra_output <- tryCatch({
        A_mat <- do.call(rbind, lapply(detector_names, function(nm)
            as.numeric(sensitivities[[nm]])))
        b_vec <- as.numeric(readings[detector_names])
        diag <- gee_fit(A_mat, b_vec, x0 = NULL, family = family,
                         corstr = corstr, regularization = regularization,
                         max_iterations = max_iterations,
                         tolerance = tolerance)
        list(alpha = diag$alpha, phi = diag$phi,
             family = diag$family, corstr = diag$corstr,
             robust_se = diag$robust_se, naive_se = diag$naive_se,
             spectrum_uncert_robust = diag$robust_se,
             pearson_chi2 = diag$pearson_chi2,
             gee_converged = diag$converged)
    }, error = function(e) NULL)

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_gee, family = family,
                                         corstr = corstr,
                                         regularization = regularization,
                                         max_iterations = max_iterations,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "GEE",
        extra_output = extra_output,
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
