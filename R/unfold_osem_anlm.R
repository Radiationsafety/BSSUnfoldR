#' OSEM-ANLM unfolding (ordered-subset EM with asymptotic non-local means)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_osem_anlm.py}.
#' Port of the OSEM-ANLM algorithm of Jamaati et al. (2026), "Enhanced sparse
#' view CT reconstruction using ordered subset expectation maximization and
#' asymptotic non-local means algorithms", Scientific Reports, adapted to
#' neutron spectrum unfolding from Bonner sphere readings.
#'
#' The algorithm alternates ordered-subset expectation maximisation (OSEM)
#' updates with the asymptotic non-local means (ANLM) filter applied to the
#' intermediate spectrum after every subset update (article pseudo-code steps
#' 4-5).  The ANLM filter is a two-stage non-local means (NLM) filter:
#' stage 1 applies the uniform parameter \eqn{h_1(i) = 0.5 \sigma}{h1(i) =
#' 0.5 * sigma} (the first filter applies h1 = 0.5 sigma uniformly), stage 2
#' applies the point-wise parameter \eqn{h_2(i) = \sigma_2(i) = \sqrt{\sum_j
#' w(i,j)^2 \sigma^2}}{h2(i) = sigma_2(i) = sqrt(sum_j w(i,j)^2 sigma^2)}
#' (article eq. 6), i.e. the noise standard deviation smoothed by the initial
#' NLM weights of the first stage.
#'
#' For one-dimensional spectra the 2D NLM windows of the article become index
#' windows on the energy grid: the search window \code{search_window} locates
#' similar neighbouring bins, and the similarity (patch) window
#' \code{similarity_window} with a Gaussian kernel of spread \code{alpha}
#' weights the squared-L2 patch distance
#' \deqn{d(i, j) = \sum_k G_\alpha(k) (x(i + k) - x(j + k))^2.}{d(i,j) =
#' sum_k G_alpha(k) (x(i+k) - x(j+k))^2,}
#' so the NLM weights are \eqn{w(i,j) = \exp(-d(i,j)/h^2)/z(i)}{w(i,j) =
#' exp(-d(i,j)/h^2) / z(i)} (article eq. 5 with the normalisation factor
#' \code{z(i)}).  With the article's optimal settings \code{N = 11},
#' \code{nu = 3} the filter preserves spectral structure (peaks and edges)
#' while suppressing the streak-like oscillations produced by sparse
#' (few-detector) data.
#'
#' The filter parameter \code{h} corresponds to the noise level \code{sigma}
#' of the reconstructed spectrum.  When \code{h = NULL} the noise level is
#' estimated automatically at every ANLM application from the current iterate
#' with a robust median-absolute-deviation estimator on the second
#' differences (\code{\link{estimate_noise_1d}}), which makes the method
#' scale-free for spectra spanning several orders of magnitude.  Because
#' Bonner-sphere spectra span several orders of magnitude (unlike CT images
#' in Hounsfield units) and the EM noise amplitude scales with the local
#' fluence, the filter by default operates on the logarithm of the spectrum
#' (\code{log_space = TRUE}); pass \code{log_space = FALSE} to reproduce the
#' raw-unit filtering of the original CT formulation.
#'
#' @param x Numeric input signal (length n).
#' @return Estimated noise standard deviation (strictly positive; a tiny
#'   relative floor is enforced so the value can safely be used as an NLM
#'   filter parameter).  For fewer than 3 samples Python returns
#'   \code{0.0}, which the caller then re-floors.
#' @keywords internal
#' @examples
#' \dontrun{
#' estimate_noise_1d(c(1, 1.2, 0.9, 1.1, 1.0))
#' }
estimate_noise_1d <- function(x) {
    x <- as.numeric(x)
    n <- length(x)
    if (n < 3L) return(0)
    d2 <- x[3:n] - 2 * x[2:(n - 1)] + x[seq_len(n - 2)]
    mad <- stats::median(abs(d2))
    sigma <- mad / (0.6745 * sqrt(6))
    floor_ <- 1e-12 * max(abs(x)) + 1e-300
    max(sigma, floor_)
}

# Map shifted indices back into [0, n) by symmetric reflection (Python
# _reflect_indices).  Returns a 0-based integer matrix with n rows and one
# column per offset.
.reflect_indices <- function(offsets, n) {
    if (n == 1L) return(matrix(0L, nrow = 1L, ncol = length(offsets)))
    idx <- outer(seq_len(n) - 1L, offsets, "+")
    period <- 2 * (n - 1L)
    idx <- abs(idx) %% period
    out <- ifelse(idx >= n, period - idx, idx)
    dim(out) <- dim(idx)
    out
}

#' Two-stage asymptotic non-local means filter for 1D spectra
#'
#' R port of \code{anlm_filter_1d} in
#' \code{bssunfold/src/bssunfold/core/unfold_osem_anlm.py}.
#'
#' \enumerate{
#'   \item Patch (similarity-window) distances \code{d(i, j)} are computed
#'         with a Gaussian kernel of spread \code{alpha} over the similarity
#'         window; indices outside the signal are reflected at the borders.
#'   \item Stage 1 applies NLM with the uniform parameter
#'         \code{h1 = 0.5 * sigma}, producing a lightly denoised intermediate
#'         spectrum and the "initial" normalised weights \code{w1(i, j)}.
#'   \item Stage 2 applies NLM with the point-wise parameter of article eq. 6,
#'         \code{h2(i) = sigma * sqrt(sum_j w1(i, j)^2)}, to the stage-1
#'         output, incrementally reducing the noise while preserving
#'         structure.
#' }
#'
#' By default (\code{log_space = TRUE}) the filter operates on the logarithm
#' of the spectrum, making the automatic \code{h} estimate scale-free (the
#' filtered value becomes a weighted geometric mean, which also preserves
#' non-negativity).  Set \code{log_space = FALSE} to filter in raw units
#' exactly as the original CT formulation of the article.
#'
#' @param x Numeric input spectrum (length n), typically non-negative.
#' @param h Optional noise level \code{sigma} used by both stages (in log
#'   units when \code{log_space = TRUE}).  Default \code{NULL} = estimate
#'   automatically with \code{\link{estimate_noise_1d}}.
#' @param search_window Positive integer; size \code{N} of the search window
#'   around each bin (article optimum: 11).  Even values are rounded down to
#'   the preceding odd size.
#' @param similarity_window Positive integer; size \code{nu} of the
#'   Gaussian-weighted similarity (patch) window (article optimum: 3).  Even
#'   values are rounded down to the preceding odd size.
#' @param alpha Positive numeric; spread of the Gaussian kernel over the
#'   similarity window.
#' @param log_space Logical; filter the logarithm of the spectrum instead of
#'   the raw values.  Default \code{TRUE}.
#' @return Filtered spectrum (numeric vector of length n).  Non-negative when
#'   \code{log_space = TRUE}; reduces to the identity for
#'   \code{search_window == 1}.
#' @keywords internal
#' @examples
#' \dontrun{
#' anlm_filter_1d(c(1, 1.4, 0.2, 1.3, 1.1), h = 0.1)
#' }
anlm_filter_1d <- function(x, h = NULL, search_window = 11L,
                           similarity_window = 3L, alpha = 1.0,
                           log_space = TRUE) {
    x <- as.numeric(x)
    n <- length(x)
    if (n == 0L) stop("input spectrum must be non-empty")
    search_window <- as.integer(search_window)
    similarity_window <- as.integer(similarity_window)
    if (search_window < 1L) {
        stop("search_window must be >= 1, got ", search_window)
    }
    if (similarity_window < 1L) {
        stop("similarity_window must be >= 1, got ", similarity_window)
    }
    if (!(alpha > 0)) stop("alpha must be positive, got ", alpha)
    if (!is.null(h) && !(as.numeric(h) > 0)) {
        stop("h must be a positive noise level, got ", h)
    }

    if (n == 1L) return(x)
    if (search_window %/% 2L == 0L) {
        # No neighbours in the search window: the NLM weights collapse to
        # the self bin, so the filter is the exact identity.
        return(x)
    }

    sigma <- if (!is.null(h)) as.numeric(h) else estimate_noise_1d(x)
    # Strictly positive floor: estimate_noise_1d returns 0.0 for fewer than
    # 3 samples, and a zero filter parameter would produce 0/0 weights.
    sigma <- max(sigma, 1e-12 * max(abs(x)) + 1e-300)

    work <- if (isTRUE(log_space)) {
        # Log domain: add a tiny relative floor so zero bins stay finite.
        log(x + 1e-12 * max(x) + 1e-300)
    } else {
        x
    }

    search_r <- search_window %/% 2L
    half_v <- similarity_window %/% 2L
    offsets <- seq(-half_v, half_v, by = 1L)
    gauss <- exp(-0.5 * (offsets / alpha)^2)
    gauss <- gauss / sum(gauss)

    # Reflected patch columns: P[i, t] = x[reflect(i + offsets[t])]
    refl <- .reflect_indices(offsets, n)
    nw <- length(offsets)
    patches <- matrix(x[refl + 1L], nrow = n, ncol = nw)

    .window_bounds <- function(i) {
        c(max(1L, i - search_r), min(n, i + search_r))
    }
    # One NLM pass with per-bin filter parameters `sigmas` (Python
    # _nlm_pass).
    nlm_pass <- function(signal, sigmas) {
        out <- numeric(n)
        sig_patches <- matrix(signal[refl + 1L], nrow = n, ncol = nw)
        for (i in seq_len(n)) {
            wb <- .window_bounds(i)
            lo <- wb[1L]; hi <- wb[2L]
            diff <- sweep(sig_patches[lo:hi, , drop = FALSE], 2L,
                          sig_patches[i, ], "-")
            # Gaussian-weighted patch distances (W_i,)
            dist <- as.numeric((diff^2) %*% gauss)
            # sqrt form avoids h^2 underflow for extremely small h; the
            # clip keeps the squared ratio finite when sigmas[i] is a
            # denormal (the weight is exp(-inf) = 0 either way).
            ratio <- pmin(sqrt(dist) / sigmas[i], 1e150)
            weights <- exp(-(ratio^2))
            z <- sum(weights)  # >= 1: the self-bin weight is exp(0) = 1
            out[i] <- sum(weights * signal[lo:hi]) / z
        }
        out
    }

    # Stage 1: uniform h1 = 0.5 * sigma (article, ANLM filter section).
    h1 <- 0.5 * sigma
    intermediate <- nlm_pass(work, rep(h1, n))
    w2_sq <- numeric(n)
    for (i in seq_len(n)) {
        wb <- .window_bounds(i)
        lo <- wb[1L]; hi <- wb[2L]
        diff <- sweep(patches[lo:hi, , drop = FALSE], 2L, patches[i, ], "-")
        dist <- as.numeric((diff^2) %*% gauss)
        ratio <- pmin(sqrt(dist) / h1, 1e150)
        weights <- exp(-(ratio^2))
        w1 <- weights / sum(weights)
        # Article eq. (6): h2(i) = sigma_2(i) = sqrt(sum_j w(i,j)^2 sigma^2)
        w2_sq[i] <- sum(w1^2)
    }
    h2 <- sigma * sqrt(w2_sq)

    # Stage 2: point-wise h2(i) applied to the stage-1 output.
    filtered <- nlm_pass(intermediate, h2)
    if (isTRUE(log_space)) exp(filtered) else filtered
}

#' Solve the unfolding problem with OSEM-ANLM
#'
#' Core solver mirroring \code{solve_osem_anlm} in
#' \code{bssunfold/src/bssunfold/core/unfold_osem_anlm.py}.
#'
#' OSEM update (article eq. 3 / pseudo-code step 4):
#' \deqn{x^{n+1} = x^n \frac{A_m^T (b_m / (A_m x^n + \epsilon))}
#'   {A_m^T 1 + \epsilon}}{
#'   x^(n+1) = x^n * A_m' (b_m / (A_m x^n + eps)) / (A_m' 1 + eps)}
#' followed by the ANLM filter (pseudo-code step 5,
#' \eqn{f^{*(n+1,b)} = \mathrm{ANLMFilter}(\mu^{*(n+1,b)})}{f* =
#' ANLMFilter(mu*)}), applied after every subset update when
#' \code{anlm_mode = "subset"} (default, per the article pseudo-code).  With
#' \code{anlm_mode = "post"} the plain OSEM solution is produced first and
#' the ANLM filter is applied once at the end ("OSEM reconstruction followed
#' by ANLM regularization" in the article abstract).  With
#' \code{n_subsets = 1} the OSEM update reduces to standard MLEM.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum guess (length n).
#' @param max_iterations Positive integer; maximum number of full iterations
#'   (sweeps over all subsets).  Default 50.
#' @param n_subsets Positive integer; number of ordered subsets over the
#'   detector readings.  Default 1 (standard MLEM with per-iteration ANLM).
#'   Must not exceed the number of detectors.
#' @param tolerance Positive numeric; relative change tolerance for early
#'   stopping.  Default 1e-6.
#' @param h Optional ANLM noise level; \code{NULL} = estimated from each
#'   intermediate spectrum.  Default \code{NULL}.
#' @param search_window ANLM search window \code{N}.  Default 11.
#' @param similarity_window ANLM similarity (patch) window \code{nu}.
#'   Default 3.
#' @param alpha Gaussian kernel spread over the similarity window.
#'   Default 1.0.
#' @param anlm_mode \code{'subset'} (default) or \code{'post'}; see above.
#' @param log_space Logical; apply the ANLM filter to the logarithm of the
#'   spectrum.  Default \code{TRUE}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_osem_anlm(A, b, c(0, 1, 1), max_iterations = 10)
solve_osem_anlm <- function(A, b, x0, max_iterations = 50L, n_subsets = 1L,
                            tolerance = 1e-6, h = NULL, search_window = 11L,
                            similarity_window = 3L, alpha = 1.0,
                            anlm_mode = "subset", log_space = TRUE) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    x0 <- as.numeric(x0)
    m <- nrow(A)

    n_subsets <- as.integer(n_subsets)
    if (n_subsets < 1L) stop("n_subsets must be >= 1")
    if (n_subsets > m) {
        stop("n_subsets (", n_subsets, ") must not exceed the number of detectors (",
             m, ")")
    }
    if (!anlm_mode %in% c("subset", "post")) {
        stop("anlm_mode must be one of 'subset', 'post', got '", anlm_mode, "'")
    }
    if (!is.null(h) && !(as.numeric(h) > 0)) {
        stop("h must be a positive noise level, got ", h)
    }
    search_window <- as.integer(search_window)
    similarity_window <- as.integer(similarity_window)
    if (search_window < 1L) {
        stop("search_window must be >= 1, got ", search_window)
    }
    if (similarity_window < 1L) {
        stop("similarity_window must be >= 1, got ", similarity_window)
    }
    if (!(alpha > 0)) stop("alpha must be positive, got ", alpha)

    eps <- 1e-11
    # np.array_split(np.arange(m), n_subsets): the first m %% n_subsets
    # subsets get ceil(m / n_subsets) rows, the rest floor(m / n_subsets).
    base <- m %/% n_subsets
    extra <- m %% n_subsets
    subset_rows <- vector("list", n_subsets)
    start <- 1L
    for (s in seq_len(n_subsets)) {
        len <- base + if (s <= extra) 1L else 0L
        subset_rows[[s]] <- if (len > 0L) start:(start + len - 1L) else integer(0)
        start <- start + len
    }

    x <- pmax(x0, 0)
    converged <- FALSE
    iterations <- 0L

    anlm <- function(current) {
        anlm_filter_1d(current, h = h, search_window = search_window,
                       similarity_window = similarity_window, alpha = alpha,
                       log_space = log_space)
    }

    for (it in seq_len(as.integer(max_iterations))) {
        iterations <- it
        x_old <- x

        for (idx in subset_rows) {
            if (length(idx) == 0L) next
            A_sub <- A[idx, , drop = FALSE]
            b_sub <- as.numeric(b[idx])
            norm <- as.numeric(colSums(A_sub))
            ratio <- b_sub / (as.numeric(A_sub %*% x) + eps)
            correction <- as.numeric(t(A_sub) %*% ratio)
            x <- pmax(x * correction / (norm + eps), 0)
            if (anlm_mode == "subset") {
                # Article pseudo-code step 5: ANLM after every subset update.
                x <- anlm(x)
            }
        }

        rel <- sqrt(sum((x - x_old)^2)) / (sqrt(sum(x_old^2)) + eps)
        if (rel < tolerance) {
            converged <- TRUE
            break
        }
    }

    if (anlm_mode == "post") {
        # Article abstract: "OSEM reconstruction followed by ANLM
        # regularization".
        x <- anlm(x)
    }

    list(spectrum = as.numeric(x), iterations = iterations,
         converged = converged)
}

#' OSEM-ANLM unfolding (unified workflow wrapper)
#'
#' Thin wrapper around \code{\link{solve_osem_anlm}} for the unified
#' workflow, mirroring \code{unfold_osem_anlm} in
#' \code{bssunfold/src/bssunfold/core/unfold_osem_anlm.py}.
#'
#' Ordered-subset expectation maximisation with asymptotic non-local means
#' regularization (Jamaati et al. 2026), adapted to Bonner sphere spectra:
#' the ANLM filter is applied to the intermediate spectrum after every OSEM
#' subset update (\code{anlm_mode = "subset"}) or once to the OSEM result
#' (\code{anlm_mode = "post"}).
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_osem_anlm
#' @param calculate_errors Logical; if \code{TRUE}, run Monte-Carlo
#'   uncertainty estimation.  Default \code{FALSE}.
#' @param noise_level Numeric; relative Gaussian noise level for Monte-Carlo.
#'   Default 0.01.
#' @param n_montecarlo Integer; number of Monte-Carlo samples.  Default 100.
#' @param random_state Optional integer seed for Monte-Carlo.
#' @param max_neutron_energy Optional numeric energy cutoff in MeV.
#'   Default \code{NULL} = no cutoff.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_osem_anlm <- function(detector_names, n_energy_bins, E_MeV,
                             sensitivities, cc_icrp116, save_result_callback,
                             readings, initial_spectrum = NULL,
                             max_iterations = 50L, n_subsets = 1L,
                             tolerance = 1e-6, h = NULL,
                             search_window = 11L, similarity_window = 3L,
                             alpha = 1.0, anlm_mode = "subset",
                             log_space = TRUE,
                             calculate_errors = FALSE,
                             noise_level = 0.01,
                             n_montecarlo = 100L,
                             save_result = FALSE, random_state = NULL,
                             max_neutron_energy = NULL) {
    x0_default <- rep(1.0, n_energy_bins)
    x0_default[1] <- 0.0

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_osem_anlm,
                                        max_iterations = max_iterations,
                                        n_subsets = n_subsets,
                                        tolerance = tolerance,
                                        h = h,
                                        search_window = search_window,
                                        similarity_window = similarity_window,
                                        alpha = alpha,
                                        anlm_mode = anlm_mode,
                                        log_space = log_space),
        solve_kwargs = list(),
        method_name = "OSEM-ANLM",
        extra_output = list(n_subsets = as.integer(n_subsets),
                            h = if (is.null(h)) NULL else as.numeric(h),
                            search_window = as.integer(search_window),
                            similarity_window = as.integer(similarity_window),
                            alpha = as.numeric(alpha),
                            anlm_mode = as.character(anlm_mode),
                            log_space = as.logical(log_space)),
        calculate_errors = calculate_errors,
        noise_level = noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
