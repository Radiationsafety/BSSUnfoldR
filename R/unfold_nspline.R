#' N-spline unfolding (Islamgulov & Lartsev, 2008)
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_nspline.py}.
#' The spectrum is parameterised as a "neutron spline"
#' \eqn{N_k(E) = \exp(a_k + q_k \ln E + r_k E)} on the segment
#' \eqn{E_k \le E \le E_{k+1}}, with \eqn{C^0/C^1} continuity imposed at the
#' interior knots through the constraint matrix \eqn{D}.  The spectrum is
#' recovered by minimising the directed divergence
#' \eqn{H = \sum_i p_{N,i}\ln(p_{N,i}/p_i) - p_{N,i} + p_i} between measured
#' and calculated normalised activations with the flux-conserving MIRD
#' gradient iteration, re-fitting (smoothing) the N-spline after *every*
#' iteration.
#'
#' @section Knots:
#' \code{knots = NULL} builds the automatic log-uniform *segment* grid from
#' \code{min(E)} to \code{max(E)} with \code{min(12, max(4, n %/% 4))}
#' segments, exactly like the Python implementation.  A character value
#' selects a preset from the internal \code{NSPLINE_KNOT_PRESETS} table, a
#' numeric value
#' is read as a full (endpoint inclusive) knot sequence and is clipped to the
#' energy range and extended to span the whole grid.  Note that this is
#' \emph{not} the interior-knot vector returned by \code{\link{auto_knots}}.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n); \code{NULL} means flat.
#' @param E_MeV Numeric energy grid (length n), strictly positive.
#' @param knots Optional knot specification: \code{NULL} (auto log-uniform
#'   grid), a preset name or a numeric knot sequence (see Sections).
#' @param max_iterations Positive integer; iteration budget. Default 500.
#'   When left at its default the reference budget of 200 iterations is used.
#' @param tolerance Positive numeric; relative stopping tolerance on the
#'   directed divergence \eqn{H}. Default 1e-6; when left at its default the
#'   reference tolerance 1e-3 is used.
#' @param dmu_start Numeric; conservative initial step factor,
#'   \code{dmu = dmu_start / sup|R - Rbar|}. Default 0.1.
#' @param continuity Character; \code{"C0C1"} (default), \code{"C0"} or
#'   \code{"none"}.
#' @param sigma_rel Optional numeric relative measurement uncertainties
#'   (length m) used by the \eqn{H}-target stopping criterion and the
#'   \code{nev} statistic; \code{NULL} means 0.1 for every detector.
#' @param smoothing Logical; re-fit the N-spline after every iteration.
#'   Default TRUE (the paper's procedure).
#' @param n_segments Optional integer; number of spline segments for the auto
#'   knot grid (\code{NULL} = adaptive).
#' @return A list \code{list(spectrum, iterations, converged, stop_reason,
#'   H, H_history, H_target, nev, nev_limit, acceptable, Qr,
#'   relative_residuals, fluence, mean_energy, knots, knots_source,
#'   continuity, params)}.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_nspline(A, b, rep(1, 60), E, max_iterations = 50L)
solve_nspline <- function(A, b, x0, E_MeV, knots = NULL,
                            max_iterations = 500L, tolerance = 1e-6,
                            dmu_start = 0.1, continuity = "C0C1",
                            sigma_rel = NULL, smoothing = TRUE,
                            n_segments = NULL) {
    maxiter <- if (missing(max_iterations)) .nspline_MAX_ITERATIONS else
        max(as.integer(max_iterations)[1L], 1L)
    tol <- if (missing(tolerance)) .nspline_TOL else
        max(as.numeric(tolerance)[1L], .Machine$double.eps)
    .nspline_solve_full(A = A, b = b, x0 = x0, E_MeV = E_MeV, knots = knots,
                        sigma_rel = sigma_rel, continuity = continuity,
                        max_iterations = maxiter, tol = tol,
                        step_theta = dmu_start, smoothing = smoothing,
                        n_segments = n_segments)
}

# Numerical guards (mirrors of the Python module constants) --------------
.nspline_PHI_FLOOR <- 1e-300
.nspline_LOG_CLIP <- 50.0
.nspline_MAX_ITERATIONS <- 200L
.nspline_TOL <- 1e-3
.nspline_CONTINUITY <- "C0C1"
.nspline_REL_UNCERTAINTY <- 0.1

# Knot presets (MeV) as used in the paper.
NSPLINE_KNOT_PRESETS <- list(
    BARS5_channel = c(1e-10, 1.3e-7, 3.83e-7, 8e-6, 2e-5, 3e-5, 7.3e-5,
                      3.2e-3, 0.38, 0.95, 7.0, 17.0, 20.0),
    IGRIK_channel = c(1e-10, 2e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 1.5e-4,
                      3e-4, 6e-4, 6e-3, 0.27, 1.0, 2.7, 7.0, 13.0, 20.0),
    IGRIK_surface = c(1e-10, 2e-8, 1e-7, 2e-7, 3e-6, 5e-6, 2.5e-4, 0.6,
                      0.8, 1.5, 2.7, 7.0, 11.5, 14.0, 20.0),
    YAGUAR_channel = c(1e-10, 2e-8, 1e-7, 6e-7, 1e-6, 3e-6, 1e-5, 4.3e-5,
                       1.8e-4, 6.3e-4, 5e-3, 0.6, 0.8, 1.0, 2.5, 7.0, 11.0,
                       13.0, 20.0)
)

# ---- internal helpers ------------------------------------------------------

# Minimum-norm least squares, mirroring numpy.linalg.lstsq(rcond=None):
# singular values below eps * max(dim) * sigma_max are discarded.
.nspline_lstsq <- function(Mmat, rhs) {
    sv <- svd(Mmat)
    cutoff <- .Machine$double.eps * max(nrow(Mmat), ncol(Mmat)) * sv$d[1L]
    keep <- which(sv$d > cutoff)
    if (!length(keep)) return(rep(0, ncol(Mmat)))
    vec <- as.numeric(crossprod(sv$u[, keep, drop = FALSE], rhs)) / sv$d[keep]
    as.numeric(sv$v[, keep, drop = FALSE] %*% vec)
}

# Rank-revealing pivoted QR followed by the truncated-SVD solve of the
# triangular factor -- the same operation sequence as LAPACK dgelsd, which
# numpy.linalg.lstsq uses.  On the near-singular KKT system the plain
# pseudo-inverse and the QR-then-SVD paths disagree in the directions whose
# singular values sit next to the rcond cut-off, so the sequence has to be
# mirrored to reproduce the Python numbers.
.nspline_gelsd <- function(Mmat, rhs) {
    qx <- qr(Mmat, pivot = TRUE)
    Rf <- qr.R(qx)
    bq <- as.numeric(qr.qty(qx, rhs))
    sv <- svd(Rf)
    cutoff <- .Machine$double.eps * max(nrow(Mmat), ncol(Mmat)) * sv$d[1L]
    keep <- which(sv$d > cutoff)
    y <- if (length(keep)) {
        vec <- as.numeric(crossprod(sv$u[, keep, drop = FALSE], bq)) / sv$d[keep]
        as.numeric(sv$v[, keep, drop = FALSE] %*% vec)
    } else {
        rep(0, ncol(Rf))
    }
    out <- numeric(ncol(Mmat))
    out[qx$pivot] <- y
    out
}

# Automatic log-uniform knot grid spanning [min(E), max(E)]; `n_segments`
# spline segments => n_segments + 1 knots (Python auto_knots).
.nspline_auto_knots_full <- function(E_MeV, n_segments = 12L) {
    E <- as.numeric(E_MeV)
    Epos <- E[E > 0]
    if (length(Epos) < 2L)
        stop("auto_knots requires at least two positive energy points, got ",
             length(Epos))
    emin <- min(Epos); emax <- max(Epos)
    if (!is.finite(emin) || !is.finite(emax) || emin >= emax)
        stop("auto_knots requires finite min(E) < max(E), got [",
             emin, ", ", emax, "]")
    n_segments <- as.integer(n_segments)[1L]
    if (n_segments < 1L) stop("n_segments must be >= 1, got ", n_segments)
    as.numeric(10^seq(log10(emin), log10(emax), length.out = n_segments + 1L))
}

# Resolve a knot specification to a valid knot vector (Python _resolve_knots).
.nspline_resolve_knots <- function(knots, E_MeV, n_segments = NULL) {
    E <- as.numeric(E_MeV)
    Epos <- E[E > 0]
    emin <- if (length(Epos)) min(Epos) else 1e-10
    emax <- max(E)

    if (is.null(knots)) {
        if (is.null(n_segments)) {
            ns <- max(length(E), 2L)
            n_segments <- min(12L, max(4L, ns %/% 4L))
        }
        return(list(knots = .nspline_auto_knots_full(E, n_segments),
                    src = "auto"))
    }

    if (is.character(knots)) {
        key <- trimws(knots[1L])
        if (!key %in% names(NSPLINE_KNOT_PRESETS))
            stop("Unknown N-spline knot preset '", key, "'. Available: ",
                 paste(sort(names(NSPLINE_KNOT_PRESETS)), collapse = ", "))
        src <- paste0("preset:", key)
        kn <- as.numeric(NSPLINE_KNOT_PRESETS[[key]])
    } else {
        src <- "user"
        kn <- as.numeric(knots)
    }
    if (length(kn) < 2L) stop("N-spline needs at least 2 knots, got ", length(kn))
    if (any(diff(kn) <= 0)) stop("N-spline knots must be strictly increasing")

    # Clip to the grid range and extend the outer knots so the spline domain
    # spans the whole energy grid.
    kn <- sort(unique(pmin(pmax(kn, emin), emax)))
    if (kn[1L] > emin) kn[1L] <- emin
    if (kn[length(kn)] < emax) kn[length(kn)] <- emax
    if (length(kn) < 2L) kn <- c(emin, emax)
    list(knots = as.numeric(kn), src = src)
}

# Segment index (1-based, 1..M) for each energy point; mirrors
# np.searchsorted(knots, E, side = "right") - 1 clipped to [0, M - 1].
.nspline_segment_indices <- function(E, knots) {
    M <- length(knots) - 1L
    idx <- vapply(as.numeric(E), function(v) sum(knots <= v), numeric(1))
    as.integer(pmin(pmax(idx - 1, 0), M - 1L)) + 1L
}

# Continuity matrix D of the paper (rows = 0, M-1 or 2(M-1), cols = 3M).
.nspline_continuity_matrix <- function(knots, continuity = "C0C1") {
    kn <- as.numeric(knots)
    M <- length(kn) - 1L
    if (M < 1L) stop("knots must contain at least 2 values")
    cont <- toupper(gsub("[[:space:]]", "", as.character(continuity)[1L]))
    if (!cont %in% c("C0C1", "C0", "NONE"))
        stop("continuity must be one of 'C0C1', 'C0', 'none', got '",
             continuity, "'")
    n_int <- M - 1L
    if (identical(cont, "NONE") || n_int == 0L)
        return(matrix(0, 0L, 3L * M))

    rows_c1 <- identical(cont, "C0C1")
    n_rows <- if (rows_c1) 2L * n_int else n_int
    D <- matrix(0, n_rows, 3L * M)
    for (k in seq_len(n_int)) {          # k = 1..n_int <-> Python k = 0..n_int-1
        Ek <- kn[k + 1L]
        uk <- log(Ek)
        D[k, k] <- -1.0
        D[k, k + 1L] <- 1.0
        D[k, M + k] <- -uk
        D[k, M + k + 1L] <- uk
        D[k, 2L * M + k] <- -Ek
        D[k, 2L * M + k + 1L] <- Ek
        if (rows_c1) {
            row <- n_int + k
            D[row, M + k] <- -1.0
            D[row, M + k + 1L] <- 1.0
            D[row, 2L * M + k] <- -Ek
            D[row, 2L * M + k + 1L] <- Ek
        }
    }
    D
}

# Evaluate the N-spline exp(a_k + q_k ln E + r_k E).
.nspline_eval <- function(E, a, q, r, knots) {
    E <- as.numeric(E)
    if (any(E <= 0)) stop("nspline_eval requires strictly positive energies")
    seg <- .nspline_segment_indices(E, knots)
    as.numeric(exp(as.numeric(a)[seg] + as.numeric(q)[seg] * log(E) +
                   as.numeric(r)[seg] * E))
}

# Directed divergence H = sum(pN log(pN/p) - pN + p) >= 0.
.nspline_divergence <- function(p_calc, p_meas) {
    pN <- pmax(as.numeric(p_calc), .nspline_PHI_FLOOR)
    p <- pmax(as.numeric(p_meas), .nspline_PHI_FLOOR)
    sum(pN * log(pN / p) - pN + p)
}

# Pointwise N-spline approximation (Eqs. 6-7): constrained weighted
# log-domain least squares solved through the KKT system.
.nspline_fit <- function(E, phi, knots = NULL, rel_err = NULL,
                         continuity = "C0C1", n_segments = NULL) {
    E <- as.numeric(E)
    phi <- as.numeric(phi)
    if (length(E) != length(phi))
        stop("E and phi length mismatch: ", length(E), " vs ", length(phi))
    if (any(E <= 0)) stop("fit_nspline requires strictly positive energies")
    if (length(phi) < 3L) stop("fit_nspline requires at least 3 spectrum points")

    res <- .nspline_resolve_knots(knots, E, n_segments)
    kn <- res$knots
    src <- res$src
    M <- length(kn) - 1L
    n <- length(E)

    seg <- .nspline_segment_indices(E, kn)
    u <- log(E)
    rows <- seq_len(n)
    G <- matrix(0, n, 3L * M)
    G[cbind(rows, seg)] <- 1.0
    G[cbind(rows, M + seg)] <- u
    G[cbind(rows, 2L * M + seg)] <- E

    # Floor tiny/zero bins and down-weight them so they do not drag the fit.
    phi_max <- max(phi)
    tiny <- max(.nspline_PHI_FLOOR, 1e-12 * phi_max)
    floored <- phi < tiny
    y <- log(ifelse(floored, tiny, phi))

    w <- if (is.null(rel_err)) rep(1, n) else
        1 / pmax(as.numeric(rel_err), 1e-12)
    w <- ifelse(floored, 1e-3 * w, w)

    D <- .nspline_continuity_matrix(kn, continuity)
    nc <- nrow(D)

    Gw <- G * w                       # row-wise scaling (column-major recycle)
    yw <- y * w
    H_norm <- crossprod(Gw, Gw)
    KKT <- matrix(0, 3L * M + nc, 3L * M + nc)
    KKT[seq_len(3L * M), seq_len(3L * M)] <- H_norm
    rhs <- c(as.numeric(crossprod(Gw, yw)), rep(0, nc))
    if (nc) {
        idx_p <- seq_len(3L * M)
        idx_c <- 3L * M + seq_len(nc)
        KKT[idx_p, idx_c] <- t(D)
        KKT[idx_c, idx_p] <- D
    }

    X <- .nspline_lstsq(KKT, rhs)[seq_len(3L * M)]
    N_E <- as.numeric(exp(G %*% X))
    resid <- w * (as.numeric(G %*% X) - y)
    rms <- sqrt(mean(resid^2)) / max(mean(w), .nspline_PHI_FLOOR)

    list(N_E = N_E,
         knots = kn, knots_source = src, continuity = continuity,
         a = X[seq_len(M)],
         q = X[M + seq_len(M)],
         r = 2L * M + seq_len(M),
         log_rms_residual = rms)
}

# np.trapezoid(y, x)
.nspline_trapezoid <- function(y, x) {
    y <- as.numeric(y); x <- as.numeric(x)
    if (length(y) < 2L) return(0)
    dx <- diff(x)
    sum(dx * (y[-length(y)] + y[-1L]) / 2)
}

# ---- full solver with diagnostics (Python solve_nspline_full) --------------
.nspline_solve_full <- function(A, b, x0 = NULL, E_MeV = NULL, knots = NULL,
                                sigma_rel = NULL, continuity = "C0C1",
                                max_iterations = 200L, tol = 1e-3,
                                step_theta = 0.1, smoothing = TRUE,
                                n_segments = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A); n <- ncol(A)
    if (length(b) != m)
        stop("b length (", length(b), ") does not match A rows (", m, ")")
    if (m < 1L) stop("At least one measurement is required")
    if (is.null(E_MeV)) stop("E_MeV (energy grid in MeV) is required")
    E <- as.numeric(E_MeV)
    if (length(E) != n)
        stop("E_MeV length (", length(E), ") does not match A columns (", n, ")")
    if (any(E <= 0)) stop("E_MeV must contain strictly positive energies")
    if (max_iterations < 1)
        stop("max_iterations must be >= 1, got ", max_iterations)
    if (!(step_theta > 0 && step_theta <= 1))
        stop("step_theta must be in (0, 1], got ", step_theta)
    if (!(tol > 0)) stop("tol must be positive, got ", tol)

    valid <- b > 0
    if (!any(valid)) stop("solve_nspline requires at least one positive measurement")
    A_v <- A[valid, , drop = FALSE]
    b_v <- b[valid]

    sigma_v <- if (is.null(sigma_rel)) rep(.nspline_REL_UNCERTAINTY, length(b_v))
        else pmax(as.numeric(sigma_rel)[which(valid)], 1e-12)

    kres <- .nspline_resolve_knots(knots, E, n_segments)
    kn <- kres$knots; knot_src <- kres$src

    p <- b_v / sum(b_v)
    H_target <- 0.5 * sum(p * sigma_v^2)

    if (is.null(x0)) {
        x <- rep(1, n)
    } else {
        x <- as.numeric(x0)
        if (length(x) != n)
            stop("x0 length (", length(x), ") does not match A columns (", n, ")")
        x[!is.finite(x)] <- 0
        x <- pmax(x, 0)
    }
    if (sum(x) <= 0) x <- rep(1, n)

    Qc0 <- as.numeric(A_v %*% x)
    scale <- sum(b_v) / max(sum(Qc0), .nspline_PHI_FLOOR)
    x <- pmax(x * scale, .nspline_PHI_FLOOR)

    # Pointwise relative errors for the per-iteration smoothing fits.
    sens <- colSums(A_v)
    sens_max <- if (length(sens)) max(sens) else 0
    smooth_rel_err <- if (isTRUE(smoothing) && sens_max > 0)
        sqrt(pmin(pmax(sens_max / pmax(sens, .nspline_PHI_FLOOR), 1), 1e12))
    else NULL

    fit_info <- NULL
    if (isTRUE(smoothing)) {
        fl <- .nspline_fit(E, x, knots = kn, rel_err = smooth_rel_err,
                           continuity = continuity)
        x <- pmax(fl$N_E, .nspline_PHI_FLOOR)
        x <- x * (sum(b_v) / max(sum(as.numeric(A_v %*% x)), .nspline_PHI_FLOOR))
        fit_info <- fl
    }

    b_total <- sum(b_v)
    eps_scale <- 1e-12 * max(b_total, .nspline_PHI_FLOOR)

    # Pin the activation scale: sum(A x) = sum(b).
    gauge <- function(xx) xx * (b_total / max(sum(as.numeric(A_v %*% xx)),
                                              .nspline_PHI_FLOOR))
    state <- function(xx) {
        xx <- pmax(xx, .nspline_PHI_FLOOR)
        Qc <- pmax(as.numeric(A_v %*% xx), eps_scale)
        pN <- Qc / max(sum(Qc), .nspline_PHI_FLOOR)
        list(x = xx, Qc = Qc, pN = pN, H = .nspline_divergence(pN, p))
    }

    st <- state(gauge(x))
    x <- st$x; Qc <- st$Qc; pN <- st$pN; H <- st$H
    H_history <- H
    converged <- FALSE
    stop_reason <- "max_iterations"
    iterations <- 0L

    if (H <= H_target) {
        converged <- TRUE
        stop_reason <- "H_target (initial)"
    }

    for (iteration in seq_len(max_iterations)) {
        iterations <- iteration

        ln_ratio <- pmax(pmin(log(pN / p), .nspline_LOG_CLIP),
                         -.nspline_LOG_CLIP)
        R <- as.numeric(crossprod(A_v, ln_ratio)) / b_total
        Rbar <- sum(x * R) / max(sum(x), .nspline_PHI_FLOOR)
        g <- R - Rbar
        g_max <- max(abs(g))
        if (!is.finite(g_max) || g_max <= 0) {
            stop_reason <- "stalled_gradient"
            iterations <- iterations - 1L
            break
        }

        mu <- step_theta / g_max
        accepted <- FALSE
        st_new <- list(x = x, Qc = Qc, pN = pN, H = H)
        for (bt in seq_len(60L)) {
            x_trial <- x * (1 - mu * g)
            if (isTRUE(smoothing))
                x_trial <- .nspline_fit(E, x_trial, knots = kn,
                                        rel_err = smooth_rel_err,
                                        continuity = continuity)$N_E
            st_new <- state(gauge(x_trial))
            if (is.finite(st_new$H) &&
                st_new$H <= H + 1e-4 * max(H, .nspline_PHI_FLOOR)) {
                accepted <- TRUE
                break
            }
            mu <- mu * 0.5
        }
        if (!accepted) {
            stop_reason <- "no_further_reduction"
            iterations <- iterations - 1L
            break
        }

        H_prev <- H
        x <- st_new$x; Qc <- st_new$Qc; pN <- st_new$pN; H <- st_new$H
        H_history <- c(H_history, H)

        if (H <= H_target) {
            converged <- TRUE
            stop_reason <- "H_target"
            break
        }
        if (abs(H_prev - H) <= tol * max(H_prev, .nspline_PHI_FLOOR)) {
            converged <- TRUE
            stop_reason <- "relative_change"
            break
        }
    }

    # Acceptability statistic: nev = RMS((Qr - Q)/dQ), nev <= 1 + 2/sqrt(N).
    Qr_full <- as.numeric(A %*% x)
    rel_res <- numeric(m)
    denom <- pmax(sigma_v * b_v, .nspline_PHI_FLOOR)
    rel_res[valid] <- (Qr_full[valid] - b_v) / denom
    cnt <- sum(valid)
    div <- if (cnt > 1) cnt - 1 else cnt
    nev <- sqrt(sum(rel_res[valid]^2) / max(div, 1))
    nev_limit <- 1 + 2 / sqrt(cnt)
    acceptable <- nev <= nev_limit

    fluence <- .nspline_trapezoid(x, E)
    mean_energy <- if (fluence > 0)
        .nspline_trapezoid(E * x, E) / fluence else NaN

    if (is.null(fit_info)) {
        M <- length(kn) - 1L
        seg <- .nspline_segment_indices(E, kn)
        rows <- seq_len(n)
        G <- matrix(0, n, 3L * M)
        G[cbind(rows, seg)] <- 1.0
        G[cbind(rows, M + seg)] <- log(E)
        G[cbind(rows, 2L * M + seg)] <- E
        Xl <- .nspline_lstsq(G, log(pmax(x, .nspline_PHI_FLOOR)))
        params <- list(a = Xl[seq_len(M)], q = Xl[M + seq_len(M)],
                       r = 2L * M + seq_len(M))
    } else {
        params <- list(a = fit_info$a, q = fit_info$q, r = fit_info$r)
    }

    list(spectrum = as.numeric(x), iterations = as.integer(iterations),
         converged = converged, stop_reason = stop_reason,
         H = H, H_history = H_history, H_target = H_target,
         nev = nev, nev_limit = nev_limit, acceptable = acceptable,
         Qr = Qr_full, relative_residuals = rel_res,
         fluence = fluence, mean_energy = mean_energy,
         knots = kn, knots_source = knot_src,
         continuity = continuity, params = params)
}

#' Auto-generate log-uniform interior knots
#'
#' @param E_MeV Numeric energy grid.
#' @param n_interior Integer; number of interior knots. Default
#'   \code{NULL} = \code{max(4, length(E_MeV)/10)}.
#' @return Numeric vector of strictly increasing interior knot energies.
#' @export
#' @examples
#' E <- 10^seq(-9, 1, length.out = 60)
#' auto_knots(E)
auto_knots <- function(E_MeV, n_interior = NULL) {
    E_MeV <- as.numeric(E_MeV)
    e_min <- min(E_MeV); e_max <- max(E_MeV)
    if (is.null(n_interior)) n_interior <- max(4L, length(E_MeV) %/% 10L)
    10^seq(log10(e_min), log10(e_max), length.out = n_interior + 2L)[-c(1L, n_interior + 2L)]
}

#' Build continuity-constraint matrix for N-spline basis
#'
#' Returns the block matrix \eqn{D} of continuity constraints
#' (\eqn{C^0/C^1}) on the N-spline parameters \eqn{X = (a, q, r)^T}.
#'
#' @param n_knots Integer; number of interior knots (segments = n_knots + 1).
#' @param order Integer; 0 (C0) or 1 (C1). Default 0.
#' @return Numeric matrix.
#' @export
#' @examples
#' D <- build_continuity_matrix(4L, order = 0L)
#' dim(D)
build_continuity_matrix <- function(n_knots, order = 0L) {
    segments <- n_knots + 1L
    M <- 3L * segments
    if (order == 0L) {
        n_constraints <- segments - 1L
        D <- matrix(0.0, n_constraints, M)
        for (i in seq_len(n_constraints)) {
            D[i, (i - 1L) * 3L + 1L] <- 1
            D[i, i * 3L + 1L] <- -1
        }
        D
    } else if (order == 1L) {
        n_constraints <- 2L * (segments - 1L)
        D <- matrix(0.0, n_constraints, M)
        for (i in seq_len(segments - 1L)) {
            D[i, (i - 1L) * 3L + 1L] <- 1
            D[i, i * 3L + 1L] <- -1
            # C1: dN/dln E = q + r E. At knot i, this gives constraints on
            # q_left + r_left * E_i = q_right + r_right * E_i.
            D[segments - 1L + i, (i - 1L) * 3L + 2L] <- 1
            D[segments - 1L + i, (i - 1L) * 3L + 3L] <- 1
            D[segments - 1L + i, i * 3L + 2L] <- -1
            D[segments - 1L + i, i * 3L + 3L] <- -1
        }
        D
    } else {
        stop("order must be 0 or 1")
    }
}

#' Wrapper around \code{\link{solve_nspline}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_nspline
#' @param relative_uncertainty Numeric; relative measurement uncertainty
#'   \code{dQ/Q} used by the stopping criteria and the \code{nev} statistic.
#'   Default 0.1.
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_nspline <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              knots = NULL,
                              max_iterations = 500L, tolerance = 1e-6,
                              dmu_start = 0.1,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL,
                              continuity = .nspline_CONTINUITY,
                              relative_uncertainty = .nspline_REL_UNCERTAINTY,
                              smoothing = TRUE, n_segments = NULL) {
    maxiter <- if (missing(max_iterations)) .nspline_MAX_ITERATIONS else
        max(as.integer(max_iterations)[1L], 1L)
    tol <- if (missing(tolerance)) .nspline_TOL else
        max(as.numeric(tolerance)[1L], .Machine$double.eps)
    x0_default <- rep(1.0, n_energy_bins) / 2.0
    solver_with_E <- function(A, b, x0 = NULL, ...) {
        .nspline_solve_full(A = A, b = b,
                            x0 = if (is.null(x0)) x0_default else x0,
                            E_MeV = E_MeV, knots = knots,
                            sigma_rel = rep(as.numeric(relative_uncertainty)[1L],
                                            length(b)),
                            continuity = continuity,
                            max_iterations = maxiter, tol = tol,
                            step_theta = dmu_start, smoothing = smoothing,
                            n_segments = n_segments)
    }
    probe <- tryCatch({
        used <- detector_names[detector_names %in% names(readings)]
        .nspline_solve_full(
            A = do.call(rbind, lapply(used,
                                      function(n) as.numeric(sensitivities[[n]]))),
            b = as.numeric(readings[used]),
            x0 = x0_default, E_MeV = E_MeV, knots = knots,
            sigma_rel = rep(as.numeric(relative_uncertainty)[1L], length(used)),
            continuity = continuity, max_iterations = maxiter, tol = tol,
            step_theta = dmu_start, smoothing = smoothing,
            n_segments = n_segments)
    }, error = function(e) NULL)
    extra_output <- if (is.null(probe)) NULL else list(
        continuity = continuity,
        relative_uncertainty = as.numeric(relative_uncertainty)[1L],
        H = probe$H, H_history = probe$H_history, H_target = probe$H_target,
        nev = probe$nev, nev_limit = probe$nev_limit,
        acceptable = probe$acceptable, stop_reason = probe$stop_reason,
        fluence = probe$fluence, mean_energy = probe$mean_energy,
        knots = probe$knots, knots_source = probe$knots_source,
        Qr = probe$Qr, relative_residuals = probe$relative_residuals)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = solver_with_E,
        solve_kwargs = list(),
        method_name = "NSpline",
        extra_output = extra_output,
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
