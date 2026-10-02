#' Non-negative K-SVD unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_nnksvd.py}.
#' Two-step pipeline following Xu et al. (NIMA 2026): (1) non-negative K-SVD
#' dictionary learning with non-negative truncation of the rank-1 SVD atom
#' update, (2) sparse inversion on the equivalent detection dictionary
#' \code{M_norm = normalize(R \%*\% D)} using NNLS-top-K, OMP or NN-OMP.
#' The reconstruction objective is the Tikhonov-regularised non-negative
#' least-squares problem (Eq. 2.5/2.6 of the article), solved through the
#' augmented-matrix form so that a plain NNLS solver can be used, plus the
#' optional training-sample-driven prior constraint
#' \code{prior_wt * ||alpha - alpha_prior||^2}.
#'
#' @section Metrics:
#' The three evaluation metrics of the article (Section 2.2.3) live in
#' \code{relative_flux_error}, \code{pearson_correlation} and
#' \code{comprehensive_score}.
#'
#' @name nnksvd
NULL

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# L2-normalise every column of D (returns a copy), mirroring
# bssunfold.core.unfold_nnksvd._normalize_dictionary.
.normalize_dictionary <- function(D) {
    D <- as.matrix(D); storage.mode(D) <- "double"
    norms <- sqrt(colSums(D^2))
    norms <- ifelse(norms == 0, 1.0, norms)
    sweep(D, 2L, norms, "/")
}

# scipy.optimize.nnls equivalent: min ||A x - b|| s.t. x >= 0.
.nnksvd_nnls <- function(A, b, max_iter = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    if (ncol(A) == 0L) return(numeric(0L))
    as.numeric(.bss_nnls(A, b, maxiter = max_iter))
}

# Equivalent (column-normalised) detection dictionary M = R %*% D.
.build_equivalent_dictionary <- function(R, D, normalize = TRUE) {
    M <- as.matrix(R) %*% as.matrix(D)
    if (isTRUE(normalize)) M <- .normalize_dictionary(M)
    M
}

# ---------------------------------------------------------------------------
# Tikhonov-regularised NNLS via the augmented form (Eq. 2.5 / 2.6)
# ---------------------------------------------------------------------------

#' Tikhonov-regularised non-negative least squares
#'
#' Solves \eqn{\min \|y - M \alpha\|^2 + \lambda_{tik}\|\alpha\|^2} s.t.
#' \eqn{\alpha \ge 0} through the augmented-matrix equivalent
#' \deqn{\min \|[y;0] - [M; \sqrt{\lambda_{tik}} I] \alpha\|^2.}
#' When \code{prior_wt > 0} the training-sample-driven prior rows
#' \eqn{\sqrt{w} I} with right-hand side \eqn{\sqrt{w} \alpha_{prior}} are
#' appended as well.
#'
#' @param M_norm Numeric matrix \code{(m x p)}.
#' @param y Numeric measurement vector (length \code{m}).
#' @param lambda_tik Numeric; Tikhonov regularisation weight. Default \code{0.01}.
#' @param prior_wt Numeric; training-sample-driven prior weight. Default \code{0}.
#' @param alpha_prior Optional prior coefficient vector (length \code{p});
#'   required when \code{prior_wt > 0}.
#' @param max_iter Integer or \code{NULL}; maximum NNLS iterations.
#' @return Numeric vector (length \code{p}) of non-negative coefficients.
#' @export
#' @examples
#' set.seed(1)
#' A <- matrix(runif(3 * 5), nrow = 3)
#' b <- c(1, 2, 3)
#' x <- solve_tikhonov_nnls(A, b, lambda_tik = 0.01)
solve_tikhonov_nnls <- function(M_norm, y, lambda_tik = 0.01, prior_wt = 0.0,
                                alpha_prior = NULL, max_iter = NULL) {
    M_norm <- as.matrix(M_norm); storage.mode(M_norm) <- "double"
    y <- as.numeric(y)
    p <- ncol(M_norm)

    A_aug <- rbind(M_norm, sqrt(max(lambda_tik, 0.0)) * diag(p))
    b_aug <- c(y, numeric(p))

    if (prior_wt > 0.0) {
        if (is.null(alpha_prior))
            stop("alpha_prior must be provided when prior_wt > 0")
        alpha_prior <- as.numeric(alpha_prior)
        if (length(alpha_prior) != p) {
            stop("alpha_prior length (", length(alpha_prior),
                 ") must match number of dictionary atoms (", p, ")")
        }
        A_aug <- rbind(A_aug, sqrt(prior_wt) * diag(p))
        b_aug <- c(b_aug, sqrt(prior_wt) * alpha_prior)
    }
    .nnksvd_nnls(A_aug, b_aug, max_iter = max_iter)
}

# ---------------------------------------------------------------------------
# Non-negative OMP (NN-OMP) sparse coding
# ---------------------------------------------------------------------------

#' Non-negative Orthogonal Matching Pursuit
#'
#' Greedy sparse coding with a non-negativity constraint: at every step the
#' atom with the largest \emph{positive} normalised projection onto the
#' residual is appended to the support and the support coefficients are
#' re-solved as an NNLS problem.
#'
#' @param D Numeric dictionary matrix \code{(n x p)}.
#' @param y Numeric signal vector (length \code{n}).
#' @param sparsity Integer; maximum number of non-zero coefficients (K).
#' @param tolerance Numeric; early-stopping residual tolerance. Default \code{1e-6}.
#' @return Numeric vector (length \code{p}) of non-negative coefficients.
#' @export
#' @examples
#' set.seed(1)
#' A <- matrix(runif(10 * 5), nrow = 10)
#' b <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' x <- solve_nn_omp(A, b, sparsity = 3L)
solve_nn_omp <- function(D, y, sparsity, tolerance = 1e-6) {
    D <- as.matrix(D); storage.mode(D) <- "double"
    y <- as.numeric(y)
    p <- ncol(D)
    alpha <- numeric(p)
    residual <- y
    support <- integer(0L)

    # Normalise the dictionary columns for atom selection (NNLS unaffected).
    norms <- sqrt(colSums(D^2))
    norms <- ifelse(norms == 0, 1.0, norms)
    D_norm <- sweep(D, 2L, norms, "/")

    for (step in seq_len(min(as.integer(sparsity), p))) {
        correlations <- as.numeric(crossprod(D_norm, residual))
        if (length(support)) correlations[support] <- -Inf
        idx <- which.max(correlations)
        if (!is.finite(correlations[idx]) || correlations[idx] <= 0) break
        support <- c(support, idx)

        D_s <- D[, support, drop = FALSE]
        coefs <- .nnksvd_nnls(D_s, y)
        alpha[support] <- coefs
        residual <- y - as.numeric(D_s %*% coefs)

        if (sqrt(sum(residual^2)) < tolerance) break
    }
    alpha
}

# ---------------------------------------------------------------------------
# NNLS + TopK sparse coding (the article's proposed strategy)
# ---------------------------------------------------------------------------

#' NNLS with top-K atom screening
#'
#' Three-stage hierarchical sparse coding: (1) a global Tikhonov-regularised
#' NNLS coarse solution on the full dictionary, (2) screening of the \code{K}
#' atoms with the largest coarse coefficients, (3) a local Tikhonov-regularised
#' NNLS re-optimisation restricted to the screened support.
#'
#' @param M_norm Numeric matrix \code{(m x p)}.
#' @param y Numeric measurement vector (length \code{m}).
#' @param sparsity Integer; target sparsity \code{K}.
#' @param lambda_tik Numeric; Tikhonov weight. Default \code{0.01}.
#' @param prior_wt Numeric; prior weight. Default \code{0}.
#' @param alpha_prior Optional prior coefficient vector (length \code{p}).
#' @param max_iter Integer or \code{NULL}; maximum NNLS iterations.
#' @return Numeric K-sparse non-negative vector (length \code{p}).
#' @export
#' @examples
#' set.seed(1)
#' A <- matrix(runif(5 * 10), nrow = 5)
#' b <- c(1, 2, 3, 4, 5)
#' x <- solve_nnls_topk(A, b, sparsity = 3L)
solve_nnls_topk <- function(M_norm, y, sparsity, lambda_tik = 0.01,
                            prior_wt = 0.0, alpha_prior = NULL,
                            max_iter = NULL) {
    M_norm <- as.matrix(M_norm); storage.mode(M_norm) <- "double"
    y <- as.numeric(y)
    p <- ncol(M_norm)
    K <- max(0L, min(as.integer(sparsity), p))

    if (K == 0L) return(numeric(p))

    # Step 1: global NNLS coarse solution.
    alpha_full <- solve_tikhonov_nnls(M_norm, y, lambda_tik = lambda_tik,
                                      prior_wt = prior_wt,
                                      alpha_prior = alpha_prior,
                                      max_iter = max_iter)

    # Step 2: top-K atom screening (numpy argsort ascending, last K).
    if (K >= p) return(alpha_full)
    ord <- order(alpha_full)
    topk_idx <- ord[(p - K + 1L):p]

    # Step 3: local NNLS fine optimisation on the screened support.
    M_topk <- M_norm[, topk_idx, drop = FALSE]
    alpha_prior_topk <- if (!is.null(alpha_prior) && prior_wt > 0.0)
                            as.numeric(alpha_prior)[topk_idx] else NULL
    alpha_topk <- solve_tikhonov_nnls(M_topk, y, lambda_tik = lambda_tik,
                                      prior_wt = prior_wt,
                                      alpha_prior = alpha_prior_topk,
                                      max_iter = max_iter)

    alpha <- numeric(p)
    alpha[topk_idx] <- alpha_topk
    alpha
}

# ---------------------------------------------------------------------------
# Non-negative K-SVD dictionary learning
# ---------------------------------------------------------------------------

#' Non-negative K-SVD dictionary learning
#'
#' K-SVD variant with non-negativity constraints on both the dictionary atoms
#' and the sparse coefficients: after the rank-1 SVD update of an atom the
#' atom and its coefficients are non-negatively truncated and the atom is
#' re-normalised, with the norm folded into the coefficients so that
#' \code{D \%*\% alpha} stays invariant.
#'
#' @param signals Numeric training-signal matrix \code{(n x m)}, one column per
#'   training sample; the columns are clipped at zero first.
#' @param n_atoms Integer; number of dictionary atoms \code{P}.
#' @param n_iterations Integer; maximum K-SVD iterations. Default \code{80L}.
#' @param sparsity Integer; target sparsity \code{K} for the sparse coders.
#'   Default \code{2L}.
#' @param lambda_tik Numeric; Tikhonov weight for \code{nnls_topk}. Default \code{0.01}.
#' @param prior_wt Numeric; prior weight (unused while training, kept for API
#'   parity with the Python original). Default \code{0.5}.
#' @param sparse_coder Character; \code{"nnls_topk"}, \code{"omp"} or
#'   \code{"nn_omp"}. Default \code{"nnls_topk"}.
#' @param random_state Integer or \code{NULL}; seed for the random atom
#'   initialisation.
#' @param tolerance Numeric; early-stopping tolerance on the dictionary change.
#'   Default \code{1e-6}.
#' @return A list \code{list(D, alpha_prior)} where \code{D} is the learned
#'   \code{(n x p)} column-normalised non-negative dictionary and
#'   \code{alpha_prior} is the mean sparse code of the training signals.
#' @export
#' @examples
#' set.seed(7)
#' S <- matrix(runif(12 * 6), nrow = 12)
#' fit <- solve_nnksvd(S, n_atoms = 3L, n_iterations = 2L)
#' dim(fit$D)
solve_nnksvd <- function(signals, n_atoms, n_iterations = 80L, sparsity = 2L,
                         lambda_tik = 0.01, prior_wt = 0.5,
                         sparse_coder = "nnls_topk", random_state = NULL,
                         tolerance = 1e-6) {
    if (!sparse_coder %in% c("nnls_topk", "omp", "nn_omp")) {
        stop("Unknown sparse_coder '", sparse_coder,
             "'. Expected 'nnls_topk', 'omp' or 'nn_omp'.", call. = FALSE)
    }

    signals <- as.matrix(signals); storage.mode(signals) <- "double"
    signals <- pmax(signals, 0.0)
    n <- nrow(signals); m <- ncol(signals)

    if (!is.null(random_state)) set.seed(as.integer(random_state))
    n_atoms <- max(1L, min(as.integer(n_atoms), m))
    # Initialise the dictionary with random training samples (non-negative).
    idx <- sample.int(m, size = n_atoms, replace = FALSE)
    D <- signals[, idx, drop = FALSE]
    # Replace an all-zero chosen sample with a smooth bump.
    zero_cols <- which(sqrt(colSums(D^2)) == 0)
    for (zc in zero_cols) {
        bump <- pmax(stats::rnorm(n), 0.0)
        if (sqrt(sum(bump^2)) == 0) bump <- rep(1, n)
        D[, zc] <- bump
    }
    D <- .normalize_dictionary(D)

    coefficients <- matrix(0.0, nrow = n_atoms, ncol = m)

    for (it in seq_len(as.integer(n_iterations))) {
        D_prev <- D

        # ---- Sparse coding stage (per training sample) ----
        for (j in seq_len(m)) {
            y_j <- signals[, j]
            coefficients[, j] <- if (sparse_coder == "omp") {
                solve_omp(D, y_j, sparsity = sparsity, tolerance = tolerance)
            } else if (sparse_coder == "nn_omp") {
                solve_nn_omp(D, y_j, sparsity = sparsity, tolerance = tolerance)
            } else {
                # nnls_topk operates on the dictionary directly: there is no
                # response matrix here, so M_norm = normalise(D) = D, and the
                # prior is not yet available during training.
                solve_nnls_topk(D, y_j, sparsity = sparsity,
                                lambda_tik = lambda_tik, prior_wt = 0.0,
                                alpha_prior = NULL)
            }
        }

        # ---- Dictionary update stage ----
        for (atom in seq_len(n_atoms)) {
            used <- which(coefficients[atom, ] != 0)
            if (length(used) == 0L) {
                # Re-initialise the unused atom from a random training sample.
                j_new <- sample.int(m, 1L)
                new_atom <- pmax(signals[, j_new], 0.0)
                nrm <- sqrt(sum(new_atom^2))
                if (nrm == 0) { new_atom <- rep(1, n); nrm <- sqrt(n) }
                D[, atom] <- new_atom / nrm
                next
            }

            # Error matrix after removing the current atom's contribution.
            D_restricted <- D
            D_restricted[, atom] <- 0
            E <- signals[, used, drop = FALSE] -
                 D_restricted %*% coefficients[, used, drop = FALSE]

            # Rank-1 SVD approximation of E (classic K-SVD step).
            svdE <- svd(E, nu = 1L, nv = 1L)
            new_atom <- svdE$u[, 1]
            new_coef <- svdE$d[1] * svdE$v[, 1]

            # Non-negative truncation (the article's key modification).
            new_atom <- pmax(new_atom, 0.0)
            new_coef <- pmax(new_coef, 0.0)

            nrm <- sqrt(sum(new_atom^2))
            if (nrm > 0) {
                # Fold the norm into the coefficients so D %*% alpha is kept.
                D[, atom] <- new_atom / nrm
                coefficients[atom, used] <- new_coef * nrm
            } else {
                # The atom collapsed to zero: re-initialise from a sample.
                j_new <- sample.int(m, 1L)
                new_atom <- pmax(signals[, j_new], 0.0)
                nrm <- sqrt(sum(new_atom^2))
                if (nrm == 0) { new_atom <- rep(1, n); nrm <- sqrt(n) }
                D[, atom] <- new_atom / nrm
                # Recompute this atom's coefficient per signal via NNLS.
                for (u_idx in used) {
                    coefficients[atom, u_idx] <-
                        .nnksvd_nnls(matrix(D[, atom], ncol = 1L),
                                     signals[, u_idx])[1L]
                }
            }
        }

        # ---- Convergence check on the dictionary change ----
        if (sqrt(sum((D - D_prev)^2)) <
                tolerance * max(1.0, sqrt(sum(D_prev^2)))) break
    }

    # Final non-negative safeguard and re-normalisation.
    D <- .normalize_dictionary(pmax(D, 0.0))
    # Training-sample-driven prior: mean sparse code across training signals.
    list(D = D, alpha_prior = as.numeric(rowMeans(coefficients)))
}

# ---------------------------------------------------------------------------
# Top-level NN-KSVD unfolding solver
# ---------------------------------------------------------------------------

#' Non-negative K-SVD neutron-spectrum unfolding (low level)
#'
#' R port of \code{bssunfold.core.unfold_nnksvd.solve_nnksvd_unfold}. Two
#' operating modes: a pre-learned \code{dictionary} (used as-is, training
#' prior disabled) or online K-SVD on \code{training_signals}, synthesising
#' log-spaced Gaussian bumps on the \code{E_MeV} grid (plus the normalised
#' initial guess) when no signals are supplied. Sparse coding then runs on the
#' equivalent detection dictionary \code{M_norm = normalize(A \%*\% D)} and the
#' spectrum is reconstructed as \code{D \%*\% (alpha * ||A D_[k]||)}, i.e. the
#' column-normalisation of \code{M} is undone so that \code{A \%*\% phi ~ b},
#' followed by the optional scale alignment used by \code{solve_cs}.
#'
#' @param A Numeric response matrix \code{(m x n)}.
#' @param b Numeric measurement vector (length \code{m}).
#' @param x0 Optional initial spectrum (length \code{n}); seeds the default
#'   training signals.
#' @param n_atoms Integer; number of dictionary atoms. Default \code{15L}.
#' @param sparsity Integer; target sparsity K. Default \code{2L}.
#' @param dictionary Optional pre-learned non-negative dictionary
#'   \code{(n x p)}; bypasses online K-SVD. Default \code{NULL}.
#' @param training_signals Optional training signals \code{(n x m)}.
#'   Default \code{NULL}.
#' @param n_dictionary_iterations Integer; K-SVD iterations. Default \code{80L}.
#' @param lambda_tik Numeric; Tikhonov weight. Default \code{0.01}.
#' @param prior_wt Numeric; training-sample prior weight. Default \code{0.5}.
#' @param sparse_coder Character; \code{"nnls_topk"}, \code{"omp"} or
#'   \code{"nn_omp"}. Default \code{"nnls_topk"}.
#' @param random_state Integer or \code{NULL}; K-SVD seed (\code{NULL} = 42).
#' @param tolerance Numeric; convergence tolerance. Default \code{1e-6}.
#' @param n_nnls_iter Integer or \code{NULL}; maximum NNLS iterations.
#' @param E_MeV Numeric energy grid (length \code{n}); used to place the
#'   default log-spaced Gaussian training signals.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' set.seed(7)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_nnksvd_unfold(A, b, rep(0.5, 60), n_atoms = 8L,
#'                          n_dictionary_iterations = 2)
solve_nnksvd_unfold <- function(A, b, x0 = NULL, n_atoms = 15L, sparsity = 2L,
                                dictionary = NULL, training_signals = NULL,
                                n_dictionary_iterations = 80L,
                                lambda_tik = 0.01, prior_wt = 0.5,
                                sparse_coder = "nnls_topk",
                                random_state = NULL, tolerance = 1e-6,
                                n_nnls_iter = NULL, E_MeV = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)

    if (is.null(random_state)) random_state <- 42L

    # ---- Dictionary acquisition ----
    if (!is.null(dictionary)) {
        D <- as.matrix(dictionary); storage.mode(D) <- "double"
        if (nrow(D) != n) {
            stop("Dictionary first dimension (", nrow(D),
                 ") must match the number of energy bins (", n, ").",
                 call. = FALSE)
        }
        D <- .normalize_dictionary(pmax(D, 0.0))
        alpha_prior <- NULL
    } else {
        if (!is.null(training_signals)) {
            signals <- as.matrix(training_signals); storage.mode(signals) <- "double"
            if (nrow(signals) != n) {
                stop("Training signals first dimension (", nrow(signals),
                     ") must match the number of energy bins (", n, ").",
                     call. = FALSE)
            }
            signals <- pmax(signals, 0.0)
        } else {
            # Synthesize training signals for online K-SVD.
            if (!is.null(x0) && any(x0 != 0)) {
                base <- pmax(as.numeric(x0), 0)
                nrm <- sqrt(sum(base^2))
                base <- if (nrm > 0) base / nrm else rep(1 / sqrt(n), n)
            } else {
                base <- rep(1 / sqrt(n), n)
            }
            n_basis <- max(as.integer(n_atoms) * 2L, 8L)
            n_basis <- min(n_basis, n)
            nE <- length(E_MeV)
            if (!is.null(E_MeV) && nE > 0 && any(as.numeric(E_MeV) > 0)) {
                log_E <- log10(pmax(as.numeric(E_MeV), 1e-15))
                centers <- seq(log_E[1], log_E[nE], length.out = n_basis)
                width <- (log_E[nE] - log_E[1]) / max(n_basis * 1.5, 1.0)
                signals <- matrix(0.0, nrow = n, ncol = n_basis + 1L)
                for (i in seq_len(n_basis)) {
                    col <- exp(-((log_E - centers[i])^2) / (2 * width^2))
                    nrm <- sqrt(sum(col^2))
                    if (nrm > 0) col <- col / nrm
                    signals[, i] <- col
                }
            } else {
                t <- seq(0, 1, length.out = n)
                signals <- matrix(0.0, nrow = n, ncol = n_basis + 1L)
                for (i in seq_len(n_basis)) {
                    center <- (i - 1) / (n_basis + 1)
                    col <- exp(-((t - center)^2) / (2 * (1.0 / n_basis)^2))
                    nrm <- sqrt(sum(col^2))
                    if (nrm > 0) col <- col / nrm
                    signals[, i] <- col
                }
            }
            signals[, ncol(signals)] <- base
        }

        fit <- solve_nnksvd(signals, n_atoms = n_atoms,
                            n_iterations = n_dictionary_iterations,
                            sparsity = sparsity, lambda_tik = lambda_tik,
                            prior_wt = prior_wt, sparse_coder = sparse_coder,
                            random_state = random_state, tolerance = tolerance)
        D <- fit$D
        alpha_prior <- fit$alpha_prior
    }

    # ---- Sparse inversion on the equivalent detection dictionary ----
    M <- .build_equivalent_dictionary(A, D, normalize = TRUE)

    # With a pre-learned dictionary the training prior is unavailable, so it is
    # disabled silently (mirrors the Python effective_prior_wt).
    effective_prior_wt <- if (is.null(alpha_prior)) 0.0 else prior_wt

    alpha <- if (sparse_coder == "nnls_topk") {
        solve_nnls_topk(M, b, sparsity = sparsity, lambda_tik = lambda_tik,
                        prior_wt = effective_prior_wt, alpha_prior = alpha_prior,
                        max_iter = n_nnls_iter)
    } else if (sparse_coder == "omp") {
        solve_omp(M, b, sparsity = sparsity, tolerance = tolerance)
    } else if (sparse_coder == "nn_omp") {
        solve_nn_omp(M, b, sparsity = sparsity, tolerance = tolerance)
    } else {
        stop("Unknown sparse_coder '", sparse_coder,
             "'. Expected 'nnls_topk', 'omp' or 'nn_omp'.", call. = FALSE)
    }

    # ---- Spectrum reconstruction ----
    # M_norm[, k] = (A D[, k]) / ||A D[, k]||, so alpha lives in the normalised
    # basis: undo the normalisation by rescaling each atom by ||A D[, k]|| so
    # that A %*% phi ~ b.
    AD <- A %*% D
    atom_norms <- sqrt(colSums(AD^2))
    atom_norms <- ifelse(atom_norms == 0, 1.0, atom_norms)
    phi <- pmax(as.numeric(D %*% (alpha * atom_norms)), 0.0)

    # Optional scale alignment (mirrors unfold_cs.solve_cs).
    computed <- as.numeric(A %*% phi)
    if (sqrt(sum(computed^2)) > 0 && sqrt(sum(b^2)) > 0) {
        scale <- (sum(b * computed)) / (sum(computed * computed) + 1e-12)
        if (scale > 0) phi <- phi * scale
    }

    residual <- sqrt(sum((as.numeric(A %*% phi) - b)^2))
    converged <- residual < tolerance * max(1.0, sqrt(sum(b^2)))

    list(spectrum = phi,
         iterations = as.integer(n_dictionary_iterations),
         converged = as.logical(converged),
         dictionary = D,
         alpha = alpha,
         alpha_prior = alpha_prior,
         M = M)
}

# ---------------------------------------------------------------------------
# Detector-level wrapper
# ---------------------------------------------------------------------------

#' Wrapper around \code{\link{solve_nnksvd_unfold}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @param ln_steps Numeric per-bin \eqn{d(\ln E)} widths. Accepted for API
#'   parity with the Python original; the R workflow derives the widths from
#'   \code{E_MeV} inside \code{\link{run_unfolding}}.
#' @param n_atoms Integer; number of dictionary atoms. Default \code{15L}.
#' @param sparsity Integer; target sparsity K. Default \code{2L}.
#' @param dictionary Optional pre-learned non-negative dictionary \code{(n x p)}.
#' @param training_signals Optional K-SVD training signals \code{(n x m)}.
#' @param n_dictionary_iterations Integer; K-SVD iterations. Default \code{80L}.
#' @param lambda_tik Numeric; Tikhonov weight. Default \code{0.01}.
#' @param prior_wt Numeric; training-sample prior weight. Default \code{0.5}.
#' @param sparse_coder Character; \code{"nnls_topk"}, \code{"omp"} or
#'   \code{"nn_omp"}. Default \code{"nnls_topk"}.
#' @param tolerance Numeric; convergence tolerance. Default \code{1e-6}.
#' @param n_nnls_iter Integer or \code{NULL}; maximum NNLS iterations.
#' @param reading_uncertainties,reading_covariance,noise_model,measurement_time
#'   Monte-Carlo options accepted for API parity with the Python original.
#' @return A result list as produced by \code{\link{run_unfolding}} with extra
#'   keys \code{n_atoms}, \code{sparsity}, \code{sparse_coder},
#'   \code{lambda_tik}, \code{prior_wt}.
#' @export
#' @examples
#' rf <- RF_PTB()
#' det <- Detector$new(response_function = rf,
#'                     detector_names = names(rf)[-1])
#' \dontrun{
#' det$unfold_nnksvd(setNames(rep(1, length(det$detector_names)),
#'                            det$detector_names))
#' }
unfold_nnksvd <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings, ln_steps = NULL,
                          initial_spectrum = NULL,
                          n_atoms = 15L,
                          sparsity = 2L,
                          dictionary = NULL,
                          training_signals = NULL,
                          n_dictionary_iterations = 80L,
                          lambda_tik = 0.01,
                          prior_wt = 0.5,
                          sparse_coder = "nnls_topk",
                          calculate_errors = FALSE,
                          noise_level = 0.01,
                          n_montecarlo = 100L,
                          save_result = FALSE,
                          random_state = NULL,
                          tolerance = 1e-6,
                          n_nnls_iter = NULL,
                          reading_uncertainties = NULL,
                          reading_covariance = NULL,
                          noise_model = "gaussian",
                          measurement_time = NULL,
                          max_neutron_energy = NULL) {
    x0_default <- rep(0, n_energy_bins)

    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_nnksvd_unfold,
                                        n_atoms = n_atoms,
                                        sparsity = sparsity,
                                        dictionary = dictionary,
                                        training_signals = training_signals,
                                        n_dictionary_iterations = n_dictionary_iterations,
                                        lambda_tik = lambda_tik,
                                        prior_wt = prior_wt,
                                        sparse_coder = sparse_coder,
                                        random_state = random_state,
                                        tolerance = tolerance,
                                        n_nnls_iter = n_nnls_iter,
                                        E_MeV = E_MeV),
        solve_kwargs = list(),
        method_name = "NNKSVD",
        extra_output = list(n_atoms = n_atoms,
                            sparsity = sparsity,
                            sparse_coder = sparse_coder,
                            lambda_tik = lambda_tik,
                            prior_wt = prior_wt),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
