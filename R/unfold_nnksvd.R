#' Non-negative K-SVD unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_nnksvd.py}.
#' Two-step pipeline: (1) non-negative K-SVD dictionary learning, (2)
#' sparse inversion via NNLS-top-K / NN-OMP / Tikhonov-NNLS.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial spectrum (used as training-signal prior).
#' @param n_atoms Integer; number of dictionary atoms. Default \code{NULL} = n/4.
#' @param sparsity Integer; sparsity for sparse coding. Default 2.
#' @param lambda_tik Numeric; Tikhonov regularization for NNLS. Default 0.01.
#' @param prior_wt Numeric; training-sample prior weight. Default 0.5.
#' @param n_dictionary_iterations Integer; K-SVD iterations. Default 80.
#' @param method Character: \code{"nnls_topk"}, \code{"nn_omp"}, \code{"omp"}.
#'   Default \code{"nnls_topk"}.
#' @param E_MeV Optional energy grid for log-spaced training signals.
#' @param random_state Optional seed.
#' @param tolerance Numeric; early-stopping tolerance. Default 1e-6.
#' @return A list \code{list(spectrum, iterations, converged, dictionary, alpha)}.
#' @export
#' @examples
#' set.seed(7)
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' r <- solve_nnksvd(A, b, rep(0.5, 60), n_atoms = 8L,
#'                    n_dictionary_iterations = 3)
solve_nnksvd <- function(A, b, x0 = NULL, n_atoms = NULL, sparsity = 2L,
                          lambda_tik = 0.01, prior_wt = 0.5,
                          n_dictionary_iterations = 80L,
                          method = "nnls_topk", E_MeV = NULL,
                          random_state = NULL, tolerance = 1e-6) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    m <- nrow(A); n <- ncol(A)
    if (is.null(n_atoms)) n_atoms <- max(4L, n %/% 4L)
    n_atoms <- as.integer(n_atoms)
    sparsity <- as.integer(sparsity)
    if (!is.null(random_state)) set.seed(as.integer(random_state))

    # Build training signals for online K-SVD.  The Python implementation
    # places log-spaced Gaussian bumps on the log-energy grid so the
    # dictionary atoms span the full spectral range.  When no E_MeV is
    # supplied we fall back to a uniform normalised index.
    if (!is.null(x0) && any(x0 != 0)) {
        base <- pmax(as.numeric(x0), 0)
        nrm <- sqrt(sum(base^2))
        base <- if (nrm > 0) base / nrm else rep(1 / sqrt(n), n)
    } else {
        base <- rep(1 / sqrt(n), n)
    }
    n_basis <- max(n_atoms * 2L, 8L)
    n_basis <- min(n_basis, n)

    if (!is.null(E_MeV) && any(E_MeV > 0)) {
        log_E <- log10(pmax(as.numeric(E_MeV), 1e-15))
        centers <- seq(log_E[1], log_E[length(log_E)], length.out = n_basis)
        width <- (log_E[length(log_E)] - log_E[1]) / max(n_basis * 1.5, 1.0)
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
            center <- i / (n_basis + 1)
            col <- exp(-((t - center)^2) / (2 * (1.0 / n_basis)^2))
            nrm <- sqrt(sum(col^2))
            if (nrm > 0) col <- col / nrm
            signals[, i] <- col
        }
    }
    signals[, ncol(signals)] <- base
    signals <- pmax(signals, 0.0)

    # Non-negative K-SVD dictionary learning: D is (n_features x n_atoms)
    learn <- .nnksvd_learn(signals, n_atoms, n_dictionary_iterations,
                            sparsity, lambda_tik, method, tolerance)
    D <- learn$D                     # n x p
    alpha_prior <- learn$alpha_prior # length p

    # Equivalent (normalized) detection dictionary M = (A %*% D), columns
    # L2-normalized to remove amplitude differences between atoms.
    AD <- A %*% D                    # m x p
    atom_norms <- sqrt(colSums(AD^2))
    atom_norms <- ifelse(atom_norms == 0, 1.0, atom_norms)
    M <- sweep(AD, 2L, atom_norms, "/")

    # Sparse inversion on the equivalent detection dictionary.
    if (method == "nnls_topk") {
        alpha <- solve_nnls_topk(M, b, K = sparsity, lambda_tik = lambda_tik,
                                  prior_wt = prior_wt,
                                  alpha_prior = alpha_prior)
    } else if (method == "nn_omp") {
        alpha <- .nn_omp(M, b, sparsity, tolerance)
    } else if (method == "omp") {
        alpha <- solve_omp(M, b, sparsity, tolerance)
    } else {
        stop("Unknown method '", method, "'. Expected 'nnls_topk', 'nn_omp' or 'omp'.")
    }

    # Spectrum reconstruction.  Compensate for the M-normalization by
    # rescaling each atom's contribution by ||A %*% D[, k]||.
    phi <- D %*% (alpha * atom_norms)
    phi <- pmax(as.numeric(phi), 0.0)

    # Optional scale alignment: make the predicted counts match the
    # measured counts in total magnitude (mirrors unfold_cs.solve_cs).
    computed <- as.numeric(A %*% phi)
    if (sqrt(sum(computed^2)) > 0 && sqrt(sum(b^2)) > 0) {
        scale <- sum(b * computed) / (sum(computed * computed) + 1e-12)
        if (scale > 0) phi <- phi * scale
    }

    residual <- sqrt(sum((A %*% phi - b)^2))
    converged <- residual < tolerance * max(1.0, sqrt(sum(b^2)))

    list(spectrum = phi,
         iterations = as.integer(n_dictionary_iterations),
         converged = converged,
         dictionary = D,
         alpha = alpha,
         alpha_prior = alpha_prior)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
.normalize_dictionary <- function(D) {
    D <- as.matrix(D); storage.mode(D) <- "double"
    norms <- sqrt(colSums(D^2))
    norms <- ifelse(norms == 0, 1.0, norms)
    sweep(D, 2L, norms, "/")
}

.nnksvd_learn <- function(signals, n_atoms, n_iterations, sparsity,
                           lambda_tik, sparse_coder, tolerance) {
    signals <- as.matrix(signals); storage.mode(signals) <- "double"
    n <- nrow(signals); m <- ncol(signals)
    n_atoms <- max(1L, min(as.integer(n_atoms), m))
    # Initialize dictionary with random training samples (non-negative).
    idx <- sample.int(m, size = n_atoms, replace = FALSE)
    D <- signals[, idx, drop = FALSE]
    # Replace all-zero columns with smooth bumps.
    col_norms <- sqrt(colSums(D^2))
    for (zc in which(col_norms == 0)) {
        bump <- pmax(rnorm(n), 0.0)
        if (sqrt(sum(bump^2)) == 0) bump <- rep(1, n)
        D[, zc] <- bump
    }
    D <- .normalize_dictionary(D)

    coefficients <- matrix(0.0, nrow = n_atoms, ncol = m)

    for (it in seq_len(n_iterations)) {
        D_prev <- D
        # Sparse coding stage (per training sample).
        for (j in seq_len(m)) {
            y_j <- signals[, j]
            if (sparse_coder == "omp") {
                coefficients[, j] <- solve_omp(D, y_j, sparsity, tolerance)
            } else if (sparse_coder == "nn_omp") {
                coefficients[, j] <- .nn_omp(D, y_j, sparsity, tolerance)
            } else {  # nnls_topk
                coefficients[, j] <- solve_nnls_topk(D, y_j, K = sparsity,
                                                      lambda_tik = lambda_tik,
                                                      prior_wt = 0.0,
                                                      alpha_prior = NULL)
            }
        }
        # Dictionary update stage.
        for (atom in seq_len(n_atoms)) {
            used <- which(coefficients[atom, ] != 0)
            if (length(used) == 0L) {
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
            # Rank-1 SVD approximation of E.
            svdE <- svd(E, nu = 1L, nv = 1L)
            new_atom <- svdE$u[, 1]
            new_coef <- svdE$d[1] * svdE$v[, 1]
            # Non-negative truncation (article's key modification).
            new_atom <- pmax(new_atom, 0.0)
            new_coef <- pmax(new_coef, 0.0)
            nrm <- sqrt(sum(new_atom^2))
            if (nrm > 0) {
                D[, atom] <- new_atom / nrm
                coefficients[atom, used] <- new_coef * nrm
            } else {
                j_new <- sample.int(m, 1L)
                new_atom <- pmax(signals[, j_new], 0.0)
                nrm <- sqrt(sum(new_atom^2))
                if (nrm == 0) { new_atom <- rep(1, n); nrm <- sqrt(n) }
                D[, atom] <- new_atom / nrm
                # Recompute coefficient for this atom via NNLS per signal.
                for (u_idx in used) {
                    coefs_u <- tryCatch(as.numeric(lsei::nnls(
                        matrix(D[, atom], ncol = 1L), signals[, u_idx])$x),
                        error = function(e) 0.0)
                    coefficients[atom, u_idx] <- coefs_u[1L]
                }
            }
        }
        # Convergence check on the dictionary change.
        if (sqrt(sum((D - D_prev)^2)) <
            tolerance * max(1.0, sqrt(sum(D_prev^2)))) break
    }
    # Final non-negative safeguard and re-normalization.
    D <- pmax(D, 0.0)
    D <- .normalize_dictionary(D)
    # Training-sample-driven prior: mean sparse code across training signals.
    alpha_prior <- rowMeans(coefficients)
    list(D = D, alpha_prior = alpha_prior)
}

.nn_omp <- function(D, y, sparsity, tolerance = 1e-6) {
    # D is (n x p), y is (n,). Returns (p,) non-negative sparse vector.
    D <- as.matrix(D); storage.mode(D) <- "double"
    y <- as.numeric(y)
    p <- ncol(D)
    alpha <- numeric(p)
    residual <- y
    support <- integer(0L)
    norms <- sqrt(colSums(D^2))
    norms <- ifelse(norms == 0, 1.0, norms)
    D_norm <- sweep(D, 2L, norms, "/")
    for (i in seq_len(min(sparsity, p))) {
        correlations <- as.numeric(t(D_norm) %*% residual)
        if (length(support) > 0L) correlations[support] <- -Inf
        idx <- which.max(correlations)
        if (!is.finite(correlations[idx]) || correlations[idx] <= 0) break
        support <- c(support, idx)
        D_s <- D[, support, drop = FALSE]
        coefs <- tryCatch(as.numeric(lsei::nnls(D_s, y)$x),
                          error = function(e)
                              pmax(qr.solve(D_s, y, tol = 1e-12), 0))
        alpha[support] <- coefs
        residual <- y - as.numeric(D_s %*% coefs)
        if (sqrt(sum(residual^2)) < tolerance) break
    }
    alpha
}

#' Tikhonov-regularized NNLS
#'
#' Solves \eqn{\min || y - A x ||^2 + \lambda ||x||^2} s.t. \eqn{x \ge 0}
#' via the augmented-matrix form.
#'
#' @param A Numeric matrix (m x n).
#' @param b Numeric vector (length m).
#' @param lambda Numeric; Tikhonov regularization. Default 0.01.
#' @param prior_wt Numeric; optional prior weight. Default 0.
#' @param alpha_prior Optional prior coefficient vector (length n).
#' @return Numeric vector (length n).
#' @export
#' @examples
#' set.seed(1)
#' A <- matrix(runif(3 * 5), nrow = 3)
#' b <- c(1, 2, 3)
#' x <- solve_tikhonov_nnls(A, b, lambda = 0.01)
solve_tikhonov_nnls <- function(A, b, lambda = 0.01, prior_wt = 0.0,
                                  alpha_prior = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    A_aug <- rbind(A, sqrt(max(lambda, 0.0)) * diag(n))
    b_aug <- c(b, rep(0, n))
    if (prior_wt > 0.0) {
        if (is.null(alpha_prior))
            stop("alpha_prior must be provided when prior_wt > 0")
        alpha_prior <- as.numeric(alpha_prior)
        if (length(alpha_prior) != n)
            stop("alpha_prior length must match number of columns of A")
        A_aug <- rbind(A_aug, sqrt(prior_wt) * diag(n))
        b_aug <- c(b_aug, sqrt(prior_wt) * alpha_prior)
    }
    tryCatch(as.numeric(lsei::nnls(A_aug, b_aug)$x),
             error = function(e) rep(0, n))
}

#' Non-negative OMP
#'
#' @param A Numeric dictionary matrix (n x k).
#' @param b Numeric signal vector (length n).
#' @param sparsity Integer; maximum non-zero entries.
#' @return Numeric vector (length k).
#' @export
#' @examples
#' set.seed(1)
#' A <- matrix(runif(10 * 5), nrow = 10)
#' b <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' x <- solve_nn_omp(A, b, sparsity = 3L)
solve_nn_omp <- function(A, b, sparsity) {
    .nn_omp(as.matrix(A), as.numeric(b), as.integer(sparsity))
}

#' NNLS with top-K screening
#'
#' Three-stage hierarchical sparse coding:
#' 1. global NNLS coarse solution, 2. top-K atom screening,
#' 3. local NNLS fine optimization on the screened support.
#'
#' @param A Numeric matrix (m x k).
#' @param b Numeric vector (length m).
#' @param K Integer; number of top atoms to keep.
#' @param lambda_tik Numeric; Tikhonov regularization. Default 0.01.
#' @param prior_wt Numeric; prior weight. Default 0.
#' @param alpha_prior Optional prior coefficient vector (length k).
#' @return Numeric vector (length k).
#' @export
#' @examples
#' set.seed(1)
#' A <- matrix(runif(5 * 10), nrow = 5)
#' b <- c(1, 2, 3, 4, 5)
#' x <- solve_nnls_topk(A, b, K = 3L)
solve_nnls_topk <- function(A, b, K, lambda_tik = 0.01, prior_wt = 0.0,
                             alpha_prior = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    k <- ncol(A)
    K <- max(0L, min(as.integer(K), k))
    if (K == 0L) return(numeric(k))
    # Step 1: global NNLS coarse solution (Tikhonov-regularized).
    alpha_full <- solve_tikhonov_nnls(A, b, lambda = lambda_tik,
                                       prior_wt = prior_wt,
                                       alpha_prior = alpha_prior)
    # Step 2: top-K atom screening (largest coefficients).
    topk <- order(abs(alpha_full), decreasing = TRUE)[seq_len(K)]
    # Step 3: local NNLS fine optimization on the screened support.
    A_k <- A[, topk, drop = FALSE]
    coefs_k <- tryCatch(as.numeric(lsei::nnls(A_k, b)$x),
                        error = function(e) rep(0, K))
    alpha <- numeric(k); alpha[topk] <- coefs_k
    alpha
}

#' Wrapper around \code{\link{solve_nnksvd}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_nnksvd
#' @export
unfold_nnksvd <- function(detector_names, n_energy_bins, E_MeV,
                             sensitivities, cc_icrp116, save_result_callback,
                             readings, initial_spectrum = NULL,
                             n_atoms = NULL, sparsity = 2L,
                             lambda_tik = 0.01, prior_wt = 0.5,
                             n_dictionary_iterations = 80L,
                             method = "nnls_topk",
                             random_state = NULL,
                             tolerance = 1e-6,
                             calculate_errors = FALSE,
                             noise_level = 0.01, n_montecarlo = 100L,
                             save_result = FALSE,
                             max_neutron_energy = NULL) {
    x0_default <- rep(0.5, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_nnksvd,
                                         n_atoms = n_atoms,
                                         sparsity = sparsity,
                                         lambda_tik = lambda_tik,
                                         prior_wt = prior_wt,
                                         n_dictionary_iterations = n_dictionary_iterations,
                                         method = method,
                                         E_MeV = E_MeV,
                                         random_state = random_state,
                                         tolerance = tolerance),
        solve_kwargs = list(),
        method_name = "NN-KSVD",
        extra_output = list(n_atoms = if (is.null(n_atoms))
                                            max(4L, n_energy_bins %/% 4L)
                                         else as.integer(n_atoms),
                            sparsity = sparsity,
                            method = method),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
