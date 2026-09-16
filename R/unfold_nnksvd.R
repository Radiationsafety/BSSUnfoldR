#' Non-negative K-SVD unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_nnksvd.py}.
#' Two-step pipeline: (1) non-negative K-SVD dictionary learning, (2)
#' sparse inversion via NNLS-top-K / NN-OMP / Tikhonov-NNLS.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial spectrum (used as training signal).
#' @param n_atoms Integer; number of dictionary atoms. Default \code{NULL} = n/5.
#' @param sparsity Integer; sparsity for NN-OMP. Default 5.
#' @param lambda_tik Numeric; Tikhonov regularization for NNLS. Default 0.01.
#' @param n_dictionary_iterations Integer; K-SVD iterations. Default 10.
#' @param method Character: \code{"nnls_topk"}, \code{"nn_omp"}, \code{"tikhonov_nnls"}. Default \code{"tikhonov_nnls"}.
#' @param random_state Optional seed.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @examples
#' A <- matrix(runif(5 * 60), nrow = 5)
#' b <- as.numeric(A %*% rep(0.5, 60))
#' set.seed(7)
#' r <- solve_nnksvd(A, b, rep(0.5, 60), n_atoms = 10L,
#'                    n_dictionary_iterations = 3)
solve_nnksvd <- function(A, b, x0 = NULL, n_atoms = NULL, sparsity = 5L,
                          lambda_tik = 0.01, n_dictionary_iterations = 10L,
                          method = "tikhonov_nnls", random_state = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    if (is.null(n_atoms)) n_atoms <- max(4L, n %/% 5L)
    if (!is.null(random_state)) set.seed(as.integer(random_state))
    # Training signals: use x0 (or flat) as the only training sample
    if (is.null(x0)) x0 <- rep(mean(b) / max(mean(rowSums(A)), 1e-10), n)
    signals <- matrix(pmax(as.numeric(x0), 1e-10), ncol = 1L)
    # Non-negative K-SVD dictionary learning
    D <- .nnksvd_learn(signals, n_atoms, n_dictionary_iterations, sparsity)
    AD <- A %*% D  # m x n_atoms
    if (method == "nnls_topk") {
        # Global NNLS, then top-K screening, then local NNLS
        coefs <- tryCatch(as.numeric(lsei::nnls(AD, b)$x),
                          error = function(e) rep(0, n_atoms))
        K <- min(sparsity, n_atoms)
        topk <- order(abs(coefs), decreasing = TRUE)[seq_len(K)]
        AD_k <- AD[, topk, drop = FALSE]
        coefs_k <- tryCatch(as.numeric(lsei::nnls(AD_k, b)$x),
                             error = function(e) rep(0, K))
        alpha <- numeric(n_atoms); alpha[topk] <- coefs_k
    } else if (method == "nn_omp") {
        alpha <- .nn_omp(AD, b, sparsity)
    } else {
        # Tikhonov-NNLS via augmented matrix
        Aw <- rbind(AD, sqrt(lambda_tik) * diag(n_atoms))
        bw <- c(b, rep(0, n_atoms))
        alpha <- tryCatch(as.numeric(lsei::nnls(Aw, bw)$x),
                          error = function(e) rep(0, n_atoms))
    }
    spectrum <- as.numeric(D %*% alpha)
    list(spectrum = pmax(spectrum, 0),
         iterations = as.integer(n_dictionary_iterations),
         converged = TRUE, dictionary = D, alpha = alpha)
}

.nnksvd_learn <- function(signals, n_atoms, n_iterations, sparsity) {
    n <- nrow(signals); m <- ncol(signals)
    n_atoms <- min(n_atoms, n)
    idx <- sample.int(n, size = n_atoms, replace = FALSE)
    D <- signals[idx, , drop = FALSE]
    # Ensure non-negative
    D <- pmax(D, 0)
    norms <- sqrt(rowSums(D^2)); norms <- ifelse(norms == 0, 1, norms)
    D <- D / norms
    for (it in seq_len(n_iterations)) {
        coefs <- matrix(0.0, nrow = n_atoms, ncol = m)
        for (j in seq_len(m)) {
            coefs[, j] <- solve_omp(t(D), signals[, j], sparsity)
        }
        for (atom in seq_len(n_atoms)) {
            used <- which(coefs[atom, ] != 0)
            if (length(used) == 0L) next
            D_restricted <- D
            D_restricted[atom, ] <- 0
            E <- signals[, used, drop = FALSE] -
                 t(D_restricted) %*% coefs[, used, drop = FALSE]
            # Wait, D is n_atoms x n (transposed from Python convention).
            # Let me fix the orientation.
            # Actually in R, D is (n_atoms x n), so D[atom,] is a row vector
            # of length n. signals is (n x m). E should be (n x length(used)).
            # D_restricted %*% coefs[, used] gives (n_atoms x length(used)).
            # We want signals - t(D_restricted) %*% coefs = (n x m).
            # That's not right because D_restricted is (n_atoms x n).
            # Let me fix: E = signals[, used] - t(D_restricted) %*% coefs[, used]
            # But D_restricted is (n_atoms x n), so t(D_restricted) is (n x n_atoms),
            # and coefs[, used] is (n_atoms x length(used)), so t(D_restricted) %*%
            # coefs[, used] is (n x length(used)). ✓
            E <- signals[, used, drop = FALSE] -
                 t(D_restricted) %*% coefs[, used, drop = FALSE]
            svdE <- svd(E, nu = min(dim(E)), nv = min(dim(E)))
            new_atom <- svdE$u[, 1]
            new_atom <- pmax(new_atom, 0)  # non-negative truncation
            norm_na <- sqrt(sum(new_atom^2))
            if (norm_na > 0) new_atom <- new_atom / norm_na
            D[atom, ] <- new_atom
            coefs[atom, used] <- svdE$d[1] * svdE$v[, 1]
        }
        # Re-normalize
        norms <- sqrt(rowSums(D^2)); norms <- ifelse(norms == 0, 1, norms)
        D <- D / norms
    }
    D  # (n_atoms x n)
}

.nn_omp <- function(D, y, sparsity) {
    # D is (m x n_atoms), y is (m,). Returns (n_atoms,) non-negative sparse vector.
    k <- ncol(D)
    alpha <- numeric(k)
    residual <- y
    support <- integer(0L)
    norms <- sqrt(colSums(D^2)); norms <- ifelse(norms == 0, 1, norms)
    D_norm <- sweep(D, 2L, norms, "/")
    for (i in seq_len(min(sparsity, k))) {
        correlations <- as.numeric(t(D_norm) %*% residual)
        correlations[support] <- -1
        idx <- which.max(correlations)
        if (correlations[idx] <= 0) break
        support <- c(support, idx)
        D_s <- D[, support, drop = FALSE]
        coefs <- tryCatch(as.numeric(lsei::nnls(D_s, y)$x),
                          error = function(e) qr.solve(D_s, y))
        residual <- y - as.numeric(D_s %*% coefs)
    }
    if (length(support) > 0L) {
        D_s <- D[, support, drop = FALSE]
        coefs <- tryCatch(as.numeric(lsei::nnls(D_s, y)$x),
                          error = function(e) qr.solve(D_s, y))
        alpha[support] <- pmax(coefs, 0)
    }
    alpha
}

#' Tikhonov-regularized NNLS
#'
#' @param A Numeric matrix (m x n).
#' @param b Numeric vector (length m).
#' @param lambda Numeric; Tikhonov regularization. Default 0.01.
#' @return Numeric vector (length n).
#' @export
#' @examples
#' A <- matrix(runif(3 * 5), nrow = 3)
#' b <- c(1, 2, 3)
#' x <- solve_tikhonov_nnls(A, b, lambda = 0.01)
solve_tikhonov_nnls <- function(A, b, lambda = 0.01) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    Aw <- rbind(A, sqrt(lambda) * diag(n))
    bw <- c(b, rep(0, n))
    tryCatch(as.numeric(lsei::nnls(Aw, bw)$x),
             error = function(e) as.numeric(qr.solve(A, b)))
}

#' Non-negative OMP
#'
#' @param A Numeric dictionary matrix (m x k).
#' @param b Numeric signal vector (length m).
#' @param sparsity Integer; maximum non-zero entries.
#' @return Numeric vector (length k).
#' @export
#' @examples
#' A <- matrix(runif(5 * 10), nrow = 5)
#' b <- c(1, 2, 3, 4, 5)
#' x <- solve_nn_omp(A, b, sparsity = 3L)
solve_nn_omp <- function(A, b, sparsity) {
    .nn_omp(as.matrix(A), as.numeric(b), as.integer(sparsity))
}

#' NNLS with top-K screening
#'
#' @param A Numeric matrix (m x k).
#' @param b Numeric vector (length m).
#' @param K Integer; number of top atoms to keep.
#' @return Numeric vector (length k).
#' @export
#' @examples
#' A <- matrix(runif(5 * 10), nrow = 5)
#' b <- c(1, 2, 3, 4, 5)
#' x <- solve_nnls_topk(A, b, K = 3L)
solve_nnls_topk <- function(A, b, K) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    k <- ncol(A)
    K <- min(as.integer(K), k)
    coefs <- tryCatch(as.numeric(lsei::nnls(A, b)$x),
                      error = function(e) rep(0, k))
    topk <- order(abs(coefs), decreasing = TRUE)[seq_len(K)]
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
                             n_atoms = NULL, sparsity = 5L,
                             lambda_tik = 0.01,
                             n_dictionary_iterations = 10L,
                             method = "tikhonov_nnls",
                             random_state = NULL,
                             calculate_errors = FALSE,
                             noise_level = 0.01, n_montecarlo = 100L,
                             save_result = FALSE) {
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
                                         n_dictionary_iterations = n_dictionary_iterations,
                                         method = method,
                                         random_state = random_state),
        solve_kwargs = list(),
        method_name = "NN-KSVD",
        extra_output = list(n_atoms = if (is.null(n_atoms))
                                            max(4L, n_energy_bins %/% 5L)
                                         else as.integer(n_atoms),
                            sparsity = sparsity,
                            method = method),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
