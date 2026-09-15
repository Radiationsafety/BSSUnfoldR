#' Sparse / dense matrix utilities used by BSSUnfoldR
#'
#' Helper constructors for finite-difference derivative matrices and the
#' Tikhonov regularization operator \eqn{L}. These mirror
#' \code{core/_matrix_utils.py} of the Python \code{bssunfold} package.
#'
#' @name matrix-utils
#' @rdname matrix-utils
NULL

#' @rdname matrix-utils
#'
#' @description
#' \code{create_derivative_matrix} returns a sparse (CsparseMatrix, via
#' \pkg{Matrix}) finite-difference operator. Order 1 gives a first-difference
#' matrix of shape \eqn{(n-1) x n}; order 2 gives a second-difference matrix
#' of shape \eqn{(n-2) x n} (stencil \code{[1, -2, 1]}).
#'
#' @param n Integer; length of the spectrum the operator acts on.
#' @param order Integer; 1 or 2 for first or second derivative.
#' @return A \code{\link[Matrix]{CsparseMatrix-class}} object.
#' @export
#' @examples
#' L1 <- create_derivative_matrix(6, 1)
#' L2 <- create_derivative_matrix(6, 2)
#' dim(L1); dim(L2)
create_derivative_matrix <- function(n, order) {
    n <- as.integer(n); order <- as.integer(order)
    if (n < 1L) stop("'n' must be a positive integer")
    if (order == 1L) {
        if (n < 2L) stop("n must be >= 2 for order=1")
        # First-difference matrix of shape (n-1) x n with stencil [-1, 1]
        i <- c(seq_len(n - 1L), seq_len(n - 1L))
        j <- c(seq_len(n - 1L), seq_len(n - 1L) + 1L)
        x <- c(rep(-1, n - 1L), rep(1, n - 1L))
        Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n - 1L, n),
                             dimnames = list(NULL, NULL))
    } else if (order == 2L) {
        if (n < 3L) stop("n must be >= 3 for order=2")
        # Second-difference matrix of shape (n-2) x n with stencil [1, -2, 1].
        # Built as an explicit sparseMatrix so it works for any n >= 3
        # (bandSparse rejects small m because max(k) = 2 but m - 1 = 0 for n=3).
        rows <- n - 2L
        col_starts <- seq_len(rows)
        i <- rep(col_starts, each = 3L)
        j <- c(rbind(col_starts, col_starts + 1L, col_starts + 2L))
        x <- rep(c(1, -2, 1), rows)
        Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(rows, n),
                             dimnames = list(NULL, NULL))
    } else {
        stop("Unsupported derivative order: ", order, ". Use 1 or 2.")
    }
}

#' @rdname matrix-utils
#'
#' @description
#' \code{make_regularization_operator} returns the dense \code{L} operator
#' used by Tikhonov-regularized solvers.
#'
#' @param smoothness_order Integer; 0 (identity), 1 or 2.
#' @param identity_for_zero Logical; if \code{TRUE} (default) returns
#'   \code{diag(n)} for order 0; if \code{FALSE} returns \code{NULL} so callers
#'   can skip the regularization term entirely.
#' @return A dense numeric matrix of shape \eqn{(n-k) x n} (or
#'   \eqn{n x n} for order 0), or \code{NULL} for order 0 when
#'   \code{identity_for_zero = FALSE}.
#' @export
#' @examples
#' L0 <- make_regularization_operator(6, 0)            # identity
#' L1 <- make_regularization_operator(6, 1)            # first-difference
#' L2 <- make_regularization_operator(6, 2)            # second-difference
#' dim(L1); dim(L2)
make_regularization_operator <- function(n, smoothness_order,
                                         identity_for_zero = TRUE) {
    n <- as.integer(n); smoothness_order <- as.integer(smoothness_order)
    if (smoothness_order == 0L) {
        if (identity_for_zero) return(diag(n))
        return(NULL)
    }
    if (!(smoothness_order %in% c(1L, 2L))) {
        stop("Unsupported smoothness_order: ", smoothness_order,
             ". Use 0, 1 or 2.")
    }
    as.matrix(create_derivative_matrix(n, smoothness_order))
}

#' @rdname matrix-utils
#'
#' @description
#' \code{build_tikhonov_system} solves the regularized least-squares system
#' \eqn{(A^T A + alpha * L^T L) x = A^T b}{(A'A + alpha * L'L) x = A'b}.
#'
#' @param A Numeric response matrix.
#' @param b Numeric measurement vector.
#' @param alpha Numeric regularization parameter.
#' @param L Numeric regularization matrix (or \code{NULL} for identity).
#' @return Numeric vector \code{x}; if the solve fails, returns \code{NULL}.
#' @export
build_tikhonov_system <- function(A, b, alpha, L = NULL) {
    A <- as.matrix(A); b <- as.numeric(b)
    if (is.null(L)) L <- diag(ncol(A))
    L <- as.matrix(L)
    tryCatch({
        P <- crossprod(A) + alpha * crossprod(L)
        x <- solve(P, crossprod(A, b))
        pmax(as.numeric(x), 0)
    }, error = function(e) NULL)
}

#' @rdname matrix-utils
#'
#' @description
#' \code{compute_svd_components} returns \eqn{U, s, V^T, s^2}{U, s, V', s^2}
#' for use by GCV / L-curve computations.
#'
#' @return A list with components \code{U}, \code{s}, \code{Vt}, \code{s_sq}.
#' @export
compute_svd_components <- function(A) {
    A <- as.matrix(A)
    sv <- svd(A, nu = min(dim(A)), nv = min(dim(A)))
    list(
        U = sv$u,
        s = sv$d,
        Vt = t(sv$v),
        s_sq = sv$d^2
    )
}

#' @rdname matrix-utils
#'
#' @description
#' \code{compute_log_steps} computes log10 bin-width steps for an energy grid
#' using edge differences and central differences for interior points.
#'
#' @param E_MeV Numeric energy grid (length n).
#' @return Numeric vector of length n.
#' @export
compute_log_steps <- function(E_MeV) {
    E_MeV <- as.numeric(E_MeV)
    n <- length(E_MeV)
    log_e <- log10(E_MeV + 1e-15)
    out <- numeric(n)
    if (n == 1L) {
        out[1L] <- 1.0
    } else if (n == 2L) {
        out[1L] <- log_e[2L] - log_e[1L]
        out[2L] <- log_e[2L] - log_e[1L]
    } else {
        out[1L] <- log_e[2L] - log_e[1L]
        out[n]  <- log_e[n] - log_e[n - 1L]
        if (n > 3L) {
            out[2:(n - 1L)] <- (log_e[3:n] - log_e[1:(n - 2L)]) / 2
        } else {
            # n == 3 -> only interior point 2
            out[2L] <- (log_e[3L] - log_e[1L]) / 2
        }
    }
    out
}
