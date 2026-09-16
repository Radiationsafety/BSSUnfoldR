#' Statistical Regularization (Turchin's method) unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_statreg.py}. Implements
#' Turchin's method of statistical regularization:
#' \deqn{\hat{\phi} = \arg\min \{ 1/2 \|\Sigma^{-1/2}(A\phi - b)\|^2 + 1/2 \alpha \|D_2 \phi\|^2 \}}
#' where \eqn{D_2} is the second-order finite-difference operator. The
#' regularization parameter \eqn{\alpha} is selected automatically via the
#' L-curve (maximum curvature) or set by the user.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Ignored; accepted for API compatibility.
#' @param E_MeV Optional energy grid (used only for log-energy penalty scaling
#'   in the Python version; here kept for API compatibility).
#' @param unfoldermethod Character; \code{"EmpiricalBayes"} (L-curve, default)
#'   or \code{"User"} (fixed alpha).
#' @param regularization Optional numeric; alpha for \code{"User"} method.
#' @param basis_name Character; ignored (API compatibility).
#' @param boundary Character; ignored (API compatibility).
#' @param derivative_degree Integer; derivative order for penalty (only 2
#'   implemented). Default 2.
#' @return A list \code{list(spectrum, iterations = 0L, converged = TRUE)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_statreg(A, b, NULL, E_MeV = c(1e-9, 1e-6, 1e-3))
solve_statreg <- function(A, b, x0 = NULL, E_MeV = NULL,
                            unfoldermethod = "EmpiricalBayes",
                            regularization = NULL,
                            basis_name = "CubicSplines",
                            boundary = NULL,
                            derivative_degree = 2L) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n_ene <- ncol(A)
    if (any(b < 0)) stop("STREG requires strictly positive measurements")
    if (any(b == 0)) {
        keep <- b > 0
        A <- A[keep, , drop = FALSE]
        b <- b[keep]
        if (length(b) == 0L) {
            stop("STREG requires strictly positive measurements")
        }
    }
    L <- as.matrix(create_derivative_matrix(n_ene, 2L))
    sigma <- pmax(b * 0.05, 1e-300)
    sigma_inv <- 1.0 / sigma
    A_tilde <- A * matrix(sigma_inv, nrow = nrow(A), ncol = ncol(A), byrow = FALSE)
    b_tilde <- b * sigma_inv

    if (unfoldermethod == "User") {
        alpha <- if (is.null(regularization)) 1e-4 else as.numeric(regularization)
    } else if (unfoldermethod == "EmpiricalBayes") {
        alpha <- .statreg_lcurve(A_tilde, b_tilde, L)
    } else {
        stop("Unknown method: ", unfoldermethod)
    }
    ATA <- crossprod(A_tilde)
    ATb <- as.numeric(crossprod(A_tilde, b_tilde))
    LTL <- crossprod(L)
    x <- tryCatch(as.numeric(qr.solve(ATA + alpha * LTL, ATb)),
                  error = function(e) {
                      qr.solve(ATA + alpha * LTL, ATb, tol = 1e-8)
                  })
    list(spectrum = pmax(as.numeric(x), 0.0), iterations = 0L,
         converged = TRUE)
}

.statreg_lcurve <- function(A_tilde, b_tilde, L,
                             n_alphas = 50L,
                             alpha_range = c(1e-8, 1e3)) {
    alphas <- 10^seq(log10(alpha_range[1L]), log10(alpha_range[2L]),
                     length.out = n_alphas)
    ATA <- crossprod(A_tilde)
    ATb <- as.numeric(crossprod(A_tilde, b_tilde))
    LTL <- crossprod(L)
    residuals <- numeric(n_alphas)
    norms <- numeric(n_alphas)
    for (i in seq_along(alphas)) {
        alpha <- alphas[i]
        x <- tryCatch({
            x_int <- qr.solve(ATA + alpha * LTL, ATb)
            pmax(as.numeric(x_int), 0)
        }, error = function(e) NULL)
        if (is.null(x)) { next }
        residuals[i] <- sqrt(sum((as.numeric(A_tilde %*% x) - b_tilde)^2))
        norms[i] <- sqrt(sum(as.numeric(L %*% x)^2))
    }
    valid <- residuals > 0 & norms > 0
    if (sum(valid) < 3L) return(1.0)
    log_res <- log(pmax(residuals[valid], 1e-300))
    log_norm <- log(pmax(norms[valid], 1e-300))
    p1 <- c(log_res[1L], log_norm[1L])
    p2 <- c(log_res[length(log_res)], log_norm[length(log_norm)])
    v <- p2 - p1
    edge <- sqrt(sum(v^2))
    if (edge < 1e-300) return(alphas[floor(length(alphas) / 2)])
    distances <- abs(v[1L] * (p1[2L] - log_norm) -
                       v[2L] * (p1[1L] - log_res)) / edge
    alphas[valid][which.max(distances)]
}

#' Wrapper around \code{\link{solve_statreg}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_statreg
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_statreg <- function(detector_names, n_energy_bins, E_MeV,
                              sensitivities, cc_icrp116, save_result_callback,
                              readings, initial_spectrum = NULL,
                              unfoldermethod = "EmpiricalBayes",
                              regularization = NULL,
                              basis_name = "CubicSplines",
                              boundary = NULL,
                              derivative_degree = 2L,
                              calculate_errors = FALSE,
                              noise_level = 0.01, n_montecarlo = 100L,
                              save_result = FALSE, random_state = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_statreg,
                                         E_MeV = E_MeV,
                                         unfoldermethod = unfoldermethod,
                                         regularization = regularization,
                                         basis_name = basis_name,
                                         boundary = boundary,
                                         derivative_degree = derivative_degree),
        solve_kwargs = list(),
        method_name = "StatReg",
        extra_output = list(
            unfoldermethod = unfoldermethod,
            basis_name = basis_name,
            derivative_degree = derivative_degree
        ),
        calculate_errors = calculate_errors,
        noise_level = noise_level, n_montecarlo = n_montecarlo,
        random_state = random_state, save_result = save_result
    )
}
