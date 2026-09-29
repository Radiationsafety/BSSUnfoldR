#' Projected gradient descent unfolding
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_pgd.py}.
#' Solves \deqn{\min_x \frac12\|Ax-b\|^2 + \frac{\text{reg}}{2}\|x\|^2 \quad
#' \text{s.t. } x \in C} where \eqn{C} is the nonnegative orthant, a box, or the
#' simplex \eqn{\{x \ge 0, \text{sum}(x) = F\}}, by projected gradient descent.
#' The step defaults to \eqn{1/L}{1/L} with \eqn{L = \|A\|_2^2 + \text{reg}}.
#' A Lagrange-duality gap certificate is appended to the result.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Initial spectrum (length n).
#' @param max_iterations Positive integer; default 1000L.
#' @param tolerance Positive numeric relative change tolerance; default 1e-6.
#' @param step_size Optional numeric fixed gradient step; when \code{NULL}
#'   (default) set to \code{1 / L}.
#' @param regularization Numeric Tikhonov (L2) strength; default 0.
#' @param constraint Character; \code{"nonnegative"}, \code{"box"} or
#'   \code{"simplex"}; default \code{"nonnegative"}.
#' @param total_fluence Numeric simplex level for \code{constraint = "simplex"}.
#' @param x_max Numeric upper bound for the \code{"box"} constraint; default \code{Inf}.
#' @param backtracking Logical; use Armijo backtracking when the nominal step
#'   fails; default \code{FALSE}.
#' @return A list \code{list(spectrum, iterations, converged)}.
#' @export
#' @keywords internal
solve_pgd <- function(A, b, x0, max_iterations = 1000L, tolerance = 1e-6,
                      step_size = NULL, regularization = 0.0,
                      constraint = "nonnegative", total_fluence = NULL,
                      x_max = Inf, backtracking = FALSE) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)

    x <- .pgd_project_onto_set(as.numeric(x0), constraint, total_fluence, x_max)

    L <- .bss_spectral_norm(A)^2 + max(as.numeric(regularization), 0.0)
    if (L <= 0) return(list(spectrum = as.numeric(x), iterations = 0L,
                            converged = FALSE))
    t_step <- if (is.null(step_size)) 1.0 / L else as.numeric(step_size)

    objective <- function(z) {
        r <- as.numeric(A %*% z) - b
        0.5 * sum(r * r) + 0.5 * regularization * sum(z * z)
    }

    converged <- FALSE
    iterations <- 0L
    for (k in seq_len(max_iterations)) {
        gradient <- as.numeric(t(A) %*% (as.numeric(A %*% x) - b))
        if (regularization) gradient <- gradient + regularization * x

        x_new <- x - t_step * gradient
        if (isTRUE(backtracking)) {
            direction <- x_new - x
            f0 <- objective(x)
            slope <- sum(gradient * direction)
            trial <- 1.0
            while (trial > 1e-12 &&
                   (objective(x + trial * direction) >
                        f0 + 1e-4 * trial * slope)) {
                trial <- trial * 0.5
            }
            x_new <- x + trial * direction
        }
        x_new <- .pgd_project_onto_set(x_new, constraint, total_fluence, x_max)

        rel_change <- sqrt(sum((x_new - x)^2)) / max(sqrt(sum(x^2)), 1e-30)
        x <- x_new
        iterations <- k
        if (rel_change < tolerance) {
            converged <- TRUE
            break
        }
    }

    list(spectrum = as.numeric(x), iterations = iterations, converged = converged)
}

# Euclidean projection onto the constraint set (nonnegative | box | simplex).
.pgd_project_onto_set <- function(x, constraint = "nonnegative",
                                  total_fluence = NULL, x_max = Inf) {
    constraint <- tolower(as.character(constraint))
    if (constraint == "nonnegative") return(pmax(x, 0.0))
    if (constraint == "box") return(pmin(pmax(x, 0.0), x_max))
    if (constraint == "simplex") {
        if (is.null(total_fluence) || as.numeric(total_fluence) <= 0) {
            stop("simplex projection requires total_fluence > 0")
        }
        return(.pgd_project_simplex(x, as.numeric(total_fluence)))
    }
    stop("Unknown constraint '", constraint,
         "'; expected 'nonnegative', 'box' or 'simplex'")
}

# Euclidean projection onto {x >= 0, sum x = total} (Duchi et al. algorithm).
.pgd_project_simplex <- function(v, total) {
    v <- as.numeric(v)
    n <- length(v)
    u <- sort(v, decreasing = TRUE)
    css <- cumsum(u)
    cond <- u * seq_len(n) + (total - css) > 0
    rho <- max(which(cond))
    theta <- (css[rho] - total) / (rho + 1.0)
    pmax(v - theta, 0.0)
}

# KKT/Lagrange-duality gap certificate for the NNLS problem, mirroring
# core/_dual_diagnostics.py::nnls_duality_gap.
.bss_nnls_duality_gap <- function(A, b, x) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    x <- pmax(as.numeric(x), 0.0)
    G <- crossprod(A)
    c <- as.numeric(crossprod(A, b))
    primal_value <- 0.5 * sum((as.numeric(A %*% x) - b)^2)
    lam <- pmax(as.numeric(crossprod(A, as.numeric(A %*% x) - b)), 0.0)
    # Moore-Penrose pseudo-inverse of G = A^T A via the (thin) SVD of A:
    # the singular values of G are d^2, so pinv(G) = V diag(1/d^2) V'.
    sv <- svd(A, nu = 0L, nv = min(dim(A)))
    V <- sv$v
    d2 <- sv$d^2
    tol <- max(dim(A)) * max(d2) * .Machine$double.eps
    dinv2 <- ifelse(d2 > tol, 1 / d2, 0)
    G_pinv <- V %*% diag(dinv2, nrow = length(dinv2)) %*% t(V)
    dual_value <- 0.5 * sum(b * b) -
        0.5 * as.numeric(t(c + lam) %*% (G_pinv %*% (c + lam)))
    gap <- primal_value - dual_value
    floor <- max(primal_value, 1e-8 * sum(b * b), 1e-30)
    list(primal_value = primal_value, dual_value = dual_value,
         duality_gap = max(gap, 0.0), relative_gap = max(gap, 0.0) / floor)
}

#' @rdname solve_pgd
#' @inheritParams run_unfolding
#' @param initial_spectrum Optional initial spectrum guess.
#' @param constraint Character; \code{"nonnegative"}, \code{"box"} or \code{"simplex"}.
#' @param total_fluence Numeric simplex level (required for \code{"simplex"}).
#' @param x_max Numeric upper bound for the box constraint.
#' @param backtracking Logical; use Armijo backtracking.
#' @param calculate_errors Logical; run Monte-Carlo uncertainty.
#' @param noise_level Numeric; MC noise level.
#' @param n_montecarlo Integer; MC samples.
#' @param save_result Logical; call the save callback.
#' @param random_state Optional integer seed.
#' @param max_neutron_energy Optional energy cutoff in MeV.
#' @return A result list produced by \code{\link{run_unfolding}} with an
#'   additional \code{duality_gap} entry.
#' @export
unfold_pgd <- function(detector_names, n_energy_bins, E_MeV,
                       sensitivities, cc_icrp116, save_result_callback,
                       readings, initial_spectrum = NULL,
                       max_iterations = 1000L, tolerance = 1e-6,
                       regularization = 0.0, constraint = "nonnegative",
                       total_fluence = NULL, x_max = Inf,
                       backtracking = FALSE, calculate_errors = FALSE,
                       noise_level = 0.01, n_montecarlo = 100L,
                       save_result = FALSE, random_state = NULL,
                       max_neutron_energy = NULL) {
    if (constraint == "simplex" && is.null(total_fluence)) {
        stop("constraint='simplex' requires total_fluence")
    }
    x0_default <- rep(0.0, n_energy_bins)

    result <- run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_pgd,
                                        max_iterations = max_iterations,
                                        tolerance = tolerance,
                                        regularization = regularization,
                                        constraint = constraint,
                                        total_fluence = total_fluence,
                                        x_max = x_max,
                                        backtracking = backtracking),
        solve_kwargs = list(),
        method_name = "Projected Gradient Descent",
        extra_output = list(constraint = constraint,
                            total_fluence = total_fluence),
        calculate_errors = calculate_errors,
        noise_level = noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)

    est <- .bss_select_system(detector_names, readings, sensitivities)
    result$duality_gap <- .bss_nnls_duality_gap(est$A, est$b, result$spectrum)
    result
}
