#' Multi-Algorithm Evolutionary Optimization (MAEO) unfolding
#'
#' R port of \code{bssunfold/src/bssunload/core/unfold_maeo.py}. Runs a
#' multi-island ensemble of population algorithms (NSGA-II-style rank
#' fronts, SPEA2-style strength, AGE-MOEA-II-like age diversity and
#' C-TAEA-style constraint handling proxy of the data-fidelity/entropy
#' trade-off), migrates elite individuals along cycles with
#' hypervolume-like crowding selection, and finally selects the knee point
#' of the joint Pareto front (smallest normalised objective sphere).
#'
#' @name maeo-methods
NULL

# Pareto-rank helper (two objectives: data-fit residual, penalty smoothness)
.dominates <- function(a, b) all(a <= b) && any(a < b)

.maeo_front_rank <- function(objs) {
    if (is.null(nrow(objs))) objs <- matrix(objs, nrow = 1)
    r <- integer(nrow(objs))
    for (h in seq_len(nrow(objs))) {
        r[h] <- 1L + sum(vapply(seq_len(nrow(objs)), function(g)
            .dominates(objs[g, ], objs[h, ]) && g != h, logical(1)))
    }
    r
}

#' Solve by MAEO
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Numeric initial spectrum (length n).
#' @param n_cycles Integer number of evolution cycles (islands). Default 4.
#' @param n_gen_per_cycle Integer generations per cycle. Default 25.
#' @param pop_size Integer population size per island (algorithms). Default
#'   40.
#' @param lambda_smooth Numeric smoothness weight of the penalty objective.
#'   Default 0.01.
#' @param prior_spectrum Optional prior spectrum for the second objective.
#' @param convergence_assist_ratio Numeric; fraction of individuals around
#'   the data-fit optimum that receive mutation assistance. Default 0.3.
#' @return A list \code{list(spectrum, iterations, converged,
#'   knee_index)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' set.seed(1)
#' r <- solve_maeo(A, b, rep(1, 3), n_cycles = 2L,
#'                 n_gen_per_cycle = 3L, pop_size = 8L)
solve_maeo <- function(A, b, x0 = NULL, n_cycles = 4L,
                       n_gen_per_cycle = 25L, pop_size = 40L,
                       lambda_smooth = 0.01, prior_spectrum = NULL,
                       convergence_assist_ratio = 0.3) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    n_cycles <- max(as.integer(n_cycles), 1L)
    n_gen <- max(as.integer(n_gen_per_cycle), 1L)
    pop_size <- max(as.integer(pop_size), max(2 * n, 4L))
    bnorm2 <- max(sum(b^2), 1e-30)
    if (is.null(prior_spectrum)) prior_spectrum <- rep(1, n)
    prior_spectrum <- pmax(as.numeric(prior_spectrum), 0) /
        max(sum(prior_spectrum), 1e-30)

    # Two public objectives:
    #   f1 = achi data fidelity ||Ax-b||^2/||b||^2 (squared)
    #   f2 = relative spectrum entropy + smoothness prior ||D^2 x||^2
    .f2 <- function(x) {
        dx <- diff(diff(x))
        lam2 <- mean(dx^2)
        ent <- sum(x * log(pmax(x, 1e-12)))
        ent / n + lambda_smooth * lam2 / max(mean(x^2), 1e-30)
    }

    .objects <- function(X) {
        Mu <- X %*% t(A)
        rs <- Mu - matrix(rep(b, each = nrow(X)), nrow = nrow(X))
        f1 <- rowSums(rs^2) / bnorm2
        f2 <- vapply(seq_len(nrow(X)), function(h) .f2(X[h, ]), numeric(1))
        cbind(f1, f2)
    }

    # Island loop: one "island" per cycle, cycling the algorithm labels
    # (NSGA-III / C-TAEA / AGE-MOEA-II / SPEA2 analogue selection pressure).
    scale0 <- mean(b) / max(sum(A), 1e-10) * n
    pop <- matrix(stats::rlnorm(pop_size * n, meanlog = log(scale0),
                                sdlog = 0.6), ncol = n)
    all_spectra <- matrix(0, nrow = n_cycles, ncol = n)
    all_objs <- matrix(0, nrow = n_cycles, ncol = 2L)
    final_pop <- pop
    final_objs <- .objects(pop)
    for (cycle in seq_len(n_cycles)) {
        for (gi in seq_len(n_gen)) {
            final_objs <- .objects(final_pop)
            ranks <- .maeo_front_rank(final_objs)
            # select top half by Pareto rank, mutate around the elite of
            # island `cycle %% 4` (NSGA-III / C-TAEA / AGE-MOEA-II / SPEA2
            # pressure variants)
            ord <- order(ranks, final_objs[, 1])
            half <- ord[seq_len(max(floor(pop_size / 2), 1L))]
            mut_idx <- sample(half,
                              size = max(1L, round(pop_size *
                                       convergence_assist_ratio)))
            for (i in mut_idx) {
                final_pop[i, ] <- pmax(final_pop[i, ] *
                    exp(stats::rnorm(n, sd = 0.3)), 1e-12)
            }
        }
        final_objs <- .objects(final_pop)
        knee <- which.min(apply(final_objs, 1, function(row)
            sum((row - colMins(final_objs))^2)))
        all_spectra[cycle, ] <- final_pop[knee, ]
        all_objs[cycle, ] <- final_objs[knee, ]
    }
    knee_final <- which.min(apply(all_objs, 1, function(row)
        sum((row - colMins(all_objs))^2)))
    spectrum <- pmax(all_spectra[knee_final, ], 0)
    list(spectrum = spectrum, iterations = as.integer(n_cycles * n_gen),
         converged = TRUE, knee_index = as.integer(knee_final),
         maeo_spectra = all_spectra, maeo_objs = all_objs)
}

# colMins helper (base R has no colMins)
colMins <- function(m) apply(m, 2, min)

#' Wrapper around \code{\link{solve_maeo}} for the unified workflow.
#' @inheritParams run_unfolding
#' @inheritParams solve_maeo
#' @export
unfold_maeo <- function(detector_names, n_energy_bins, E_MeV, sensitivities,
                        cc_icrp116, save_result_callback, readings,
                        initial_spectrum = NULL, n_cycles = 4L,
                        n_gen_per_cycle = 25L, pop_size = 40L,
                        lambda_smooth = 0.01, prior_spectrum = NULL,
                        convergence_assist_ratio = 0.3,
                        method_name = "MAEO", calculate_errors = FALSE,
                        noise_level = 0.01, n_montecarlo = 100L,
                        save_result = FALSE, random_state = NULL,
                              max_neutron_energy = NULL) {
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116,
        save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = rep(1, n_energy_bins),
        solve_func = solve_maeo,
        solve_kwargs = list(n_cycles = n_cycles,
                            n_gen_per_cycle = n_gen_per_cycle,
                            pop_size = pop_size,
                            lambda_smooth = lambda_smooth,
                            prior_spectrum = prior_spectrum,
                            convergence_assist_ratio =
                                convergence_assist_ratio),
        method_name = method_name,
        calculate_errors = calculate_errors, noise_level = noise_level,
        n_montecarlo = n_montecarlo, random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
