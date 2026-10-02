# Gnowee engine: R port of bssunfold/core/unfold_gnowee.py and _gnowee.py.
# The public roxygen topic for solve_gnowee lives just above its definition.

# ---- boundary helpers (ports of GnoweeHeuristics.simple_bounds /
#      rejection_bounds) ---------------------------------------------------
.gnowee_simple_bounds <- function(child, lb, ub) pmin(pmax(child, lb), ub)

.gnowee_rejection_bounds <- function(parent, child, step_size, lb, ub,
                                     max_reductions = 5L) {
    parent <- as.numeric(parent)
    child <- as.numeric(child)
    step <- as.numeric(step_size)
    oob <- (child < lb) | (child > ub)
    for (i in seq_len(max_reductions)) {
        if (!any(oob)) break
        step[oob] <- step[oob] * 0.5
        child[oob] <- child[oob] - step[oob]
        oob <- (child < lb) | (child > ub)
    }
    # Fall back to the parent for anything still out of bounds.
    child[oob] <- parent[oob]
    child
}

# ---- symmetric Levy stable sampler (Mantegna), port of Sampling.levy -----
.gnowee_levy <- function(nc, nr = 0L, alpha = 1.5, gam = 1.0, n = 1L) {
    invalpha <- 1 / alpha
    sigx <- ((gamma(1 + alpha) * sin(pi * alpha / 2)) /
             (gamma((1 + alpha) / 2) * alpha * 2^((alpha - 1) / 2)))^invalpha
    shape <- if (nr != 0L) c(n, nr, nc) else c(n, nc)
    g1 <- matrix(stats::rnorm(prod(shape)), nrow = n)
    g2 <- matrix(stats::rnorm(prod(shape)), nrow = n)
    v <- sigx * g1 / (abs(g2)^invalpha)
    kappa <- (alpha * gamma((alpha + 1) / (2 * alpha))) / gamma(invalpha) *
        ((alpha * gamma((alpha + 1) / 2)) /
         (gamma(1 + alpha) * sin(pi * alpha / 2)))^invalpha
    # Mantegna's polynomial fit for the temperature parameter c(alpha),
    # evaluated in Horner form (polyval(p, alpha) == p[0]*a^(k-1)+...+p[k-1]).
    p <- c(-17.7767, 113.3855, -281.5879, 337.5439, -193.5494, 44.8754)
    cval <- 0
    for (coef in p) cval <- cval * alpha + coef
    w <- ((kappa - 1) * exp(-abs(v) / cval) + 1) * v
    z <- if (n > 1L) colSums(w) / n^invalpha else w[1L, , drop = TRUE]
    z <- gam^invalpha * z
    if (nr != 0L) matrix(as.numeric(z), nrow = nr, byrow = TRUE)
    else matrix(as.numeric(z), nrow = 1L)[1L, ]
}

# ---- scale-consistent objective factory -----------------------------------
.gnowee_build_fitness <- function(A, b, alpha, norm_, L,
                                  smoothness_weight, entropy_weight) {
    denom <- as.numeric(crossprod(b))
    if (denom <= 0) denom <- 1
    A_fro <- sqrt(sum(A^2))
    if (A_fro <= 0) A_fro <- 1
    x_scale <- sqrt(denom) / A_fro
    x_scale2 <- x_scale * x_scale
    function(y) {
        x <- exp(as.numeric(y))
        residual <- as.numeric(A %*% x) - b
        value <- as.numeric(crossprod(residual)) / denom
        if (alpha > 0) {
            if (norm_ == 2L) {
                value <- value + alpha * as.numeric(crossprod(x)) / x_scale2
            } else if (norm_ == 1L) {
                value <- value + alpha * sum(abs(x)) / x_scale
            }
        }
        if (!is.null(L) && smoothness_weight > 0) {
            Lx <- as.numeric(L %*% x)
            value <- value +
                smoothness_weight * as.numeric(crossprod(Lx)) / x_scale2
        }
        if (entropy_weight > 0) {
            total <- sum(x)
            if (total > 0) {
                p <- x / total
                logp <- log(pmax(p, 1e-300))
                value <- value - entropy_weight * as.numeric(crossprod(p, logp))
            }
        }
        value
    }
}

# ---- initial population sampler -------------------------------------------
.gnowee_initialize <- function(num_samples, method, lb, ub) {
    n_dim <- length(lb)
    if (identical(method, "random")) {
        unit <- matrix(stats::runif(num_samples * n_dim),
                       nrow = num_samples, ncol = n_dim)
    } else {
        # Crude Latin-hypercube: a random permutation per dimension, jittered
        # into the [0, 1] strata.  Matches _gnowee.py's documented fallback.
        unit <- matrix(0, nrow = num_samples, ncol = n_dim)
        for (j in seq_len(n_dim)) {
            perm <- sample.int(num_samples)
            unit[, j] <- (as.numeric(perm) - 1 + stats::runif(num_samples)) /
                num_samples
        }
    }
    lb + (ub - lb) * unit
}

# ---- Gnowee core ----------------------------------------------------------
.gnowee_run <- function(lb, ub, objective, s, seed_solution) {
    n_dim <- length(lb)
    init_num <- max(s$population * 2L, n_dim * 10L)
    init_vars <- .gnowee_initialize(init_num, s$init_sampling, lb, ub)
    if (!is.null(seed_solution)) {
        init_vars[1L, ] <- .gnowee_simple_bounds(seed_solution, lb, ub)
    }
    pop <- lapply(seq_len(nrow(init_vars)), function(i) {
        v <- as.numeric(init_vars[i, ])
        list(variables = v, fitness = objective(v),
             change_count = 0L, stall_count = 0L)
    })
    ord <- order(vapply(pop, function(p) p$fitness, numeric(1)))
    pop <- pop[ord]
    if (length(pop) > s$population) {
        pop <- pop[seq_len(s$population)]
    } else {
        s$population <- length(pop)
    }

    # Timeline: single improvement record (generation, evaluations, fitness,
    # design).  Implemented as a small mutable list-of-lists.
    timeline <- list(list(generation = 0L,
                          evaluations = length(pop),
                          fitness = pop[[1]]$fitness,
                          design = pop[[1]]$variables))

    # --- heuristics (each returns a list of children, sometimes parents) ---
    cont_levy <- function(pop) {
        np <- length(pop)
        if (np == 0L) return(list(children = list(), used = integer(0)))
        n_take <- min(max(1L, as.integer(s$frac_levy * s$population)), np)
        idx <- sample.int(np, size = n_take, replace = FALSE) - 1L
        step <- .gnowee_levy(length(pop[[1]]$variables), n_take,
                              alpha = s$alpha, gam = s$gam, n = s$n)
        children <- vector("list", n_take)
        used <- integer(n_take)
        for (k in seq_len(n_take)) {
            pi_ <- idx[k] + 1L          # 1-based index into pop
            base <- pop[[pi_]]$variables
            ss <- as.numeric(step[k, ]) / s$scaling_factor
            child <- base + ss
            children[[k]] <- .gnowee_rejection_bounds(base, child, ss, lb, ub)
            # adopted_parents is 0-based in _gnowee.py (rng.choice indices);
            # population_update maps child i (1-based) -> pop[[adopted[i]+1]].
            used[k] <- idx[k]
        }
        list(children = children, used = used)
    }

    crossover <- function(pop) {
        np <- length(pop)
        n_take <- as.integer(s$frac_elite * np)
        if (n_take == 0L || np < 2L) return(list(children = list(),
                                                 used = integer(0)))
        golden <- (1 + sqrt(5)) / 2
        children <- list(); used <- integer(0)
        for (i in seq_len(n_take)) {
            r <- sample.int(np, 1L)
            attempts <- 0L
            while (r == i && attempts < 10L) {
                r <- sample.int(np, 1L); attempts <- attempts + 1L
            }
            if (r == i) next
            # crossover's adopted_parents are discarded by run_gnowee (see
            # _gnowee.py), but keep them 0-based for consistency.
            used <- c(used, i - 1L)
            base <- pop[[r]]$variables
            dx <- abs(pop[[i]]$variables - base) / golden
            children[[length(children) + 1L]] <-
                .gnowee_simple_bounds(base + dx, lb, ub)
        }
        list(children = children, used = used)
    }

    scatter_search <- function(pop) {
        np <- length(pop)
        n_take <- as.integer(s$frac_elite * np)
        if (n_take == 0L || np < 2L) return(list(children = list(),
                                                 used = integer(0)))
        children <- list(); used <- integer(0)
        # All indices follow _gnowee.py's 0-based convention; population_update
        # maps child i (1-based) -> parent pop[[used[i] + 1]].
        for (i0 in seq_len(n_take)) {
            i <- i0 - 1L
            j <- sample.int(np, 1L) - 1L
            attempts <- 0L
            while ((j == i || j %in% used) && attempts < 10L) {
                j <- sample.int(np, 1L) - 1L; attempts <- attempts + 1L
            }
            if (j == i || j %in% used) next
            used <- c(used, i)
            xi <- pop[[i + 1L]]$variables; xj <- pop[[j + 1L]]$variables
            d <- (xj - xi) / 2
            a_ <- if (i < j) 1 else -1
            beta <- (abs(j - i) - 1) / max(np - 2L, 1L)
            c1 <- xi - d * (1 + a_ * beta)
            c2 <- xi + d * (1 - a_ * beta)
            r <- stats::runif(length(xi))
            children[[length(children) + 1L]] <-
                .gnowee_simple_bounds(c1 + (c2 - c1) * r, lb, ub)
        }
        list(children = children, used = used)
    }

    mutate <- function(pop) {
        np <- length(pop)
        if (np == 0L) return(list())
        dim_ <- length(pop[[1]]$variables)
        pop_arr <- do.call(rbind, lapply(pop, function(p) p$variables))
        perm1 <- sample.int(np); perm2 <- sample.int(np)
        r <- stats::runif(1)
        k <- matrix(stats::runif(np * dim_), nrow = np, ncol = dim_) >
            (s$frac_mutation * stats::runif(1))
        diff <- pop_arr[perm1, , drop = FALSE] - pop_arr[perm2, , drop = FALSE]
        children <- pop_arr + r * diff * k
        lapply(seq_len(np), function(i)
            .gnowee_simple_bounds(children[i, ], lb, ub))
    }

    # --- population update (elitism + MH fallback + stall restart) ---------
    # Mirrors GnoweeHeuristics.population_update exactly, including the
    # timeline bookkeeping: the FIRST update of a run appends an unconditional
    # Event(generation = 1, evaluations = feval), i.e. the evaluation counter
    # restarts from that update's own cost and does NOT include the initial
    # population.  Later updates either append an improvement event or add
    # their feval to the last event.
    population_update <- function(pop, children, adopted = NULL,
                                  mh_frac = 0, random_parents = FALSE) {
        n_parents <- length(pop); n_children <- length(children)
        if (n_children == 0L) return(list(pop = pop, feval = 0L))
        feval <- 0L
        for (i in seq_len(n_children)) {
            child_vars <- as.numeric(children[[i]])
            fnew <- objective(child_vars); feval <- feval + 1L
            j <- if (random_parents) {
                sample.int(n_parents, 1L)
            } else if (!is.null(adopted) && length(adopted) == n_children) {
                adopted[i] + 1L
            } else if (i <= n_parents) {
                i
            } else {
                sample.int(n_parents, 1L)
            }
            if (fnew < pop[[j]]$fitness) {
                pop[[j]]$fitness <- fnew
                pop[[j]]$variables <- child_vars
                pop[[j]]$change_count <- pop[[j]]$change_count + 1L
                pop[[j]]$stall_count <- 0L
                if (pop[[j]]$change_count >= 25L &&
                    (j - 1L) >= as.integer(s$population * s$frac_elite)) {
                    nv <- .gnowee_initialize(1L, "random", lb, ub)[1L, ]
                    pop[[j]]$variables <- as.numeric(nv)
                    pop[[j]]$fitness <- objective(nv); feval <- feval + 1L
                    pop[[j]]$change_count <- 0L
                }
            } else {
                pop[[j]]$stall_count <- pop[[j]]$stall_count + 1L
                if (pop[[j]]$stall_count > 50000L && j != 1L) {
                    nv <- .gnowee_initialize(1L, "random", lb, ub)[1L, ]
                    pop[[j]]$variables <- as.numeric(nv)
                    pop[[j]]$fitness <- objective(nv); feval <- feval + 1L
                    pop[[j]]$change_count <- 0L; pop[[j]]$stall_count <- 0L
                }
                if (mh_frac > 0 && stats::runif(1) < mh_frac) {
                    r2 <- sample.int(n_parents, 1L)
                    if (fnew < pop[[r2]]$fitness) {
                        pop[[r2]]$fitness <- fnew
                        pop[[r2]]$variables <- child_vars
                        pop[[r2]]$change_count <- pop[[r2]]$change_count + 1L
                        pop[[r2]]$stall_count <- pop[[r2]]$stall_count + 1L
                    }
                }
            }
        }
        ord <- order(vapply(pop, function(p) p$fitness, numeric(1)))
        pop <- pop[ord]
        # Timeline: identical to _gnowee.py.
        if (length(timeline) < 2L) {
            timeline[[length(timeline) + 1L]] <<- list(
                generation = 1L, evaluations = feval,
                fitness = pop[[1]]$fitness, design = pop[[1]]$variables)
        } else {
            tl <- timeline[[length(timeline)]]
            if (pop[[1]]$fitness < tl$fitness &&
                abs((tl$fitness - pop[[1]]$fitness) /
                    max(abs(pop[[1]]$fitness), 1e-300)) > s$conv_tol) {
                timeline[[length(timeline) + 1L]] <<- list(
                    generation = tl$generation,
                    evaluations = tl$evaluations + feval,
                    fitness = pop[[1]]$fitness,
                    design = pop[[1]]$variables)
            } else {
                # generation is bumped once per loop by the caller
                timeline[[length(timeline)]]$evaluations <<- tl$evaluations +
                    feval
            }
        }
        list(pop = pop, feval = feval)
    }

    converge <- FALSE
    while (!converge) {
        lv <- cont_levy(pop)
        if (length(lv$children)) {
            res <- population_update(pop, lv$children, adopted = lv$used,
                                     mh_frac = 0.2, random_parents = TRUE)
            pop <- res$pop
        }
        cr <- crossover(pop)
        if (length(cr$children)) {
            # Python passes no adopted_parents for crossover
            res <- population_update(pop, cr$children)
            pop <- res$pop
        }
        sc <- scatter_search(pop)
        if (length(sc$children)) {
            res <- population_update(pop, sc$children, adopted = sc$used)
            pop <- res$pop
        }
        mu <- mutate(pop)
        if (length(mu)) {
            res <- population_update(pop, mu)
            pop <- res$pop
        }

        tl <- timeline[[length(timeline)]]
        gen <- tl$generation + 1L
        evals <- tl$evaluations

        if (evals > s$stall_limit && length(timeline) >= 2L) {
            if (evals > timeline[[length(timeline) - 1L]]$evaluations +
                s$stall_limit) converge <- TRUE
        }
        if (gen > s$max_gens) converge <- TRUE
        if (evals > s$max_fevals) converge <- TRUE
        if (s$optimum == 0) {
            if (pop[[1]]$fitness < s$opt_conv_tol) converge <- TRUE
        } else if (abs((pop[[1]]$fitness - s$optimum) / s$optimum) <=
                   s$opt_conv_tol) {
            converge <- TRUE
        } else if (pop[[1]]$fitness < s$optimum) {
            converge <- TRUE
        }
        timeline[[length(timeline)]]$generation <- gen
    }

    list(best_y = pop[[1]]$variables,
         best_f = pop[[1]]$fitness,
         generations = timeline[[length(timeline)]]$generation,
         evaluations = timeline[[length(timeline)]]$evaluations)
}

# ---- public solver --------------------------------------------------------
#' Gnowee-based unfolding method for neutron spectrum reconstruction
#'
#' R port of \code{bssunfold/src/bssunfold/core/unfold_gnowee.py} (and its
#' engine \code{core/_gnowee.py}).  Gnowee is a hybrid metaheuristic
#' optimiser (Bevins & Parsons, UC Berkeley / Slaybaugh Lab) that combines
#' Levy flights (Cuckoo Search), golden-ratio crossover, scatter search and
#' DE-style mutation in an elitist population with Metropolis-Hastings
#' acceptance and stall-driven restarts.
#'
#' The unfolding problem is posed in log space \eqn{y = \log(x)}{y = log(x)}
#' so positivity is automatic.  The population is seeded with a Landweber
#' warm start (or the user \code{initial_spectrum}) and the search is bounded
#' to \code{log(seed) +/- half_range} decades.  The objective is the
#' scale-consistent sum of a relative L2 residual, a Tikhonov term, a
#' second-difference smoothness term and (optionally) a negative Shannon
#' entropy.  The port is faithful to \code{_gnowee.py}: it runs the same
#' heuristics and population update and returns the best point found.
#'
#' @param A Numeric response matrix (m x n).
#' @param b Numeric measurement vector (length m).
#' @param x0 Optional initial spectrum (length n).  If \code{NULL} or all-zero
#'   a Landweber warm start is used to seed the population.  Default
#'   \code{NULL}.
#' @param population Positive integer; population size.  Default 25.
#' @param max_gens Positive integer; generation cap.  Default 200.
#' @param max_fevals Positive integer; fitness-evaluation cap.  Default 5000.
#' @param stall_limit Positive integer; evaluations without a timeline
#'   improvement before termination.  Default 200.
#' @param conv_tol Positive numeric; relative improvement required to extend
#'   the timeline.  Default 1e-6.
#' @param opt_conv_tol Positive numeric; absolute fitness-convergence
#'   tolerance.  Default 1e-2.
#' @param frac_elite Numeric in 0-to-1 range; elite fraction.  Default 0.2.
#' @param frac_levy Numeric in 0-to-1 range; Levy-flight fraction.  Default 1.0.
#' @param frac_mutation Numeric in 0-to-1 range; mutation discovery probability.
#'   Default 0.2.
#' @param alpha_levy Numeric; Levy exponent in (0.3, 1.99).  Default 1.5.
#' @param gamma_levy Numeric; Levy scale.  Default 1.0.
#' @param n_levy Positive integer; independent Levy samples.  Default 1.
#' @param scaling_factor Numeric; Levy step-length divisor.  Default 10.
#' @param init_sampling Character; \code{"lhc"} (default) or \code{"random"}.
#' @param regularization Numeric; Tikhonov weight.  Default 1e-2.
#' @param norm Integer 1 or 2; regularisation norm.  Default 2.
#' @param smoothness_order Integer 0, 1 or 2; derivative penalty order.
#'   Default 2.
#' @param smoothness_weight Numeric; smoothness-term weight.  Default 1.0.
#' @param entropy_weight Numeric; negative-entropy weight (0 disables).
#'   Default 0.
#' @param half_range Numeric; half-width of the log-space box in decades.
#'   Default 2.
#' @param random_state Optional integer RNG seed.
#' @param verbose Logical; print progress.  Default \code{FALSE}.
#' @return A list \code{list(spectrum, iterations, converged, diagnostics)}.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05,
#'               0.10, 0.8, 0.10,
#'               0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' r <- solve_gnowee(A, b, rep(0, 3), max_gens = 5L, max_fevals = 200L)
#' length(r$spectrum)
solve_gnowee <- function(A, b, x0 = NULL,
                         population = 25L, max_gens = 200L,
                         max_fevals = 5000L, stall_limit = 200L,
                         conv_tol = 1e-6, opt_conv_tol = 1e-2,
                         frac_elite = 0.2, frac_levy = 1.0,
                         frac_mutation = 0.2, alpha_levy = 1.5,
                         gamma_levy = 1.0, n_levy = 1L,
                         scaling_factor = 10.0, init_sampling = "lhc",
                         regularization = 1e-2, norm = 2L,
                         smoothness_order = 2L, smoothness_weight = 1.0,
                         entropy_weight = 0.0, half_range = 2.0,
                         random_state = NULL, verbose = FALSE) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    n <- ncol(A)
    norm_ <- as.integer(norm)
    smoothness_order <- as.integer(smoothness_order)
    if (!(norm_ %in% c(1L, 2L)))
        stop("Unsupported norm type: ", norm, ". Use 1 or 2.")
    if (!(smoothness_order %in% c(0L, 1L, 2L)))
        stop("Unsupported smoothness order: ", smoothness_order,
             ". Use 0, 1 or 2.")
    isam <- tolower(as.character(init_sampling))
    if (!(isam %in% c("lhc", "lhs", "random")))
        stop("Unsupported init_sampling: ", init_sampling,
             ". Use 'lhc' or 'random'.")
    isam <- if (isam %in% c("lhc", "lhs")) "lhc" else "random"
    if (!is.null(random_state)) set.seed(as.integer(random_state))

    L <- if (smoothness_order %in% c(1L, 2L))
        as.matrix(create_derivative_matrix(n, smoothness_order)) else NULL
    objective <- .gnowee_build_fitness(A, b, regularization, norm_, L,
                                       smoothness_weight, entropy_weight)

    # _build_seed: user initial if it has any positive mass, else Landweber.
    seed_from_x0 <- function(xv) {
        xv <- as.numeric(xv)
        pmax(xv, 1e-12)
    }
    if (!is.null(x0) && any(as.numeric(x0) > 0)) {
        seed <- seed_from_x0(x0)
    } else {
        lw <- solve_landweber(A, b, rep(0, n), max_iterations = 500L)
        seed <- pmax(as.numeric(lw$spectrum), 1e-12)
    }
    y0 <- log(pmax(as.numeric(seed), 1e-300))
    span <- half_range * log(10)
    lb <- y0 - span
    ub <- y0 + span

    s <- list(population = as.integer(population),
              init_sampling = isam,
              frac_mutation = as.numeric(frac_mutation),
              frac_elite = as.numeric(frac_elite),
              frac_levy = as.numeric(frac_levy),
              alpha = as.numeric(alpha_levy),
              gam = as.numeric(gamma_levy),
              n = as.integer(n_levy),
              scaling_factor = as.numeric(scaling_factor),
              max_gens = as.integer(max_gens),
              max_fevals = as.integer(max_fevals),
              conv_tol = as.numeric(conv_tol),
              stall_limit = as.integer(stall_limit),
              opt_conv_tol = as.numeric(opt_conv_tol),
              optimum = 0,
              penalty = 0,
              verbose = isTRUE(verbose))

    out <- .gnowee_run(lb = lb, ub = ub, objective = objective, s = s,
                       seed_solution = y0)
    spectrum <- pmax(exp(out$best_y), 0)
    n_evals <- out$evaluations
    hit_feval_cap <- n_evals >= s$max_fevals
    hit_gen_cap <- out$generations >= s$max_gens
    converged <- !(hit_feval_cap || hit_gen_cap)
    diagnostics <- list(best_fitness = out$best_f,
                        generations = out$generations,
                        evaluations = n_evals)
    list(spectrum = as.numeric(spectrum), iterations = as.integer(n_evals),
         converged = converged, diagnostics = diagnostics)
}

#' Wrapper around \code{\link{solve_gnowee}} for the unified workflow.
#'
#' @inheritParams run_unfolding
#' @inheritParams solve_gnowee
#' @return A result list as produced by \code{\link{run_unfolding}}.
#' @export
unfold_gnowee <- function(detector_names, n_energy_bins, E_MeV,
                          sensitivities, cc_icrp116, save_result_callback,
                          readings,
                          initial_spectrum = NULL,
                          population = 25L, max_gens = 200L,
                          max_fevals = 5000L, stall_limit = 200L,
                          conv_tol = 1e-6, opt_conv_tol = 1e-2,
                          frac_elite = 0.2, frac_levy = 1.0,
                          frac_mutation = 0.2, alpha_levy = 1.5,
                          gamma_levy = 1.0, n_levy = 1L,
                          scaling_factor = 10.0, init_sampling = "lhc",
                          regularization = 1e-2, norm = 2L,
                          smoothness_order = 2L, smoothness_weight = 1.0,
                          entropy_weight = 0.0, half_range = 2.0,
                          calculate_errors = FALSE,
                          noise_level = 0.01,
                          n_montecarlo = 100L,
                          save_result = FALSE, random_state = NULL,
                          verbose = FALSE,
                          max_neutron_energy = NULL) {
    x0_default <- rep(0.0, n_energy_bins)
    run_unfolding(
        detector_names = detector_names, n_energy_bins = n_energy_bins,
        E_MeV = E_MeV, sensitivities = sensitivities,
        cc_icrp116 = cc_icrp116, save_result_callback = save_result_callback,
        readings = readings, initial_spectrum = initial_spectrum,
        default_initial = x0_default,
        solve_func = make_solve_wrapper(solve_gnowee,
                                         population = population,
                                         max_gens = max_gens,
                                         max_fevals = max_fevals,
                                         stall_limit = stall_limit,
                                         conv_tol = conv_tol,
                                         opt_conv_tol = opt_conv_tol,
                                         frac_elite = frac_elite,
                                         frac_levy = frac_levy,
                                         frac_mutation = frac_mutation,
                                         alpha_levy = alpha_levy,
                                         gamma_levy = gamma_levy,
                                         n_levy = n_levy,
                                         scaling_factor = scaling_factor,
                                         init_sampling = init_sampling,
                                         regularization = regularization,
                                         norm = norm,
                                         smoothness_order = smoothness_order,
                                         smoothness_weight = smoothness_weight,
                                         entropy_weight = entropy_weight,
                                         half_range = half_range,
                                         random_state = random_state,
                                         verbose = verbose),
        solve_kwargs = list(),
        method_name = "Gnowee",
        extra_output = list(population = population,
                            max_gens = max_gens,
                            max_fevals = max_fevals,
                            regularization = regularization,
                            norm = norm,
                            smoothness_order = smoothness_order),
        calculate_errors = calculate_errors,
        noise_level = noise_level,
        n_montecarlo = n_montecarlo,
        random_state = random_state,
        save_result = save_result,
        max_neutron_energy = max_neutron_energy)
}
