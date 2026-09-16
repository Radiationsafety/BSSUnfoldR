#' Benchmark unfolding methods on synthetic spectra
#'
#' R port of
#' \code{bssunload/utils/comparison.py:benchmark_unfold_methods}.
#' Generates a family of synthetic spectra (a 1/E component plus lognormal
#' bumps plus a flat component), forward-projects them onto the supplied
#' detector set, unfolds with each requested method and scores the results
#' against the ground truth. Recognized method names mirror the registry
#' of \code{\link{solve_combined}}.
#'
#' @name benchmark-methods
NULL

# Synthetic mixed spectrum: thermal-spot / epithermal-1/E / flat mixture
.make_probe_spec <- function(E_MeV, amp_1e, amp_bump, center, width,
                             flat) {
    E <- as.numeric(E_MeV)
    lx <- log10(E)
    spec <- amp_1e * (0.5 / E) +
        amp_bump * exp(-((lx - center)^2) / width^2) + flat
    pmax(spec, 1e-12)
}

.build_probe_readings <- function(detector, n_probes, seed) {
    n_probes <- max(as.integer(n_probes), 1L)
    E <- as.numeric(detector$E_MeV)
    n_bins <- length(E)
    set.seed(as.integer(seed))
    centers <- seq(-3.5, 2.5, length.out = n_probes)
    out <- vector("list", n_probes)
    probes <- matrix(0, nrow = n_probes, ncol = n_bins)
    lx <- log10(E)
    for (i in seq_len(n_probes)) {
        spec <- (0.35 + 0.35 * runif(1)) / E +
            exp(-((lx - centers[i])^2) / max(0.5, runif(1))) +
            rep(runif(1) * 0.1, n_bins)
        spec <- spec / sum(spec)
        probes[i, ] <- spec
        out[[i]] <- sapply(detector$detector_names, function(n)
            sum(detector$sensitivities[[n]] * spec))
    }
    list(readings = out, probes = probes)
}

#' Benchmark selected unfolding methods on synthetic probes
#'
#' @param detector Detector R6 object (uses \code{sensitivities},
#'   \code{detector_names}, \code{E_MeV}).
#' @param methods Named list where each entry is the parameter list for
#'   the solver \code{method} given as the list name.
#' @param n_probes Integer number of synthetic spectra. Default 4.
#' @param seed Integer RNG seed. Default 1234.
#' @return A list \code{list(results)} with per-method entries
#'   \code{method}, \code{mean_residual}, \code{mean_integral_diff},
#'   \code{mean_kl}, \code{per_probe}.
#' @export
#' @examples
#' \donttest{
#' det <- Detector$new(response_function = "PTB",
#'                     detector_names = c("0in", "3in", "5in", "8in", "12in"))
#' bm <- benchmark_unfold_methods(det,
#'     methods = list(landweber = list(max_iterations = 20L)), n_probes = 2L)
#' bm$results[[1]]$mean_residual}
benchmark_unfold_methods <- function(detector,
                                     methods = list(
                                         landweber = list(max_iterations =
                                                              30L)),
                                     n_probes = 4L, seed = 1234L) {
    n_probes <- max(as.integer(n_probes), 1L)
    n_bins <- length(as.numeric(detector$E_MeV))
    pb <- .build_probe_readings(detector, n_probes, seed)
    registry <- .combined_solver_registry()
    results <- lapply(names(methods), function(mname) {
        if (!mname %in% names(registry)) {
            stop("Unknown benchmark method: ", mname,
                 ". See solve_combined for the registry of method names.")
        }
        solver <- get(registry[[mname]], mode = "function")
        params <- methods[[mname]]
        sys <- .build_system(pb$readings[[1]], detector$detector_names,
                             detector$sensitivities)
        res_list <- vector("list", n_probes)
        for (i in seq_len(n_probes)) {
            sys_i <- .build_system(pb$readings[[i]], detector$detector_names,
                                   detector$sensitivities)
            r <- tryCatch(do.call(solver, c(list(A = sys_i$A,
                                                b = sys_i$b),
                                             params,
                                             list(x0 = rep(1, n_bins)))),
                          error = function(e) conditionMessage(e))
            spec <- if (is.list(r)) as.numeric(r$spectrum) else NULL
            if (is.null(spec)) {
                res_list[[i]] <- NULL
                next
            }
            truth <- pb$probes[i, ]
            spec <- pmax(spec, 0)
            res_list[[i]] <- list(
                residual_norm =
                    sqrt(sum((as.numeric(sys_i$A %*% spec) - sys_i$b)^2)),
                relative_integral_diff = abs(sum(spec) - sum(truth)) /
                    max(sum(truth), 1e-30),
                kl_divergence = sum((spec / max(sum(spec), 1e-30)) *
                    log(pmax(spec / max(sum(spec), 1e-30) , 1e-30) /
                            pmax(truth / max(sum(truth), 1e-30), 1e-30))
                )
            )
        }
        ok <- Filter(Negate(is.null), res_list)
        # extract nested numeric statistics safely when the method has completed
        get_scalar <- function(field) {
            if (length(ok) == 0L) return(NA_real_)
            vv <- vapply(ok, function(x) x[[field]], numeric(1))
            vv[!is.finite(vv)] <- 1e6
            mean(vv)
        }
        list(method = mname,
             mean_residual = get_scalar("residual_norm"),
             mean_integral_diff = get_scalar("relative_integral_diff"),
             mean_kl = get_scalar("kl_divergence"),
             per_probe = res_list,
             method_ok = length(ok) / max(n_probes, 1L))
    })
    names(results) <- names(methods)
    list(results = results)
}
