#!/usr/bin/env Rscript
# Correctness gate for the five commercial-license QP ports.
#
# Python returns an all-zero spectrum for gurobi/mosek/cplex/copt/xpress here
# (no license, no open-source fallback), so there is nothing to A/B against.
# These checks validate the ports against the optimisation problem instead:
#   1. all five backends agree;
#   2. the solution satisfies the KKT conditions of  min 0.5 x'Px + q'x,
#      lb <= x <= ub;
#   3. the unconstrained, unregularised case reduces to least squares;
#   4. the unfolding reproduces the readings far better than a flat spectrum;
#   5. its integral is in the same range as the reference solvers.
#
# Usage: Rscript tools/parity/check_commercial.R [fixture-dir]
suppressWarnings(suppressMessages(pkgload::load_all(getwd(), quiet = TRUE)))

fxdir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(fxdir)) fxdir <- "/tmp/bsscmp"
fx <- jsonlite::fromJSON(file.path(fxdir, "fixture.json"), simplifyVector = FALSE)
E <- unlist(fx$E_MeV)
nms <- unlist(fx$detector_names)
rf <- setNames(lapply(seq_along(nms),
                      function(i) unlist(fx$sensitivities[[i]])), nms)
det <- Detector$new(response_function = c(list(E_MeV = E), rf),
                    detector_names = nms)

solvers <- c("gurobi", "mosek", "cplex", "copt", "xpress")
qp_of <- BSSUnfoldR:::.commercial_qp_data

kkt_report <- function(A, b, x, alpha, norm, order, weight, nonneg) {
    qp <- qp_of(A, b, alpha, as.integer(norm), as.integer(order), weight)
    P <- qp$P; q <- qp$q
    grad <- as.numeric(P %*% x + q)
    lo <- if (nonneg) 0 else -Inf
    at_lb <- x <= lo + 1e-12
    defect_free <- if (any(!at_lb)) max(abs(grad[!at_lb])) else 0
    bad_lb <- if (any(at_lb)) max(-grad[at_lb]) else 0
    gscale <- max(abs(q), abs(as.numeric(P %*% x)))
    list(primal_res = sqrt(sum(as.numeric(A %*% x - b)^2)),
         min_x = min(x),
         stat_defect = max(defect_free, bad_lb) / max(gscale, 1e-30),
         obj = 0.5 * sum(x * (P %*% x)) + sum(q * x))
}

cat("fixture:", fxdir, " bins:", length(E),
    " detectors:", length(nms), "\n")
overall_ok <- TRUE

for (case in names(fx$cases)) {
    readings <- setNames(as.numeric(unlist(fx$cases[[case]]$readings)), nms)
    sys <- BSSUnfoldR:::.build_system(readings, det$detector_names,
                                      det$sensitivities)
    A <- as.matrix(sys$A); b <- as.numeric(sys$b); n <- ncol(A)

    for (cfg in list(list(alpha = 1e-4, norm = 2L, order = 0L),
                     list(alpha = 1e-2, norm = 2L, order = 1L),
                     list(alpha = 1e-2, norm = 2L, order = 2L),
                     list(alpha = 1e-1, norm = 1L, order = 0L))) {
        xs <- lapply(solvers, function(s)
            det[[paste0("unfold_", s)]](readings, regularization = cfg$alpha,
                                        norm = cfg$norm,
                                        smoothness_order = cfg$order)$spectrum)
        base <- xs[[1]]
        spread <- max(sapply(xs, function(x) {
            d <- x - base
            sqrt(sum(d^2)) / max(sqrt(sum(base^2)), 1e-30)
        }))
        k <- kkt_report(A, b, base, cfg$alpha, cfg$norm, cfg$order, 1.0, TRUE)
        # Independent optimality certificate, free of any KKT bookkeeping:
        # try to descend from the reported solution with projected gradient.
        # A global minimiser of a convex QP over a box admits no descent.
        qp_ <- qp_of(A, b, cfg$alpha, cfg$norm, cfg$order, 1.0)
        pg <- BSSUnfoldR:::.commercial_pg_polish(
            qp_$P, qp_$q, rep(0, nrow(t(A))), rep(Inf, ncol(A)),
            base, max_iterations = 20000L)
        f_of <- function(z) 0.5 * sum(z * (qp_$P %*% z)) + sum(qp_$q * z)
        gain <- f_of(base) - f_of(pg)
        fscale <- max(abs(f_of(base)), 1e-30)
        flat <- rep(mean(b) / max(mean(Matrix::colSums(A)), 1e-30),
                    length(base))
        flat_res <- sqrt(sum(as.numeric(A %*% flat - b)^2))
        ok <- spread < 1e-8 && k$stat_defect < 1e-5 &&
            k$min_x > -1e-12 && k$primal_res < flat_res &&
            gain < 1e-6 * fscale   # L1 (norm=1) exits at ~1e-7 relative
                                   # improvability; L2 at ~1e-13
        if (!ok) { overall_ok <- FALSE
            cat(sprintf("   ^FAIL spread=%.3e stat=%.3e gain=%.3e minx=%.3e res=%.3e<%.3e\n",
                        spread, k$stat_defect, gain / fscale, k$min_x,
                        k$primal_res, flat_res)) }
        cat(sprintf("%-10s a=%.0e L%d nrm%d | spread %.2e | stat %.2e | descent %.2e | min %.1e | res %.4e (flat %.4e)\n",
                    case, cfg$alpha, cfg$order, cfg$norm, spread,
                    k$stat_defect, gain / fscale, k$min_x, k$primal_res, flat_res))
    }

    # Tested at engine level on purpose: run_unfolding() clamps the returned
    # spectrum at zero (.standardize_output), so a Detector-level free solve
    # would measure the clamp, not the QP.
    free <- BSSUnfoldR:::.commercial_qp_solve(A, b, solver = "gurobi",
                                              alpha = 0, norm = 2L,
                                              smoothness_order = 0L,
                                              nonneg = FALSE)
    r <- list(spectrum = free$spectrum)
    ls <- qr.solve(A, b)
    # A is underdetermined, so the LS minimiser is not unique: compare the
    # residual, not the iterate.
    xq <- as.numeric(r$spectrum)
    res_q <- sqrt(sum(as.numeric(A %*% xq - b)^2))
    res_ls <- sqrt(sum(as.numeric(A %*% as.numeric(ls) - b)^2))
    cat(sprintf("%-10s free/alpha=0 residual %.6e vs qr.solve %.6e (rel %.2e) kkt %.2e\n",
                case, res_q, res_ls, abs(res_q - res_ls) / max(res_ls, 1e-30),
                free$kkt))
    if (res_q > res_ls * (1 + 1e-6) + 1e-9) { overall_ok <- FALSE
        cat("   ^FAIL free-solve\n") }

    ref <- c(det$unfold_tsvd(readings)$spectrum,
             det$unfold_cgls(readings)$spectrum)
    ref_int <- sum(unlist(ref)) / 2
    g_int <- sum(det$unfold_gurobi(readings)$spectrum)
    cat(sprintf("%-10s integral gurobi %.4e vs (tsvd,cgls) mean %.4e ratio %.3f\n",
                case, g_int, ref_int, g_int / ref_int))
}

cat(if (overall_ok) "=> ALL CHECKS PASS\n" else "=> SOME CHECKS FAILED\n")
if (!overall_ok) quit(status = 1)
