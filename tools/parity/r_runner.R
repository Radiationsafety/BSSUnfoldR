#!/usr/bin/env Rscript
# Mirror of tools/parity/py_runner.py: run every R Detector unfold_* method on
# the shared fixture and write results-per-method JSON for the comparer.

pkg_root <- function() {
    a <- grep("^--root=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(a)) sub("^--root=", "", a) else getwd()
}

suppressWarnings(suppressMessages(pkgload::load_all(pkg_root(), quiet = TRUE)))
suppressMessages(library(jsonlite))

out_dir <- if (length(commandArgs(trailingOnly = TRUE)) >= 1)
    commandArgs(trailingOnly = TRUE)[1] else
    file.path(Sys.getenv("PARITY_OUT", "/tmp/bsscmp"), "r")
seed <- as.integer(Sys.getenv("PARITY_SEED", "42"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fx_path <- file.path(Sys.getenv("PARITY_OUT", "/tmp/bsscmp"), "fixture.json")
fx <- jsonlite::fromJSON(fx_path, simplifyVector = FALSE)

rf <- c(list(E_MeV = unlist(fx$E_MeV)),
        setNames(lapply(seq_along(fx$detector_names),
                        function(i) unlist(fx$sensitivities[[i]])),
                 unlist(fx$detector_names)))

det <- Detector$new(response_function = rf,
                    detector_names = unlist(fx$detector_names))

methods <- sort(grep("^unfold_", names(Detector$public_methods), value = TRUE))

numeric_keys <- c("residual_norm", "iterations", "converged", "chi2",
                  "reduced_chi2", "n_iterations")

run_case <- function(case_name) {
    readings <- setNames(as.numeric(unlist(fx$cases[[case_name]]$readings)),
                         unlist(fx$detector_names))
    res <- list()
    for (m in methods) {
        fn <- det[[m]]
        rec <- list()
        base_kwargs <- list(readings = readings)
        try_with <- function(kwargs)
            tryCatch(do.call(fn, kwargs), error = function(e) e)
        out <- try_with(c(base_kwargs, list(random_state = seed)))
        if (inherits(out, "error") &&
            grepl("random_state", conditionMessage(out), fixed = TRUE)) {
            out <- try_with(base_kwargs)
        }
        if (inherits(out, "error")) {
            rec$error <- conditionMessage(out)
            out <- NULL
        }
        if (!is.null(out)) {
            spec <- out$spectrum
            if (is.null(spec) || !is.numeric(spec) || length(spec) == 0) {
                rec$error <- "no numeric spectrum"
            } else if (any(!is.finite(spec))) {
                rec$error <- "non-finite spectrum"
                rec$spectrum <- as.numeric(spec)
            } else {
                rec$spectrum <- as.numeric(spec)
                for (k in numeric_keys)
                    if (!is.null(out[[k]]) && length(out[[k]]) == 1)
                        rec[[k]] <- out[[k]]
                if (!is.null(out$doserates))
                    rec$doserates <- as.list(out$doserates)
            }
        }
        res[[m]] <- rec
    }
    res
}

for (case in names(fx$cases)) {
    r <- run_case(case)
    jsonlite::write_json(r, file.path(out_dir, paste0(case, ".json")),
                         digits = 17, auto_unbox = TRUE, pretty = FALSE)
    ok <- sum(vapply(r, function(v) !is.null(v$spectrum), logical(1)))
    cat(case, ":", ok, "/", length(methods), "methods produced a spectrum\n",
        file = stderr())
}
