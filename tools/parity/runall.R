#!/usr/bin/env Rscript
# Run every Detector unfold_* method on a fixture and list the ones that error.
fxdir <- if (length(commandArgs(trailingOnly = TRUE)) >= 1)
    commandArgs(trailingOnly = TRUE)[1] else "/tmp/bsscmp"
Sys.setenv(PARITY_OUT = fxdir)
suppressWarnings(suppressMessages(
    pkgload::load_all("/Users/spiralfractal/Work/code/BSSUnfoldR", quiet = TRUE)))
suppressMessages(library(jsonlite))

fx <- jsonlite::fromJSON(file.path(fxdir, "fixture.json"), simplifyVector = FALSE)
rf <- c(list(E_MeV = unlist(fx$E_MeV)),
        setNames(lapply(seq_along(fx$detector_names),
                        function(i) unlist(fx$sensitivities[[i]])),
                 unlist(fx$detector_names)))
det <- Detector$new(response_function = rf,
                    detector_names = unlist(fx$detector_names))
methods <- sort(grep("^unfold_", names(Detector$public_methods), value = TRUE))
nbin <- length(unlist(fx$E_MeV))

for (case in names(fx$cases)) {
    readings <- setNames(as.numeric(unlist(fx$cases[[case]]$readings)),
                         unlist(fx$detector_names))
    ok <- 0L; bad <- list()
    for (m in methods) {
        fn <- det[[m]]
        r <- tryCatch({
            out <- tryCatch(fn(readings = readings, random_state = 42L),
                            error = function(e) {
                                if (grepl("random_state", conditionMessage(e), fixed = TRUE))
                                    fn(readings = readings) else stop(e)
                            })
            sp <- out$spectrum
            if (is.null(sp) || !is.numeric(sp) || length(sp) != nbin)
                stop(sprintf("bad spectrum (numeric=%s len=%d want=%d)",
                             is.numeric(sp), length(sp), nbin))
            if (any(!is.finite(sp))) stop("non-finite spectrum")
            NULL
        }, error = function(e) conditionMessage(e))
        if (is.null(r)) ok <- ok + 1L
        else bad[[length(bad) + 1L]] <- c(m, gsub("\n", " | ", r))
    }
    cat(sprintf("== %s / %s : methods %d  ran %d  errored %d\n",
                basename(fxdir), case, length(methods), ok, length(bad)))
    for (b in bad) cat(sprintf("   %-30s %s\n", b[1], substr(b[2], 1, 120)))
}
