#!/usr/bin/env Rscript
# Build a shared parity fixture from the R package's own response-function
# data, so that Python and R are fed byte-identical numbers. Any difference in
# the unfolded spectra is then attributable to the method, not to the data.
#
# Parameterised by environment variables so the harness can check the same
# method against grids / detector sets it was not tuned on:
#   PARITY_RF     response function name (PTB, GSF, LANL, JINR, ...)  [PTB]
#   PARITY_DETS   comma-separated detector names to keep              [all]
#   PARITY_EMAX   drop bins above this energy in MeV, non-uniform cut [unset]
#   PARITY_NBIN   subsample the grid to at most N bins                [unset]

suppressWarnings(suppressMessages(pkgload::load_all(getwd(), quiet = TRUE)))

library(jsonlite)

rf_name <- Sys.getenv("PARITY_RF", "PTB")
keep <- strsplit(Sys.getenv("PARITY_DETS", ""), ",", fixed = TRUE)[[1]]
keep <- trimws(keep)
keep <- keep[nzchar(keep)]
emax <- Sys.getenv("PARITY_EMAX", "")
nbin <- Sys.getenv("PARITY_NBIN", "")

det <- Detector$new(response_function = rf_name)
if (length(keep)) {
    missing <- setdiff(keep, names(det$sensitivities))
    if (length(missing)) {
        stop("unknown detectors for ", rf_name, ": ", paste(missing, collapse = ", "))
    }
    det$set_detector_names(keep)
}

E <- det$E_MeV
ln <- det$ln_steps

# Detector stores sensitivities already weighted by ln_steps; the fixture
# carries the RAW response functions and each language re-applies the
# weighting itself.
raw_sens <- lapply(det$detector_names,
                   function(n) det$sensitivities[[n]] / ln)
names(raw_sens) <- det$detector_names

# Non-uniform grid stress test: drop every 3rd bin, or cap the energy.
if (nzchar(emax)) {
    idx <- which(E <= as.numeric(emax))
    E <- E[idx]
    ln <- ln[idx]
    raw_sens <- lapply(raw_sens, function(s) s[idx])
}
if (nzchar(nbin)) {
    step <- ceiling(length(E) / as.numeric(nbin))
    idx <- seq(1L, length(E), by = step)
    E <- E[idx]
    ln <- ln[idx]
    raw_sens <- lapply(raw_sens, function(s) s[idx])
}

mk_phi <- function(kind, E) {
    switch(kind,
        moderated = 0.5 / E +
            exp(-((log10(E) + 7)^2) / 2) +
            exp(-((log10(E))^2) / 0.4),
        fission = E * exp(-E / 1.3),
        evap = exp(-sqrt(E) / 0.35),
        split = {
            p <- rep(0, length(E))
            p[seq(1L, length(E), length.out = 3L)] <- c(1, 40, 0.2)
            p
        },
        stop("unknown spectrum kind: ", kind))
}

mk_case <- function(kind, E, raw_sens, ln, det_names, target) {
    phi <- mk_phi(kind, E)
    phi <- phi / sum(phi)
    resp <- function(n) sum(raw_sens[[n]] * ln * phi)
    phi <- phi * (target / resp(det_names[1]))
    list(
        truth = unname(phi),
        readings = unname(vapply(det_names, resp, numeric(1)))
    )
}

det_names <- det$detector_names
kinds <- c("moderated", "fission", "evap", "split")
targets <- c(moderated = 1000, fission = 800, evap = 500, split = 600)
cases <- setNames(lapply(kinds, function(k)
    mk_case(k, E, raw_sens, ln, det_names, targets[[k]])), kinds)

fx <- list(
    E_MeV = E,
    detector_names = det_names,
    sensitivities = unname(raw_sens),
    cases = cases
)

outfile <- file.path(Sys.getenv("PARITY_OUT", "/tmp/bsscmp"), "fixture.json")
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
write_json(fx, outfile, digits = 17, auto_unbox = TRUE, pretty = FALSE)
cat("wrote", outfile, "rf =", rf_name, "nbins =", length(E),
    "ndet =", length(det_names), "\n")
