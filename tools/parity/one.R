#!/usr/bin/env Rscript
# Print one method's spectrum (and key scalars) from the R package as JSON.
# Usage: Rscript tools/parity/one.R <method> [case]
suppressWarnings(suppressMessages(pkgload::load_all(getwd(), quiet = TRUE)))
suppressMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
method <- args[1]
case <- if (length(args) >= 2) args[2] else "moderated"

fx <- fromJSON(file.path(Sys.getenv("PARITY_OUT", "/tmp/bsscmp"),
                         "fixture.json"), simplifyVector = FALSE)
E <- unlist(fx$E_MeV)
nms <- unlist(fx$detector_names)
rf <- setNames(lapply(seq_along(nms),
                      function(i) unlist(fx$sensitivities[[i]])), nms)
readings <- setNames(as.numeric(unlist(fx$cases[[case]]$readings)), nms)

det <- Detector$new(response_function = c(list(E_MeV = E), rf),
                    detector_names = nms)

extra <- list()
kvs <- args[-(1:2)]
for (kv in kvs) {
    p <- strsplit(kv, "=", fixed = TRUE)[[1]]
    extra[[p[1]]] <- suppressWarnings(as.numeric(p[2]))
    if (is.na(extra[[p[1]]])) extra[[p[1]]] <- p[2]
}

out <- tryCatch(do.call(det[[method]], c(list(readings), extra)),
                error = function(e) list(error = conditionMessage(e)))
if (!is.null(out$error)) {
    write_json(list(error = out$error), stdout(), auto_unbox = TRUE)
} else {
    write_json(list(spectrum = out$spectrum,
                    residual_norm = out$residual_norm,
                    iterations = out$iterations,
                    doserates = out$doserates),
               stdout(), digits = 17, auto_unbox = TRUE)
}
