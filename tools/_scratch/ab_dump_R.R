# Dump R wrapper spectra for all 4 parametric optimizers, per fixture/case,
# to JSON so a single comparison pass can compute rel-L2 vs Python.
suppressWarnings(suppressMessages(pkgload::load_all(".", quiet = TRUE)))
library(jsonlite)
dirs <- c(main="/tmp/bsscmp", gsf="/tmp/alt_gsf",
          coarse="/tmp/alt_coarse", jinr="/tmp/alt_jinr")
opts <- c("lmfit","cvxpy","qpsolvers","combined")
out <- list()
for (rf in names(dirs)) {
  fx <- fromJSON(file.path(dirs[[rf]], "fixture.json"), simplifyVector = FALSE)
  E <- unlist(fx$E_MeV); nms <- unlist(fx$detector_names)
  raw <- setNames(lapply(seq_along(nms), function(i) unlist(fx$sensitivities[[i]])), nms)
  det <- Detector$new(response_function=c(list(E_MeV=E), raw), detector_names=nms)
  for (case in names(fx$cases)) {
    rd <- setNames(as.numeric(unlist(fx$cases[[case]]$readings)), nms)
    for (opt in opts) {
      r <- tryCatch(det$unfold_parametric(rd, optimizer=opt),
                    error=function(e) list(err_=conditionMessage(e)))
      key <- paste(rf, case, opt, sep="/")
      out[[key]] <- if (!is.null(r$err_)) list(err=r$err_)
                    else list(spectrum=r$spectrum, iterations=r$iterations,
                              converged=r$converged)
    }
  }
}
write_json(out, "/tmp/ab_R.json", digits=17, auto_unbox=TRUE)
cat("R dump keys:", length(out), "\n")
