# Resolve why det$unfold_parametric(optimizer=) differs from calling the
# raw engine directly with the fixture A/b. (Diagnostic only.)
suppressWarnings(suppressMessages(pkgload::load_all(".", quiet = TRUE)))
fx <- jsonlite::fromJSON("/tmp/bsscmp/fixture.json", simplifyVector = FALSE)
E <- unlist(fx$E_MeV); nms <- unlist(fx$detector_names)
raw <- setNames(lapply(seq_along(nms), function(i) unlist(fx$sensitivities[[i]])), nms)
rd  <- setNames(as.numeric(unlist(fx$cases$fission$readings)), nms)
det <- Detector$new(response_function = c(list(E_MeV = E), raw), detector_names = nms)
ln <- compute_log_steps(E) * log(10)
sys <- BSSUnfoldR:::.build_system(rd, nms, raw)
A <- sys$A; b <- sys$b

# 1) raw engine, default args
e1 <- solve_parametric_cvxpy(A, b, E, ln)
# 2) raw engine, alpha=1e-4 (what the wrapper uses)
e2 <- solve_parametric_cvxpy(A, b, E, ln, alpha = 1e-4)
# 3) wrapper
w  <- det$unfold_parametric(rd, optimizer = "cvxpy")

cat(sprintf("raw default  sum=%.6f\n", sum(e1$spectrum)))
cat(sprintf("raw a=1e-4   sum=%.6f\n", sum(e2$spectrum)))
cat(sprintf("wrapper      sum=%.6f\n", sum(w$spectrum)))
cat(sprintf("rel(raw-def, raw-1e-4)=%.3e\n",
  max(abs(e1$spectrum-e2$spectrum))/(max(e2$spectrum)+1e-30)))
cat(sprintf("rel(raw-1e-4, wrap)=%.3e\n",
  max(abs(e2$spectrum-w$spectrum))/(max(w$spectrum)+1e-30)))
cat("raw params:"); print(round(e2$params,6))
cat("wrap params:"); print(round(w$params,6))
