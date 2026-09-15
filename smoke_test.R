# End-to-end smoke test for BSSUnfoldR
library(BSSUnfoldR)

cat("=== Constants ===\n")
icrp116 <- ICRP116_COEFF_EFFECTIVE_DOSE()
cat("ICRP116: length(E_MeV) =", length(icrp116$E_MeV),
    "; geometries:", paste(names(icrp116), collapse = ", "), "\n")
ptb <- RF_PTB()
cat("RF_PTB detectors:", paste(names(ptb), collapse = ", "), "\n\n")

cat("=== Detector construction ===\n")
det <- Detector$new(response_function = "PTB",
                    detector_names = c("0in", "3in", "5in", "8in", "12in"))
cat("n_energy_bins =", det$n_energy_bins, "\n")
cat("energy grid first/last:", head(det$E_MeV, 2), " ... ", tail(det$E_MeV, 2), "\n")
cat("sensitivities available:", paste(names(det$sensitivities), collapse = ", "), "\n\n")

cat("=== Synthetic ground-truth spectrum ===\n")
# Use a simple Watt fission spectrum shape on the energy grid of the detector
E <- det$E_MeV
# A crude 1/E shape with a thermal peak around 1e-7 MeV and a fission bump around 1 MeV
phi_true <- 0.5 / E + exp(-((log10(E) + 7)^2) / 2) + exp(-((log10(E))^2) / 0.4)
phi_true <- phi_true / sum(phi_true)
# Scale to count rate of "0in" = 1000 counts/s
target_0in <- 1000
scale <- target_0in / sum(det$sensitivities[["0in"]] * phi_true)
phi_true <- phi_true * scale
cat("Spectrum integral:", sum(phi_true), "\n")

cat("\n=== Synthesize readings from 5 spheres ===\n")
readings <- sapply(det$detector_names, function(n) {
    sum(det$sensitivities[[n]] * phi_true)
})
print(readings)

cat("\n=== Test each unfolding method ===\n")
detectors_used <- det$detector_names
E_MeV <- det$E_MeV
sens <- det$sensitivities
cc116 <- det$cc_icrp116
cb <- function(out) det$save_result(out)

methods_to_test <- list(
    list(name = "MLEM", fn = unfold_mlem, args = list(max_iterations = 200L)),
    list(name = "Landweber", fn = unfold_landweber, args = list(max_iterations = 500L)),
    list(name = "GRAVEL", fn = unfold_gravel, args = list(max_iterations = 200L)),
    list(name = "SAND-II", fn = unfold_sandii, args = list(max_iterations = 100L)),
    list(name = "BUNKI", fn = unfold_bunki, args = list(max_iterations = 200L)),
    list(name = "MAXED", fn = unfold_maxed, args = list(max_iterations = 200L)),
    list(name = "FERDOR", fn = unfold_ferdor, args = list(max_iterations = 50L)),
    list(name = "STAY'SL", fn = unfold_staysl, args = list()),
    list(name = "CGLS", fn = unfold_cgls, args = list(max_iterations = 50L)),
    list(name = "TSVD", fn = unfold_tsvd, args = list(method = "gcv"))
)

results <- list()
for (m in methods_to_test) {
    cat(sprintf("\n--- %s ---\n", m$name))
    args <- c(list(
        detector_names = detectors_used, n_energy_bins = det$n_energy_bins,
        E_MeV = E_MeV, sensitivities = sens, cc_icrp116 = cc116,
        save_result_callback = cb, readings = readings
    ), m$args)
    r <- do.call(m$fn, args)
    cat("method =", r$method, "\n")
    cat("residual_norm =", format(r$residual_norm, digits = 4), "\n")
    cat("spectrum length:", length(r$spectrum),
        "; integral:", format(sum(r$spectrum), digits = 4), "\n")
    if (!is.null(r$iterations)) cat("iterations:", r$iterations,
                                    "; converged:", r$converged, "\n")
    cat("dose rates AP/ISO:", format(r$doserates["AP"], digits = 4),
        "/", format(r$doserates["ISO"], digits = 4), " pSv/s\n")
    results[[m$name]] <- r
}

cat("\n=== Compare to ground truth ===\n")
# Simple comparison: ratio of integrals + cosine similarity
phi_true_norm <- phi_true / sqrt(sum(phi_true^2))
for (n in names(results)) {
    phi_est <- results[[n]]$spectrum
    phi_est_norm <- phi_est / sqrt(sum(phi_est^2))
    cos_sim <- sum(phi_true_norm * phi_est_norm)
    cat(sprintf("%-12s cosine_sim = %.4f  integral_ratio = %.4f\n",
                n, cos_sim, sum(phi_est) / sum(phi_true)))
}

cat("\n=== Detector R6 convenience wrappers ===\n")
det$history <- list()
r1 <- det$unfold_mlem(readings, max_iterations = 100L, save_result = TRUE)
r2 <- det$unfold_gravel(readings, max_iterations = 100L, save_result = TRUE)
r3 <- det$unfold_maxed(readings, max_iterations = 100L, save_result = TRUE)
cat("history length after 3 saves:", length(det$history), "\n")
cat("history methods:", sapply(det$history, function(h) h$method), "\n")

cat("\nAll tests completed successfully.\n")
