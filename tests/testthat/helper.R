# Shared fixtures used across the BSSUnfoldR test suite.
#
# A small (3 detector x 3 energy bin) synthetic system is reused by most
# tests so that they run in well under a second on CI. The PTB-based
# 60-bin fixtures are kept in their own helper files for the slow tests.

# ---- tiny 3x3 system ----

tiny_A <- matrix(c(0.90, 0.05, 0.05,
                   0.10, 0.80, 0.10,
                   0.30, 0.30, 0.40), nrow = 3, byrow = TRUE)
tiny_b <- c(1.0, 0.6, 0.4)
tiny_E  <- c(1e-9, 1e-6, 1e-3)
tiny_x0_default <- rep(0.5, 3)

# ---- PTB-based system for integration tests ----

make_ptb_detector <- function() {
    Detector$new(
        response_function = "PTB",
        detector_names = c("0in", "3in", "5in", "8in", "12in")
    )
}

make_ptb_readings <- function(det = make_ptb_detector(),
                              target_0in = 1000) {
    E <- det$E_MeV
    phi <- 0.5 / E +
        exp(-((log10(E) + 7)^2) / 2) +
        exp(-((log10(E))^2) / 0.4)
    phi <- phi / sum(phi)
    scale <- target_0in / sum(det$sensitivities[["0in"]] * phi)
    phi <- phi * scale
    sapply(det$detector_names, function(n)
        sum(det$sensitivities[[n]] * phi))
}
