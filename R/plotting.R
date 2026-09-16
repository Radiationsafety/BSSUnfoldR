#' Plotting, comparison and conversion utilities (ggplot2)
#'
#' R port of the \code{bssunload/utils/plotting.py},
#' \code{comparison.py} and \code{converters.py} utility modules.
#' All plot functions return a \code{ggplot} object; viewing/printing it
#' is up to the caller. \code{ggplot2} is a suggested dependency: every
#' function calls \code{require_ggplot2()} which stops with a helpful
#' message when it is not installed.
#'
#' @name plotting-utils
NULL

#' Require the ggplot2 namespace
require_ggplot2 <- function() {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Plotting functions require the suggested package 'ggplot2'. ",
             "Install it with install.packages('ggplot2').")
    }
    invisible(TRUE)
}

#' Extract a tidy data frame from an unfold result
#'
#' @param result A result list with \code{energy} and \code{spectrum} elements.
#' @keywords internal
.to_df <- function(result) {
    data.frame(E_MeV = as.numeric(result$energy),
               spectrum = as.numeric(result$spectrum))
}

#' Plot an unfolded spectrum
#'
#' @param result A result list (\code{energy}, \code{spectrum}).
#' @param logx Logical; log-scale energy axis. Default TRUE.
#' @param logy Logical; log-scale fluence axis. Default TRUE.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' A <- matrix(c(0.9, 0.05, 0.05, 0.10, 0.8, 0.10, 0.30, 0.30, 0.40),
#'             nrow = 3, byrow = TRUE)
#' b <- c(1, 0.6, 0.4)
#' det_names <- c("d1", "d2", "d3")
#' sens <- lapply(setNames(det_names, det_names), function(n) A[which(det_names == n), ])
#' E <- c(1e-9, 1e-6, 1e-3)
#' res <- unfold_mlem(det_names, 3L, E, sens, NULL, NULL, c(d1=1, d2=.6, d3=.4))
#' \donttest{if (requireNamespace("ggplot2", quietly = TRUE)) { print(plot_spectrum(res)) }}
plot_spectrum <- function(result, logx = TRUE, logy = TRUE) {
    require_ggplot2()
    df <- cbind(.to_df(result),
                data.frame(method = as.character(result$method)))
    p <- ggplot2::ggplot(df,
                         ggplot2::aes(x = .data$E_MeV, y = .data$spectrum,
                                      colour = .data$method)) +
        ggplot2::geom_point() +
        ggplot2::labs(x = "Energy [MeV]", y = "Fluence (per leth.)",
                      colour = "Method")
    if (isTRUE(logx)) p <- p + ggplot2::scale_x_log10()
    if (isTRUE(logy)) p <- p + ggplot2::scale_y_log10(
        limits = c(NA, NA))
    p
}

#' Plot measurement residuals of an unfold result
#'
#' @param result A result list (\code{effective_readings},
#'   \code{residual}).
#' @param detector_names Optional names of the detectors (used for
#'   x-labelling).
#' @return A \code{ggplot} object.
#' @export
plot_residuals <- function(result, detector_names = NULL) {
    require_ggplot2()
    names <- if (!is.null(detector_names)) detector_names
             else names(result$residual)
    if (is.null(names)) names <- seq_along(result$residual)
    df <- data.frame(detector = names, residual = as.numeric(result$residual))
    ggplot2::ggplot(df, ggplot2::aes(x = .data$detector,
                    y = .data$residual)) +
        ggplot2::geom_col() +
        ggplot2::labs(x = "Detector", y = "Residual (measured - fitted)",
                      title = "Residuals of the unfolding")
}

#' Overlay comparison of two or more unfold results
#'
#' @param results Named list of result lists; each needs
#'   \code{energy}, \code{spectrum}.
#' @return A \code{ggplot} object.
#' @export
plot_comparison <- function(results) {
    require_ggplot2()
    if (length(results) == 0L) stop("results must be a non-empty list.")
    if (is.null(names(results))) {
        names(results) <- paste0("method_", seq_along(results))
    }
    df <- do.call(rbind, lapply(names(results), function(k) {
        r <- results[[k]]
        data.frame(E_MeV = as.numeric(r$energy),
                   spectrum = as.numeric(r$spectrum), method = k)
    }))
    ggplot2::ggplot(df, ggplot2::aes(x = .data$E_MeV, y = .data$spectrum,
                                     colour = .data$method)) +
        ggplot2::geom_line() +
        ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
        ggplot2::labs(x = "Energy [MeV]", y = "Fluence per unit lethargy",
                      title = "Comparison of unfolded spectra")
}

# ---- comparison metrics -----------------------------------------------------

#' Compare, divergence and related metrics for two spectra
#'
#' @param x,y Numeric spectra (same length).
#' @return A named list of scalar metrics: \code{relative_integral_diff},
#'   \code{chi_squared}, \code{max_relative_deviation},
#'   \code{ratio_mean}, \code{ratio_median}, \code{kl_divergence}
#'   (asymmetric, base e), \code{euclidean_distance},
#'   \code{cosine_similarity}, \code{fluence_mean_energy_x},
#'   \code{fluence_mean_energy_y}.
#' @export
#' @examples
#' x <- c(1, 0.5, 0.25); y <- c(0.9, 0.55, 0.2)
#' compare_spectra(x, y)
compare_spectra <- function(x, y) {
    x <- as.numeric(x); y <- as.numeric(y)
    if (length(x) != length(y) || length(x) == 0L) {
        stop("compare_spectra: spectra must have equal non-zero length.")
    }
    sx <- sum(x); sy <- sum(y)
    rel_diff <- abs(sx - sy) / max(abs(sy), 1e-30)
    chi2 <- sum((x - y)^2 / pmax(y^2, 1e-30))
    maxdev <- max(abs(x - y) / pmax(abs(y), 1e-30))
    ratio <- x / pmax(y, 1e-30)
    p <- x / max(sx, 1e-30)
    q <- y / max(sy, 1e-30)
    kl <- sum(p * log(p / pmax(q, 1e-30)))
    euc <- sqrt(sum((x - y)^2))
    cossim <- sum(x * y) / max(sqrt(sum(x^2) * sum(y^2)), 1e-30)
    fuex <- sum(x * seq_along(x)) / max(sx, 1e-30)
    fuey <- sum(y * seq_along(x)) / max(sy, 1e-30)
    list(relative_integral_diff = rel_diff, chi_squared = chi2,
         max_relative_deviation = maxdev,
         ratio_mean = mean(ratio), ratio_median = stats::median(ratio),
         kl_divergence = kl, euclidean_distance = euc,
         cosine_similarity = cossim,
         fluence_mean_energy_x = fuex, fluence_mean_energy_y = fuey)
}

# ---- converters -------------------------------------------------------------

#' Convert an unfold result to a data frame
#' @param result A result list.
#' @param sig_figs Optional integer; round numeric columns to that many
#'   significant figures.
#' @return A \code{data.frame} with columns \code{E_MeV},
#'   \code{spectrum}, \code{doserate} (when available).
#' @export
convert_to_dataframe <- function(result, sig_figs = NULL) {
    df <- .to_df(result)
    if (!is.null(result$doserates) && length(result$doserates) ==
        nrow(df)) {
        df$doserate <- as.numeric(result$doserates)
    }
    if (!is.null(sig_figs)) {
        sig_figs <- as.integer(sig_figs)
        df[] <- lapply(df, function(cl)
            if (is.numeric(cl)) round(cl, sig_figs) else cl)
    }
    df
}

#' Convert an unfold result to a plain named list of vectors
#' @inheritParams convert_to_dataframe
#' @return A list \code{list(energy, spectrum, doserates)}.
#' @export
convert_to_dict <- function(result) {
    out <- list(energy = as.numeric(result$energy),
                spectrum = as.numeric(result$spectrum))
    if (!is.null(result$doserates)) {
        out$doserates <- as.numeric(result$doserates)
    }
    out
}

#' Discretize an analytic spectrum on an energy grid (helper for input
#' generation)
#'
#' @param f Function of energy (MeV) returning a non-negative spectral
#'   density.
#' @param E_MeV Energy grid.
#' @return Numeric vector of bin-integrated fluence approximated by the
#'   midpoint of consecutive grid edges (log grid assumed).
#' @export
discretize_spectrum <- function(f, E_MeV) {
    E_MeV <- sort(as.numeric(E_MeV))
    edgesL <- log(E_MeV)
    lore <- exp(edgesL - step / 2)
    hire <- exp(edgesL + step / 2)
    vapply(seq_len(length(E_MeV)), function(i) {
        stats::integrate(f, lore[i], hire[i])$value
    }, numeric(1))
}
