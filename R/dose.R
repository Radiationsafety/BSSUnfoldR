#' Dose-rate calculation
#'
#' Computes dose rates from an unfolded neutron spectrum using one of the
#' built-in ICRP-74 / ICRP-116 / NRB99-2009 conversion-coefficient datasets.
#' Mirrors \code{core/dose_calculation.py} of the Python \code{bssunfold}
#' package.
#'
#' @param spectrum Numeric vector with the unfolded neutron spectrum (length
#'   equal to the energy grid of the conversion coefficients, or shorter;
#'   missing tail is zero-padded).
#' @param cc_icrp116 Named list with an \code{"E_MeV"} entry plus one numeric
#'   vector per geometry/quantity. If \code{NULL} (default), the built-in
#'   ICRP-116 effective-dose dataset is used.
#' @param dlnE Numeric; uniform log-E step used for integration. Default 0.2.
#' @return Named numeric vector with one dose rate per geometry/quantity, in
#'   pico-Sievert per second (pSv/s).
#' @export
#' @examples
#' cc <- ICRP116_COEFF_EFFECTIVE_DOSE()
#' spec <- rep(1.0, length(cc$E_MeV))
#' doses <- calculate_dose_rates(spec, cc)
#' head(doses)
calculate_dose_rates <- function(spectrum, cc_icrp116 = NULL, dlnE = 0.2) {
    if (is.null(cc_icrp116)) {
        cc_icrp116 <- ICRP116_COEFF_EFFECTIVE_DOSE()
    }
    if (!is.list(cc_icrp116) || length(cc_icrp116) == 0L) {
        return(numeric(0L))
    }
    if (!is.numeric(spectrum)) stop("'spectrum' must be a numeric vector")
    spectrum <- as.numeric(spectrum)
    n_spec <- length(spectrum)
    if (n_spec == 0L) return(numeric(0L))

    ln10 <- log(10.0) * dlnE
    geoms <- setdiff(names(cc_icrp116), "E_MeV")
    if (length(geoms) == 0L) return(numeric(0L))

    # Stack coefficient vectors into a matrix (n_geoms x n_spec)
    cc_mat <- matrix(0.0, nrow = length(geoms), ncol = n_spec)
    for (idx in seq_along(geoms)) {
        k_arr <- as.numeric(cc_icrp116[[geoms[idx]]])
        min_len <- min(length(k_arr), n_spec)
        if (min_len > 0L) {
            cc_mat[idx, seq_len(min_len)] <- k_arr[seq_len(min_len)]
        }
    }
    doses <- as.numeric(cc_mat %*% spectrum) * ln10
    names(doses) <- geoms
    doses
}
