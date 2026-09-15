#' Access conversion-coefficient and response-function constants
#'
#' The BSSUnfoldR package ships a private environment variable
#' \code{.BSSUNFOLD_CONSTANTS} (defined in \code{R/constants_data.R}) with the
#' same numeric values as the Python \code{bssunfold} package's
#' \code{constants.py}. These accessors expose the individual datasets to the
#' user and also build the dose-coefficient registry used by
#' \code{\link{calculate_dose_rates}}.
#'
#' @name constants
#' @rdname constants
#'
#' @return Each accessor returns a named list with an \code{"E_MeV"} entry
#'   (energy grid in MeV) plus one numeric vector per geometry or detector.
#'   \itemize{
#'     \item \code{ICRP116_COEFF_EFFECTIVE_DOSE} -- geometries: AP, PA, LLAT, RLAT, ROT, ISO.
#'     \item \code{ICRP74_COEFF_EFFECTIVE_DOSE} -- geometries: AP, PA, RLAT, ROT, ISO.
#'     \item \code{NRB99_2009_COEFF_EFFECTIVE_DOSE} -- geometries: AP, ISO.
#'     \item \code{ICRP74_COEFF_OPERATIONAL_QUANTITIES} -- quantities: ADE, PDE0, PDE45, PDE60, PDE75.
#'     \item \code{RF_GSF}, \code{RF_PTB}, \code{RF_LANL}, \code{RF_JINR}, \code{RF_FERMILAB}, \code{RF_EURADOS}, \code{RF_IHEP} -- Bonner-sphere response functions keyed by sphere diameter (e.g. \code{"0in"}, \code{"3in"}, \code{"12in"}, ...).
#'   }
#'
#' @examples
#' cc <- ICRP116_COEFF_EFFECTIVE_DOSE()
#' str(cc[1:3])
#' rf <- RF_PTB()
#' cat("PTB detectors:", paste(names(rf), collapse = ", "), "\n")
NULL

#' @rdname constants
#' @export
ICRP116_COEFF_EFFECTIVE_DOSE <- function() {
    .BSSUNFOLD_CONSTANTS$ICRP116_COEFF_EFFECTIVE_DOSE
}

#' @rdname constants
#' @export
ICRP74_COEFF_EFFECTIVE_DOSE <- function() {
    .BSSUNFOLD_CONSTANTS$ICRP74_COEFF_EFFECTIVE_DOSE
}

#' @rdname constants
#' @export
NRB99_2009_COEFF_EFFECTIVE_DOSE <- function() {
    .BSSUNFOLD_CONSTANTS$NRB99_2009_COEFF_EFFECTIVE_DOSE
}

#' @rdname constants
#' @export
ICRP74_COEFF_OPERATIONAL_QUANTITIES <- function() {
    .BSSUNFOLD_CONSTANTS$ICRP74_COEFF_OPERATIONAL_QUANTITIES
}

#' @rdname constants
#' @export
RF_GSF <- function() .BSSUNFOLD_CONSTANTS$RF_GSF

#' @rdname constants
#' @export
RF_PTB <- function() .BSSUNFOLD_CONSTANTS$RF_PTB

#' @rdname constants
#' @export
RF_LANL <- function() .BSSUNFOLD_CONSTANTS$RF_LANL

#' @rdname constants
#' @export
RF_JINR <- function() .BSSUNFOLD_CONSTANTS$RF_JINR

#' @rdname constants
#' @export
RF_FERMILAB <- function() .BSSUNFOLD_CONSTANTS$RF_FERMILAB

#' @rdname constants
#' @export
RF_EURADOS <- function() .BSSUNFOLD_CONSTANTS$RF_EURADOS

#' @rdname constants
#' @export
RF_IHEP <- function() .BSSUNFOLD_CONSTANTS$RF_IHEP

#' Registry of available dose conversion coefficient datasets
#'
#' Returns a named list of the four built-in dose-coefficient datasets
#' (ICRP-116 effective dose, ICRP-74 effective dose, NRB99-2009 effective dose,
#' and ICRP-74 operational quantities). Used internally by
#' \code{\link{get_coefficients}} and \code{\link{calculate_dose_rates}}.
#'
#' @return A named list with entries \code{ICRP116}, \code{ICRP74_effective},
#'   \code{NRB99_2009_effective}, \code{ICRP74_operational}.
#' @export
#' @examples
#' names(DOSE_COEFFICIENTS_REGISTRY())
DOSE_COEFFICIENTS_REGISTRY <- function() {
    list(
        ICRP116 = ICRP116_COEFF_EFFECTIVE_DOSE(),
        ICRP74_effective = ICRP74_COEFF_EFFECTIVE_DOSE(),
        NRB99_2009_effective = NRB99_2009_COEFF_EFFECTIVE_DOSE(),
        ICRP74_operational = ICRP74_COEFF_OPERATIONAL_QUANTITIES()
    )
}

#' Look up a dose conversion coefficient dataset by name
#'
#' @param name Character scalar naming a registry entry: \code{"ICRP116"},
#'   \code{"ICRP74_effective"}, \code{"NRB99_2009_effective"}, or
#'   \code{"ICRP74_operational"}.
#' @return A named list with an \code{"E_MeV"} entry plus one numeric vector
#'   per geometry/quantity.
#' @export
#' @examples
#' cc <- get_coefficients("ICRP74_effective")
#' cat("Geometries:", paste(names(cc), collapse = ", "), "\n")
get_coefficients <- function(name = c("ICRP116", "ICRP74_effective",
                                      "NRB99_2009_effective", "ICRP74_operational")) {
    if (length(name) != 1L || is.na(name) || !nzchar(name)) {
        stop("'name' must be a single non-empty character string.")
    }
    reg <- DOSE_COEFFICIENTS_REGISTRY()
    if (!(name %in% names(reg))) {
        stop("Unknown dose coefficient name: '", name, "'. Available options: ",
             paste(names(reg), collapse = ", "))
    }
    reg[[name]]
}

#' Interpolate conversion coefficients to a target energy grid
#'
#' Linear interpolation in energy of any of the built-in conversion-coefficient
#' datasets. Values outside the source energy range are replaced with
#' \code{fill_value} (default 0).
#'
#' @param cc A named list (such as returned by \code{\link{get_coefficients}})
#'   containing an \code{"E_MeV"} entry and one or more numeric vectors.
#' @param E_target Numeric vector giving the target energy grid in MeV.
#' @param fill_value Numeric scalar; value to assign for energies outside the
#'   source grid range. Default 0.
#' @return A named list with the same keys as \code{cc} but with the geometry
#'   vectors replaced by interpolated values on \code{E_target}.
#' @export
#' @examples
#' cc <- get_coefficients("NRB99_2009_effective")
#' E_det <- 10^seq(-9, 2, length.out = 50)
#' cc_interp <- interpolate_coefficients(cc, E_det)
interpolate_coefficients <- function(cc, E_target, fill_value = 0.0) {
    if (!is.list(cc) || is.null(cc$E_MeV)) {
        stop("'cc' must be a list with an 'E_MeV' entry.")
    }
    if (!is.numeric(E_target)) {
        stop("'E_target' must be a numeric vector.")
    }
    E_source <- as.numeric(cc$E_MeV)
    E_target <- as.numeric(E_target)
    out <- list(E_MeV = E_target)
    for (key in names(cc)) {
        if (key == "E_MeV") next
        values <- as.numeric(cc[[key]])
        n_src <- length(values)
        n_tgt <- length(E_target)
        if (n_src == 0L || n_tgt == 0L) {
            out[[key]] <- numeric(0)
            next
        }
        interp <- approx(x = E_source, y = values, xout = E_target,
                         method = "linear", rule = 2)$y
        below <- E_target < E_source[1L]
        above <- E_target > E_source[n_src]
        interp[below] <- fill_value
        interp[above] <- fill_value
        out[[key]] <- interp
    }
    out
}
