#' Input validation utilities for BSSUnfoldR
#'
#' Helper functions that validate detector readings, energy grids, spectra,
#' and the response matrix / measurement vector pair used by every solver.
#' These mirror the corresponding validators in the Python \code{bssunfold}
#' package (\code{utils/validators.py}).
#'
#' @name validators
#' @rdname validators
NULL

#' @rdname validators
#'
#' @description
#' \code{validate_readings} filters the readings list to the names present in
#' \code{detector_names} and rejects negative, infinite or NaN values.
#'
#' @param readings Named numeric vector of detector readings.
#' @param detector_names Character vector of allowed detector names.
#' @param allow_zero Logical; if \code{FALSE} zero readings are rejected.
#'   Default \code{TRUE}.
#' @return A named numeric vector with the validated readings.
#' @export
#'
#' @examples
#' readings <- c("0in" = 100, "3in" = 200, "12in" = 50)
#' validate_readings(readings, c("0in", "3in", "12in"))
validate_readings <- function(readings, detector_names, allow_zero = TRUE) {
    if (!is.numeric(readings)) {
        stop("'readings' must be a named numeric vector, got ",
             class(readings)[1L])
    }
    if (length(readings) == 0L) {
        stop("'readings' must be a non-empty numeric vector.")
    }
    if (is.null(names(readings))) {
        stop("'readings' must be a named vector.")
    }
    valid <- stats::setNames(numeric(0L), character(0L))
    for (det in detector_names) {
        if (det %in% names(readings)) {
            v <- as.numeric(readings[[det]])
            if (is.na(v)) stop("Reading '", det, "' is NA")
            if (is.infinite(v)) stop("Reading '", det, "' is infinite")
            if (v < 0) stop("Reading '", det, "' is negative: ", v)
            if (v == 0 && !allow_zero) {
                stop("Reading '", det, "' is zero, which is not allowed")
            }
            valid[det] <- v
        }
    }
    if (length(valid) == 0L) {
        stop("No valid detector readings provided. Available detectors: ",
             paste(detector_names, collapse = ", "))
    }
    valid
}

#' @rdname validators
#'
#' @description
#' \code{validate_energy_grid} checks that the energy grid is a strictly
#' increasing positive numeric vector of length at least \code{min_points}.
#'
#' @param E_MeV Numeric vector of energy grid points in MeV.
#' @param min_points Integer; minimum required length. Default 2.
#' @param Emin Optional numeric lower bound. Default \code{NULL} (no bound).
#' @param Emax Optional numeric upper bound. Default \code{NULL} (no bound).
#' @return A numeric vector (the validated energy grid).
#' @export
#'
#' @examples
#' E <- 10^seq(-9, 1, length.out = 10)
#' validate_energy_grid(E)
validate_energy_grid <- function(E_MeV, min_points = 2L, Emin = NULL, Emax = NULL) {
    E_MeV <- as.numeric(E_MeV)
    if (any(is.na(E_MeV))) stop("Energy grid contains NA values")
    if (length(E_MeV) < min_points) {
        stop("Energy grid must have at least ", min_points, " points, got ",
             length(E_MeV))
    }
    if (any(E_MeV <= 0)) stop("All energy values must be positive")
    if (any(diff(E_MeV) <= 0)) stop("Energy grid must be strictly increasing")
    if (!is.null(Emin) && E_MeV[1L] < Emin) {
        stop("Minimum energy ", E_MeV[1L], " is below allowed minimum ", Emin)
    }
    if (!is.null(Emax) && E_MeV[length(E_MeV)] > Emax) {
        stop("Maximum energy ", E_MeV[length(E_MeV)],
             " is above allowed maximum ", Emax)
    }
    E_MeV
}

#' @rdname validators
#'
#' @description
#' \code{validate_spectrum} checks that a spectrum vector is finite, the right
#' length, and (by default) non-negative.
#'
#' @param spectrum Numeric vector of spectrum values.
#' @param E_MeV Numeric vector of energy grid points.
#' @param allow_negative Logical; if \code{TRUE} negative values are accepted.
#'   Default \code{FALSE}.
#' @return The validated spectrum as a numeric vector.
#' @export
#'
#' @examples
#' E <- 10^seq(-9, 1, length.out = 10)
#' spec <- rep(0.5, 10)
#' validate_spectrum(spec, E)
validate_spectrum <- function(spectrum, E_MeV, allow_negative = FALSE) {
    spectrum <- as.numeric(spectrum)
    if (any(is.na(spectrum))) stop("Spectrum contains NA values")
    if (any(is.infinite(spectrum))) stop("Spectrum contains infinite values")
    if (length(spectrum) != length(E_MeV)) {
        stop("Spectrum length (", length(spectrum),
             ") must match energy grid length (", length(E_MeV), ")")
    }
    if (!allow_negative && any(spectrum < 0)) {
        n_neg <- sum(spectrum < 0)
        stop("Spectrum contains ", n_neg, " negative values. ",
             "Set allow_negative=TRUE to allow negative values.")
    }
    spectrum
}

#' @rdname validators
#'
#' @description
#' \code{validate_response_matrix} checks the dimensions and finiteness of the
#' response matrix \code{A} and measurement vector \code{b}.
#'
#' @param A Numeric matrix (m x n).
#' @param b Numeric vector of length m.
#' @param check_rank Logical; if \code{TRUE} warns when A is rank-deficient.
#'   Default \code{FALSE}.
#' @return A list \code{list(A = A, b = b)} with validated values.
#' @export
#'
#' @examples
#' A <- matrix(c(1.0, 0.5, 0.2, 0.4), nrow = 2)
#' b <- c(1, 0.6)
#' validate_response_matrix(A, b)
validate_response_matrix <- function(A, b, check_rank = FALSE) {
    A <- as.matrix(A)
    b <- as.numeric(b)
    if (!is.numeric(A)) stop("'A' must be a numeric matrix")
    if (nrow(A) == 0L || ncol(A) == 0L) stop("Response matrix A is empty")
    if (any(is.na(A)) || any(is.infinite(A))) {
        stop("Response matrix A contains NA/Inf values")
    }
    if (length(b) == 0L) stop("Measurement vector b is empty")
    if (any(is.na(b)) || any(is.infinite(b))) {
        stop("Measurement vector b contains NA/Inf values")
    }
    if (nrow(A) != length(b)) {
        stop("Number of rows in A (", nrow(A),
             ") must match length of b (", length(b), ")")
    }
    if (check_rank) {
        r <- qr(A)$rank
        if (r < min(dim(A))) {
            warning("Response matrix is rank-deficient: rank=", r,
                    ", shape=(", nrow(A), ", ", ncol(A), ")")
        }
    }
    list(A = A, b = b)
}

#' @rdname validators
#'
#' @description
#' \code{validate_system} is the convenience wrapper used by every iterative
#' solver; it checks shapes, dimension compatibility and (optionally) the
#' positivity of \code{max_iterations} and \code{tolerance}.
#'
#' @param A Numeric matrix (m x n).
#' @param b Numeric vector of length m.
#' @param x0 Optional initial guess (length n). Default \code{NULL}.
#' @param max_iterations Optional positive integer. Default \code{NULL}.
#' @param tolerance Optional non-negative numeric. Default \code{NULL}.
#' @return A list \code{list(A = A, b = b, x0 = x0)} with validated values.
#' @export
#'
#' @examples
#' A <- matrix(c(1.0, 0.5, 0.2, 0.4), nrow = 2)
#' b <- c(1, 0.6)
#' x0 <- c(0.5, 0.5)
#' v <- validate_system(A, b, x0 = x0, max_iterations = 100L, tolerance = 1e-6)
validate_system <- function(A, b, x0 = NULL, max_iterations = NULL,
                           tolerance = NULL) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    b <- as.numeric(b)
    if (!is.numeric(A)) stop("Response matrix A must be numeric")
    if (nrow(A) == 0L || ncol(A) == 0L) stop("Response matrix A is empty")
    if (length(b) == 0L) stop("Measurement vector b is empty")
    if (nrow(A) != length(b)) {
        stop("Row count of A (", nrow(A),
             ") must match length of b (", length(b), ")")
    }
    if (!is.null(x0)) {
        x0 <- as.numeric(x0)
        if (length(x0) != ncol(A)) {
            stop("Length of x0 (", length(x0),
                 ") must match column count of A (", ncol(A), ")")
        }
    }
    if (!is.null(max_iterations)) {
        if (!is.numeric(max_iterations) || length(max_iterations) != 1L ||
            max_iterations <= 0 || max_iterations != as.integer(max_iterations)) {
            stop("max_iterations must be a positive integer, got ",
                 max_iterations)
        }
    }
    if (!is.null(tolerance)) {
        if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance < 0) {
            stop("tolerance must be a non-negative number, got ", tolerance)
        }
    }
    list(A = A, b = b, x0 = x0)
}
