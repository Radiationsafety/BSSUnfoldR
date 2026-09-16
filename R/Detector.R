#' Detector R6 class for Bonner-sphere neutron spectrum unfolding
#'
#' The \code{Detector} class bundles together everything a user needs to unfold
#' a neutron spectrum from Bonner-sphere measurements:
#'
#' \itemize{
#'   \item an energy grid (MeV);
#'   \item a sensitivity (response function) for each detector sphere;
#'   \item ICRP-116 conversion coefficients for dose-rate calculation;
#'   \item a result-history list (for retrospective comparison);
#'   \item wrappers around all the classic unfolding algorithms provided by
#'         this package (\code{\link{unfold_mlem}}, \code{\link{unfold_gravel}},
#'         \code{\link{unfold_maxed}}, \code{\link{unfold_sandii}},
#'         \code{\link{unfold_bunki}}, \code{\link{unfold_ferdor}},
#'         \code{\link{unfold_staysl}}, \code{\link{unfold_landweber}},
#'         \code{\link{unfold_cgls}}, \code{\link{unfold_tsvd}}).
#' }
#'
#' Construction options cover the common cases:
#'
#' \describe{
#'   \item{\code{response_function = "PTB"}}{use the built-in PTB response
#'     function (\code{\link{RF_PTB}}); other options: \code{"GSF"},
#'     \code{"LANL"}, \code{"JINR"}, \code{"FERMILAB"}, \code{"EURADOS"},
#'     \code{"IHEP"}.}
#'   \item{\code{E_MeV = NULL}}{use the energy grid shipped with the chosen
#'     response function.}
#'   \item{\code{cc_icrp116 = NULL}}{use the built-in
#'     \code{\link{ICRP116_COEFF_EFFECTIVE_DOSE}} dataset interpolated to
#'     \code{E_MeV}.}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{add_response_function(name, sensitivities)}}{add a new named
#'     response function to the detector.}
#'   \item{\code{set_detector_names(...)}}{set or update the ordered list of
#'     detector names that should be used when calling the unfold methods.}
#'   \item{\code{unfold_mlem(readings, ...)}}{call \code{\link{unfold_mlem}}
#'     with the detector's settings. The same pattern repeats for every
#'     other unfold method (see the Examples section).}
#'   \item{\code{save_result(output)}}{append a result to the history list.}
#'   \item{\code{history}}{accessor returning the list of saved results.}
#' }
#'
#' @param E_MeV Optional numeric energy grid. Default \code{NULL} = take the
#'   energy grid of the chosen response function.
#' @param response_function Either a character name (one of \code{"GSF"},
#'   \code{"PTB"}, \code{"LANL"}, \code{"JINR"}, \code{"FERMILAB"},
#'   \code{"EURADOS"}, \code{"IHEP"}) or a named list with an \code{"E_MeV"}
#'   entry plus one numeric vector per detector.
#' @param cc_icrp116 Optional named list of ICRP-116 conversion coefficients.
#'   Default \code{NULL} = use the built-in dataset interpolated to
#'   \code{E_MeV}.
#' @param detector_names Optional character vector of detector names. If
#'   \code{NULL}, all non-\code{E_MeV} keys of the response function are used.
#'
#' @return An R6 \code{Detector} object.
#' @export
#' @examples
#' det <- Detector$new(response_function = "PTB",
#'                     detector_names = c("0in", "3in", "5in", "8in", "12in"))
#' cat("Energy bins:", det$n_energy_bins, "\n")
#' cat("Detectors:", paste(det$detector_names, collapse = ", "), "\n")
Detector <- R6::R6Class(
    classname = "Detector",
    public = list(
        #' @field E_MeV Numeric energy grid (MeV).
        E_MeV = numeric(0L),
        #' @field n_energy_bins Integer; number of energy bins.
        n_energy_bins = 0L,
        #' @field sensitivities Named list of numeric sensitivity vectors.
        sensitivities = list(),
        #' @field detector_names Character vector of detector names to use.
        detector_names = character(0L),
        #' @field cc_icrp116 Named list of ICRP-116 conversion coefficients.
        cc_icrp116 = list(),
        #' @field history List of saved unfolding results.
        history = list(),

        #' @description
        #' Create a new Detector object. See \code{\link{Detector}} for the
        #' parameter documentation.
        initialize = function(E_MeV = NULL,
                              response_function = "PTB",
                              cc_icrp116 = NULL,
                              detector_names = NULL) {
            rf <- if (is.character(response_function) &&
                         length(response_function) == 1L) {
                fn_name <- paste0("RF_", response_function)
                if (!exists(fn_name, mode = "function")) {
                    stop("Unknown built-in response function '",
                         response_function,
                         "'. Available: GSF, PTB, LANL, JINR, FERMILAB, ",
                         "EURADOS, IHEP.")
                }
                get(fn_name, mode = "function")()
            } else if (is.list(response_function)) {
                response_function
            } else {
                stop("'response_function' must be a string name or a list.")
            }
            if (is.null(rf$E_MeV)) {
                stop("Response function must have an 'E_MeV' entry.")
            }
            if (is.null(E_MeV)) {
                private$.E_MeV <- as.numeric(rf$E_MeV)
            } else {
                private$.E_MeV <- validate_energy_grid(as.numeric(E_MeV))
            }
            self$E_MeV <- private$.E_MeV
            self$n_energy_bins <- length(private$.E_MeV)

            # Sensitivities: interpolate every detector vector to E_MeV
            sens <- list()
            for (key in setdiff(names(rf), "E_MeV")) {
                arr <- as.numeric(rf[[key]])
                if (length(arr) == length(rf$E_MeV)) {
                    sens[[key]] <- arr
                } else {
                    sens[[key]] <- approx(x = rf$E_MeV, y = arr,
                                          xout = private$.E_MeV,
                                          rule = 2)$y
                }
            }
            self$sensitivities <- sens
            if (is.null(detector_names)) {
                self$detector_names <- names(sens)
            } else {
                self$detector_names <- detector_names
            }
            # ICRP-116 conversion coefficients interpolated to E_MeV
            if (is.null(cc_icrp116)) {
                self$cc_icrp116 <- interpolate_coefficients(
                    ICRP116_COEFF_EFFECTIVE_DOSE(), self$E_MeV)
            } else {
                if (!is.null(cc_icrp116$E_MeV) &&
                    !identical(as.numeric(cc_icrp116$E_MeV),
                               private$.E_MeV)) {
                    self$cc_icrp116 <- interpolate_coefficients(
                        cc_icrp116, private$.E_MeV)
                } else {
                    self$cc_icrp116 <- cc_icrp116
                }
            }
            invisible(self)
        },

        #' @description Add or overwrite a named sensitivity entry.
        #' @param name Character detector name.
        #' @param sensitivity Numeric sensitivity vector (length
        #'   \code{n_energy_bins}).
        add_response_function = function(name, sensitivity) {
            if (!is.character(name) || length(name) != 1L) {
                stop("'name' must be a single character string.")
            }
            sensitivity <- as.numeric(sensitivity)
            if (length(sensitivity) != self$n_energy_bins) {
                stop("sensitivity length (", length(sensitivity),
                     ") must match n_energy_bins (", self$n_energy_bins, ")")
            }
            self$sensitivities[[name]] <- sensitivity
            if (!(name %in% self$detector_names)) {
                self$detector_names <- c(self$detector_names, name)
            }
            invisible(self)
        },

        #' @description Set the ordered list of detector names.
        #' @param detector_names Character vector of detector names; must be a
        #'   subset of \code{names(self$sensitivities)}.
        set_detector_names = function(detector_names) {
            if (!is.character(detector_names) || length(detector_names) == 0L) {
                stop("'detector_names' must be a non-empty character vector.")
            }
            missing <- setdiff(detector_names, names(self$sensitivities))
            if (length(missing) > 0L) {
                stop("The following detector names are not present in the ",
                     "sensitivities list: ",
                     paste(missing, collapse = ", "))
            }
            self$detector_names <- detector_names
            invisible(self)
        },

        #' @description Append a result list to the history.
        #' @param output A list (as produced by \code{\link{run_unfolding}}).
        save_result = function(output) {
            self$history <- c(self$history, list(output))
            invisible(self)
        },

        # --- Unfold methods ---

        #' @description Unfold using MLEM. See \code{\link{unfold_mlem}}.
        #' @param readings Named numeric vector.
        #' @param ... Extra arguments forwarded to \code{\link{unfold_mlem}}.
        unfold_mlem = function(readings, ...) {
            unfold_mlem(self$detector_names, self$n_energy_bins,
                        self$E_MeV, self$sensitivities, self$cc_icrp116,
                        function(out) self$save_result(out),
                        readings, ...)
        },

        #' @description Unfold using GRAVEL. See \code{\link{unfold_gravel}}.
        unfold_gravel = function(readings, ...) {
            unfold_gravel(self$detector_names, self$n_energy_bins,
                          self$E_MeV, self$sensitivities, self$cc_icrp116,
                          function(out) self$save_result(out),
                          readings, ...)
        },

        #' @description Unfold using MAXED. See \code{\link{unfold_maxed}}.
        unfold_maxed = function(readings, ...) {
            unfold_maxed(self$detector_names, self$n_energy_bins,
                         self$E_MeV, self$sensitivities, self$cc_icrp116,
                         function(out) self$save_result(out),
                         readings, ...)
        },

        #' @description Unfold using SAND-II. See \code{\link{unfold_sandii}}.
        unfold_sandii = function(readings, ...) {
            unfold_sandii(self$detector_names, self$n_energy_bins,
                          self$E_MeV, self$sensitivities, self$cc_icrp116,
                          function(out) self$save_result(out),
                          readings, ...)
        },

        #' @description Unfold using BUNKI. See \code{\link{unfold_bunki}}.
        unfold_bunki = function(readings, ...) {
            unfold_bunki(self$detector_names, self$n_energy_bins,
                         self$E_MeV, self$sensitivities, self$cc_icrp116,
                         function(out) self$save_result(out),
                         readings, ...)
        },

        #' @description Unfold using FERDOR. See \code{\link{unfold_ferdor}}.
        unfold_ferdor = function(readings, ...) {
            unfold_ferdor(self$detector_names, self$n_energy_bins,
                          self$E_MeV, self$sensitivities, self$cc_icrp116,
                          function(out) self$save_result(out),
                          readings, ...)
        },

        #' @description Unfold using STAY'SL. See \code{\link{unfold_staysl}}.
        unfold_staysl = function(readings, ...) {
            unfold_staysl(self$detector_names, self$n_energy_bins,
                           self$E_MeV, self$sensitivities, self$cc_icrp116,
                           function(out) self$save_result(out),
                           readings, ...)
        },

        #' @description Unfold using Landweber. See \code{\link{unfold_landweber}}.
        unfold_landweber = function(readings, ...) {
            unfold_landweber(self$detector_names, self$n_energy_bins,
                             self$E_MeV, self$sensitivities, self$cc_icrp116,
                             function(out) self$save_result(out),
                             readings, ...)
        },

        #' @description Unfold using CGLS. See \code{\link{unfold_cgls}}.
        unfold_cgls = function(readings, ...) {
            unfold_cgls(self$detector_names, self$n_energy_bins,
                         self$E_MeV, self$sensitivities, self$cc_icrp116,
                         function(out) self$save_result(out),
                         readings, ...)
        },

        #' @description Unfold using TSVD. See \code{\link{unfold_tsvd}}.
        unfold_tsvd = function(readings, ...) {
            unfold_tsvd(self$detector_names, self$n_energy_bins,
                        self$E_MeV, self$sensitivities, self$cc_icrp116,
                        function(out) self$save_result(out),
                        readings, ...)
        },

        # ------------------------------------------------------------------
        # Second batch of 10 algorithms (added in v0.1.1)
        # ------------------------------------------------------------------

        #' @description Unfold using OSEM. See \code{\link{unfold_osem}}.
        unfold_osem = function(readings, ...) {
            unfold_osem(self$detector_names, self$n_energy_bins,
                         self$E_MeV, self$sensitivities, self$cc_icrp116,
                         function(out) self$save_result(out),
                         readings, ...)
        },

        #' @description Unfold using MAP-EM. See \code{\link{unfold_mapem}}.
        unfold_mapem = function(readings, ...) {
            unfold_mapem(self$detector_names, self$n_energy_bins,
                         self$E_MeV, self$sensitivities, self$cc_icrp116,
                         function(out) self$save_result(out),
                         readings, ...)
        },

        #' @description Unfold using BSREM. See \code{\link{unfold_bsrem}}.
        unfold_bsrem = function(readings, ...) {
            unfold_bsrem(self$detector_names, self$n_energy_bins,
                          self$E_MeV, self$sensitivities, self$cc_icrp116,
                          function(out) self$save_result(out),
                          readings, ...)
        },

        #' @description Unfold using SART. See \code{\link{unfold_sart}}.
        unfold_sart = function(readings, ...) {
            unfold_sart(self$detector_names, self$n_energy_bins,
                         self$E_MeV, self$sensitivities, self$cc_icrp116,
                         function(out) self$save_result(out),
                         readings, ...)
        },

        #' @description Unfold using Kaczmarz. See \code{\link{unfold_kaczmarz}}.
        unfold_kaczmarz = function(readings, ...) {
            unfold_kaczmarz(self$detector_names, self$n_energy_bins,
                              self$E_MeV, self$sensitivities, self$cc_icrp116,
                              function(out) self$save_result(out),
                              readings, ...)
        },

        #' @description Unfold using Randomized Kaczmarz. See \code{\link{unfold_randomized_kaczmarz}}.
        unfold_randomized_kaczmarz = function(readings, ...) {
            unfold_randomized_kaczmarz(self$detector_names, self$n_energy_bins,
                                         self$E_MeV, self$sensitivities,
                                         self$cc_icrp116,
                                         function(out) self$save_result(out),
                                         readings, ...)
        },

        #' @description Unfold using Lanczos. See \code{\link{unfold_lanczos}}.
        unfold_lanczos = function(readings, ...) {
            unfold_lanczos(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold using Tikhonov-Legendre. See \code{\link{unfold_tikhonov_legendre}}.
        unfold_tikhonov_legendre = function(readings, ...) {
            unfold_tikhonov_legendre(self$detector_names, self$n_energy_bins,
                                        self$E_MeV, self$sensitivities,
                                        self$cc_icrp116,
                                        function(out) self$save_result(out),
                                        readings, ...)
        },

        #' @description Unfold using ReBUNKI. See \code{\link{unfold_rebunki}}.
        unfold_rebunki = function(readings, ...) {
            unfold_rebunki(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold using Doroshenko. See \code{\link{unfold_doroshenko}}.
        unfold_doroshenko = function(readings, ...) {
            unfold_doroshenko(self$detector_names, self$n_energy_bins,
                                self$E_MeV, self$sensitivities,
                                self$cc_icrp116,
                                function(out) self$save_result(out),
                                readings, ...)
        },

        # ------------------------------------------------------------------
        # Third batch of 10 algorithms (added in v0.1.2)
        # ------------------------------------------------------------------

        #' @description Unfold using Bayes (D'Agostini). See \code{\link{unfold_bayes}}.
        unfold_bayes = function(readings, ...) {
            unfold_bayes(self$detector_names, self$n_energy_bins,
                          self$E_MeV, self$sensitivities, self$cc_icrp116,
                          function(out) self$save_result(out),
                          readings, ...)
        },

        #' @description Unfold using Directed divergence. See \code{\link{unfold_directed_divergence}}.
        unfold_directed_divergence = function(readings, ...) {
            unfold_directed_divergence(self$detector_names, self$n_energy_bins,
                                          self$E_MeV, self$sensitivities,
                                          self$cc_icrp116,
                                          function(out) self$save_result(out),
                                          readings, ...)
        },

        #' @description Unfold using Express. See \code{\link{unfold_express}}.
        unfold_express = function(readings, ...) {
            unfold_express(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold using Iterative refinement. See \code{\link{unfold_iterative_refinement}}.
        unfold_iterative_refinement = function(readings, ...) {
            unfold_iterative_refinement(self$detector_names, self$n_energy_bins,
                                           self$E_MeV, self$sensitivities,
                                           self$cc_icrp116,
                                           function(out) self$save_result(out),
                                           readings, ...)
        },

        #' @description Unfold using BUNKI-UT. See \code{\link{unfold_bunkiut}}.
        unfold_bunkiut = function(readings, ...) {
            unfold_bunkiut(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold using MLEM-STOP. See \code{\link{unfold_mlem_stop}}.
        unfold_mlem_stop = function(readings, ...) {
            unfold_mlem_stop(self$detector_names, self$n_energy_bins,
                              self$E_MeV, self$sensitivities, self$cc_icrp116,
                              function(out) self$save_result(out),
                              readings, ...)
        },

        #' @description Unfold using StatReg. See \code{\link{unfold_statreg}}.
        unfold_statreg = function(readings, ...) {
            unfold_statreg(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold using Tikhonov-TV. See \code{\link{unfold_tikhonov_tv}}.
        unfold_tikhonov_tv = function(readings, ...) {
            unfold_tikhonov_tv(self$detector_names, self$n_energy_bins,
                                self$E_MeV, self$sensitivities, self$cc_icrp116,
                                function(out) self$save_result(out),
                                readings, ...)
        },

        #' @description Unfold using GKS. See \code{\link{unfold_gks}}.
        unfold_gks = function(readings, ...) {
            unfold_gks(self$detector_names, self$n_energy_bins,
                        self$E_MeV, self$sensitivities, self$cc_icrp116,
                        function(out) self$save_result(out),
                        readings, ...)
        },

        #' @description Unfold using Crystal Ball. See \code{\link{unfold_crystal_ball}}.
        unfold_crystal_ball = function(readings, ...) {
            unfold_crystal_ball(self$detector_names, self$n_energy_bins,
                                  self$E_MeV, self$sensitivities,
                                  self$cc_icrp116,
                                  function(out) self$save_result(out),
                                  readings, ...)
        },

        # ------------------------------------------------------------------
        # Fourth batch of algorithms (added in v0.1.3)
        # ------------------------------------------------------------------

        #' @description Unfold using IMAXED. See \code{\link{unfold_imaxed}}.
        unfold_imaxed = function(readings, ...) {
            unfold_imaxed(self$detector_names, self$n_energy_bins,
                           self$E_MeV, self$sensitivities, self$cc_icrp116,
                           function(out) self$save_result(out),
                           readings, ...)
        },

        #' @description Unfold using AMAXED. See \code{\link{unfold_amaxed}}.
        unfold_amaxed = function(readings, ...) {
            unfold_amaxed(self$detector_names, self$n_energy_bins,
                           self$E_MeV, self$sensitivities, self$cc_icrp116,
                           function(out) self$save_result(out),
                           readings, ...)
        },

        #' @description Unfold using FISTA. See \code{\link{unfold_fista}}.
        unfold_fista = function(readings, ...) {
            unfold_fista(self$detector_names, self$n_energy_bins,
                          self$E_MeV, self$sensitivities, self$cc_icrp116,
                          function(out) self$save_result(out),
                          readings, ...)
        },

        #' @description Unfold using Bayes-spline. See \code{\link{unfold_bayes_spline_regularization}}.
        unfold_bayes_spline_regularization = function(readings, ...) {
            unfold_bayes_spline_regularization(self$detector_names,
                                                  self$n_energy_bins,
                                                  self$E_MeV, self$sensitivities,
                                                  self$cc_icrp116,
                                                  function(out) self$save_result(out),
                                                  readings, ...)
        },

        #' @description Unfold using NSDUAZ. See \code{\link{unfold_nsduaz}}.
        unfold_nsduaz = function(readings, ...) {
            unfold_nsduaz(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold using MLEM-BS. See \code{\link{unfold_mlem_bs}}.
        unfold_mlem_bs = function(readings, ...) {
            unfold_mlem_bs(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold using NSpline. See \code{\link{unfold_nspline}}.
        unfold_nspline = function(readings, ...) {
            unfold_nspline(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold using MCMC. See \code{\link{unfold_mcmc}}.
        unfold_mcmc = function(readings, ...) {
            unfold_mcmc(self$detector_names, self$n_energy_bins,
                         self$E_MeV, self$sensitivities, self$cc_icrp116,
                         function(out) self$save_result(out),
                         readings, ...)
        },

        #' @description Unfold using Reconst. See \code{\link{unfold_reconst}}.
        unfold_reconst = function(readings, ...) {
            unfold_reconst(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold using Ensemble. See \code{\link{unfold_ensemble}}.
        unfold_ensemble = function(readings, ...) {
            unfold_ensemble(self$detector_names, self$n_energy_bins,
                              self$E_MeV, self$sensitivities, self$cc_icrp116,
                              function(out) self$save_result(out),
                              readings, ...)
        },

        #' @description Unfold using Cascade. See \code{\link{unfold_cascade}}.
        unfold_cascade = function(readings, ...) {
            unfold_cascade(self$detector_names, self$n_energy_bins,
                              self$E_MeV, self$sensitivities, self$cc_icrp116,
                              function(out) self$save_result(out),
                              readings, ...)
        },

        #' @description Unfold using Composite. See \code{\link{unfold_composite}}.
        unfold_composite = function(readings, ...) {
            unfold_composite(self$detector_names, self$n_energy_bins,
                              self$E_MeV, self$sensitivities, self$cc_icrp116,
                              function(out) self$save_result(out),
                              readings, ...)
        }
    ),
    private = list(
        .E_MeV = numeric(0L)
    )
)
