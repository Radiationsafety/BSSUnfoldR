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
        },

        # ------------------------------------------------------------------
        # Fifth batch of algorithms (added in v0.1.4)
        # ------------------------------------------------------------------

        #' @description Unfold using Scipy direct method. See \code{\link{unfold_scipy_direct_method}}.
        unfold_scipy_direct_method = function(readings, ...) {
            unfold_scipy_direct_method(self$detector_names, self$n_energy_bins,
                                          self$E_MeV, self$sensitivities,
                                          self$cc_icrp116,
                                          function(out) self$save_result(out),
                                          readings, ...)
        },

        #' @description Unfold using EKI. See \code{\link{unfold_eki}}.
        unfold_eki = function(readings, ...) {
            unfold_eki(self$detector_names, self$n_energy_bins,
                        self$E_MeV, self$sensitivities, self$cc_icrp116,
                        function(out) self$save_result(out),
                        readings, ...)
        },

        #' @description Unfold using RFSP-JUL. See \code{\link{unfold_rfsp_jul}}.
        unfold_rfsp_jul = function(readings, ...) {
            unfold_rfsp_jul(self$detector_names, self$n_energy_bins,
                              self$E_MeV, self$sensitivities, self$cc_icrp116,
                              function(out) self$save_result(out),
                              readings, ...)
        },

        #' @description Unfold using FRUIT-like. See \code{\link{unfold_fruit_like}}.
        unfold_fruit_like = function(readings, ...) {
            unfold_fruit_like(self$detector_names, self$n_energy_bins,
                                 self$E_MeV, self$sensitivities, self$cc_icrp116,
                                 function(out) self$save_result(out),
                                 readings, ...)
        },

        #' @description Unfold using AMAXED-Reg. See \code{\link{unfold_amaxed_regularization}}.
        unfold_amaxed_regularization = function(readings, ...) {
            unfold_amaxed_regularization(self$detector_names,
                                             self$n_energy_bins,
                                             self$E_MeV, self$sensitivities,
                                             self$cc_icrp116,
                                             function(out) self$save_result(out),
                                             readings, ...)
        },

        #' @description Unfold using CS. See \code{\link{unfold_cs}}.
        unfold_cs = function(readings, ...) {
            unfold_cs(self$detector_names, self$n_energy_bins,
                        self$E_MeV, self$sensitivities, self$cc_icrp116,
                        function(out) self$save_result(out),
                        readings, ...)
        },

        #' @description Unfold using Bayesian parametric. See \code{\link{unfold_bayesian_parametric}}.
        unfold_bayesian_parametric = function(readings, ...) {
            unfold_bayesian_parametric(self$detector_names,
                                          self$n_energy_bins,
                                          self$E_MeV, self$sensitivities,
                                          self$cc_icrp116,
                                          function(out) self$save_result(out),
                                          readings, ...)
        },

        #' @description Unfold using Hybrid parametric. See \code{\link{unfold_hybrid_parametric}}.
        unfold_hybrid_parametric = function(readings, ...) {
            unfold_hybrid_parametric(self$detector_names,
                                         self$n_energy_bins,
                                         self$E_MeV, self$sensitivities,
                                         self$cc_icrp116,
                                         function(out) self$save_result(out),
                                         readings, ...)
        },

        #' @description Unfold using Hybrid GMRES. See \code{\link{unfold_hybrid_gmres}}.
        unfold_hybrid_gmres = function(readings, ...) {
            unfold_hybrid_gmres(self$detector_names, self$n_energy_bins,
                                   self$E_MeV, self$sensitivities,
                                   self$cc_icrp116,
                                   function(out) self$save_result(out),
                                   readings, ...)
        },

        #' @description Unfold using Binned. See \code{\link{unfold_binned}}.
        unfold_binned = function(readings, ...) {
            unfold_binned(self$detector_names, self$n_energy_bins,
                             self$E_MeV, self$sensitivities, self$cc_icrp116,
                             function(out) self$save_result(out),
                             readings, ...)
        },

        # ------------------------------------------------------------------
        # Sixth batch of algorithms (added in v0.1.5)
        # ------------------------------------------------------------------

        #' @description Unfold using Parametric (FRUIT). See \code{\link{unfold_parametric}}.
        unfold_parametric = function(readings, ...) {
            unfold_parametric(self$detector_names, self$n_energy_bins,
                                 self$E_MeV, self$sensitivities, self$cc_icrp116,
                                 function(out) self$save_result(out),
                                 readings, ...)
        },

        #' @description Unfold using Parametric2 (BON95). See \code{\link{unfold_parametric2}}.
        unfold_parametric2 = function(readings, ...) {
            unfold_parametric2(self$detector_names, self$n_energy_bins,
                                  self$E_MeV, self$sensitivities, self$cc_icrp116,
                                  function(out) self$save_result(out),
                                  readings, ...)
        },

        #' @description Unfold using EPIC. See \code{\link{unfold_epic}}.
        unfold_epic = function(readings, ...) {
            unfold_epic(self$detector_names, self$n_energy_bins,
                          self$E_MeV, self$sensitivities, self$cc_icrp116,
                          function(out) self$save_result(out),
                          readings, ...)
        },

        #' @description Unfold using NN-KSVD. See \code{\link{unfold_nnksvd}}.
        unfold_nnksvd = function(readings, ...) {
            unfold_nnksvd(self$detector_names, self$n_energy_bins,
                             self$E_MeV, self$sensitivities, self$cc_icrp116,
                             function(out) self$save_result(out),
                             readings, ...)
        },

        #' @description Unfold using Genetic. See \code{\link{unfold_genetic}}.
        unfold_genetic = function(readings, ...) {
            unfold_genetic(self$detector_names, self$n_energy_bins,
                              self$E_MeV, self$sensitivities, self$cc_icrp116,
                              function(out) self$save_result(out),
                              readings, ...)
        },

        #' @description Unfold using Mystic. See \code{\link{unfold_mystic}}.
        unfold_mystic = function(readings, ...) {
            unfold_mystic(self$detector_names, self$n_energy_bins,
                              self$E_MeV, self$sensitivities, self$cc_icrp116,
                              function(out) self$save_result(out),
                              readings, ...)
        },

        #' @description Unfold using QUBO. See \code{\link{unfold_qubo}}.
        unfold_qubo = function(readings, ...) {
            unfold_qubo(self$detector_names, self$n_energy_bins,
                          self$E_MeV, self$sensitivities, self$cc_icrp116,
                          function(out) self$save_result(out),
                          readings, ...)
        },

        #' @description Unfold using LMfit. See \code{\link{unfold_lmfit}}.
        unfold_lmfit = function(readings, ...) {
            unfold_lmfit(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold using QPsolvers. See \code{\link{unfold_qpsolvers}}.
        unfold_qpsolvers = function(readings, ...) {
            unfold_qpsolvers(self$detector_names, self$n_energy_bins,
                                 self$E_MeV, self$sensitivities, self$cc_icrp116,
                                 function(out) self$save_result(out),
                                 readings, ...)
        },

        #' @description Unfold using CVXPY. See \code{\link{unfold_cvxpy}}.
        unfold_cvxpy = function(readings, ...) {
            unfold_cvxpy(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        # ------------------------------------------------------------------
        # Seventh batch of algorithms (added in v0.2.0)
        # ------------------------------------------------------------------

        #' @description Unfold by a named-method pipeline. See \code{\link{unfold_combined}}.
        unfold_combined = function(readings, ...) {
            unfold_combined(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold using P-spline REML. See \code{\link{unfold_pspline_reml}}.
        unfold_pspline_reml = function(readings, ...) {
            unfold_pspline_reml(self$detector_names, self$n_energy_bins,
                                self$E_MeV, self$sensitivities,
                                self$cc_icrp116,
                                function(out) self$save_result(out),
                                readings, ...)
        },

        #' @description Unfold with AMG/preconditioned Krylov. See \code{\link{unfold_amg}}.
        unfold_amg = function(readings, ...) {
            unfold_amg(self$detector_names, self$n_energy_bins,
                       self$E_MeV, self$sensitivities, self$cc_icrp116,
                       function(out) self$save_result(out),
                       readings, ...)
        },

        #' @description Unfold using SSR (sign-simplicity regression). See \code{\link{unfold_ssr}}.
        unfold_ssr = function(readings, ...) {
            unfold_ssr(self$detector_names, self$n_energy_bins,
                       self$E_MeV, self$sensitivities, self$cc_icrp116,
                       function(out) self$save_result(out),
                       readings, ...)
        },

        #' @description Unfold using the Mystic hybrid (DE + L-BFGS-B). See \code{\link{unfold_mystic_hybrid}}.
        unfold_mystic_hybrid = function(readings, ...) {
            unfold_mystic_hybrid(self$detector_names, self$n_energy_bins,
                                 self$E_MeV, self$sensitivities,
                                 self$cc_icrp116,
                                 function(out) self$save_result(out),
                                 readings, ...)
        },

        #' @description Unfold with the SMT-style exact solver. See \code{\link{unfold_smt}}.
        unfold_smt = function(readings, ...) {
            unfold_smt(self$detector_names, self$n_energy_bins,
                       self$E_MeV, self$sensitivities, self$cc_icrp116,
                       function(out) self$save_result(out),
                       readings, ...)
        },

        #' @description Unfold with the SCIP-style QP engine. See \code{\link{unfold_scip}}.
        unfold_scip = function(readings, ...) {
            unfold_scip(self$detector_names, self$n_energy_bins,
                        self$E_MeV, self$sensitivities, self$cc_icrp116,
                        function(out) self$save_result(out),
                        readings, ...)
        },

        #' @description Unfold with the CPLEX-style QP engine. See \code{\link{unfold_docplex}}.
        unfold_docplex = function(readings, ...) {
            unfold_docplex(self$detector_names, self$n_energy_bins,
                           self$E_MeV, self$sensitivities, self$cc_icrp116,
                           function(out) self$save_result(out),
                           readings, ...)
        },

        #' @description Unfold with Poisson-likelihood inference. See \code{\link{unfold_zfit}}.
        unfold_zfit = function(readings, ...) {
            unfold_zfit(self$detector_names, self$n_energy_bins,
                        self$E_MeV, self$sensitivities, self$cc_icrp116,
                        function(out) self$save_result(out),
                        readings, ...)
        },

        #' @description Unfold with MAEO (multi-algorithm evolution). See \code{\link{unfold_maeo}}.
        unfold_maeo = function(readings, ...) {
            unfold_maeo(self$detector_names, self$n_energy_bins,
                        self$E_MeV, self$sensitivities, self$cc_icrp116,
                        function(out) self$save_result(out),
                        readings, ...)
        },

        #' @description Unfold with operator-based MLEM. See \code{\link{unfold_mlem_odl}}.
        unfold_mlem_odl = function(readings, ...) {
            unfold_mlem_odl(self$detector_names, self$n_energy_bins,
                            self$E_MeV, self$sensitivities, self$cc_icrp116,
                            function(out) self$save_result(out),
                            readings, ...)
        },

        #' @description Unfold with PDHG (L2 + TV). See \code{\link{unfold_pdhg}}.
        unfold_pdhg = function(readings, ...) {
            unfold_pdhg(self$detector_names, self$n_energy_bins,
                        self$E_MeV, self$sensitivities, self$cc_icrp116,
                        function(out) self$save_result(out),
                        readings, ...)
        },

        #' @description Unfold with Douglas-Rachford splitting. See \code{\link{unfold_douglas_rachford}}.
        unfold_douglas_rachford = function(readings, ...) {
            unfold_douglas_rachford(self$detector_names,
                                    self$n_energy_bins,
                                    self$E_MeV, self$sensitivities,
                                    self$cc_icrp116,
                                    function(out) self$save_result(out),
                                    readings, ...)
        },

        #' @description Unfold and attach an interpretation report. See \code{\link{unfold_interpret}}.
        unfold_interpret = function(readings, ...) {
            unfold_interpret(self$detector_names, self$n_energy_bins,
                             self$E_MeV, self$sensitivities, self$cc_icrp116,
                             function(out) self$save_result(out),
                             readings, ...)
        },

        #' @description Unfold with the cvxpy-SQP parametric variant. See \code{\link{unfold_parametric_cvxpy}}.
        unfold_parametric_cvxpy = function(readings, ...) {
            unfold_parametric_cvxpy(self$detector_names, self$n_energy_bins,
                                    self$E_MeV, self$sensitivities,
                                    self$cc_icrp116,
                                    function(out) self$save_result(out),
                                    readings, ...)
        },

        #' @description Unfold with the qpsolvers-SQP parametric variant. See \code{\link{unfold_parametric_qpsolvers}}.
        unfold_parametric_qpsolvers = function(readings, ...) {
            unfold_parametric_qpsolvers(self$detector_names,
                                        self$n_energy_bins,
                                        self$E_MeV, self$sensitivities,
                                        self$cc_icrp116,
                                        function(out) self$save_result(out),
                                        readings, ...)
        },

        #' @description Unfold with the combined parametric variant. See \code{\link{unfold_parametric_combined}}.
        unfold_parametric_combined = function(readings, ...) {
            unfold_parametric_combined(self$detector_names,
                                       self$n_energy_bins,
                                       self$E_MeV, self$sensitivities,
                                       self$cc_icrp116,
                                       function(out) self$save_result(out),
                                       readings, ...)
        },

        # ------------------------------------------------------------------
        # Result management / dose coefficients / comparison
        # ------------------------------------------------------------------

        #' @description Retrieve a saved result by history index.
        #' @param index Integer index into \code{self\$history} (1 = oldest).
        #'   Default \code{NULL} = the most recent entry.
        get_result = function(index = NULL) {
            if (length(self$history) == 0L) return(NULL)
            if (is.null(index)) index <- length(self$history)
            index <- as.integer(index)
            if (index <= 0L || index > length(self$history)) return(NULL)
            self$history[[index]]
        },

        #' @description List descriptions of all saved results
        #'   (method and residual norm per entry).
        list_results = function() {
            out <- lapply(self$history, function(h)
                list(method = if (!is.null(h$method)) h$method else "?",
                     residual_norm = if (!is.null(h$residual_norm))
                         h$residual_norm else NA_real_))
            out
        },

        #' @description Clear the saved-result history.
        clear_results = function() {
            self$history <- list()
            invisible(self)
        },

        #' @description Switch the dose conversion coefficient dataset
#'   after construction (Python \code{Detector.set_dose_coefficients}).
        #' @param cc_type One of \code{"ICRP116"},
        #'   \code{"ICRP74_effective"}, \code{"ICRP74_operational"},
        #'   \code{"NRB99_2009_effective"} or a named list of coefficients
        #'   with an \code{E_MeV} entry.
        set_dose_coefficients = function(cc_type) {
            if (is.list(cc_type)) {
                cc <- cc_type
            } else {
                key <- tolower(cc_type)
                factory <- switch(key,
                    icrp116 = ICRP116_COEFF_EFFECTIVE_DOSE,
                    icrp74_effective = ICRP74_COEFF_EFFECTIVE_DOSE,
                    icrp74_operational =
                        ICRP74_COEFF_OPERATIONAL_QUANTITIES,
                    nrb99_2009_effective =
                        NRB99_2009_COEFF_EFFECTIVE_DOSE,
                    stop("Unknown dose coefficient set: ", cc_type))
                cc <- factory()
            }
            self$cc_icrp116 <- interpolate_coefficients(cc, self$E_MeV)
            invisible(self)
        },

        #' @description Compare two unfolding results (or a result and a
        #'   reference spectrum) with standard metrics.
        #' @param reference Either a result list with a \code{spectrum}
        #'   entry, or an explicit numeric spectrum of length
        #'   \code{n_energy_bins}.
        #' @return A list with relative-integral difference, chi-squared,
        #'   ratio statistics, and dose-rate relative difference if
        #'   available.
        compare = function(reference) {
            ref_spec <- if (is.list(reference) &&
                                !is.null(reference$spectrum))
                            as.numeric(reference$spectrum)
                        else as.numeric(reference)
            if (length(ref_spec) != self$n_energy_bins) {
                stop("Reference spectrum length mismatch.")
            }
            last <- if (length(self$history)) {
                self$history[[length(self$history)]]
            } else {
                stop("No unfolded result stored in history to compare.")
            }
            xy <- last$spectrum
            integrate_rel <- if (length(ref_spec) == length(xy))
                abs(sum(xy) - sum(ref_spec)) / max(abs(sum(ref_spec)),
                    1e-30)
            chi2 <- sum((xy - ref_spec)^2 / pmax(ref_spec^2, 1e-30))
            list(method = last$method,
                 relative_integral_difference = integrate_rel,
                 chi_squared = chi2,
                 max_relative_deviation =
                     max(abs(xy - ref_spec) / pmax(abs(ref_spec), 1e-30)),
                 reference = ref_spec)
        },

        # --- Plot methods (ggplot2) -----------------------------------------

        #' @description Plot the response (sensitivity) functions;
        #'   Python \code{Detector.plot_response_functions}.
        #' \code{requireNamespace("ggplot2")} is required.
        plot_response_functions = function() {
            require_ggplot2()
            df <- do.call(rbind, lapply(names(self$sensitivities), function(k)
                data.frame(E_MeV = self$E_MeV,
                           sensitivity = self$sensitivities[[k]],
                           detector = k)))
            ggplot2::ggplot(df, ggplot2::aes(x = .data$E_MeV,
                    y = .data$sensitivity, color = .data$detector)) +
                ggplot2::geom_line() + ggplot2::scale_x_log10() +
                ggplot2::scale_y_log10() +
                ggplot2::labs(x = "Energy [MeV]", y = "Response",
                              title = "Bonner sphere response functions")
        },

        #' @description Plot the most recent unfolded spectrum with
        #'   uncertainty bands (Monte-Carlo when available); Python
        #'   \code{Detector.plot_with_uncertainty}.
        #' @param result Optional result list; defaults to the latest in
        #'   \code{self\$history}.
        plot_with_uncertainty = function(result = NULL) {
            require_ggplot2()
            if (is.null(result)) {
                if (!length(self$history)) {
                    stop("No unfolded result to plot; unfold first.")
                }
                result <- self$history[[length(self$history)]]
            }
            e <- result$energy; s <- result$spectrum
            if (!is.null(result$montecarlo_std)) {
                lo <- pmax(result$spectrum - 2 * result$montecarlo_std, 0)
                up <- result$spectrum + 2 * result$montecarlo_std
            } else if (!is.null(result$montecarlo_samples)) {
                ms <- result$montecarlo_samples
                if (!is.null(result$montecarlo_mean)) {
                    sd <- result$montecarlo_std
                    lo <- pmax(result$spectrum - 2 * sd, 0)
                    up <- result$spectrum + 2 * sd
                } else {
                    lo <- result$spectrum; up <- result$spectrum
                }
            } else {
                lo <- result$spectrum; up <- result$spectrum
            }
            if (is.null(lo)) lo <- result$spectrum
            df <- data.frame(E_MeV = e, spectrum = s, lo = lo, up = up)
            G <- ggplot2::ggplot(df, ggplot2::aes(x = .data$E_MeV,
                    y = .data$spectrum)) +
                ggplot2::geom_point() + ggplot2::geom_line() +
                ggplot2::geom_ribbon(
                    ggplot2::aes(ymin = .data$lo, ymax = .data$up),
                    alpha = 0.3) +
                ggplot2::scale_x_log10() + ggplot2::scale_y_log10(
                    limits = c(NA, NA)) +
                ggplot2::labs(x = "Energy [MeV]",
                              y = "Fluence per unit lethargy",
                              title = "Unfolded neutron spectrum")
            G
        }
    ),
    active = list(
        #' @field n_detectors Integer; number of entry response functions
        #'   currently registered.
        n_detectors = function() length(self$detector_names)
    ),
    private = list(
        .E_MeV = numeric(0L)
    )
)
