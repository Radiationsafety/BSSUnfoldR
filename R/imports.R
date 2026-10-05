#' Internal package hooks.
#'
#' Registers the stats import used across the solvers and silences R CMD
#' check NOTEs for ggplot2's \code{.data} pronoun, which is only available
#' when the Suggests-only ggplot2 package is loaded.
#' @importFrom stats runif
NULL

utils::globalVariables(".data")
