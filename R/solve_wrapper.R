#' Build a closure wrapper around a solver function
#'
#' \code{make_solve_wrapper} returns a function that accepts
#' \code{(A, b, ..., x0)} and forwards them to \code{solve_func}, merging
#' in the supplied \code{fixed_params}. It is the R analogue of the Python
#' \code{bssunfold/core/_base_unfolder.py:make_solve_wrapper}.
#'
#' @param solve_func Function with signature \code{solve_func(A, b, x0, ...)}
#'   returning a numeric spectrum or a list
#'   \code{list(spectrum, iterations, converged)}.
#' @param ... Extra named arguments that will always be forwarded to
#'   \code{solve_func}.
#' @return A function with signature \code{(A, b, x0 = NULL, ...)}.
#' @export
#' @keywords internal
#' @examples
#' toy <- function(A, b, x0 = NULL, tol = 1e-6) {
#'     list(spectrum = as.numeric(solve(A, b)), iterations = 1L, converged = TRUE)
#' }
#' wrapped <- make_solve_wrapper(toy, tol = 1e-3)
#' A <- matrix(c(1.0, 0.5, 0.2, 0.4), nrow = 2)
#' b <- c(1, 0.6)
#' wrapped(A, b, x0 = c(0.5, 0.5))
make_solve_wrapper <- function(solve_func, ...) {
    fixed <- list(...)
    function(A, b, x0 = NULL, ...) {
        do.call(solve_func, c(list(A = A, b = b, x0 = x0),
                              fixed, list(...)))
    }
}
