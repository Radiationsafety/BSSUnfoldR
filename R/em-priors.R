#' Nearest-neighbour priors for penalised EM unfolding
#'
#' R port of \code{bssunfold/core/_em_priors.py}. The prior reads
#' \deqn{V(f) = beta * sum_r sum_{s in NN(r)} w_{r,s} * phi0(f_r, f_s)}
#' and the gradient used by the one-step-late EM updates is
#' \deqn{grad_r V(f) = beta * sum_{s in NN(r)} w_{r,s} * phi1(f_r, f_s)}
#' with unit weights \eqn{w_{r,s} = 1} along the energy axis (zero-padded at
#' the boundaries). Three priors are supported: \code{quadratic},
#' \code{logcosh} and \code{relative_difference}.
#'
#' @param x Numeric spectrum vector (length n).
#' @param prior Character scalar: \code{"quadratic"}, \code{"logcosh"} or
#'   \code{"relative_difference"}.
#' @param beta Numeric prior weight. Default 1e-3.
#' @param delta Numeric width parameter / additive floor. Default 1.0.
#' @param gamma Numeric edge-preservation parameter for the relative-difference
#'   prior. Default 1.0.
#'
#' @return \code{prior_gradient} returns a numeric vector of length n;
#'   \code{prior_value} returns a numeric scalar.
#' @name em-priors
#' @rdname em-priors
NULL

.em_priors_neighbours <- function(x) {
    n <- length(x)
    left <- numeric(n)
    right <- numeric(n)
    if (n > 1L) {
        left[2:n]   <- x[1:(n - 1L)]
        right[1:(n - 1L)] <- x[2:n]
    }
    list(left = left, right = right)
}

.em_priors_phi1 <- function(fr, fs, prior, delta, gamma) {
    if (prior == "quadratic") {
        return((fr - fs) / delta)
    }
    if (prior == "logcosh") {
        return(tanh((fr - fs) / delta))
    }
    if (prior == "relative_difference") {
        absd <- abs(fr - fs)
        denom <- gamma * absd + fr + fs + delta
        (fr - fs) * (gamma * absd + 3.0 * fs + fr + 2.0 * delta) / (denom^2)
    } else {
        stop("Unknown prior '", prior,
             "'. Choose from 'quadratic', 'logcosh', 'relative_difference'.")
    }
}

.em_priors_phi0 <- function(fr, fs, prior, delta, gamma) {
    if (prior == "quadratic") {
        return(0.25 * ((fr - fs) / delta)^2)
    }
    if (prior == "logcosh") {
        t <- (fr - fs) / delta
        t_abs <- abs(t)
        return(t_abs + log1p(exp(-2.0 * t_abs)) - log(2.0))
    }
    if (prior == "relative_difference") {
        absd <- abs(fr - fs)
        return((fr - fs)^2 / (fr + fs + gamma * absd + delta))
    }
    stop("Unknown prior '", prior,
         "'. Choose from 'quadratic', 'logcosh', 'relative_difference'.")
}

#' @rdname em-priors
#' @export
prior_gradient <- function(x, prior = "quadratic", beta = 1e-3,
                           delta = 1.0, gamma = 1.0) {
    prior <- tolower(as.character(prior))
    if (!(prior %in% c("quadratic", "logcosh", "relative_difference"))) {
        stop("Unknown prior '", prior,
             "'. Choose from 'quadratic', 'logcosh', 'relative_difference'.")
    }
    x <- as.numeric(x)
    nb <- .em_priors_neighbours(x)
    grad <- .em_priors_phi1(x, nb$left, prior, delta, gamma) +
             .em_priors_phi1(x, nb$right, prior, delta, gamma)
    beta * grad
}

#' @rdname em-priors
#' @export
prior_value <- function(x, prior = "quadratic", beta = 1e-3,
                        delta = 1.0, gamma = 1.0) {
    prior <- tolower(as.character(prior))
    if (!(prior %in% c("quadratic", "logcosh", "relative_difference"))) {
        stop("Unknown prior '", prior,
             "'. Choose from 'quadratic', 'logcosh', 'relative_difference'.")
    }
    x <- as.numeric(x)
    nb <- .em_priors_neighbours(x)
    val <- .em_priors_phi0(x, nb$left, prior, delta, gamma) +
           .em_priors_phi0(x, nb$right, prior, delta, gamma)
    as.numeric(beta * sum(val))
}
