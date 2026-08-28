# ============================================================================
# 4PL Reference Trajectory Functions (fixed asymptotes at 0, 1)
# ============================================================================

#' Time-derivative of the 4PL reference trajectory
#'
#' Derivative with respect to time of the 4-parameter logistic (4PL) curve,
#' with asymptotes fixed at 0 and 1.
#'
#' @param times numeric vector of time points.
#' @param params numeric vector of length 2: `c(inflection, slope)`.
#'
#' @returns numeric vector of derivative values, one per element of `times`.
#' @export
time_d_4PL <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  
  exp_term <- exp(-slope * (times - inflection))
  denom <- (1 + exp_term)^2
  slope * exp_term / denom
}

#' 4PL reference trajectory
#'
#' Evaluate the 4-parameter logistic (4PL) curve, with asymptotes fixed at 0
#' and 1.
#'
#' @inheritParams time_d_4PL
#'
#' @returns numeric vector of 4PL values, one per element of `times`.
#' @export
function_4PL <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  
  # Compute the 4PL model values given times, slope, and inflection point.
  1 / (1 + exp(-slope * (times - inflection)))
}

#' Jacobian of the 4PL reference trajectory
#'
#' @inheritParams time_d_4PL
#'
#' @returns numeric matrix with `length(times)` rows and 2 columns (partial
#'   derivatives with respect to `inflection` and `slope`).
#' @export
jacobian_4PL <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  
  # Compute the Jacobian of the 4PL model with respect to the slope and inflection point parameters.
  exp_term <- exp(-slope * (times - inflection))
  denom <- (1 + exp_term)^2
  d_inflection <- -slope * exp_term / denom
  d_slope <- (times - inflection) * exp_term / denom
  cbind(d_inflection, d_slope)
}

# ============================================================================
# 4PL Reference Trajectory Functions (asymptotes as parameters)
# ============================================================================

#' Time-derivative of the 4PL reference trajectory with free asymptotes
#'
#' Derivative with respect to time of the 4-parameter logistic (4PL) curve,
#' where the upper and lower asymptotes are also estimated parameters.
#'
#' @param times numeric vector of time points.
#' @param params numeric vector of length 4:
#'   `c(inflection, slope, upper_asymptote, lower_asymptote)`.
#'
#' @returns numeric vector of derivative values, one per element of `times`.
#' @export
time_d_4PL_asympt <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  upper_asymptote <- params[3]
  lower_asymptote <- params[4]
  
  exp_term <- exp(-slope * (times - inflection))
  denom <- (1 + exp_term)^2
  (upper_asymptote - lower_asymptote) * slope * exp_term / denom
}

#' 4PL reference trajectory with free asymptotes
#'
#' @inheritParams time_d_4PL_asympt
#'
#' @returns numeric vector of 4PL values, one per element of `times`.
#' @export
function_4PL_asympt <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  upper_asymptote <- params[3]
  lower_asymptote <- params[4]
  
  # Compute the 4PL model values given times, slope, and inflection point.
  (upper_asymptote - lower_asymptote) / (1 + exp(-slope * (times - inflection)))
}

#' Jacobian of the 4PL reference trajectory with free asymptotes
#'
#' @inheritParams time_d_4PL_asympt
#'
#' @returns numeric matrix with `length(times)` rows and 4 columns (partial
#'   derivatives with respect to `inflection`, `slope`, `upper_asymptote`, and
#'   `lower_asymptote`).
#' @export
jacobian_4PL_asympt <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  upper_asymptote <- params[3]
  lower_asymptote <- params[4]
  
  # Compute the Jacobian of the 4PL model with respect to the slope and
  # inflection point parameters.
  exp_term <- exp(-slope * (times - inflection))
  denom <- (1 + exp_term)^2
  d_inflection <- -slope * exp_term / denom
  d_slope <- (times - inflection) * exp_term / denom
  cbind(d_inflection, d_slope, 1 / (1 + exp(-slope * (times - inflection))), -1 / (1 + exp(-slope * (times - inflection))))
}
