# ============================================================================
# 4PL Reference Trajectory Functions (fixed asymptotes at 0, 1)
# ============================================================================


# Derivative of the 4PL curve with respect to time, given the slope and inflection point parameters.
time_d_4PL <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  
  exp_term <- exp(-slope * (times - inflection))
  denom <- (1 + exp_term)^2
  slope * exp_term / denom
}

function_4PL <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  
  # Compute the 4PL model values given times, slope, and inflection point.
  1 / (1 + exp(-slope * (times - inflection)))
}

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

# Derivative of the 4PL curve with respect to time, given the slope and inflection point parameters.
time_d_4PL_asympt <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  upper_asymptote <- params[3]
  lower_asymptote <- params[4]
  
  exp_term <- exp(-slope * (times - inflection))
  denom <- (1 + exp_term)^2
  (upper_asymptote - lower_asymptote) * slope * exp_term / denom
}

function_4PL_asympt <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  upper_asymptote <- params[3]
  lower_asymptote <- params[4]
  
  # Compute the 4PL model values given times, slope, and inflection point.
  (upper_asymptote - lower_asymptote) / (1 + exp(-slope * (times - inflection)))
}

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
