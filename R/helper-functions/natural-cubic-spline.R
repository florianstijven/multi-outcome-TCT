# ============================================================================
# Natural Cubic Spline Interpolation Reference Trajectory Functions
# ============================================================================

build_nc_spline_basis <- function(times, knots, boundary_knots, params = NULL, ...) {
  if (!is.numeric(times)) {
    stop("times must be a numeric vector.")
  }
  if (length(times) == 0) {
    stop("times must not be empty.")
  }
  
  if (!is.numeric(boundary_knots) || length(boundary_knots) != 2) {
    stop("boundary knots must be a numeric vector of length 2.")
  }
  
  basis <- splines::ns(
    x = times,
    knots = knots,
    Boundary.knots = boundary_knots,
    intercept = TRUE
  )
  
  if (!is.null(params) && length(params) != ncol(basis)) {
    stop(
      "Length of params (",
      length(params),
      ") does not match spline basis dimension (",
      ncol(basis),
      ")."
    )
  }
  
  basis
}

time_d_nc_spline <- function(times, knots, boundary_knots, params) {
  if (!is.numeric(params)) {
    stop("params must be a numeric vector.")
  }
  if (!is.numeric(times)) {
    stop("times must be a numeric vector.")
  }
  
  # Central finite difference for d/dt f(t), with adaptive step size.
  eps <- sqrt(.Machine$double.eps) * pmax(1, abs(times))
  f_plus <- function_nc_spline(
    times + eps,
    knots = knots,
    boundary_knots = boundary_knots,
    params = params
  )
  f_minus <- function_nc_spline(
    times - eps,
    knots = knots,
    boundary_knots = boundary_knots,
    params = params
  )
  
  as.numeric((f_plus - f_minus) / (2 * eps))
}

function_nc_spline <- function(times, params, knots, boundary_knots) {
  if (!is.numeric(params)) {
    stop("params must be a numeric vector.")
  }
  basis <- build_nc_spline_basis(
    times = times,
    knots = knots,
    boundary_knots = boundary_knots,
    params = params
  )
  as.numeric(basis %*% params)
}

jacobian_nc_spline <- function(times, knots, boundary_knots, params) {
  if (!is.numeric(params)) {
    stop("params must be a numeric vector.")
  }
  # The model is linear in coefficients, so the Jacobian is the spline basis.
  build_nc_spline_basis(
    times = times,
    knots = knots,
    boundary_knots = boundary_knots,
    params = params
  )
}