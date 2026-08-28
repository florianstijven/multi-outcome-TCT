# ============================================================================
# Natural Cubic Spline Interpolation Reference Trajectory Functions
# ============================================================================

#' Natural cubic spline basis matrix
#'
#' Builds a natural cubic spline basis via [splines::ns()], with an
#' intercept column, at the given `knots` and `boundary_knots`.
#'
#' @param times numeric vector of time points at which to evaluate the basis.
#' @param knots numeric vector of interior knot locations.
#' @param boundary_knots numeric vector of length 2, the boundary knots.
#' @param params optional numeric vector of spline coefficients, used only to
#'   validate that its length matches the basis dimension.
#' @param ... unused; retained for interface consistency.
#'
#' @returns numeric matrix, the spline basis evaluated at `times`.
#' @export
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

#' Time-derivative of the natural cubic spline reference trajectory
#'
#' Computed via a central finite difference with an adaptive step size, since
#' the spline basis derivative is not evaluated in closed form here.
#'
#' @inheritParams build_nc_spline_basis
#' @param params numeric vector of spline coefficients.
#'
#' @returns numeric vector of derivative values, one per element of `times`.
#' @export
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

#' Natural cubic spline reference trajectory
#'
#' @inheritParams time_d_nc_spline
#'
#' @returns numeric vector of fitted values, one per element of `times`.
#' @export
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

#' Jacobian of the natural cubic spline reference trajectory
#'
#' The model is linear in coefficients, so the Jacobian is simply the spline
#' basis matrix.
#'
#' @inheritParams time_d_nc_spline
#'
#' @returns numeric matrix, the spline basis (see [build_nc_spline_basis()]).
#' @export
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
