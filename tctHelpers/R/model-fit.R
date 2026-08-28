# ============================================================================
# `model_fit` S3 class: the result of fitting a working model to first-stage
# estimates via (possibly stratified) generalized least squares.
# ============================================================================

#' Construct a model_fit object
#'
#' Internal constructor for the result of fitting a working model to
#' first-stage estimates via generalized least squares.
#'
#' @param coefficients named or unnamed numeric vector of estimated
#'   parameters.
#' @param vcov numeric matrix, the estimated covariance matrix of
#'   `coefficients`.
#' @param criterion numeric, the value of the GLS criterion function at the
#'   optimum, summed across strata.
#' @param optim list of raw [stats::optim()] outputs, one per stratum.
#' @param working_model an object of class `model`.
#' @param null_model logical, whether the fit is under the null hypothesis of
#'   no treatment effect.
#' @param data list with elements `m_tilde` and `Sigma`, the first-stage
#'   estimates used to fit the model.
#'
#' @returns an object of class `model_fit`.
#' @keywords internal
new_model_fit <- function(coefficients,
                          vcov,
                          criterion,
                          optim,
                          working_model,
                          null_model,
                          data) {
  structure(
    list(
      coefficients = coefficients,
      vcov = vcov,
      criterion = criterion,
      optim = optim,
      working_model = working_model,
      null_model = null_model,
      data = data
    ),
    class = "model_fit"
  )
}

#' Sandwich covariance matrix for a (possibly working) GLS fit
#'
#' Computes the M-estimation sandwich covariance matrix
#' `(J' W J)^+ (J' W Sigma W J) (J' W J)^+`, with `W = Sigma_working^{-1}`
#' and `^+` the Moore-Penrose pseudo-inverse (used because the Jacobian can be
#' rank-deficient, e.g. for saturated working models).
#' When `Sigma_working` equals `Sigma` (proper, not working-independence, GLS),
#' this reduces to the usual efficient GLS covariance `(J' Sigma^{-1} J)^+`.
#'
#' @param jacobian numeric matrix, the Jacobian of the mean function
#'   evaluated at the estimated parameters.
#' @param Sigma numeric matrix, the (true) first-stage covariance estimate.
#' @param Sigma_working numeric matrix, the covariance matrix actually used
#'   as the GLS weight (the identity matrix when `ols = TRUE`).
#'
#' @returns numeric matrix, the estimated covariance matrix of the parameter
#'   estimates.
#' @keywords internal
gls_vcov <- function(jacobian, Sigma, Sigma_working) {
  jacobian <- as.matrix(jacobian)
  W <- MASS::ginv(Sigma_working)
  bread <- MASS::ginv(t(jacobian) %*% W %*% jacobian)
  meat <- t(jacobian) %*% W %*% Sigma %*% W %*% jacobian
  bread %*% meat %*% bread
}

#' Estimated coefficients from a model fit
#'
#' @param object an object of class `model_fit`.
#' @param ... unused; retained for S3 method consistency.
#'
#' @returns numeric vector of estimated parameters.
#' @export
coef.model_fit <- function(object, ...) {
  object$coefficients
}

#' Estimated covariance matrix from a model fit
#'
#' @inheritParams coef.model_fit
#'
#' @returns numeric matrix, the estimated covariance matrix of the
#'   coefficients.
#' @export
vcov.model_fit <- function(object, ...) {
  object$vcov
}

#' Print a model fit
#'
#' @param x an object of class `model_fit`.
#' @param ... unused; retained for S3 method consistency.
#'
#' @returns `x`, invisibly.
#' @export
print.model_fit <- function(x, ...) {
  cat(if (x$null_model) "Model fit (null model)\n" else "Model fit (full model)\n")
  cat("Coefficients:\n")
  print(stats::coef(x))
  cat("\nGLS criterion value:", x$criterion, "\n")
  invisible(x)
}

#' Summarize a model fit
#'
#' Computes standard errors and Wald z-statistics/p-values for each
#' coefficient, based on [vcov.model_fit()].
#'
#' @inheritParams coef.model_fit
#'
#' @returns an object of class `summary.model_fit`, a list with elements
#'   `coefficients` (a matrix with columns `Estimate`, `Std. Error`,
#'   `z value`, and `Pr(>|z|)`), `criterion`, and `null_model`.
#' @export
summary.model_fit <- function(object, ...) {
  est <- stats::coef(object)
  se <- sqrt(diag(stats::vcov(object)))
  z_value <- est / se
  p_value <- 2 * stats::pnorm(-abs(z_value))

  coefficients <- cbind(
    Estimate = est,
    `Std. Error` = se,
    `z value` = z_value,
    `Pr(>|z|)` = p_value
  )

  structure(
    list(
      coefficients = coefficients,
      criterion = object$criterion,
      null_model = object$null_model
    ),
    class = "summary.model_fit"
  )
}

#' Print a model fit summary
#'
#' @param x an object of class `summary.model_fit`.
#' @param ... unused; retained for S3 method consistency.
#'
#' @returns `x`, invisibly.
#' @export
print.summary.model_fit <- function(x, ...) {
  cat(if (x$null_model) "Model fit (null model)\n" else "Model fit (full model)\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nGLS criterion value:", x$criterion, "\n")
  invisible(x)
}
