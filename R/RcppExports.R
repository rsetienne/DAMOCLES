# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @useDynLib DAMOCLES
NULL

DAMOCLES_integrate_odeint <- function(ry, times, M, atol, rtol, stepper) {
    .Call('_DAMOCLES_DAMOCLES_integrate_odeint', PACKAGE = 'DAMOCLES', ry, times, M, atol, rtol, stepper)
}

