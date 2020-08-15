#' @examples
#' 
#' context("DAMOCLES")

test_that("DAMOCLES_ML", {
  data(NWPrimates_data)
  out <- DAMOCLES_ML(
    phy = NWPrimates_data[[1]],
    pa = NWPrimates_data[[2]],
    initparsopt = c(0.01,1.8),
    idparsopt = c(1,2),
    parsfix = NULL,
    idparsfix = NULL,
    pars2 = c(1E-3,1E-4,1E-5,1000),
    pchoice = 0,
    methode = 'analytical',
    optimmethod = 'subplex')
  testthat::expect_equal(out$loglik, -36.87957, tolerance = 1E-5)
  out2 <- DAMOCLES_ML(
    phy = NWPrimates_data[[1]],
    pa = NWPrimates_data[[2]],
    initparsopt = c(0.01,1.8),
    idparsopt = c(1,2),
    parsfix = NULL,
    idparsfix = NULL,
    pars2 = c(1E-3,1E-4,1E-5,1000),
    pchoice = 0,
    methode = 'Matrix',
    optimmethod = 'subplex')
  testthat::expect_equal(out,out2,tolerance = 1E-5)
  out3 <- DAMOCLES_ML(
    phy = NWPrimates_data[[1]],
    pa = NWPrimates_data[[2]],
    initparsopt = c(0.01,1.8),
    idparsopt = c(1,2),
    parsfix = NULL,
    idparsfix = NULL,
    pars2 = c(1E-3,1E-4,1E-5,1000),
    pchoice = 0,
    methode = 'expm',
    optimmethod = 'subplex')
  testthat::expect_equal(out,out3,tolerance = 1E-5)
})

test_that("DAMOCLES_bootstrap", {
  data(NWPrimates_data)
  set.seed(42)
  out <- DAMOCLES_bootstrap(
    phy = NWPrimates_data[[1]],
    pa = NWPrimates_data[[2]],
    initparsopt = c(0.01,1.8),
    idparsopt = c(1,2),
    idparsfix = NULL,
    parsfix = NULL,
    pars2 = c(1E-3,1E-4,1E-5,1000),
    pchoice = 1,
    runs = 2,
    estimate_pars = TRUE,
    conf.int = 0.95)
  outcome1 <- c(-38.03012,-38.15837)
  outcome2 <- c(-37.35113,-31.88564)
  expected_outcome <- out$null_community_data$loglik.DAMOCLES +
                      (all(abs(out$null_community_data$loglik.DAMOCLES - outcome1) < 1E-5)) * (outcome1 - out$null_community_data$loglik.DAMOCLES) +
                      (all(abs(out$null_community_data$loglik.DAMOCLES - outcome2) < 1E-5)) * (outcome2 - out$null_community_data$loglik.DAMOCLES)
  testthat::expect_equal(out$null_community_data$loglik.DAMOCLES,
                         expected_outcome,
                         tolerance = 1E-5)
})