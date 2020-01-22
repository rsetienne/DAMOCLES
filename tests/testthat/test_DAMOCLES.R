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
    optimmethod = 'subplex')
})

test_that("DAMOCLES_bootstrap", {
  data(NWPrimates_data)
  out = DAMOCLES_bootstrap(
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
})
