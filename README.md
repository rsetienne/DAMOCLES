# DAMOCLES

<!-- badges: start -->

![CRAN/METACRAN Version](https://img.shields.io/cran/v/DAMOCLES) ![GitHub R package version](https://img.shields.io/github/r-package/v/rsetienne/DAMOCLES) [![](http://cranlogs.r-pkg.org/badges/grand-total/DAMOCLES)](https://CRAN.R-project.org/package=DAMOCLES) [![](http://cranlogs.r-pkg.org/badges/DAMOCLES)](https://CRAN.R-project.org/package=DAMOCLES)

<!-- badges: end -->

| Branch | [![GitHub Actions logo](man/figures/github_actions_logo.png)](https://github.com/features/actions) | [![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io) |
|----|----|----|
| `master` | [![R-CMD-check](https://github.com/rsetienne/DAMOCLES/actions/workflows/R-CMD-check.yaml/badge.svg?branch=master)](https://github.com/rsetienne/DAMOCLES/actions/workflows/R-CMD-check.yaml) | [![Codecov test coverage](https://codecov.io/gh/rsetienne/damocles/branch/master/graph/badge.svg)](https://app.codecov.io/gh/rsetienne/damocles?branch=master) |
| `develop` | [![R-CMD-check](https://github.com/rsetienne/DAMOCLES/actions/workflows/R-CMD-check.yaml/badge.svg?branch=develop)](https://github.com/rsetienne/DAMOCLES/actions/workflows/R-CMD-check.yaml) | [![Codecov test coverage](https://codecov.io/gh/rsetienne/damocles/branch/develop/graph/badge.svg)](https://app.codecov.io/gh/rsetienne/damocles?branch=develop) |

Dynamic Assembly Model of Colonization, Local Extinction and Speciation in `R`.

This is a development version before the official release on CRAN.

## Installing DAMOCLES

The DAMOCLES package has a stable version on CRAN and a development version on GitHub.

### From CRAN

From within R, do:

```         
install.packages("DAMOCLES")
```

### From GitHub

Install DAMOCLES from this GitHub repository by running:

``` r
install.packages("remotes")
remotes::install_github("rsetienne/DAMOCLES")
```

## Using DAMOCLES as a package dependency

### From CRAN

To your DESCRIPTION file, add `DAMOCLES` as any normal package.

If your package directly uses `DAMOCLES`:

```         
Imports:
  DAMOCLES
```

If your package uses `DAMOCLES` in its peripherals (e.g. vignettes and tests):

```         
Suggests:
  DAMOCLES
```

### From GitHub

```         
Remotes:
  rsetienne/DAMOCLES
```

## `git` branching model

-   `master`: build should always pass. [\@rsetienne](https://github.com/rsetienne) has control over `develop` to `master` merges.
-   `develop`: merge of topic branches, merge with `master` by [\@rsetienne](https://github.com/rsetienne) iff build passes. \## Contributors

DAMOCLES was originally developed by Rampal S. Etienne and Alex Pigot

Additionally there are others working on expanding DAMOCLES at the [TECE lab](https://github.com/tece-lab), University of Groningen.

## References

Pigot, A. L., & Etienne, R. S. (2015). A new dynamic null model for phylogenetic community structure. Ecology Letters, 18(2), 153-163. <https://doi.org/10.1111/ele.12395>
