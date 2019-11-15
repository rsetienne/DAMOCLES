\name{DAMOCLES-package}
\alias{DAMOCLES-package}
\alias{DAMOCLES}
\docType{package}
\title{Dynamic Assembly Model Of Colonization, Local Extinction and Speciation}
\description{ Simulation and likelihood methods for a dynamical community assembly that accounts for phylogenetic history.\cr\cr
New in version 1.1:\cr
- Function added to do bootstrap\cr\cr
New in version 2.0:\cr
- Computes likelihood for a model where colonization and local extinction rates depend on binary trait value.\cr
- Computes likelihood for a model where colonization and local extinction rates depend on geographical trait.\cr
- Allows use of subplex as optimization algorithm.\cr
- Uses matrix exponentiation to solve ODE system.\cr\cr}
\details{                   
\tabular{ll}{
Package: \tab DAMOCLES\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2016-06-10\cr
License: \tab GPL 2.0\cr
}
DAMOCLES_loglik computes the likelihood of presence-absence data in a local community given a set of parameters and a phylogeny under an immigration-extinction model

DAMOCLES_ML finds the parameters that maximizes the likelihood computed by DAMOCLES_loglik.

DAMOCLES_sim simulates presence-absence data for a given phylogeny

DAMOCLES_bootstrap computes the maximum likelihood estimates of colonisation and local extinction rate for a given phylogeny and presence-absence data under the DAMOCLES model.
   These rate estimates are used to simulate null communities under the DAMOCLES model.
   Standardized effect size of mean nearest taxon distance (mntd), mean phylogentic distance (mpd) and loglikelihood are calculated. 
   For comparison, standardised effect sizes are also calculated relative to a "Random-Draw" null model i.e. presence absence randomised across tips.  
}

\author{Rampal S. Etienne & Alex L. Pigot\cr
Maintainer: Rampal S. Etienne (r.s.etienne@rug.nl)
}
\references{
Pigot, A.L. & R.S. Etienne (2015). A new dynamic null model for phylogenetic community structure. Ecology Letters 18: 153-163.}
\keyword{ models }
\examples{
DAMOCLES_ML()
}
