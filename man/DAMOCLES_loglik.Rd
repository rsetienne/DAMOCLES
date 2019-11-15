\name{DAMOCLES_loglik}
\alias{DAMOCLES_loglik}
\title{Likelihood for DAMOCLES model}
\description{
Computes likelihood for the presence-absence data of species in a local community for a given phylogeny of species in the region.
}
\usage{DAMOCLES_loglik(
   phy,
   pa,
   pars,
   pchoice = 0,
   edgeTList = NULL
   )
}
\arguments{
  \item{phy}{phylogeny in phylo format}
  \item{pa}{presence-absence table with the first column the species labels and the second column the presence (1) or absence (0) of the species}
  \item{pars}{Vector of model parameters: \cr
  \code{pars[1]} corresponds to mu (extinction rate in local community) \cr
  \code{pars[2]} corresponds to gamma_0 in formula  gamma(t) = gamma_0/(1 + gamma_1 * t) where gamma(t) is immigration rate into local community) \cr
  \code{pars[3]} corresponds to gamma_1 in formula  gamma(t) = gamma_0/(1 + gamma_1 * t) where gamma(t) is immigration rate into local community)
  }
  \item{pchoice}{sets the p-value to optimize: \cr
    pchoice == 0 corresponds to the sum of p_0f + p_1f \cr
    pchoice == 1 corresponds to p_0f \cr
    pchoice == 2 corresponds to p_1f \cr
  }
  \item{edgeTList}{list of edge lengths that need to be succesively pruned; if not specified, it will computed using compute_edgeTList}
}
\value{The loglikelihood}
\references{
Pigot, A.L. & R.S. Etienne (2015). A new dynamic null model for phylogenetic community structure. Ecology Letters 18: 153-163.
}
\author{Rampal S. Etienne}
\seealso{
 \code{\link{DAMOCLES_ML}}
 \code{\link{DAMOCLES_sim}}
}
\examples{
  #TEST IT WORKS
  library(ape)
  phy = rcoal(100)
  pars = c(0.5,0.1,0.1)
  pa = rbinom(100,c(0,1),0.5)
  pa = matrix(c(phy$tip.label,pa),nrow = length(phy$tip.label),ncol = 2)

  # - without a root edge
  loglik = DAMOCLES_loglik(phy,pa,pars)
  loglik

  # - with a root edge
  phy$root.edge = 2
  loglik = DAMOCLES_loglik(phy,pa,pars)
  loglik
}
\keyword{ models }
