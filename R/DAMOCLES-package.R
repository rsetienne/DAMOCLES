

#' Internal DAMOCLES Functions
#' 
#' Internal DAMOCLES functions
#' 
#' These are not to be called by the user.
#' 
#' @aliases compute_edgeTList DAMOCLES_loglik_rhs DAMOCLES_integrate_ODE
#' DAMOCLES_loglik_choosepar DAMOCLES_simplex DAMOCLES_bin_trait_sim
#' DAMOCLES_all_loglik DAMOCLES_all_loglik_rh DAMOCLES_all_loglik_old
#' DAMOCLES_all_loglik_rhs DAMOCLES_all_integrate_ODE
#' DAMOCLES_all_integrate_ODE_old DAMOCLES_all_loglik_choosepar DAMOCLES_all_ML
#' DAMOCLES_all_M DAMOCLES_all_node DAMOCLES_mats DAMOCLES_r
#' DAMOCLES_check_Mlist DAMOCLES_check_edgeTList
#' HERACLES_find_trait_stait_combinations o.dectobin combsUse
#' HERACLES_extractCI ses.mn_trait_d ses.m_trait_d mtd mntd
#' HERACLES_speciation_sim HERACLES_assign_state integral_peak2
#' @keywords internal
NULL





#' Dynamic Assembly Model Of Colonization, Local Extinction and Speciation
#' 
#' Simulation and likelihood methods for a dynamical community assembly that
#' accounts for phylogenetic history.\cr\cr New in version 1.1:\cr - Function
#' added to do bootstrap\cr\cr New in version 2.0:\cr - Computes likelihood for
#' a model where colonization and local extinction rates depend on binary trait
#' value.\cr - Computes likelihood for a model where colonization and local
#' extinction rates depend on geographical trait.\cr - Allows use of subplex as
#' optimization algorithm.\cr - Uses matrix exponentiation to solve ODE
#' system.\cr\cr
#' 
#' \tabular{ll}{ Package: \tab DAMOCLES\cr Type: \tab Package\cr Version: \tab
#' 2.0\cr Date: \tab 2016-06-10\cr License: \tab GPL 2.0\cr } DAMOCLES_loglik
#' computes the likelihood of presence-absence data in a local community given
#' a set of parameters and a phylogeny under an immigration-extinction model
#' 
#' DAMOCLES_ML finds the parameters that maximizes the likelihood computed by
#' DAMOCLES_loglik.
#' 
#' DAMOCLES_sim simulates presence-absence data for a given phylogeny
#' 
#' DAMOCLES_bootstrap computes the maximum likelihood estimates of colonisation
#' and local extinction rate for a given phylogeny and presence-absence data
#' under the DAMOCLES model.  These rate estimates are used to simulate null
#' communities under the DAMOCLES model.  Standardized effect size of mean
#' nearest taxon distance (mntd), mean phylogentic distance (mpd) and
#' loglikelihood are calculated.  For comparison, standardised effect sizes are
#' also calculated relative to a "Random-Draw" null model i.e. presence absence
#' randomised across tips.
#' 
#' @name DAMOCLES-package
#' @aliases DAMOCLES-package DAMOCLES
#' @docType package
#' @author Rampal S. Etienne & Alex L. Pigot\cr Maintainer: Rampal S. Etienne
#' (r.s.etienne@@rug.nl)
#' @references Pigot, A.L. & R.S. Etienne (2015). A new dynamic null model for
#' phylogenetic community structure. Ecology Letters 18: 153-163.
#' @keywords models
#' @examples
#' 
#' DAMOCLES_ML()
#' 
NULL





#' Dated phylogenetic tree of the New World Primates in nexus format and
#' presence-absence matrix for species in Manu
#' 
#' A list with two elements. \cr.  \code{phy} is a dated molecular phylogeny
#' for 94 species of New World Primates extracted from the maximum likelihood
#' tree (AUTOsoft dated) of Springer et al. (2012). 1 time unit = 100 million
#' years. \cr \code{pa} is the presence-absence matrix of NW Primates in Manu
#' from Solari et al. (2006). The first column indicate the species tip labels
#' and the second column indicates presence (1) and absence (0).
#' 
#' 
#' @name NWPrimates_data
#' @docType data
#' @format A list with two elements. The first element (\code{phy}) is the
#' primate phylogeny in nexus format. The second element (\code{pa}) is the
#' presence-absence matrix with 94 rows and 2 columns.
#' @seealso \code{\link{DAMOCLES_sim}}, \code{\link{DAMOCLES_ML}},
#' \code{\link{DAMOCLES_loglik}}
#' @source Solari, S., Pacheco, V., Luna, L., Velazco, P.M. & Patterson, B.D.
#' 2006 Mammals of the manu biosphere reserve. Fieldiana Zoology 110, 13-22.\cr
#' Springer, M.S., Meredith, R.W., Gatesy, J., Emerling, C.A., Park, J.,
#' Rabosky, D.L., Stadler, T., Steiner, C., Ryder, O.A., Janecka, J.E., et al.
#' 2012 Macroevolutionary dynamics and historical biogeography of primate
#' diversification inferred from a species supermatrix. Plos One 7. (doi:ARTN
#' e49521 DOI 10.1371/journal.pone.0049521).
#' @keywords datasets
NULL



