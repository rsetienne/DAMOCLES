#' Simulating DAMOCLES
#' 
#' Simulates DAMOCLES
#' 
#' 
#' @param phy phylogeny in phylo format
#' @param gamma_0 initial per lineage rate of immigration (gamma)
#' @param mu per lineage rate of local extinction
#' @param root.state geographic state of ancestor i.e. present (1) or absent(0)
#' @param keepExtinct whether to retain data for extinct lineages
#' @param stateOnly whether to return only the presence or absence states of the extant species
#' @return Either a table containing the following
#' columns: The first two columns contain the edge list of the phylogenetic tree, and the following two
#' contain when each edge starts and ends, and the fifth column if lineage is extant.
#' The sixth column contains the presence (1) or absence (0) of the species, which is identical
#' to the output if stateOnly == TRUE, when keepExtinct == FALSE and stateOnly == FALSE. Column
#' seven and eight contain information on respectively whetever the edge is a tip and its tip number.
#' @author Bouwe R. Reijenga
#' @seealso \code{\link{DAMOCLES_sim}}
#' @references Pigot, A.L. & R.S. Etienne (2015). A new dynamic null model for
#' phylogenetic community structure. Ecology Letters 18: 153-163.
#' @keywords models
#' @examples
#' 
#' #create random phylogeny
#' library(ape)
#' phy = ape::rcoal(10)
#' 		
#' #run DAMOCLES		
#' patable = DAMOCLES_sim_fast(
#'   phy,
#'   gamma_0 = 1.5,
#'   mu = 0,
#'   root.state = 1,
#'   keepExtinct = FALSE,
#'   stateOnly = FALSE
#'   )
#' 
#' #show presence/absence on the tree
#' patable$col = rep("black",dim(patable)[1])
#' patable$col[which(patable$state == 1)] = "red"
#' plot(phy,tip.col = patable$col)
#' 
#' @export DAMOCLES_sim_fast
DAMOCLES_sim_fast = function(phy, gamma_0, mu, root.state, keepExtinct = FALSE, stateOnly = FALSE){
  # first construct the presence absence table from the phylogeny
  dn = ape::dist.nodes(phy)
  ntips = length(phy$tip.label)
  nbranch = 2 * ntips - 2
  patable = data.frame(p = phy$edge[, 1], d = phy$edge[, 2], 
                       age.start = rep(NA, nbranch), age.end = rep(NA, nbranch), 
                       extant = rep(0, nbranch), state = rep(NA, nbranch))
  patable$age.start = dn[patable$p, ntips + 1]
  patable$age.end = dn[patable$d, ntips + 1]
  patable$age.end[which(patable$d < length(phy$tip.label) + 
                          1)] = max(patable$age.end)
  ce = which(patable$p == min(patable$p))
  patable$extant[ce] = 1
  patable$state = rep(0, dim(patable)[1])
  patable$tips[which(patable$d <= ntips)] = 1
  patable$y[which(patable$tips == 1)] = patable$d[which(patable$tips == 1)]
  # sample which species has which root.state
  patable$state[ce] = sample(c(0, root.state))
  # simulte
  patable$state <- DAMOCLES_sim_cpp(pa = patable, gamma = gamma_0, mu = mu)
  # return the full table or only the extant states
  if(stateOnly == TRUE){
    patable <- patable$state[which(patable$tips == 1)]
  }else{
    if(keepExtinct == TRUE){
      patable$extant[which(patable$tips == 1)]
    }else{
      patable <- patable[which(patable$tips == 1),]
      patable$extant <- 1
    }
  }
  return(patable)
}