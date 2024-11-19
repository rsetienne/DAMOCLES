compute_edgeTList = function(tree)
{
  tree$node.label = NULL
  nNode = ape::Nnode(tree)
  edgeTList = vector("list", nNode)
  ca = max(ape::branching.times(tree))
  for (i in 1:nNode)
  {
    ntips = ape::Ntip(tree)
    if(ntips > 2)
    {
      nodeDepths = ape::dist.nodes(tree)[(ntips + 1):(2 * ntips - 1), ntips + 1]
      w = which(nodeDepths == max(nodeDepths))[1]
      toDrop1 = caper::clade.members(as.numeric(names(nodeDepths[w])), tree, tip.labels = TRUE, include.nodes=TRUE)$tips
      toDrop2 = c(which(tree$tip.label == toDrop1[1]), which(tree$tip.label == toDrop1[2]))
      edgeT = ca - c(ape::dist.nodes(tree)[toDrop2, ntips + 1], nodeDepths[w])
    } else
    {
      nodeDescendants = caper::clade.members.list(tree, tip.labels = TRUE)
      sisterNodes = which(sapply(nodeDescendants, length) == 2)
      sisterDepths = nodeDepths[sisterNodes]
      w = which(sisterDepths == max(sisterDepths))[1]
      toDrop1 = nodeDescendants[[sisterNodes[w]]]
      toDrop2 = c(which(tree$tip.label == toDrop1[1]), which(tree$tip.label == toDrop1[2]))
      edgeT = ca - c(ape::dist.nodes(tree)[toDrop2, ntips + 1], sisterDepths[w])
    }
    names(edgeT)[1:2] = toDrop1
    edgeTList[[i]] = edgeT
    tree = suppressWarnings(ape::drop.tip(tree, toDrop1, trim.internal = FALSE))
    tree$tip.label[which(tree$tip.label == "NA")] = paste("p",i - 1, sep = "")
  }
  return(edgeTList)
}

DAMOCLES_all_M = function(
   pars,
   model = 0
   )
{  
   if(model <= 0)
   {
      mu = pars[1]
      ga = pars[2]
      M = matrix(0, nrow = 2, ncol = 2)
      M[1,1] = -ga
      M[2,2] = -mu
      M[1,2] = ga
      M[2,1] = mu
   }
   if(model == 1)
   {
      mu0 = pars[1]
      gam0 = pars[2]
      mu1 = pars[3]
      gam1 = pars[4]
      q0t1 = pars[5]
      q1t0 = pars[6]
      M = matrix(0, nrow = 4, ncol = 4)
      M[1,1] = -(gam0 + q0t1)
      M[2,2] = -(gam1 + q1t0)
      M[3,3] = -(mu0 + q0t1)
      M[4,4] = -(mu1 + q1t0)
      M[1,2] = q0t1
      M[2,1] = q1t0
      M[1,3] = gam0
      M[3,1] = mu0
      M[2,4] = gam1
      M[4,2] = mu1
      M[3,4] = q0t1
      M[4,3] = q1t0
   }
   if(model == 2)
   {
      mu1 = pars[1]
      gam1 = pars[2]
      mu2 = pars[3]
      gam2 = pars[4]
      q0t2 = pars[5]
      q2t0 = pars[6]
      q1t2 = pars[7]
      q2t1 = pars[8]
      M = matrix(0, nrow = 5, ncol = 5)
      M[1,1] = -q0t2
      M[2,2] = -(gam1 + q1t2)
      M[3,3] = -(gam2 + q2t0 + q2t1)
      M[4,4] = -(mu1 + q1t2)
      M[5,5] = -(mu2 + q2t0 + q2t1)
      M[1,3] = q0t2
      M[2,3] = q1t2
      M[2,4] = gam1
      M[3,1] = q2t0
      M[3,2] = q2t1
      M[3,5] = gam2
      M[4,2] = mu1
      M[4,5] = q1t2
      M[5,1] = q2t0
      M[5,3] = mu2
      M[5,4] = q2t1
      if(sum(is.nan(M)) > 0 | sum(M == Inf) > 0) { print(M); utils::flush.console()}
   }
   if(model == 0.1 | model == 0.2 | model == 2.2)
   {
      mu0 = (model != 2.2) * pars[1] + (model == 2.2) * 10^10 * max(pars)
      mu1 = pars[1]
      gam0 = (model != 2.2) * pars[2]
      gam1 = pars[2]
      mu2 = (model != 2.2) * pars[1] + (model == 2.2) * pars[3]
      gam2 = (model != 2.2) * pars[2] + (model == 2.2) * pars[4]
      q0t2 = pars[3 + 2 * (model == 2.2)]
      q2t0 = pars[4 + 2 * (model == 2.2)]
      q1t2 = pars[5 + 2 * (model == 2.2)]
      q2t1 = pars[6 + 2 * (model == 2.2)]
      M = matrix(0, nrow = 6, ncol = 6)
      M[1,1] = -(gam0 + q0t2)
      M[2,2] = -(gam1 + q1t2)
      M[3,3] = -(gam2 + q2t0 + q2t1)
      M[4,4] = -(mu0 + q0t2)
      M[5,5] = -(mu1 + q1t2)
      M[6,6] = -(mu2 + q2t0 + q2t1)
      M[1,3] = q0t2
      M[1,4] = gam0
      M[2,3] = q1t2
      M[2,5] = gam1
      M[3,1] = q2t0
      M[3,2] = q2t1
      M[3,6] = gam2
      M[4,1] = mu0
      M[4,6] = q0t2
      M[5,2] = mu1
      M[5,6] = q1t2
      M[6,3] = mu2
      M[6,4] = q2t0
      M[6,5] = q2t1
   }
   return(M)
}

DAMOCLES_all_loglik_rhs = function(
   t,
   p,
   parscaMmodel
   )
{
   pars = parscaMmodel[[1]]
   ca = parscaMmodel[[2]]
   M = parscaMmodel[[3]]
   model = parscaMmodel[[4]]
   if(model == -1)
   {
      mu = pars[1]
      ga = pars[2]/(1 + pars[3] * (ca - t))
      dp = c(ga * (p[2] - p[1]),mu * (p[1] - p[2]))
   } else {
      dp = M %*% p
   }
   return(list(dp))
}        

DAMOCLES_check_Mlist = function(Mlist,pars,model,methode = 'analytical')
{
  if(is.null(Mlist))
  {
    M = DAMOCLES_all_M(pars,model)
    if(methode == 'analytical')
    {
      eigenstuff = eigen(M) 
      eigs = eigenstuff$values
      S = eigenstuff$vectors
      invS = solve(S)
      rm(eigenstuff)
      Mlist = list(M = M, S = S, invS = invS, eigs = eigs)
    } else
    {
      Mlist = list(M = M, S = NULL, invS = NULL, eigs = NULL)
    }
  }
  return(Mlist)
}

DAMOCLES_all_integrate_ODE = function(
   Mlist = NULL,
   pars,
   p,
   tt,
   ca,
   methode = 'analytical',
   model,
   numvar
   )
{
   Mlist = DAMOCLES_check_Mlist(Mlist,pars,model,methode)
   if ((methode == 'analytical' | methode == 'expm' | methode == 'Matrix') & model != -1)
   # TAKE ANALYTICAL SOLUTION
   { 
      difft = tt[2] - tt[1]
      if(model == 0)
      { 
         mu = pars[1]
         ga = pars[2]
         difft = tt[2] - tt[1]
         p0f = mu * p[1] + ga * p[2] + ga * (p[1] - p[2]) * exp(-difft * (ga + mu))
         p1f = mu * p[1] + ga * p[2] - mu * (p[1] - p[2]) * exp(-difft * (ga + mu))
         p = 1/(ga + mu) * c(p0f,p1f)
      } else if (methode == 'expm') p <- expm::expAtv(A = Mlist$M, v = p, t = difft)[[1]]
        else if (methode == 'Matrix') p <- Matrix::expm(Mlist$M * difft) %*% p
        else if (methode == 'analytical') {
         p = (Mlist$S %*% (diag(exp(Mlist$eigs * difft))) %*% Mlist$invS) %*% p
         if(prod(abs(Im(p)) > 1E-3 * abs(Re(p))))
         {
            warning('Numerical problems in matrix exponentiation detected. Results may be unreliable. Rerun with analytical = FALSE.\n')
         }
         p = Re(p)
      }
   } else {
     # SOLVE ODE NUMERICALLY
     if (startsWith(methode, "odeint::")) {
       y <- DAMOCLES_integrate_odeint(p, tt, Mlist$M, atol = 1E-16, rtol = 1E-10, stepper = methode)
       p <- y
     } else {  
       y = deSolve::ode(p,tt,DAMOCLES_all_loglik_rhs,list(pars,ca,Mlist$M,model),rtol = 1E-10,atol = 1E-16, method = methode)
       p = y[2,2:(1 + numvar)]
     }
   }
   return(p)
}

DAMOCLES_all_node = function(
   pc,
   r,
   rmat,
   model,
   numvar
   )
{   
   pnew = rep(0,numvar)
   if(model <= 0)
   {
     pnew = rep(0,2)
     pnew[1] = pc[1,1] * pc[2,1] #p0
     pnew[2] = 1/2 * (pc[1,2] * pc[2,1] + pc[2,2] * pc[1,1]) #p1
   } 
   if(model == 1)
   {
     pnew[1] = (1 - r[1]) * pc[1,1] * pc[2,1] + 1/2 * r[1] * (pc[1,1] * pc[2,2] + pc[1,2] * pc[2,1]) #p00
     pnew[2] = (1 - r[2]) * pc[1,2] * pc[2,2] + 1/2 * r[2] * (pc[1,2] * pc[2,1] + pc[1,1] * pc[2,2]) #p01
     pnew[3] = 1/2 * (1 - r[1]) * (pc[1,3] * pc[2,1] + pc[1,1] * pc[2,3]) + 1/4 * r[1] * (pc[1,3] * pc[2,2] + pc[1,1] * pc[2,4] + pc[1,4] * pc[2,1] + pc[1,2] * pc[2,3]) #p10
     pnew[4] = 1/2 * (1 - r[2]) * (pc[1,4] * pc[2,2] + pc[1,2] * pc[2,4]) + 1/4 * r[2] * (pc[1,3] * pc[2,2] + pc[1,1] * pc[2,4] + pc[1,4] * pc[2,1] + pc[1,2] * pc[2,3]) #p11
   }
   if(model == 2)
   {
     pnew[1] = (1 - r[1]) * pc[1,1] * pc[2,1] + 1/2 * r[1] * (pc[1,2] * pc[2,1] + pc[1,1] * pc[2,2]) #p00
     pnew[2] = (1 - r[2]) * pc[1,2] * pc[2,2] + 1/2 * r[2] * (pc[1,2] * pc[2,1] + pc[1,1] * pc[2,2]) #p01
     pnew[3] = (1 - r[3]) * pc[1,3] * pc[2,3] + 1/2 * r[3] * (1 - r[4]) * (pc[1,2] * pc[2,1] + pc[1,1] * pc[2,2]) + 1/4 * r[3] * r[4] * (pc[1,3] * pc[2,1] + pc[1,1] * pc[2,3] + pc[1,3] * pc[2,2] + pc[1,2] * pc[2,3]) #p02
     pnew[4] = 1/2 * (1 - r[2]) * (pc[1,4] * pc[2,2] + pc[1,2] * pc[2,4]) + 1/2 * r[2] * (pc[1,1] * pc[2,4] + pc[1,4] * pc[2,1]) #p11
     pnew[5] = 1/2 * (1 - r[3]) * (pc[1,5] * pc[2,3] + pc[1,3] * pc[2,5]) + 1/2 * r[3] * (1 - r[4]) * (pc[1,1] * pc[2,4] + pc[1,4] * pc[2,1]) + 1/8 * r[3] * r[4] * (pc[1,5] * pc[2,1] + pc[1,1] * pc[2,5] + pc[1,3] * pc[2,4] + pc[1,4] * pc[2,3] + pc[1,2] * pc[2,5] + pc[1,5] * pc[2,2]) #p12
   } 
   if(model == 0.1)
   {
     pnew[1] = (1 - r[1]) * pc[1,1] * pc[2,1] + 1/2 * r[1] * (pc[1,2] * pc[2,1] + pc[1,1] * pc[2,2]) #p00
     pnew[2] = (1 - r[2]) * pc[1,2] * pc[2,2] + 1/2 * r[2] * (pc[1,2] * pc[2,1] + pc[1,1] * pc[2,2]) #p01
     pnew[3] = (1 - r[3]) * pc[1,3] * pc[2,3] + 1/2 * r[3] * (1 - r[4]) * (pc[1,2] * pc[2,1] + pc[1,1] * pc[2,2]) + 1/4 * r[3] * r[4] * (pc[1,3] * pc[2,1] + pc[1,1] * pc[2,3] + pc[1,3] * pc[2,2] + pc[1,2] * pc[2,3]) #p02
     pnew[4] = 1/2 * (1 - r[1]) * (pc[1,4] * pc[2,1] + pc[1,1] * pc[2,4]) + 1/4 * r[1] * (pc[1,1] * pc[2,5] + pc[1,5] * pc[2,1] + pc[1,2] * pc[2,4] + pc[1,4] * pc[2,2]) #p10 
     pnew[5] = 1/2 * (1 - r[2]) * (pc[1,5] * pc[2,2] + pc[1,2] * pc[2,5]) + 1/4 * r[2] * (pc[1,1] * pc[2,5] + pc[1,5] * pc[2,1] + pc[1,5] * pc[2,4] + pc[1,4] * pc[2,5]) #p11
     pnew[6] = 1/2 * (1 - r[3]) * (pc[1,6] * pc[2,3] + pc[1,3] * pc[2,6]) + 1/4 * r[3] * (1 - r[4]) * (pc[1,1] * pc[2,5] + pc[1,5] * pc[2,1] + pc[1,4] * pc[2,2] + pc[1,2] * pc[2,4]) + 1/8 * r[3] * r[4] * (pc[1,6] * pc[2,1] + pc[1,1] * pc[2,6] + pc[1,3] * pc[2,5] + pc[1,5] * pc[2,3] + pc[1,2] * pc[2,6] + pc[1,6] * pc[2,2] + pc[1,4] * pc[2,3] + pc[1,3] * pc[2,4]) #p12
   }
   if(model == 0.2 | model == 2.2)
   {  
     ns = 6
     pnew = rep(0,ns)
     for(i in 1:ns)
     {
        #for(j in 1:ns) for(k in 1:ns)
        #{
        #   pnew[i] = pnew[i] + rmat[[i]][j,k] * pc[1,j] * pc[2,k]
        #}
        #pnew[i] = sum(rmat[[i]] * pcmat) 
        pnew[i] = sum(pc[1,1 + rmat[[i]]@i] * pc[2,1 + rmat[[i]]@j] * rmat[[i]]@x)
     }
   } 
   return(pnew)
}   

DAMOCLES_r = function(
   r,
   model
)
{
   ns = 6
   rmat = list()
   if(model == 0.2 | model == 2.2)
   {
      rmat[[1]] = Matrix::sparseMatrix(i = c(1,1),j = c(1,2),x = c(1 - r[1], r[1]/2), symmetric = TRUE, dims = c(ns,ns))
      rmat[[2]] = Matrix::sparseMatrix(i = c(2,1),j = c(2,2),x = c(1 - r[2], r[2]/2), symmetric = TRUE, dims = c(ns,ns))
      rmat[[3]] = Matrix::sparseMatrix(i = c(3,2,1,1),j = c(3,3,3,2),x = c(1 - r[3],r[3] * r[4]/4,r[3] * r[4]/4,r[3] * (1 - r[4])/2), symmetric = TRUE, dims = c(ns,ns))
      if(model == 0.2)
      {
         rmat[[4]] = Matrix::sparseMatrix(i = c(1,1,2),j = c(4,5,4),x = c((1 - r[1])/2,r[1]/4,r[1]/4), symmetric = TRUE, dims = c(ns,ns))
         rmat[[5]] = Matrix::sparseMatrix(i = c(2,4,1),j = c(5,5,5),x = c((1 - r[2])/2,r[2]/4,r[2]/4), symmetric = TRUE, dims = c(ns,ns))
         rmat[[6]] = Matrix::sparseMatrix(i = c(3,2,3,1,3,1,2),j = c(6,6,5,6,4,5,4),x = c((1 - r[3])/2,r[3] * r[4]/8, r[3] * r[4]/8,r[3] * r[4]/8,r[3] * r[4]/8,r[3] * (1 - r[4])/2,r[3] * (1 - r[4])/4), symmetric = TRUE, dims = c(ns,ns))
      }
      if(model == 2.2)
      {
         rmat[[4]] = Matrix::sparseMatrix(i = 1,j = 1,x = 0, dims = c(ns,ns))
         rmat[[5]] = Matrix::sparseMatrix(i = c(2,1),j = c(5,5),x = c((1 - r[2])/2,r[2]/2), symmetric = TRUE, dims = c(ns,ns))
         rmat[[6]] = Matrix::sparseMatrix(i = c(3,2,3,1,1),j = c(6,6,5,6,5),x = c((1 - r[3])/2,r[3] * r[4]/8,r[3] * r[4]/8,r[3] * r[4]/8,r[3] * (1 - r[4])/2), symmetric = TRUE, dims = c(ns,ns))
      }
      for(i in 1:length(rmat))
      {
         rmat[[i]] = methods::as(rmat[[i]], "dgTMatrix")
      }
   }
   return(rmat)
}

DAMOCLES_check_edgeTList = function(phy,edgeTList)
{
  if(is.null(edgeTList))
  {
     edgeTList = compute_edgeTList(phy)
  }
  return(edgeTList)
}

#' Likelihood for DAMOCLES model
#' 
#' Computes likelihood for the presence-absence data of species in a local
#' community for a given phylogeny of species in the region.
#' 
#' 
#' @param phy phylogeny in phylo format
#' @param pa presence-absence table with the first column the species labels
#' and the second column the presence (1) or absence (0) of the species
#' @param pars Vector of model parameters:\cr
#' \code{pars[1]} corresponds to mu (extinction rate in local community)\cr
#' \code{pars[2]} corresponds to gamma_0 in formula
#' gamma(t) = gamma_0/(1 + gamma_1 * t) where gamma(t) is immigration rate
#' into local community)\cr 
#' \code{pars[3]} corresponds to
#' gamma_1 in formula gamma(t) = gamma_0/(1 + gamma_1 * t) where gamma(t) is
#' immigration rate into local community)
#' @param pchoice sets the p-value to optimize:\cr
#' pchoice == 0 corresponds to
#' the sum of p_0f + p_1f\cr
#' pchoice == 1 corresponds to p_0f\cr
#' pchoice == 2 corresponds to p_1f\cr
#' @param edgeTList list of edge lengths that need to be succesively pruned; if
#' not specified, it will computed using compute_edgeTList
#' @param methode method used to solve the ODE. Either 'analytical' for the analytical
#' solution, 'Matrix' for matrix exponentiation using package Matrix or 'expm' using
#' package 'expm' or any of the numerical solvers, used in deSolve.
#' @param model model used. Default is 0 (standard null model). Other options are 1 (binary traits)
#' 2 (trinary environmental trait) or 3 (diversity-dependent colonization - beta version)
#' @param Mlist list of M matrices that can be specified when methode = 'analytical'. If set
#' at NULL (default) and methode = 'analytical', Mlist will be computed.
#' @param verbose Whether intermediate output should be printed. Default is FALSE.
#' @return The loglikelihood
#' @author Rampal S. Etienne
#' @seealso \code{\link{DAMOCLES_ML}} \code{\link{DAMOCLES_sim}}
#' @references Pigot, A.L. & R.S. Etienne (2015). A new dynamic null model for
#' phylogenetic community structure. Ecology Letters 18: 153-163.
#' @keywords models
#' @examples
#' 
#'   #TEST IT WORKS
#'   library(ape)
#'   phy = ape::rcoal(100)
#'   pars = c(0.5,0.1,0.1)
#'   pa = rbinom(100,c(0,1),0.5)
#'   pa = matrix(c(phy$tip.label,pa),nrow = length(phy$tip.label),ncol = 2)
#' 
#'   # - without a root edge
#'   loglik = DAMOCLES_loglik(phy,pa,pars)
#'   loglik
#' 
#'   # - with a root edge
#'   phy$root.edge = 2
#'   loglik = DAMOCLES_loglik(phy,pa,pars)
#'   loglik
#' 
#' @export DAMOCLES_loglik
DAMOCLES_loglik <- DAMOCLES_all_loglik <- function(
   phy,
   pa,
   pars,
   pchoice = 0,
   edgeTList = NULL,
   methode = 'analytical',
   model = 0, 
   Mlist = NULL,
   verbose = FALSE
   )
{
  edgeTList = DAMOCLES_check_edgeTList(phy,edgeTList)
  Mlist = DAMOCLES_check_Mlist(Mlist,pars,model,methode)
  patrait = data.frame(as.character(pa[,1]),matrix(as.numeric(pa[,2:ncol(pa)]),ncol = ncol(pa) - 1,byrow = FALSE))
  names(patrait) = c('label',paste('pa',1:(ncol(patrait) - 1),sep = ""))
  pchoicevec = pchoice
  if(model == -1 || model == 0)
  {
      numvar = 2
      #patraittable = cbind(as.character(pa[,1]),as.character(1 - as.numeric(pa[,2])),pa[,2])
      patraittable = data.frame(as.character(patrait[,1]),1 - patrait[,2],patrait[,2])
      if(length(pchoice) == 1)
      {
         pchoicevec = as.numeric(c(is.element(pchoice,c(0,1)),is.element(pchoice,c(0,2))))
      }
      r = 0
  } else
  if(model == 1)
  {
      numvar = 4
      #patraittable = cbind(as.character(patrait[,1]),
      patraittable = data.frame(as.character(patrait[,1]),
        (patrait[,2] == 0) * (patrait[,3] == 0),
        (patrait[,2] == 0) * (patrait[,3] == 1),
        (patrait[,2] == 1) * (patrait[,3] == 0),
        (patrait[,2] == 1) * (patrait[,3] == 1))
      if(length(pchoice) == 1)
      {
        pchoicevec = as.numeric(c(rep(is.element(pchoice,c(0,1)),2),rep(is.element(pchoice,c(0,2)),2)))
      }
      r = pars[7:8]  
  } else
  if(model == 2)
  {
      numvar = 5
      #patraittable = cbind(as.character(patrait[,1]),
      patraittable = data.frame(as.character(patrait[,1]),
        (patrait[,2] == 0) * (patrait[,3] == 0),
        (patrait[,2] == 0) * (patrait[,3] == 1),
        (patrait[,2] == 0) * (patrait[,3] == 2),
        (patrait[,2] == 1) * (patrait[,3] == 1),
        (patrait[,2] == 1) * (patrait[,3] == 2))
      if(length(pchoice) == 1)
      {
        pchoicevec = as.numeric(c(rep(is.element(pchoice,c(0,1)),3),rep(is.element(pchoice,c(0,2)),2)))
      }
      r = pars[9:12]  
  } else
  if(model == 0.1 | model == 0.2 | model == 2.2)
  {
      numvar = 6
      #patraittable = cbind(as.character(patrait[,1]),
      patraittable = data.frame(as.character(patrait[,1]),
        (patrait[,2] == 0) * (patrait[,3] == 0),
        (patrait[,2] == 0) * (patrait[,3] == 1),
        (patrait[,2] == 0) * (patrait[,3] == 2),
        (patrait[,2] == 1) * (patrait[,3] == 0),
        (patrait[,2] == 1) * (patrait[,3] == 1),
        (patrait[,2] == 1) * (patrait[,3] == 2))
      if(length(pchoice) == 1)
      {
         pchoicevec = as.numeric(c(rep(is.element(pchoice,c(0,1)),3),(model != 2.2) * is.element(pchoice,c(0,2)), rep(is.element(pchoice,c(0,2)),2)))
      }
      r = pars[7:10 + 2 * (model == 2.2)] 
  }
  names(patraittable) = c('label',paste('pat',1:(ncol(patraittable) - 1),sep = ""))
  pchoicevec = pchoicevec / sum(pchoicevec)
  loglik = 0
  if(max(r) > 1 | min(pars) < 0)
  {
     loglik = -Inf
     return(loglik)
  }  
 	nNode = dim(patraittable)[1] - 1
  rmat = DAMOCLES_r(r,model)   
	for(i in 1:nNode)
	{
 	   ca = max(unlist(edgeTList))
   	 branchesTimes = edgeTList[[i]]
     pc = matrix(0,nrow = 2,ncol = numvar)
     w2 = c(0,0)
     for(j in 1:2)
     {
        tt = c(branchesTimes[j],branchesTimes[3])         
        w2[j] = which(patraittable[,1] == names(branchesTimes)[j])
        p = as.numeric(patraittable[w2[j],2:(numvar + 1)])
        p = DAMOCLES_all_integrate_ODE(Mlist,pars,p,tt,ca,methode,model,numvar)
        if(min(p) < 0 | is.nan(min(p)))
        {
          cat('Numerical problems: negative probabilities encountered. Results may not be reliable.\n')
          utils::flush.console()
        }
        p[p < 0] = 0
        p[is.nan(p)] = 0
        sump = sum(p)
        if(sump == 0)
        {
          cat('Numerical problems: all probabilities are zero or negative!\n')
          utils::flush.console()
          loglik = -Inf
          return(loglik)       
        }
        loglik = loglik + log(sump)
        p = p/sump
        patraittable[w2[j],2:(numvar + 1)] = p
        pc[j,] = p
     }   
     #print(pc)
     #print(r)
     pnew = DAMOCLES_all_node(pc,r,rmat,model,numvar)
     if(min(pnew) < 0 || is.nan(min(pnew)))
     {
       cat('Numerical problems: negative probabilities encountered. Results may not be reliable.\n')
       utils::flush.console()
       pnew[pnew < 0] = 0
       pnew[is.nan(pnew)] = 0
       if(sum(pnew) == 0)
       {
         cat('Numerical problems: all probabilities are zero or negative!\n')
         utils::flush.console()
         loglik = -Inf
         return(loglik)       
       }
     }
     sumpnew = sum(pnew)
     if(sumpnew == 0)
     {
       loglik = -Inf
       return(loglik)       
     }
     loglik = loglik + log(sumpnew)
     pnew = pnew/sumpnew
     updateTips = as.character(patraittable[,1])
     updateTips[w2[1]] = as.character(paste("p",i - 1,sep = ""))
     patraittable[,1] = updateTips
     patraittable[w2[1],2:(numvar + 1)] = pnew
     patraittable = patraittable[-w2[2],]
	}
	if(!is.null(phy$root.edge))
	{
     # SPECIFY START AND FINISH OF BRANCH j
     tt = c(0,phy$root.edge)
     #tt = c(branchesTimes[3],0)
     # READ INITIAL VALUES FOR THIS BRANCH
     p = as.numeric(patraittable[2:(numvar + 1)])
     p = DAMOCLES_all_integrate_ODE(Mlist,pars,p,tt,ca,methode,model,numvar)     
     patraittable[2:(numvar + 1)] = p
	}  
  #if(min(patraittable[2:(numvar + 1)]) < 0) { print(as.numeric(patraittable[2:(numvar + 1)])) }
	loglik = loglik + log(sum(pchoicevec * as.numeric(patraittable[2:(numvar + 1)])))
  if (verbose) { 
    s1 = sprintf('Parameters:')
    s2 = sprintf('%f ',pars)
    s3 = sprintf('\nLoglikelihood: %f',loglik)
    cat(s1,s2,s3,"\n",sep = "")
    utils::flush.console()
  }
  return(loglik)
}

DAMOCLES_all_loglik_rh = function(
   phy,
   pa,
   pars,
   pchoice = 0,
   edgeTList = NULL,
   methode = 'analytical',
   model = 0 
   )
{
   edgeTList = DAMOCLES_check_edgeTList(phy,edgeTList)
   logintegrand = function(x)
   {
       logintegrandvec = rep(0,length(x))
       for(i in 1:length(x))
       {
          pars1 = c(pars[1],x[i],pars[3],x[i],pars[-c(1:4)])
          Mlist = DAMOCLES_check_Mlist(Mlist = NULL,pars = pars1,model = model,methode = methode)
          mu = pars[2]
          sigma = pars[length(pars)]
          logintegrandvec[i] = log(stats::dgamma(x[i], shape = mu^2/sigma^2, rate = mu/sigma^2)) + DAMOCLES_all_loglik(pars = pars1, Mlist = Mlist, phy = phy,pa = pa,pchoice = pchoice,edgeTList = edgeTList,methode = methode,model = model)
          # mu is the mean of the gamma distribution, and sigma the standard deviation
          # mu = shape/rate, sigma^2 = shape/rate^2
       }
       return(logintegrandvec)
   }
   
   integrand = function(x)
   {
       integrandvec = rep(0,length(x))
       for(i in 1:length(x))
       {
          pars1 = c(pars[1],x[i],pars[-c(1,2)])
          Mlist = DAMOCLES_check_Mlist(Mlist = NULL,pars = pars1,model = model,methode = methode)
          mu = pars[2]
          sigma = pars[length(pars)]
          integrandvec[i] = stats::dgamma(x[i], shape = mu^2/sigma^2, rate = mu/sigma^2) * exp(DAMOCLES_all_loglik(pars = pars1, Mlist = Mlist, phy = phy,pa = pa,pchoice = pchoice,edgeTList = edgeTList,methode = methode,model = model))
          # mu is the mean of the gamma distribution, and sigma the standard deviation
          # mu = shape/rate, sigma^2 = shape/rate^2
       }
       return(integrandvec)
   }
   
   loglik = integral_peak2(logintegrand, xx = seq(-100,10,2), xcutoff = 2, ycutoff = 40, ymaxthreshold = 1E-12)
   #print(loglik)
   #loglik0 = log(stats::integrate(integrand, lower = 0, upper = Inf, rel.tol = 1E-12, abs.tol = 1E-12)$value)
   #print(loglik0)
   return(loglik)
}

integral_peak2 <- function(logfun, xx = seq(-100,10,2), xcutoff = 2, ycutoff = 40, ymaxthreshold = 1E-12)
{
  # 1/ determine integrand peak
  yy <- xx + logfun(exp(xx));
  yy[which(is.na(yy) | is.nan(yy))] <- -Inf;
  yymax <- max(yy);
  if(yymax == -Inf)
  {
     logQ <- -Inf;
     return(logQ);
  }
  iimax <- which(yy >= (yymax - ymaxthreshold));
  xlft <- xx[iimax[1]] - xcutoff;
  xrgt <- xx[iimax[length(iimax)]] + xcutoff;
  optfun <- function(x) x + logfun(exp(x));
  optres <- stats::optimize(f = optfun, interval = c(xlft,xrgt), maximum = TRUE, tol = 1e-10);
  xmax <- optres$maximum;
  ymax <- optres$objective;

  # 2/ determine peak width
  iilft <- which((xx < xmax) & (yy < (ymax - ycutoff)));
  if(length(iilft) == 0)
  {
     xlft <- xx[1] - xcutoff;
  } else
  {
     ilft <- iilft[length(iilft)];
     xlft <- xx[ilft];
  }
  iirgt <- which((xx > xmax) & (yy < (ymax - ycutoff)));
  if(length(iirgt) == 0)
  {
     xrgt <- xx[length(xx)] + xcutoff;
  } else
  {
     irgt <- iirgt[1];
     xrgt <- xx[irgt];
  }
  
  # 3/ compute integral
  intfun <- function(x)
  {
     y <- exp((x + logfun(exp(x))) - ymax);
     y[is.nan(y)] <- 0
     return(y)
  }
  intres <- stats::integrate(f = intfun, lower = xlft, upper = xrgt, rel.tol = 1e-10, abs.tol = 1e-10, subdivisions = 1000L, stop.on.error = FALSE);
  corrfact <- intres$value;
  logQ <- ymax + log(corrfact);
  return(logQ);
}
