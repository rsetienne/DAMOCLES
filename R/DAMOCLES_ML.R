DAMOCLES_all_loglik_choosepar <- function(trparsopt,
                                          idparsopt,
                                          trparsfix,
                                          idparsfix,
                                          idparsequal,
                                          phy,
                                          patrait,
                                          edgeTList,
                                          locatenode,
                                          pchoice,
                                          methode,
                                          model,
                                          verbose)
{
   trpars1 = rep(0,3 * (model == -1) + 2 * (model == 0) + 10 * (model == 0.1 | model == 0.2) + 8 * (model == 1) + 12 * (model == 2 | model == 2.2))
   trpars1[idparsopt] = trparsopt
   if(length(idparsfix) != 0)
   { 
       trpars1[idparsfix] = trparsfix
   }
   if(length(idparsequal) > 0)
   { 
       trpars1[idparsequal] = trpars1[idparsequal - 1 - (model == 2 | model == 2.2) - (model == 1) * (idparsequal <= 4) + (model == 2 | model == 2.2) * (idparsequal == 10)]    
   }
   if(max(trpars1) > 1 | min(trpars1) < 0 | ((model == 2 | model == 2.2) & max(trpars1[9:12]) > 0.5))
   {
       loglik = -Inf
   } else {
       pars1 = trpars1/(1 - trpars1)
       if(model < 3)
       {
         loglik = DAMOCLES_all_loglik(phy = phy,pa = patrait,pars = pars1,pchoice = pchoice,edgeTList = edgeTList,methode = methode,model = model,Mlist = NULL, verbose = verbose)
       } else
       {
         loglik = DAMOCLES_DD_loglik(phy = phy,pa = patrait,pars = pars1,pchoice = pchoice,locatenode = locatenode,methode = methode,verbose = verbose)
       }
   }
   return(loglik)
}



#' Maximization of the loglikelihood under the DAMOCLES model
#' 
#' This function computes the maximum likelihood estimates of the parameters of
#' the DAMOCLES model for a given phylogeny and presence-absence data.  It also
#' outputs the corresponding loglikelihood that can be used in model
#' comparisons.
#' 
#' The output is a dataframe containing estimated parameters and maximum
#' loglikelihood.
#' 
#' @param phy phylogeny in phylo format
#' @param pa presence-absence table.\cr The first column contains the labels of
#' the species (corresponding to the tip labels in the phylogeny.\cr The second
#' column contains the presence (1) or absence (0) of species in the local
#' community.
#' @param initparsopt The initial values of the parameters that must be
#' optimized
#' @param idparsopt The ids of the parameters that must be optimized, e.g. 1:2
#' for extinction rate, and offset of immigration rate The ids are defined as
#' follows: \cr id == 1 corresponds to mu (extinction rate) \cr id == 2
#' corresponds to gamma_0 (offset of immigration rate) \cr id == 3 corresponds
#' to gamma_1 (parameter controlling decline in immigration rate with time)
#' @param parsfix The values of the parameters that should not be optimized.
#' See idparsfix.
#' @param idparsfix The ids of the parameters that should not be optimized,
#' e.g. c(1,3) if mu and gamma_1 should not be optimized, but only gamma_0. In
#' that case idparsopt must be c(2). The default is to fix all parameters not
#' specified in idparsopt.
#' @param idparsequal The ids of the parameters that should be set equal to the 
#' first parameter of the same type.
#' @param pars2 Vector of settings: \cr
#' \code{pars2[1]} sets the relative tolerance in the parameters \cr \cr
#' \code{pars2[2]} sets the relative tolerance in the function \cr \cr
#' \code{pars2[3]} sets the absolute tolerance in the parameters \cr \cr
#' \code{pars2[4]} sets the maximum number of iterations
#' @param optimmethod Method used in optimization of the likelihood. Current
#' default is 'subplex'. Alternative is 'simplex' (default of previous version)
#' @return \item{mu}{ gives the maximum likelihood estimate of mu}
#' \item{gamma_0}{ gives the maximum likelihood estimate of gamma_0}
#' \item{gamma_1}{ gives the maximum likelihood estimate of gamma_1}
#' \item{loglik}{ gives the maximum loglikelihood}
#' \item{df}{ gives the number of estimated parameters, i.e. degrees of feedom} 
#' \item{conv}{ gives a message on convergence of optimization; conv = 0 means convergence}
#' @param pchoice sets the p-value to optimize: \cr
#' pchoice == 0 corresponds to the sum of p_0f + p_1f \cr
#' pchoice == 1 corresponds to p_0f \cr
#' pchoice == 2 corresponds to p_1f \cr
#' @param edgeTList list of edge lengths that need to be succesively pruned; if
#' not specified, it will computed using compute_edgeTList
#' @param methode method used to solve the ODE. Either 'analytical' for the analytical
#' solution, 'Matrix' for matrix exponentiation using package Matrix or 'expm' using
#' package 'expm' or any of the numerical solvers, used in deSolve.
#' @param model model used. Default is 0 (standard null model). Other options are 1 (binary traits)
#' 2 (trinary environmental trait) or 3 (diversity-dependent colonization - beta version)
#' @param num_cycles the number of cycles of opimization. If set at Inf, it will
#' do as many cycles as needed to meet the tolerance set for the target function.
#' @param verbose Whether intermediate output should be printed. Default is FALSE.
#' @author Rampal S. Etienne
#' @seealso \code{\link{DAMOCLES_loglik}} \code{\link{DAMOCLES_sim}}
#' @references Pigot, A.L. & R.S. Etienne (2015). A new dynamic null model for
#' phylogenetic community structure. Ecology Letters 18: 153-163.
#' @keywords models
#' @export DAMOCLES_ML
DAMOCLES_ML <- DAMOCLES_all_ML <- function(
   phy,
   pa,
   initparsopt,
   idparsopt = 1:length(initparsopt),
   parsfix = NULL,
   idparsfix = NULL,
   idparsequal = NULL,
   pars2 = c(1E-3,1E-4,1E-5,1000),
   optimmethod = 'simplex',
   pchoice = 0,
   edgeTList = NULL,
   methode = 'analytical',
   model = 0,
   num_cycles = 1,
   verbose = FALSE)
{
  locatenode <- NULL
  if(model < 3)
  {
    edgeTList = DAMOCLES_check_edgeTList(phy,edgeTList)
  } else
  {
    locatenode <- DAMOCLES_check_locatenode(phy,locatenode)
  }
  patrait = pa
  if(model == -1)
  {
     numpars = 3
     namepars = c("mu","gamma_0","gamma_1")
  }
  if(model == 0 | model == 0.1 | model == 0.2)
  {
     numpars = 2
     namepars = c("mu","gamma_0")
  }
  if(model == 1)
  {
     numpars = 8
     namepars = c("mu_0","gamma_0","mu_1","gamma_1","q_0_to_1","q_1_to_0","r_0","r_1")
  }
  if(model == 2 | model == 2.2)
  {
     numpars = 12
     namepars = c("mu_1","gamma_1","mu_2","gamma_2","q_0_to_2","q_2_to_0","q_1_to_2","q_2_to_1","r_0","r_1","r_2","r_3")
  } 
  if(model == 3)
  {
    numpars = 3
    namepars = c("mu","gamma_0","K")
  }
  if(!is.matrix(patrait) & !is.data.frame(patrait))
  {
    warning('The trait data are assumed to be in the order of how they appear in the phylogeny.')
    patrait = matrix(c(phy$tip.label,patrait),nrow = length(patrait),ncol = 2 + (model > 0))
  }
  out2 = -1
  idpars = sort(c(idparsopt,idparsfix,idparsequal))
  if(length(idpars) != numpars)
  {
     stop("Incorrect number of parameters specified.\n")
  } else {
    if(prod(idpars == (1:numpars)) != 1)
    {
      cat("The parameters to be optimized and fixed are incoherent.\n")
    } else {
      cat('\nYou are running model',model,"\n")
      if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
      cat("You are optimizing",optstr,"\n")
      if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
      cat("You are fixing",fixstr,"\n")
      if(model > 0 & model < 3)
      {
         if(length(namepars[idparsequal]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsequal] }
         cat("You are setting",fixstr,"equal to their counterparts.\n")
      }
      trparsopt = initparsopt/(1 + initparsopt)
      if(model == 0.1 | model == 0.2)
      {
          parsfix = c(parsfix,0.5,0.5,0.5,0.5,0,0,0,0)
          idparsfix = c(idparsfix,3:10)
      }
      trparsfix = parsfix/(1 + parsfix)
      trparsfix[parsfix == Inf] = 1
      utils::flush.console()
      initloglik = DAMOCLES_all_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsequal = idparsequal,phy = phy,patrait = patrait,edgeTList = edgeTList,locatenode = locatenode,methode = methode,pchoice = pchoice,model = model, verbose = verbose)
      cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
      utils::flush.console()
      if(initloglik == -Inf)
      {
         cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
         out = -1
         return(out)
      } 
      cat("Optimizing ...\n")
      optimpars = pars2
      utils::flush.console()                  
      out = DDD::optimizer(optimmethod = optimmethod,
                           optimpars = optimpars,
                           fun = DAMOCLES_all_loglik_choosepar,
                           trparsopt = trparsopt,
                           trparsfix = trparsfix,
                           idparsopt = idparsopt,
                           idparsfix = idparsfix,
                           idparsequal = idparsequal,
                           phy = phy,
                           patrait = patrait,
                           edgeTList = edgeTList,
                           locatenode = locatenode,
                           methode = methode,
                           pchoice = pchoice,
                           model = model,
                           verbose = verbose,
                           num_cycles = num_cycles)
      if(out$conv > 0)
      {
        cat("Optimization has not converged. Try again with different starting values.\n")
      } else {
        MLtrpars = unlist(out$par)
        MLpars = DDD::roundn(MLtrpars/(1 - MLtrpars),14)
        out$par = list(MLpars)
        MLpars1 = rep(0,numpars)
        MLpars1[idparsopt] = MLpars
        if(length(idparsfix) != 0)
        { 
          MLpars1[idparsfix] = parsfix
        }
        if(length(idparsequal) != 0)
        { 
          MLpars1[idparsequal] = MLpars1[idparsequal - 1 - (model == 2 | model == 2.2) - (model == 1) * (idparsequal <= 4) + (model == 2 | model == 2.2) * (idparsequal == 10)]       
        }
        ML = as.numeric(unlist(out$fvalues))
        if(model == -1)
        {
           out2 = data.frame(mu = MLpars1[1], gamma_0 = MLpars1[2], gamma_1 = MLpars1[3], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
           s1 = sprintf('Maximum likelihood parameter estimates: mu: %f, gamma_0: %f, gamma_1: %f',MLpars1[1],MLpars1[2],MLpars1[3])
        }
        if(model == 0 | model == 0.1 | model == 0.2)
        {
           out2 = data.frame(mu = MLpars1[1], gamma_0 = MLpars1[2], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
           s1 = sprintf('Maximum likelihood parameter estimates: mu: %f, gamma_0: %f',MLpars1[1],MLpars1[2])
        }
        if(model == 1)
        {
           out2 = data.frame(mu_0 = MLpars1[1], gamma_0 = MLpars1[2], mu_1 = MLpars1[3], gamma_1 = MLpars1[4], q_0_to_1 = MLpars1[5], q_1_to_0 = MLpars1[6], r_0 = MLpars1[7], r_1 = MLpars1[8], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
           s1 = sprintf('Maximum likelihood parameter estimates: mu_0: %f, gamma_0: %f, mu_1: %f, gamma_1: %f, q_0_to_1: %f, q_1_to_0: %f, r_0: %f, r_1: %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8])
        }
        if(model == 2 | model == 2.2)
        {
           out2 = data.frame(mu_1 = MLpars1[1], gamma_1 = MLpars1[2], mu_2 = MLpars1[3], gamma_2 = MLpars1[4], q_0_to_2 = MLpars1[5], q_2_to_0 = MLpars1[6], q_1_to_2 = MLpars1[7], q_2_to_1 = MLpars1[8], r_0 = MLpars1[9], r_1 = MLpars1[10], r_2 = MLpars1[11], r_3 = MLpars1[12], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
           s1 = sprintf('Maximum likelihood parameter estimates: mu_1: %f, gamma_1: %f, mu_2: %f, gamma_2: %f, q_0_to_2: %f, q_2_to_0: %f, q_1_to_2: %f, q_2_to_1: %f, r_0: %f, r_1: %f, r_2: %f, r_3: %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10],MLpars1[11],MLpars1[12])
        }
        s2 = sprintf('Maximum loglikelihood: %f',out$fvalues)
        cat("\n",s1,"\n",s2,"\n\n")
      }
    }  
  }
  return(invisible(out2))  
}
