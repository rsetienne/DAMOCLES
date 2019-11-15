DAMOCLES_all_loglik_choosepar = function(trparsopt,idparsopt,trparsfix,idparsfix,idparsequal,phy,patrait,edgeTList,locatenode,pchoice,methode,model)
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
         loglik = DAMOCLES_all_loglik(phy = phy,pa = patrait,pars = pars1,pchoice = pchoice,edgeTList = edgeTList,methode = methode,model = model,Mlist = NULL)
       } else
       {
         loglik = DAMOCLES_DD_loglik(phy = phy,pa = patrait,pars = pars1,pchoice = pchoice,locatenode = locatenode,methode = methode)
       }
   }
   return(loglik)
}

DAMOCLES_ML <- DAMOCLES_all_ML <- function(
   phy,
   pa,
   initparsopt,
   idparsopt = 1:length(initparsopt),
   parsfix = NULL,
   idparsfix = NULL,
   idparsequal = NULL,
   pars2 = c(1E-3,1E-4,1E-5,1000),
   optimmethod = 'subplex',
   pchoice = 0,
   edgeTList = NULL,
   locatenode = NULL,
   methode = 'analytical',
   model = 0)
{
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
  if(is.matrix(patrait) == 0)
  {
     patrait = matrix(c(phy$tip.label,patrait),nrow = length(patrait),ncol = 2 + (model > 0))
  }
  options(warn = -1)
  out2 = -1
  idpars = sort(c(idparsopt,idparsfix,idparsequal))
  if(length(idpars) != numpars)
  {
     cat("Incorrect number of parameters specified.\n")
  } else {
    if(prod(idpars == (1:numpars)) != 1)
    {
      cat("The parameters to be optimized and fixed are incoherent.\n")
    } else {
      cat('You are running model',model,"\n")
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
      flush.console()
      initloglik = DAMOCLES_all_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsequal = idparsequal,phy = phy,patrait = patrait,edgeTList = edgeTList,locatenode = locatenode,methode = methode,pchoice = pchoice,model = model)
      cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
      flush.console()
      if(initloglik == -Inf)
      {
         cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
         out = -1
         return(out)
      } 
      cat("Optimizing ...\n")
      optimpars = pars2
      flush.console()                  
      out = DDD::optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = DAMOCLES_all_loglik_choosepar,trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsequal = idparsequal,phy = phy,patrait = patrait,edgeTList = edgeTList,locatenode = locatenode,methode = methode,pchoice = pchoice,model = model)
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
