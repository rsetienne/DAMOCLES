dec2bin = function(y,ly)
{
  stopifnot(length(y) == 1, mode(y) == 'numeric')
  q1 = (y / 2) %/% 1
  r = y - q1 * 2
  res = c(r)
  while(q1 >= 1)
  {
    q2 = (q1 / 2) %/% 1
    r = q1 - q2 * 2
    q1 = q2
    res = c(r, res)
  }
  res = c(rep(0,ly - length(res)),res)
  return(res)
}

dec2binmat = function(y)
{
  numrows = 2^y
  res = matrix(0,numrows,y)
  for(i in 0:(numrows-1))
  {
    res[i + 1,] = dec2bin(i,y)
  }
  return(res)
}

bin2dec <- function(y)
{
  res <- y %*% 2^((length(y) - 1):0)
  return(as.numeric(res))
}

kimat <- function(dec2binmatk)
{
  ki <- matrix(0,dim(dec2binmatk)[1],dim(dec2binmatk)[1])
  for(i in 2:dim(dec2binmatk)[1])
  {
    locationones <- which(dec2binmatk[i,] == 1)
    for(j in 1:length(locationones))
    {
      dec2binmatki <- dec2binmatk[i,]
      dec2binmatki[locationones[j]] = 0
      j2 <- 1 + bin2dec(dec2binmatki)
      ki[i,j2] <- 1 
    }
  }
  return(ki)
}

gan <- function(ga0,K,n)
{
  return(ga0 * pmax(rep(0,length(n)),(1 - n/K)))
}

DAMOCLES_DD_loglik_rhs <- function(
  t,
  p,
  pars
)
{
  dp <- pars %*% p
  return(list(dp))
}

DAMOCLES_DD_odeintr <- function(
  initprobs,
  t0,
  t,
  totmat,
  methode
  )
{
   Sys.setenv( "PKG_CXXFLAGS"="-std=c++11" )
   sys <- '
for(int i = 0; i != N; ++i)
{
  dxdt[i] = 0;
  for (int j = 0; j != N; ++j)
  dxdt[i] += pars[i * N + j] * x[j];
}
'
   n_pars <- length(initprobs)
   compile_sys("DAMOCLES_DD_branch", sys, n_pars^2, sys_dim = n_pars,atol = 1e-16, rtol = 1e-16, method = methode)
   totvec <- t(totmat)
   dim(totvec) <- c(n_pars^2,1)
   DAMOCLES_DD_branch_set_params(totvec)
   probs <- as.numeric(DAMOCLES_DD_branch(initprobs,t,(t - t0)/10,t0)[11,2:(1 + n_pars)])
   return(probs)
}

DAMOCLES_DD_FORTRAN <- function(
  initprobs,
  t0,
  t,
  totmat,
  methode
)
{
  #system("R CMD SHLIB d:/data/ms/allocatables.f")
  #system("R CMD SHLIB d:/data/ms/DAMOCLES/DAMOCLES_DD_loglik_rhs.f")
  dyn.load(paste("d:/data/ms/DAMOCLES/DAMOCLES_DD_loglik_rhs", .Platform$dynlib.ext, sep = ""))
  n_pars <- length(initprobs)
  totvec <- t(totmat)
  dim(totvec) <- c(n_pars^2,1)
  probs <- ode(y = initprobs, parms = c(n_pars + 0.), rpar = totvec, 
             times = c(t0,t), func = "runmod", initfunc = "initmod", 
             ynames = c("SV"), dimens = n_pars^2, nout = 1, outnames = c("Sum"), 
             dllname = "DAMOCLES_DD_loglik_rhs",method = methode)[2,2:(n_pars + 1)]
  dyn.unload(paste("d:/data/ms/DAMOCLES/DAMOCLES_DD_loglik_rhs", .Platform$dynlib.ext, sep = ""))
  return(probs)
}

DAMOCLES_DD_integrate <- function(
  t0,
  t,
  initprobs,
  pars,
  direction = 'forward',
  methode = 'lsoda'
)
{
  dime <- length(initprobs)
  N <- log2(dime)
  dec2binmatk <- dec2binmat(N)
  posc <- Matrix::rowSums(dec2binmatk)
  kimin <- kimat(dec2binmatk)
  kiplus <- t(kimin)
  mu <- pars[1] * posc
  mumat <- matrix(pars[1],nrow = dime,ncol = dime)
  if(direction == 'forward')
  {
     ga <- gan(pars[2],pars[3],posc)
     gamat <- t(replicate(dime,ga))
     totmat <- kimin * gamat + kiplus * mumat
     totmat <- totmat - diag(colSums(totmat))
  } else
  {
     ga <- gan(pars[2],pars[3],posc - 1)
     gamat <- t(replicate(dime,ga))
     totmat <- kiplus * gamat + kimin * mumat
     totmat <- totmat - diag(rowSums(totmat))
  }
  if(methode == 'analytical')
  {
     difft <- abs(t - t0)
     probs <- as.vector((Matrix::expm(totmat * difft)) %*% initprobs)
  } else if(is.element(methode,c('euler','rk4','rk54','rk54_a','rk5','rk5_a','rk5_i','rk78','rk78_a','abN','abmN','bs','bsd')))
  {
     probs <- DAMOCLES_DD_odeintr(initprobs = initprobs, t0 = t0, t = t, totmat = totmat, methode = methode)
  } else if(methode == 'experimental')
  {
     probs <- DAMOCLES_DD_FORTRAN(initprobs = initprobs, t0 = t0, t = t, totmat = totmat, methode = 'lsoda')
  } else
  {
     probs <- ode(y = initprobs,times = c(t0,t), func = DAMOCLES_DD_loglik_rhs, parms = totmat, method = methode,rtol = 1E-10,atol = 1E-16)[2,2:(1 + dime)]
  }
  return(probs)
}

locate_node <- function(phy)
{
  phy$node.label <- NULL
  nNode <- Nnode(phy)
  nodeorder <- rep(0,nNode)
  oldlabels <- phy$tip.label
  phy$tip.label <- as.character(1:Ntip(phy))
  for (i in 1:nNode)
  {
    ntips <- Ntip(phy)
    if(ntips > 2)
    {
      nodeDepths <- dist.nodes(phy)[(ntips + 1):(2 * ntips - 1), ntips + 1]
      w = which(nodeDepths == max(nodeDepths))[1]
      toDrop <- clade.members(as.numeric(names(nodeDepths[w])), phy, tip.labels = TRUE, include.nodes = TRUE)$tips
    } else
    {
      nodeDescendants <- clade.members.list(phy, tip.labels = TRUE)
      sisterNodes <- which(sapply(nodeDescendants, length) == 2)
      sisterDepths <- nodeDepths[sisterNodes]
      w <- which(sisterDepths == max(sisterDepths))[1]
      toDrop <- nodeDescendants[[sisterNodes[w]]]
    }
    nodeorder[i] <- min(as.numeric(toDrop))
    phy <- suppressWarnings(drop.tip(phy, toDrop, trim.internal = FALSE))
    phy$tip.label[which(phy$tip.label == "NA")] <- paste(nodeorder[i], sep = "")
    if(nodeorder[i] < (ntips - 1))
    {
      k = nodeorder[i] + 1
      for(j in (nodeorder[i] + 2):ntips)
      {
         phy$tip.label[which(phy$tip.label == j)] <- k
         k = k + 1
      }
    }
  }
  return(nodeorder)
}

DAMOCLES_DD_node <- function(
  initprobs,
  locnode,
  direction = 'forward'
)
{
  dime <- length(initprobs)
  if(direction == 'forward')
  {
    N <- log2(dime) + 1
    dec2binmatk <- dec2binmat(N)
    probs <- rep(0,dime * 2)
    for(i in 1:(dime * 2))
    {
      sumtwostates <- sum(dec2binmatk[i,locnode:(locnode + 1)])
      if(sumtwostates == 0) # the branches are both 0
      {
        probs[i] <- initprobs[1 + bin2dec(dec2binmatk[i,-(locnode + 1)])]
      } else
      if(sumtwostates == 1) # the branches are 1 and 0
      {
        oldstate <- append(x = dec2binmatk[i,], values = 1, after = locnode + 1)
        oldstate <- oldstate[-c(locnode:(locnode + 1))]
        probs[i] <- 0.5 * initprobs[1 + bin2dec(oldstate)]
      } # twostates = c(1,1) gets prob = 0
    }
  } else # direction backward
  {
    N <- log2(dime) - 1
    dec2binmatk <- dec2binmat(N)
    probs <- rep(0,dime/2)
    for(i in 1:(dime/2))
    {
      if(dec2binmatk[i,locnode] == 0) # the branch has a 0
      {
         probs[i] <- initprobs[1 + bin2dec(append(x = dec2binmatk[i,],values = 0,after = locnode))]
      } else # the branch has a 1
      {
         probs[i] <- 0.5 * initprobs[1 + bin2dec(append(x = dec2binmatk[i,],values = 0,after = locnode))] +
                     0.5 * initprobs[1 + bin2dec(append(x = dec2binmatk[i,],values = 0,after = locnode - 1))]
      }
    }
  }  
  return(probs)
}

DAMOCLES_check_locatenode <- function(
  phy,
  locatenode
)
{
  if(is.null(locatenode))
  {
    locatenode <- locate_node(phy)
  }
  return(locatenode)
}

DAMOCLES_DD_loglik <- function(
  phy,
  pa,
  pars,
  pchoice = c(1,0),
  direction = 'backward',
  locatenode = NULL,
  methode = 'lsoda'
)
{
  if(pars[3] < sum(as.numeric(pa[,2])))
  {
    loglik = -Inf
    return(loglik)
  }
  locatenode <- DAMOCLES_check_locatenode(phy,locatenode)
  if(length(pchoice) == 1)
  {
    pchoice = (pchoice == 0) * c(0.5,0.5) + (pchoice == 1) * c(1,0) + (pchoice == 2) * c(0,1)
  }
  S <- Ntip(phy)
  pastate <- rep(0,S) 
  for(i in 1:S)
  {
    pastate[i] <- as.numeric(pa[which(pa[,1] == phy$tip.label[i]),2])
  }
  brts <- as.numeric(sort(branching.times(phy)))
  loglik <- 0
  if(brts[1] != 0)
  {
    brts <- c(0,brts)
  }
  if(direction == 'forward')
  {
    locatenode <- rev(locatenode)
    brts <- -rev(abs(brts))
    #for(state in 0:1)
    #{
       state <- which(pchoice == 1) - 1
       probs <- c(1 - state,state)
       for(i in 1:(S - 1))
       {
         probs <- DAMOCLES_DD_node(initprobs = probs, locnode = locatenode[i], direction = direction)  
         cp <- DAISIE:::checkprobs2(NULL,loglik,probs)
         loglik <- cp[[1]]
         probs <- cp[[2]]
         probs <- DAMOCLES_DD_integrate(t0 = brts[i],t = brts[i + 1], initprobs = probs, pars = pars, direction = direction, methode = methode)
         cp <- DAISIE:::checkprobs2(NULL,loglik,probs)
         loglik <- cp[[1]]
         probs <- cp[[2]]
       }
       loglik <- loglik + log(probs[1 + bin2dec(pastate)])
    #}
  } else # backward
  {
    probs <- rep(0,2^S)
    probs[1 + bin2dec(pastate)] = 1
    for(i in 1:(S - 1))
    {
      probs <- DAMOCLES_DD_integrate(t0 = brts[i],t = brts[i + 1], initprobs = probs, pars = pars, direction = direction, methode = methode)
      cp <- DAISIE:::checkprobs2(NULL,loglik,probs)
      loglik <- cp[[1]]
      probs <- cp[[2]]
      probs <- DAMOCLES_DD_node(initprobs = probs, locnode = locatenode[i], direction = direction)  
      cp <- DAISIE:::checkprobs2(NULL,loglik,probs)
      loglik <- cp[[1]]
      probs <- cp[[2]]
    }
    loglik <- loglik + log(pchoice %*% probs)
  }
  return(as.numeric(loglik))
}