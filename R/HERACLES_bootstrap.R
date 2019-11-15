HERACLES_find_trait_stait_combinations = function(pa,paTrait)
{
	Np1 = length(which(pa  ==  1 & paTrait  ==  1))
	Np2 = length(which(pa  ==  1 & paTrait  ==  2))
	Na1 = length(which(pa  ==  0 & paTrait  ==  1))
	Na2 = length(which(pa  ==  0 & paTrait  ==  2))
	N0 = length(which(paTrait  ==  0))
	N1 = length(which(paTrait  ==  1))
	N2 = length(which(paTrait  ==  2))

	traitStateComb = c(Np1,Np2,Na1,Na2,N0,N2,N1,N2)

	return(traitStateComb)
}	

######################################################################################################################################
######################################################################################################################################

o.dectobin = function(y,ly = 0)
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
   res = c(rep(0,max(ly - length(res)),0),res)
   return(res)
}

######################################################################################################################################
######################################################################################################################################

combsUse = function(nRegional,nSample = 1000000)
{ 
	nposs = 2^nRegional
	nposs = nposs - 1
	combs = DDD::sample2(0:nposs,size = nSample,replace = TRUE)
	return(combs)
}	

######################################################################################################################################
######################################################################################################################################

#this now includes the option to calculate trait based metrics of phylogenetic structure by setting traitMetric = TRUE and including traitdist, a trait distance matrix
Heracles_ImportanceSampling = function(nSamples,n,regionalSpecies,S_regional,p,pa,phy,phydist,parsDAM,Mlist,model,pchoice,samptype,edgeObj,traitMetric,traitdist)
{	
	#create a matrix to store metrics of community structure
	#create a matrix to store the loglikelihood of each community and its sampling probability  
	loglikMatrix = matrix(ncol = 1 + 2 * (samptype  ==  'binomial'),nrow = nSamples)
	metricMatrix = matrix(ncol = 5,nrow = nSamples)
	
	#I have now included the possibility of calculate two trait based metrics
	if(samptype == 'binomial')
  {
     colnames(loglikMatrix) = c("P_i","X_i","n")
  } else
  {
     colnames(loglikMatrix) = c("P_i")
  }
  colnames(metricMatrix) = c("mnpd","mpd","mntd","mtd","n")
	
	DAMOCLES.samp = matrix(rep(0,n),nrow = 1)
	colnames(DAMOCLES.samp) = phy$tip.label
	
	#I ADDED THESE IF STATEMENTS BECAUSE IF THE REGIONAL POOL IS  LARGE THEN WE CANNOT USE combsUse BECAUSE WE HAVE TOO MANY POSSIBLE CONFIGURATIONS
	#FOR LARGE REGIONAL POOLS WE WILL HAVE TO USE BINOMIAL SAMPLING
	
	if(samptype == 'uniform')
  {
		if(S_regional > 30)
    {
  		 cat("The regional species pool is large. Uniform sampling will require exploring many community configurations and may exceed memory limitations. Use binomial sampling instead")
     	 break
		} else
    {
			 combs = combsUse(S_regional,nSamples)
		}
	}
		
 	pafoc = pa
	pafoc[,2] = as.character(rep(0,n))
	
	for(i in 1:nSamples)
  {	
	#if samptype is binomial, sample a local community richness S_loc from a binomial distribution with parameters S_regional and p
	 	 if(samptype == 'binomial')
     {
       	S_loc = stats::rbinom(n = 1,size = S_regional,prob = p)
     		sam = sample(c(rep(1,S_loc),rep(0,S_regional - S_loc)))
  		 	#calculate the sampling probability of the sampled community under a binomial distribution with parameters S_regional and p
			
	  	 	#the log probability of having local community richness S_loc is
		   	loglikMatrix[i,2] = log(stats::dbinom(S_loc, size = S_regional, prob = p)) - (dbinom(0, size = S_regional, prob = p) + stats::dbinom(1, size = S_regional, prob = p))
					
			  #the number of configurations with local richness S_loc is
		  	loglikMatrix[i,3] = lgamma(S_regional + 1) - lgamma(S_loc + 1) - lgamma(S_regional - S_loc + 1)
		 } else
	   if(samptype == 'uniform')
     {
    		sam = o.dectobin(combs[i],S_regional)
      	S_loc = sum(as.numeric(pa[,2]))
     }  
		 
	   #enter the sampled local community into our pa dataframe 	
		 pafoc[regionalSpecies,2] = sam
	   DAMOCLES.samp[1,] = pafoc[,2]
		
		 #calculate phylogenetic metric for the sampled community	
		 metricMatrix[i,1] = picante::mntd(DAMOCLES.samp,phydist)
		 metricMatrix[i,2] = picante::mpd(DAMOCLES.samp,phydist)
		
		 if(traitMetric == TRUE)
     {
        #calculate trait metrics for the sampled community	
		    metricMatrix[i,3] = mntd(traitdist,pafoc)
			  metricMatrix[i,4] = mtd(traitdist,pafoc)
	   }
		
		 metricMatrix[i,5] = S_loc
		  
		 #calculate the loglikelihood of the sampled community under DAMOCLES
		 loglikMatrix[i,1] = DAMOCLES_all_loglik(phy = phy,pa = pafoc,pars = parsDAM,pchoice = pchoice,edgeTList = edgeObj,analytical = TRUE, model = model,Mlist = Mlist)
				
		#print(i)
	}	

	resultsList = list()
	resultsList[[1]] = metricMatrix
	resultsList[[2]] = loglikMatrix
	return(resultsList)
}

######################################################################################################################################
######################################################################################################################################

HERACLES_extractCI = function(loglikMatrix,metricMatrix,ci_lower = 0.025,ci_upper = 0.975,traitMetric = FALSE)
{	
    # check whether there are weights for the logliks; if not, set them to 0
    loglikMatrix = as.matrix(loglikMatrix)
    if(min(dim(loglikMatrix))  ==  1)
    {
       loglikMatrix = cbind(loglikMatrix,rep(0,length(loglikMatrix)),rep(0,length(loglikMatrix)))
    }

    #only use samples that have more than 1 species locally present and no NA in the loglikelihood
    #tokeep = which(is.na(loglikMatrix[,1])  ==  FALSE & metricMatrix[,3] > 1)
    #NEEDED TO CHANGE metricMatrix[,3] TO metricMatrix[,5] AS WE HAVE TWO MORE METRICS BEFORE THE RICHNESS VALUES
    tokeep = which(is.na(loglikMatrix[,1])  ==  FALSE & metricMatrix[,5] > 1)
    if(length(tokeep) > 0)
    {
		   loglikMatrix = loglikMatrix[tokeep,]
		   metricMatrix = metricMatrix[tokeep,]
    }

    #calculate the likelihood ratio of each community
    likMatrix = loglikMatrix[,1] - loglikMatrix[,2] + loglikMatrix[,3]		
    likMatrix = likMatrix - max(likMatrix)
    likMatrix = exp(likMatrix) 
     
    nmets = 2
    if(traitMetric == TRUE)
    {
       nmets = 4
    }
    out = list()
    for(it in 1:nmets)
    {
	    out[[it]] = list()
    	out[[it]][[1]] = Hmisc::wtd.quantile(x = metricMatrix[,it],weights = likMatrix,probs = c(ci_lower,ci_upper),normwt = T,type = 'quantile')
    	out[[it]][[2]] = Hmisc::wtd.quantile(x = metricMatrix[,it],weights = likMatrix,probs = seq(0.01,0.99,by = 0.01),normwt = T,type = 'quantile')
	  }
    return(out)	
}	

######################################################################################################################################
######################################################################################################################################

#functions to calculate metrics of trait structure. These are used in "Heracles_ImportanceSampling"

#mean nearest trait distance
mntd = function(traitdist,pa)
{
	present = which(pa[,2] == "1")
	if(length(present)>1)
  {
		traitdistfoc = traitdist[present,present]
		mntd = mean(rowMins(traitdistfoc,na.rm=TRUE))
	} else {
		mntd = NA
	}	
	return(mntd)
}	

#mean trait distance

mtd = function(traitdist,pa)
{
	present = which(pa[,2] == "1")
	if(length(present)>1)
  {
		traitdistfoc = traitdist[present,present]
		mtd = mean(rowMeans(traitdistfoc,na.rm=TRUE))
	} else {
		mtd = NA
	}	
	return(mtd)
}		

######################################################################################################################################
######################################################################################################################################
	
#the two functions below can be ignored for now. They just calculate the expected mean nearest and mean trait distance under a random draw model

ses.mn_trait_d = function(traitdist,pa,runs)
{
	mntd = rep(NA,runs)
	for(i in 1:runs)
  {
		present = which(sample(pa[,2]) == "1")
		if(length(present) > 1)
    {
			 traitdistfoc = traitdist[present,present]
			 mntd[i] = mean(rowMins(traitdistfoc,na.rm = TRUE))
		} else {
			 mntd[i] = NA
		}
	}
	return(mntd)
}
	
ses.m_trait_d = function(traitdist,pa,runs)
{
	mtd = rep(NA,runs)
	for(i in 1:runs)
  {
		present = which(sample(pa[,2]) == "1")
		if(length(present) > 1)
    {
			traitdistfoc = traitdist[present,present]
			mtd[i] = mean(rowMeans(traitdistfoc,na.rm = TRUE))
		} else {
			mtd[i] = NA
		}
	}			
	return(mtd)	
}	
