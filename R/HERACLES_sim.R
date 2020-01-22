HERACLES_sim = function(phy, mu_1,gamma_1,mu_2, gamma_2,q_0_to_2,q_2_to_0,q_1_to_2,q_2_to_1,r_0,r_1,r_2,r_3,
root.state, root.trait.state, plotit = FALSE,keepExtinct = FALSE) 
{	
	
	if(plotit == TRUE)
	{
    ape::plot.phylo(phy,main="trait 2 (green), trait 1 (red), trait 0 (blue); dark (present), light (absent)")
	}
  dn = DDD::roundn(ape::dist.nodes(phy),6)
  ntips = length(phy$tip.label)
  nbranch = 2 * ntips - 2
  patable = data.frame(p = phy$edge[, 1], d = phy$edge[, 2], 
  age.start = rep(NA, nbranch), age.end = rep(NA, nbranch), 
  extant = rep(0, nbranch), state = rep(NA, nbranch), traits = rep(NA, nbranch))
  patable$age.start = dn[patable$p, ntips + 1]
  patable$age.end = dn[patable$d, ntips + 1]
  patable$age.end[which(patable$d < length(phy$tip.label) +1)] = max(patable$age.end)
  	    	
  dedge = which(patable$p  ==  min(patable$p))
  patable$extant[dedge] = 1  
    	    	
	if(root.state == 1 & root.trait.state  == 0)
	{
     stop("Root state and root trait state are incoherent")
	}
	if(r_2 == 0 & r_3 > 0)
	{
     stop("r_3 cannot > 0 when r_2 = 0")
	}    
  patable$tips[which(patable$d<=ntips)] = 1
  patable$y[which(patable$tips == 1)] = patable$d[which(patable$tips == 1)]
    
  if(plotit == TRUE)
  {
  	 while(length(stats::na.omit(patable$y)) < length(patable$y))
     {
	    	for(i in 1:length(patable$y))
        {
	    	   if(is.na(patable$y[i]) == FALSE)
           {
		    	    focalp = patable$p[i]
		    	    focaly = patable$y[which(patable$p == focalp)]
   		    	  if(length(stats::na.omit(focaly)) == 2)
              {
	   	   			   focald = which(patable$d == focalp)
		    		     patable$y[focald] = (focaly[1]+focaly[2])/2
	   	   		  }
			     }
		    }
	    }
  }
	
	parentTrait = root.trait.state   
	
	patable = HERACLES_speciation_sim(patable,parentTrait,dedge,r_0,r_1,r_2,r_3)
	
	parentState = root.state
	
	patable = HERACLES_assign_state(patable,parentState,dedge)                 

	tstep = min(patable$age.start)
	while(tstep < max(patable$age.end))
	{
		    
 		 #numbers of species in each state 
     extant = which(patable$extant  ==  1)
     extant.id = patable$d[extant]
     pa = patable$state[extant]
     paTrait = patable$traits[extant]
       
     #transition rates
     rates = c(mu_1,mu_2,gamma_1,gamma_2,q_0_to_2,q_2_to_0,q_1_to_2,q_2_to_1)
     
     traitStateComb = HERACLES_find_trait_stait_combinations(pa,paTrait)
          
     #gillespie algorithm    
     rates.real = (rates * traitStateComb)
     totalRate = sum(rates.real)
     tstep0 = tstep
     wt = stats::rexp(n = 1,rate = as.numeric(totalRate))
     tstep1 = tstep0 + wt
     tstep = min(tstep1, min(patable$age.end[extant]))
                	
		 if(tstep < max(patable$age.end)) 
		 {
		  	if(tstep < min(patable$age.end[extant]))
			  {
			   	 event = sample(x = 1:8, size = 1, replace = F,prob = c(rates.real))
				   if(event  ==  1) #local extinction of species locally present and endemic to region 1
           {
              patable$state[extant[DDD::sample2(which(pa  ==  1 & paTrait  ==  1), 1)]] = 0
           } else
    		   if(event  ==  2) #local extinction of species locally present and in both regions
           {
              patable$state[extant[DDD::sample2(which(pa  ==  1 & paTrait  ==  2), 1)]] = 0
           } else 
		       if(event  ==  3) #immigration of species locally absent and endemic to region 1
           {
              patable$state[extant[DDD::sample2(which(pa  ==  0 & paTrait  ==  1), 1)]] = 1
           } else
           if(event  ==  4) #immigration of species locally absent and in both regions
           {
              patable$state[extant[DDD::sample2(which(pa  ==  0 & paTrait  ==  2), 1)]] = 1
           } else 
			  	 if(event  ==  5) #entry of species in region 0 into region 1
           {
              patable$traits[extant[DDD::sample2(which(paTrait  ==  0), 1)]] = 2
           } else
			     if(event  ==  6) #this may seem odd because the regional extinction rate also influences the local extinction rate but is consistent with the model
           { 
					    #reg1extinct = DDD::sample2(which(pa  ==  0 & paTrait  ==  2), 1) #extinction of species from region 1
				 	    reg1extinct = DDD::sample2(which(paTrait  ==  2), 1) #extinction of species from region 1
              patable$traits[extant[reg1extinct]] = 0
              patable$state[extant[reg1extinct]] = 0
       	   } else     								
				   if(event  ==  7) #entry of species in region 1 into region 0
           {
              patable$traits[extant[DDD::sample2(which(paTrait  ==  1), 1)]] = 2
           } else
				   if(event  ==  8) #extinction of species from region 0
           {
              patable$traits[extant[DDD::sample2(which(paTrait  ==  2), 1)]] = 1
           } 
        } else
        {
    	   	 focedge = extant[which(patable$age.end[extant]  ==  min(patable$age.end[extant]))]
		    	 for(fe in unique(focedge))
           {
              dedge = which(patable$p  ==  patable$d[fe])
					    if(length(dedge) > 0)
              {
       		     	 patable$extant[fe] = 0
        	  		 patable$extant[dedge] = 1
          			
          			 parentTrait = patable$traits[fe]      
	
			      		 patable = HERACLES_speciation_sim(patable,parentTrait,dedge,r_0,r_1,r_2,r_3)
                   
               	 parentState = patable$state[fe] 
	
						     patable = HERACLES_assign_state(patable,parentState,dedge)                 
  				    } else
              {
      		       patable$extant[fe] = 0
              }
           }
        }
	   }
  }
	patablelist = list()
	if(keepExtinct  ==  FALSE)
	{
		patable = patable[which(patable$extant  ==  1), ]
	}
	patablelist[[1]] = patable
	return(patablelist)
}


# for a species endemic to a region (either 1 or 0) speciation could conceivably be associated with two kinds of event
# first, within region speciation { 0->0,0 or 1->1,1 } occurs with probabilities 1-r_0 or 1-r_1
# second, peripatric speciation, whereby the other region is colonised and this event results in speciation
# { 0->0,1 or 1->1,0 } occurs with probabilies r_0 and r_1
 
# event 11 speciation of species in both regions
					
# there are four posibilities here
# first, the split could occur in the middle of the region so that two daughters with state 2 are produced { 2->2,2 }
# imagine a species in the stretching from central america to the atlantic forests being split by the amazon river
# both daughters would be present in amazonia but not endemic to it 
					
# we denote these possibilities with probabilities r_2 and r_3 respectively
					
# so { 2->2,2 } occurs with probability 1-r_2
# and { 2->1,2 or 2->0,2 or 2->0,1 } with probability r_2
					
# and { 2->1,2 or 2->0,2 } occurs with probability r_3
# and { 2->0,1 } occurs with probability 1-r_3

######################################################################################################################################
######################################################################################################################################

HERACLES_speciation_sim = function(patable,parentTrait,dedge,r_0,r_1,r_2,r_3)
{        
     
  
  if(parentTrait  ==  0) # event 9 speciation of species endemic to region 0
  {
 	   patable$traits[dedge] = sample(c(0, sample(c(0,  1), 1, prob = c(1 - r_0, r_0))))
	} else
  if(parentTrait  ==  1) # event 10 speciation of species endemic to region 1
  {
 	   patable$traits[dedge] = sample(c(1, sample(c(1,  0), 1, prob = c(1 - r_1, r_1))))
  } else
	if(parentTrait  ==  2) # event 11 speciation of global species	
  {
		 specmode = sample(c("within","between"),1,prob = c(1 - r_2, r_2))
     if(specmode == "within")
     {
     		patable$traits[dedge] = c(2,2)
     } else
     {
        specmode = sample(c("within","between"),1,prob = c(r_3, 1 - r_3))
        if(specmode == "within")
        {
       	   patable$traits[dedge] = sample(c(sample(0:1,1),2))
        } else
        { 
           patable$traits[dedge] = sample(0:1)
        }	
     }
	}
	return(patable)
}	
	     
######################################################################################################################################
######################################################################################################################################

HERACLES_assign_state = function(patable,parentState,dedge)
{        
	#if the root species is locally present then we have to ensure that the daughter species
	#which inherits this local presence is also present in region 1
	if(parentState > 0)
  {
		#both daughters may be present in the region so just select one of these at random
		#this also works if just a single daughter is regionally present
		InheritLocal = DDD::sample2(which(patable$traits[dedge] > 0),1)
    if(InheritLocal  ==  1)
    {
       patable$state[dedge] = c(1,0)
    } else {
   	   patable$state[dedge] = c(0,1)
		}
	} else
  {
		#if the root species is initially absent then these are both zeroes		
    patable$state[dedge] = c(0,0)
	}
	return(patable)
}		     

		