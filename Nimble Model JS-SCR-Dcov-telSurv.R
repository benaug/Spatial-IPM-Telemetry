NimModel <- nimbleCode({
  #Density covariates
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  # D.intercept <- exp(D.beta0)*cellArea
  D.intercept <- D0*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells]) #separate this component so s's do not depend on D.intercept
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  
  ##Abundance##
  lambda.y1 <- D.intercept*pi.denom #Expected starting population size
  N[1] ~ dpois(lambda.y1) #Realized starting population size
  for(g in 2:n.year){
    N[g] <- N.survive[g-1] + N.recruit[g-1] #yearly abundance
    #N.recruit and N.survive information also contained in z/z.start + z.stop
    #N.recruit has distributions assigned below, but survival distributions defined on z
  }
  N.super <- N[1] + sum(N.recruit[1:(n.year-1)]) #size of superpopulation
  
  #Recruitment
  gamma.fixed ~ dunif(0,2)
  for(g in 1:(n.year-1)){
    #gamma[g] ~ dunif(0,2) #yearly recruitment priors
    gamma[g] <- gamma.fixed #fixed recruitment rate
    ER[g] <- N[g]*gamma[g] #yearly expected recruits
    N.recruit[g] ~ dpois(ER[g]) #yearly realized recruits
  }
  
  #Individual covariates
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]]) #categorical likelihood for this cell, equivalent to zero's trick
  }
  
  #Survival (phi must have M x n.year - 1 dimension for custom updates to work)
  #without individual or year effects, use for loop to plug into phi[i,g]
  phi.fixed ~ dunif(0,1)
  for(i in 1:M){
    #phi[g] ~ dunif(0,1) #yearly survival priors
    for(g in 1:(n.year-1)){#plugging same individual phi's into each year for custom update
      phi[i,g] <- phi.fixed
    }
    #survival likelihood (bernoulli) that only sums from z.start to z.stop
    z[i,1:n.year] ~ dSurvival(phi=phi[i,1:(n.year-1)],z.start=z.start[i],z.stop=z.stop[i],z.super=z.super[i])
    #telemetry survival likelihood
    #fixes z states, 1=known alive, 0=known dead, NA=unknown
    #currently assumes censoring is uninformative
    tel.z.states[i,1:n.year] ~ dSurvivalTel(z=z[i,1:n.year],z.super=z.super[i])
  }
  
  ##Detection##
  sigma.fixed ~ dunif(0,10)
  for(g in 1:n.year){ #sigma informed by data except in years with no capture effort and no telemetry
    # sigma[g] ~ dunif(0,10) #sigma varies by year, shared across methods
    sigma[g] <- sigma.fixed #sigma fixed across years, shared across methods
  }
  #Method 1 - only loop over sampled years
  for(g in 1:n.sampled.m1){
    p01[g] ~ dunif(0,1)#p01 varies by year
    for(i in 1:M){
      pd1[i,sampled.years.m1[g],
          1:J1[sampled.years.m1[g]]] <- GetDetectionProb(s=s[i,1:2],
                                                         X=X1[sampled.years.m1[g],1:J1[sampled.years.m1[g]],1:2],
                                                         J=J1[sampled.years.m1[g]],sigma=sigma[sampled.years.m1[g]],
                                                         p0=p01[g],z=z[i,sampled.years.m1[g]],
                                                         z.super=z.super[i])
      y1[i,sampled.years.m1[g],
         1:J1[sampled.years.m1[g]]] ~ dBinomialVector(pd=pd1[i,sampled.years.m1[g],1:J1[sampled.years.m1[g]]],
                                                      K=K1D1[sampled.years.m1[g],1:J1[sampled.years.m1[g]]],
                                                      z=z[i,sampled.years.m1[g]],z.super=z.super[i])
    }
  }
  
  #Method 2 - only loop over sampled years
  for(g in 1:n.sampled.m2){
    p02[g] ~ dunif(0,1) #p02 varies by year
    for(i in 1:M){
      pd2[i,sampled.years.m2[g],
          1:J2[sampled.years.m2[g]]] <- GetDetectionProb(s=s[i,1:2],
                                                         X=X2[sampled.years.m2[g],1:J2[sampled.years.m2[g]],1:2],
                                                         J=J2[sampled.years.m2[g]],sigma=sigma[sampled.years.m2[g]],
                                                         p0=p02[g],z=z[i,sampled.years.m2[g]],
                                                         z.super=z.super[i])
      y2[i,sampled.years.m2[g],
         1:J2[sampled.years.m2[g]]] ~ dBinomialVector(pd=pd2[i,sampled.years.m2[g],1:J2[sampled.years.m2[g]]],
                                                      K=K1D2[sampled.years.m2[g],1:J2[sampled.years.m2[g]]],
                                                      z=z[i,sampled.years.m2[g]],z.super=z.super[i])
    }
  }
  #Telemetry informs activity centers and sigma
  for(i in 1:n.tel.inds){
    for(g in 1:n.tel.years[i]){
      locs[i,g,1:n.locs.ind[i,g],1:2] ~ dNormVector(s=s[tel.ID[i],1:2],sigma=sigma[tel.year[i,g]],
                                                    n.locs.ind=n.locs.ind[i,g])
    }
  }
})

#custom updates:
#1) for detected individuals: update z.start, then update z.stop
#2) for undetected individuals: update entire z vectors
#3) N.super/z.super update