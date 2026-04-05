#telemetry survival vector distribution
dSurvivalTel <- nimbleFunction(
  run = function(x = double(1), z = double(1), z.super = double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){
      for(g in 1:length(z)){
        if(x[g]!=2){ #2 codes NA (not yet collared or uninformatively censored), 1 is alive, 0 is dead
          #z invalid if conflicts with tel.
          if(x[g] != z[g]){
            logProb <- -Inf
          }
        }
      }
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

rSurvivalTel <- nimbleFunction(
  run = function(n = integer(0), z = double(1),z.super = double(0)) {
    returnType(double(1))
    return(rep(0,length(z)))
  }
)

#telemetry location vector distribution
dNormVector <- nimbleFunction(
  run = function(x = double(2), s = double(1), sigma = double(0), n.locs.ind = double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(n.locs.ind>0){
      for(i in 1:n.locs.ind){
        logProb <- logProb + dnorm(x[i,1], mean = s[1], sd = sigma, log = TRUE)
        logProb <- logProb + dnorm(x[i,2], mean = s[2], sd = sigma, log = TRUE)
      }
    }
    return(logProb) 
  })

rNormVector <- nimbleFunction(
  run = function(n = integer(0), s = double(1), sigma = double(0), n.locs.ind = double(0)) {
    returnType(double(2))
    return(matrix(0,n.locs.ind,2))
  }
)

dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0),log = integer(0)) {
    returnType(double(0))
    logProb <- log(pi.cell)
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0)) {
    returnType(double(0))
    return(0)
  }
)

#this is used to restrict likelihood evaluation to only the years relevant for survival for each individual
dSurvival <- nimbleFunction(
  run = function(x = double(1), phi = double(1), z.start = double(0), z.stop = double(0),
                 z.super = double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){
      n.year <- length(phi)+1
      #extract first and last survival event years
      surv.start <- z.start+1
      surv.stop <- z.stop+1 #count death events, first z[i,]=0
      if(surv.start <= n.year){ #if surv.start beyond last year, no survival events, logProb=0
        if(surv.stop > n.year){ #but can't survive past n.year
          surv.stop <- n.year 
        }
        for(g in surv.start:surv.stop){ #sum logprob over survival event years
          logProb <- logProb + dbinom(x[g], size = 1, p = phi[g-1], log = TRUE)
        }
      }
    }
    return(logProb)
  }
)

#make dummy random vector generator to make nimble happy
rSurvival <- nimbleFunction(
  run = function(n = integer(0),phi = double(1), z.start = double(0), z.stop = double(0), z.super = double(0)){
    returnType(double(1))
    n.year <- length(phi)
    return(rep(0,n.year))
  }
)

GetDetectionProb <- nimbleFunction(
  run = function(s = double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0), z.super=double(0)){ 
    returnType(double(1))
    if(z.super==0 | z.super==1&z==0){
      return(rep(0,J)) #skip calculation if not is superpop, or in superpop, but not alive in this year
    }
    if(z==1){ #otherwise calculate
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dBinomialVector <- nimbleFunction(
  run = function(x = double(1), pd = double(1), K = double(1), z = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z.super==0 | z.super==1&z==0){#skip calculation if not is superpop, or in superpop, but not alive in this year
      return(0)
    }else{
      logProb <- sum(dbinom(x, size = K, p = pd, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBinomialVector <- nimbleFunction(
  run = function(n = integer(0), pd = double(1), K = double(1), z = double(0), z.super = double(0)) {
    returnType(double(1))
    J <- nimDim(pd)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

#all z updates live here
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    J1 <- control$J1
    J2 <- control$J2
    y2D <- control$y2D
    n.sampled.m1 <- control$n.sampled.m1
    n.sampled.m2 <- control$n.sampled.m2
    sampled.years.m1 <- control$sampled.years.m1
    sampled.years.m2 <- control$sampled.years.m2
    z.super.ups <- control$z.super.ups
    n.year <- control$n.year
    z.obs <- control$z.obs
    z.nodes <- control$z.nodes
    y1.nodes <- control$y1.nodes
    y2.nodes <- control$y2.nodes
    pd1.nodes <- control$pd1.nodes
    pd2.nodes <- control$pd2.nodes
    tel.z.states.nodes <- control$tel.z.states.nodes
    N.nodes <- control$N.nodes
    ER.nodes <- control$ER.nodes
    N.survive.nodes <- control$N.survive.nodes
    N.recruit.nodes <- control$N.recruit.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    #precompute entry counts
    entry.counts.curr <- rep(0, n.year+1)
    for(g in 1:n.year){
      entry.counts.curr[g] <- sum(model$z.start==g & model$z.super==1)
    }
    entry.counts.curr[n.year + 1] <- sum(model$z.super==0)
    
    #1) Detected guy updates: z.start, z.stop
    # 1a) z start update (z.stop update below): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&y2D[i,1]==0){ #for detected guys, skip if observed 1st year
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        N.curr <- model$N
        N.recruit.curr <- model$N.recruit
        dets <- which(y2D[i,]>0)
        first.det <- min(dets)
        lp.start <- rep(-Inf,n.year)
        i.idx.m1 <- seq(i,M*n.sampled.m1,M) #used to reference correct y1 and pd1 nodes
        i.idx.m2 <- seq(i,M*n.sampled.m2,M) #used to reference correct y2 and pd2 nodes
        for(g in 1:first.det){ #must be recruited in year with first detection or before
          z.start.prop <- g
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.year)
          z.prop[g:first.det] <- 1 #must be alive until first detection
          if(first.det < n.year){
            z.prop[(first.det+1):n.year] <- z.curr[(first.det+1):n.year] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          
          #update N, N.recruit, N.survive. These individuals always in superpopulation
          #1) Update N
          model$N <<- N.curr - z.curr + z.prop
          #2) Update N.recruit
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd1.nodes[i.idx.m1]) #update pd1 nodes when a z changes
          model$calculate(pd2.nodes[i.idx.m2]) #update pd2 nodes when a z changes
          #get these logProbs
          lp.N1 <- model$calculate(N.nodes[1])
          lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.y1 <- model$calculate(y1.nodes[i.idx.m1])
          lp.y2 <- model$calculate(y2.nodes[i.idx.m2])
          lp.surv <- model$calculate(z.nodes[i])
          lp.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          # Add the full multinomial coefficient prior log-prob for this proposed configuration
          entry.counts.prop <- entry.counts.curr
          #z.super always 1 for detected guys
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          lp.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
          lp.start[g] <- lp.N1 + lp.N.recruit + lp.y1 + lp.y2 + lp.surv + lp.tel.z.states + lp.prior
        }
        maxlp <- max(lp.start) #deal with overflow
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        
        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr #set back to original
        
        if(model$z.start[i]!=z.start.prop){#if proposal is same as current, no need to replace anything
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.year)
          z.prop[model$z.start[i]:first.det] <- 1 #must be alive until first detection
          if(first.det < n.year){
            z.prop[(first.det+1):n.year] <- z.curr[(first.det+1):n.year] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd1.nodes[i.idx.m1]) #update pd1 nodes
          model$calculate(pd2.nodes[i.idx.m2]) #update pd2 nodes
          #update these logProbs
          model$calculate(y1.nodes[i.idx.m1])
          model$calculate(y2.nodes[i.idx.m2])
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.sampled.m1){
            g.use <- sampled.years.m1[g]
            for(j in 1:J1[g.use]){
              mvSaved["pd1",1][i,g.use,j] <<- model[["pd1"]][i,g.use,j]
            }
          }
          for(g in 1:n.sampled.m2){
            g.use <- sampled.years.m2[g]
            for(j in 1:J2[g.use]){
              mvSaved["pd2",1][i,g.use,j] <<- model[["pd2"]][i,g.use,j]
            }
          }
          #recompute entry counts
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.sampled.m1){
            g.use <- sampled.years.m1[g]
            for(j in 1:J1[g.use]){
              model[["pd1"]][i,g.use,j] <<- mvSaved["pd1",1][i,g.use,j]
            }
          }
          for(g in 1:n.sampled.m2){
            g.use <- sampled.years.m2[g]
            for(j in 1:J2[g.use]){
              model[["pd2"]][i,g.use,j] <<- mvSaved["pd2",1][i,g.use,j]
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y1.nodes[i.idx.m1])
          model$calculate(y2.nodes[i.idx.m2])
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
        }
      }
    }
    
    #1b) z stop update (z.start update above): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&y2D[i,n.year]==0){ #for detected guys, skip if observed in final year
        z.curr <- model$z[i,]
        z.stop.curr <- model$z.stop[i]
        N.curr <- model$N
        dets <- which(y2D[i,]>0)
        last.det <- max(dets)
        lp.stop <- rep(-Inf,n.year)
        i.idx.m1 <- seq(i,M*n.sampled.m1,M) #used to reference correct y1 and pd1 nodes
        i.idx.m2 <- seq(i,M*n.sampled.m2,M) #used to reference correct y2 and pd2 nodes
        for(g in (last.det):n.year){ #can't die on or before year of last detection
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.year)
          z.prop[last.det:g] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          #update N, number of recruits does not change going backwards
          model$N <<- N.curr - z.curr + z.prop
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd1.nodes[i.idx.m1]) #update pd1 nodes when a z changes
          model$calculate(pd2.nodes[i.idx.m2]) #update pd2 nodes when a z changes
          #get these logProbs
          lp.N1 <- model$calculate(N.nodes[1])
          lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.y1 <- model$calculate(y1.nodes[i.idx.m1])
          lp.y2 <- model$calculate(y2.nodes[i.idx.m2])
          lp.surv <- model$calculate(z.nodes[i])
          lp.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          #no prior term, z.stop update does not change it
          lp.stop[g] <- lp.N1 + lp.N.recruit + lp.y1 + lp.y2 + lp.surv + lp.tel.z.states
        }
        maxlp <- max(lp.stop) #deal with overflow
        prop.probs <- exp(lp.stop-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.stop.prop <- rcat(1,prop.probs)
        model$z.stop[i] <<- z.stop.curr #set back to original
        if(model$z.stop[i]!=z.stop.prop){#if proposal differs from current
          model$z.stop[i] <<- z.stop.prop
          z.prop <- rep(0,n.year)
          z.prop[last.det:model$z.stop[i]] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd1.nodes[i.idx.m1]) #update pd1 nodes when z changes
          model$calculate(pd2.nodes[i.idx.m2]) #update pd2 nodes when z changes
          #update these logProbs
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y1.nodes[i.idx.m1])
          model$calculate(y2.nodes[i.idx.m2])
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.sampled.m1){
            g.use <- sampled.years.m1[g]
            for(j in 1:J1[g.use]){
              mvSaved["pd1",1][i,g.use,j] <<- model[["pd1"]][i,g.use,j]
            }
          }
          for(g in 1:n.sampled.m2){
            g.use <- sampled.years.m2[g]
            for(j in 1:J2[g.use]){
              mvSaved["pd2",1][i,g.use,j] <<- model[["pd2"]][i,g.use,j]
            }
          }
        }else{
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.sampled.m1){
            g.use <- sampled.years.m1[g]
            for(j in 1:J1[g.use]){
              model[["pd1"]][i,g.use,j] <<- mvSaved["pd1",1][i,g.use,j]
            }
          }
          for(g in 1:n.sampled.m2){
            g.use <- sampled.years.m2[g]
            for(j in 1:J2[g.use]){
              model[["pd2"]][i,g.use,j] <<- mvSaved["pd2",1][i,g.use,j]
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y1.nodes[i.idx.m1])
          model$calculate(y2.nodes[i.idx.m2])
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
        }
      }
    }
    #2) undetected guy update. Only if in the superpopulation.
    # Metropolis-Hastings, Propose z vectors from priors
    #entry counts current after z.start update
    for(i in 1:M){
      if(z.obs[i]==0&model$z.super[i]==1){
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        z.stop.curr <- model$z.stop[i]
        i.idx.m1 <- seq(i,M*n.sampled.m1,M) #used to reference correct y1 and pd1 nodes
        i.idx.m2 <- seq(i,M*n.sampled.m2,M) #used to reference correct y2 and pd2 nodes
        #get forwards recruitment probabilities
        recruit.probs.for <- c(model$lambda.y1,model$ER)
        recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
        #get initial logProbs
        lp.initial.entry <- model$getLogProb(N.nodes[1])
        lp.initial.entry <- lp.initial.entry + model$getLogProb(N.recruit.nodes)
        lp.initial.y1 <- model$getLogProb(y1.nodes[i.idx.m1])
        lp.initial.y2 <- model$getLogProb(y2.nodes[i.idx.m2])
        lp.initial.surv <- model$getLogProb(z.nodes[i])
        lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[i])
        log.prior.curr <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr + 1)))

        #track proposal probs - survival is symmetric, but not recruitment and detection
        log.prop.for <- log.prop.back <- 0

        #simulate recruitment
        z.start.prop <- rcat(1,recruit.probs.for)
        z.prop <- rep(0,n.year)
        z.prop[z.start.prop] <- 1
        log.prop.for <- log.prop.for + log(recruit.probs.for[z.start.prop])

        #simulate survival
        if(z.start.prop < n.year){#if you don't recruit in final year
          for(g in (z.start.prop+1):n.year){
            z.prop[g] <- rbinom(1,1,model$phi[i,g-1]*z.prop[g-1])
            log.prop.for <- log.prop.for + dbinom(z.prop[g],1,model$phi[i,g-1]*z.prop[g-1],log=TRUE)
          }
        }
        z.on.prop <- which(z.prop==1)
        z.stop.prop <- max(z.on.prop)
        model$z[i,] <<- z.prop
        model$z.start[i] <<- z.start.prop
        model$z.stop[i] <<- z.stop.prop

        #update N, N.recruit, N.survive only if individual is in superpopulation
        #1) Update N
        model$N <<- model$N - z.curr + z.prop
        #2) Update N.recruit
        if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
          model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
        }
        if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
          model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
        }
        #3) Update N.survive
        model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g

        model$calculate(ER.nodes) #update ER when N updated
        model$calculate(pd1.nodes[i.idx.m1])
        model$calculate(pd2.nodes[i.idx.m2])
        #get proposed logProbs
        lp.proposed.entry <- model$calculate(N.nodes[1])
        lp.proposed.entry <- lp.proposed.entry + model$calculate(N.recruit.nodes)
        lp.proposed.y1 <- model$calculate(y1.nodes[i.idx.m1])
        lp.proposed.y2 <- model$calculate(y2.nodes[i.idx.m2])
        lp.proposed.surv <- model$calculate(z.nodes[i])
        lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[i])

        # Full multinomial coefficient prior for proposed configuration
        entry.counts.prop <- entry.counts.curr
        entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
        entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
        log.prior.prop <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))

        #get backwards proposal probs
        recruit.probs.back <- c(model$lambda.y1,model$ER)
        recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
        log.prop.back <- log.prop.back + log(recruit.probs.back[z.start.curr])
        if(z.start.curr < n.year){#if you don't recruit in final year
          for(g in (z.start.curr+1):n.year){
            log.prop.back <- log.prop.back + dbinom(z.curr[g],1,model$phi[i,g-1]*z.curr[g-1],log=TRUE)
          }
        }
        lp.initial.total <- lp.initial.entry + lp.initial.y1 + lp.initial.y2 +
          lp.initial.surv + lp.initial.tel.z.states + log.prior.curr
        lp.proposed.total <- lp.proposed.entry + lp.proposed.y1 + lp.proposed.y2 +
          lp.proposed.surv + lp.proposed.tel.z.states + log.prior.prop

        #MH step
        log_MH_ratio <- (lp.proposed.total + log.prop.back) - (lp.initial.total + log.prop.for)
        accept <- decide(log_MH_ratio)

        if(accept) {
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.sampled.m1){
            g.use <- sampled.years.m1[g]
            for(j in 1:J1[g.use]){
              mvSaved["pd1",1][i,g.use,j] <<- model[["pd1"]][i,g.use,j]
            }
          }
          for(g in 1:n.sampled.m2){
            g.use <- sampled.years.m2[g]
            for(j in 1:J2[g.use]){
              mvSaved["pd2",1][i,g.use,j] <<- model[["pd2"]][i,g.use,j]
            }
          }
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.sampled.m1){
            g.use <- sampled.years.m1[g]
            for(j in 1:J1[g.use]){
              model[["pd1"]][i,g.use,j] <<- mvSaved["pd1",1][i,g.use,j]
            }
          }
          for(g in 1:n.sampled.m2){
            g.use <- sampled.years.m2[g]
            for(j in 1:J2[g.use]){
              model[["pd2"]][i,g.use,j] <<- mvSaved["pd2",1][i,g.use,j]
            }
          }
          #set these logProbs back
          model$calculate(N.recruit.nodes)
          model$calculate(N.nodes[1])
          model$calculate(y1.nodes[i.idx.m1])
          model$calculate(y2.nodes[i.idx.m2])
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
        }
      }
    }
    #3) update z.super: Metropolis-Hastings
    #entry counts current coming out of undetected ind update
    for(up in 1:z.super.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected individual
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z.super==1)
        non.init <- length(z.on)
        pick <- rcat(1,rep(1/non.init,non.init))
        pick <- z.on[pick]
        if(z.obs[pick]==1){ #is this individual detected?
          reject <- TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
        }
        if(!reject){
          z.start.curr <- model$z.start[pick]
          z.curr <- model$z[pick,]

          #p select off guy
          log.p.select.for <- log(1/non.init)
          #log multinomial coefficient prior
          log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
          pick.idx.m1 <- seq(pick,M*n.sampled.m1,M)
          pick.idx.m2 <- seq(pick,M*n.sampled.m2,M)

          #get initial logProbs (survival logProb does not change)
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y1 <- model$getLogProb(y1.nodes[pick.idx.m1])
          lp.initial.y2 <- model$getLogProb(y2.nodes[pick.idx.m2])
          lp.initial.surv <- model$getLogProb(z.nodes[pick])
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[pick])

          # propose new N.super/z.super/z.start/z.stop
          model$N.super <<-  model$N.super - 1
          model$z.super[pick] <<- 0
          model$z.start[pick] <<- 0
          model$z.stop[pick] <<- 0
          model$z[pick,] <<- rep(0,n.year)

          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N - z.curr
          #2) Update N.recruit
          if(z.start.curr > 1){ #if wasn't in pop in year 1
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd1.nodes[pick.idx.m1]) #update pd1 nodes when z changes
          model$calculate(pd2.nodes[pick.idx.m2]) #update pd2 nodes when z changes

          #Reverse proposal probs
          recruit.probs.back <- c(model$lambda.y1, model$ER)
          recruit.probs.back <- recruit.probs.back / sum(recruit.probs.back)
          log.prop.back <- log(recruit.probs.back[z.start.curr])
          if(z.start.curr < n.year){
            for(g in (z.start.curr+1):n.year){
              log.prop.back <- log.prop.back + dbinom(z.curr[g],1,model$phi[pick,g-1]*z.curr[g-1],log=TRUE)
            }
          }

          #get proposed logProbs for N, N.recruit, and y
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y1 <- model$calculate(y1.nodes[pick.idx.m1]) #will always be 0
          lp.proposed.y2 <- model$calculate(y2.nodes[pick.idx.m2]) #will always be 0
          lp.proposed.surv <- model$calculate(z.nodes[pick]) #will always be 0
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[pick]) #will always be 0


          lp.initial.total <- lp.initial.N + lp.initial.y1 + lp.initial.y2 +
            lp.initial.N.recruit + lp.initial.surv + lp.initial.tel.z.states
          lp.proposed.total <- lp.proposed.N + lp.proposed.y1 + lp.proposed.y2 +
            lp.proposed.N.recruit + lp.proposed.surv + lp.proposed.tel.z.states

          #backwards prior and select probs
          #move from class z.start.curr in z.super==0 to class g in z.super==1
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[n.year + 1] <- entry.counts.prop[n.year + 1] + 1

          #p select on guy
          noff.back <- sum(model$z.super == 0)
          log.p.select.back <- log(1/noff.back)
          #log multinomial coefficient prior
          log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
          log.prop.for <- 0
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back) -
            (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for)

          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1] <<- model[["z.super"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            for(g in 1:n.sampled.m1){
              g.use <- sampled.years.m1[g]
              for(j in 1:J1[g.use]){
                mvSaved["pd1",1][pick,g.use,j] <<- model[["pd1"]][pick,g.use,j]
              }
            }
            for(g in 1:n.sampled.m2){
              g.use <- sampled.years.m2[g]
              for(j in 1:J2[g.use]){
                mvSaved["pd2",1][pick,g.use,j] <<- model[["pd2"]][pick,g.use,j]
              }
            }
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]] <<- mvSaved["z.super",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            for(g in 1:n.sampled.m1){
              g.use <- sampled.years.m1[g]
              for(j in 1:J1[g.use]){
                model[["pd1"]][pick,g.use,j] <<- mvSaved["pd1",1][pick,g.use,j]
              }
            }
            for(g in 1:n.sampled.m2){
              g.use <- sampled.years.m2[g]
              for(j in 1:J2[g.use]){
                model[["pd2"]][pick,g.use,j] <<- mvSaved["pd2",1][pick,g.use,j]
              }
            }
            #set these logProbs back
            model$calculate(y1.nodes[pick.idx.m1])
            model$calculate(y2.nodes[pick.idx.m2])
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            model$calculate(z.nodes[pick])
            model$calculate(tel.z.states.nodes[pick])
          }
        }
      }else{#add
        if(model$N.super[1] < M){ #cannot update if z.super maxed out. Need to raise M
          z.off <- which(model$z.super==0)
          noff.init <- length(z.off)
          pick <- rcat(1,rep(1/noff.init,noff.init)) #select one of these individuals
          pick <- z.off[pick]
          pick.idx.m1 <- seq(pick,M*n.sampled.m1,M)
          pick.idx.m2 <- seq(pick,M*n.sampled.m2,M)

          non.init <- sum(model$z.super == 1)

          #p select off guy
          log.p.select.for <- log(1/noff.init)

          #log multinomial coefficient prior
          log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))

          #get initial logProbs (survival logProb does not change)
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y1 <- model$getLogProb(y1.nodes[pick.idx.m1]) #will always be 0
          lp.initial.y2 <- model$getLogProb(y2.nodes[pick.idx.m2]) #will always be 0
          lp.initial.surv <- model$getLogProb(z.nodes[pick]) #will always be 0
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[pick]) #will always be 0

          # Propose new z.start for the new on individual
          recruit.probs.for <- c(model$lambda.y1, model$ER)
          recruit.probs.for <- recruit.probs.for / sum(recruit.probs.for)
          z.start.prop <- rcat(1, recruit.probs.for)  # propose entry cohort
          log.prop.for <- log(recruit.probs.for[z.start.prop])
          model$z.start[pick] <<- z.start.prop

          # Simulate survival path
          model$z[pick,] <<- 0 # initialize to 0
          model$z[pick, z.start.prop] <<- 1
          if(z.start.prop < n.year){
            for(g in (z.start.prop+1):n.year){
              model$z[pick, g] <<- rbinom(1, 1, model$phi[pick, g-1] * model$z[pick, g-1])
              log.prop.for <- log.prop.for + dbinom(model$z[pick, g], 1, model$phi[pick, g-1] * model$z[pick, g-1], log=TRUE)
            }
          }
          # Update z.stop
          z.on.prop <- which(model$z[pick,] == 1)
          z.stop.prop <- max(z.on.prop)
          model$z.stop[pick] <<- z.stop.prop

          #propose new N/z
          model$N.super <<-  model$N.super + 1
          model$z.super[pick] <<- 1

          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N + model$z[pick,]
          #2) Update N.recruit
          if(model$z.start[pick] > 1){ #if wasn't in pop in year 1
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year] - model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd1.nodes[pick.idx.m1]) #update pd1 nodes when z changes
          model$calculate(pd2.nodes[pick.idx.m2]) #update pd2 nodes when z changes
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y1 <- model$calculate(y1.nodes[pick.idx.m1])
          lp.proposed.y2 <- model$calculate(y2.nodes[pick.idx.m2])
          lp.proposed.surv <- model$calculate(z.nodes[pick])
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[pick])

          lp.initial.total <- lp.initial.N + lp.initial.y1 + lp.initial.y2 +
            lp.initial.N.recruit + lp.initial.surv + lp.initial.tel.z.states
          lp.proposed.total <- lp.proposed.N + lp.proposed.y1 + lp.proposed.y2 +
            lp.proposed.N.recruit + lp.proposed.surv + lp.proposed.tel.z.states

          #backwards prior and select probs
          #move from class g in z.super==0 to class g in z.super==1
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          entry.counts.prop[n.year + 1] <- entry.counts.prop[n.year + 1] - 1

          #p select on guy
          non.back <- sum(model$z.super == 1)
          log.p.select.back <- log(1/non.back)
          #log multinomial coefficient prior
          log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
          log.prop.back <- 0

          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back) -
            (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for)

          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1] <<- model[["z.super"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            for(g in 1:n.sampled.m1){
              g.use <- sampled.years.m1[g]
              for(j in 1:J1[g.use]){
                mvSaved["pd1",1][pick,g.use,j] <<- model[["pd1"]][pick,g.use,j]
              }
            }
            for(g in 1:n.sampled.m2){
              g.use <- sampled.years.m2[g]
              for(j in 1:J2[g.use]){
                mvSaved["pd2",1][pick,g.use,j] <<- model[["pd2"]][pick,g.use,j]
              }
            }
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]] <<- mvSaved["z.super",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            for(g in 1:n.sampled.m1){
              g.use <- sampled.years.m1[g]
              for(j in 1:J1[g.use]){
                model[["pd1"]][pick,g.use,j] <<- mvSaved["pd1",1][pick,g.use,j]
              }
            }
            for(g in 1:n.sampled.m2){
              g.use <- sampled.years.m2[g]
              for(j in 1:J2[g.use]){
                model[["pd2"]][pick,g.use,j] <<- mvSaved["pd2",1][pick,g.use,j]
              }
            }
            #set these logProbs back
            model$calculate(y1.nodes[pick.idx.m1])
            model$calculate(y2.nodes[pick.idx.m2])
            model$calculate(z.nodes[pick])
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            model$calculate(tel.z.states.nodes[pick])
          }
        }
      }
    }
    #copy back to mySaved to update logProbs.
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)