e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

rtruncpois <- function(n,lambda,lower=0,upper=Inf){
  p.lo <- ppois(lower-1,lambda)
  p.hi <- ppois(upper,lambda)
  u <- runif(n,min=p.lo,max=p.hi)
  qpois(u,lambda)
}

sim.JS.SCR.Dcov.telSurv <- function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,
                                    gamma=NA,n.year=NA,phi=NA,
                                    p01=NA,p02=NA,sigma=NA,X1=NA,X2=NA,buff=buff,K1=NA,K2=NA,xlim=NA,ylim=NA,res=NA,
                                    n.tel.locs=NA,col.year.pars=NA,col.protocol=NA){
  if(col.year.pars[2]<1)stop("col.year.par[2] must be 1 or larger")
  #Population dynamics
  N <- rep(NA,n.year)
  N.recruit <- N.survive <- ER <- rep(NA,n.year-1)
  #get expected N in year 1 from D.cov parameters
  cellArea <- res^2
  lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
  lambda.y1 <- sum(lambda.cell)
  N[1] <- rpois(1,lambda.y1)
  if(N[1]==0)stop("Simulated starting population size of 0")

  #recreate some Dcov things so we can pass fewer arguments into this function
  x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
  y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
  dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
  cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
  n.cells <- nrow(dSS)
  n.cells.x <- length(x.vals)
  n.cells.y <- length(y.vals)

  #Easiest to increase dimension of z as we simulate bc size not known in advance.
  z <- matrix(0,N[1],n.year)
  z[1:N[1],1] <- 1
  cov <- rnorm(N[1],0,1) #simulate ind survival covariate for 1st year guys
  phi.mat <- matrix(NA,N[1],n.year-1)
  for(g in 2:n.year){
    #Simulate recruits
    ER[g-1] <- N[g-1]*gamma[g-1]
    N.recruit[g-1] <- rpois(1,ER[g-1])
    if(N.recruit[g-1]>0){
      #add recruits to z
      z.dim.old <- nrow(z)
      z <- rbind(z,matrix(0,nrow=N.recruit[g-1],ncol=n.year))
      z[(z.dim.old+1):(z.dim.old+N.recruit[g-1]),g] <- 1
      #and to phi matrix
      phi.mat <- rbind(phi.mat,matrix(NA,nrow=N.recruit[g-1],ncol=n.year-1))
    }
    phi.mat[,g-1] <- phi[g-1]
    idx <- which(z[,g-1]==1)
    z[idx,g] <- rbinom(length(idx),1,phi.mat[idx,g-1])
    N.survive[g-1] <- sum(z[,g-1]==1&z[,g]==1)
    N[g] <- N.recruit[g-1]+N.survive[g-1]
  }

  if(any(N.recruit+N.survive!=N[2:n.year]))stop("Simulation bug")
  if(any(colSums(z)!=N))stop("Simulation bug")

  #detection
  J1 <- unlist(lapply(X1,nrow)) #extract number of traps per year
  J1.max <- max(J1)
  J2 <- unlist(lapply(X2,nrow)) #extract number of traps per year
  J2.max <- max(J2)

  #simulate activity centers - fixed through time
  N.super <- nrow(z)
  # simulate a population of activity centers
  pi.cell <- lambda.cell/sum(lambda.cell)
  s.cell <- sample(1:n.cells,N.super,prob=pi.cell,replace=TRUE)
  #distribute activity centers uniformly inside cells
  s <- matrix(NA,nrow=N.super,ncol=2)
  for(i in 1:N.super){
    s.xlim <- dSS[s.cell[i],1] + c(-res,res)/2
    s.ylim <- dSS[s.cell[i],2] + c(-res,res)/2
    s[i,1] <- runif(1,s.xlim[1],s.xlim[2])
    s[i,2] <- runif(1,s.ylim[1],s.ylim[2])
  }

  pd1 <- y1 <- array(0,dim=c(N.super,n.year,J1.max))
  pd2 <- y2 <- array(0,dim=c(N.super,n.year,J2.max))
  for(g in 1:n.year){
    #marking (collar deployment) capture process
    D1 <- e2dist(s,X1[[g]])
    pd1[,g,1:J1[g]] <- p01[g]*exp(-D1*D1/(2*sigma[g]^2))
    for(i in 1:N.super){
      if(z[i,g]==1){
        y1[i,g,1:J1[g]] <- rbinom(J1[g],K1[g],pd1[i,g,1:J1[g]])
      }
    }
    #Method 2
    D2 <- e2dist(s,X2[[g]])
    pd2[,g,1:J2[g]] <- p02[g]*exp(-D2*D2/(2*sigma[g]^2))
    for(i in 1:N.super){
      if(z[i,g]==1){
        y2[i,g,1:J2[g]] <- rbinom(J2[g],K2[g],pd2[i,g,1:J2[g]])
      }
    }
  }
  
  #telemetry data
  #deploy collars to bears captured in marking process
  mark.caps <- 1*apply(y1,c(1,2),sum)
  tel.z.states <- z*NA
  col.lifes <- vector("list",n.year)
  #observed data, not true states (because we don't know if dead)
  eligible.states <- matrix(1,N.super,n.year) #eligible based on tel.z.states collaring history, may be dead and eligible
  tel.ID.g <- vector("list",n.year)
  for(g in 1:n.year){
    collar.g <- which(mark.caps[,g]>0&eligible.states[,g]==1)
    if(length(collar.g)>0){
      for(i in collar.g){
        col.life <- rtruncpois(1,lambda=col.year.pars[1],lower=col.year.pars[2],upper=col.year.pars[3])
        end.year <- min(g+col.life-1,n.year)
        tel.z.states[i,g:end.year] <- 1
        if(col.life>1&col.protocol==1){ #if we don't replace collars on capture, make ineligible
          if(g<n.year){
            eligible.states[i,(g+1):end.year] <- 0
          }
        }
      }
      tel.ID.g[[g]] <- collar.g
    }
  }
  #switch states to observed deaths when z==0 
  tel.z.states[which(tel.z.states==1&z==0)] <- 0
  #if you observe a death, fill in 0s to the end
  for(i in 1:N.super){
    idx <- which(tel.z.states[i,]==0)
    if(length(idx)>0){
      tel.z.states[i,max(idx):n.year] <- 0
    }
  }
  #simulate telemetry locations
  if(n.tel.locs>0&sum(y1)>0){
    n.tel.years <- rowSums(tel.z.states==1,na.rm=TRUE)
    tel.ID <- which(n.tel.years>0)
    n.tel.years <- n.tel.years[tel.ID]
    n.tel.inds <- sum(n.tel.years>0)
    tel.year <- matrix(NA,n.tel.inds,n.year)
    max.n.tel.years <- max(n.tel.years)
    locs <- array(NA,dim=c(n.tel.inds,max.n.tel.years,n.tel.locs,2))
    for(i in 1:n.tel.inds){
      tel.year[i,1:n.tel.years[i]] <- which(tel.z.states[tel.ID[i],]==1)
      for(g in 1:n.tel.years[i]){
        #if adding movement, reference correct s years
        locs[i,g,,] <- c(rnorm(n.tel.locs,s[tel.ID[i],1],sigma[tel.year[i,g]]),
                         rnorm(n.tel.locs,s[tel.ID[i],2],sigma[tel.year[i,g]]))
      }
    }
    n.locs.ind <- apply(!is.na(locs[,,,1]),c(1,2),sum)
    if(dim(locs)[2]==1){
      n.locs.ind <- matrix(rowSums(n.locs.ind),ncol=1)
    }
  }else{
    print("no individuals captured, no telemetry")
    locs <- tel.ID <- tel.year <- NA
    n.tel.inds <- 0
    n.tel.years <- NA
    n.locs.ind <- NA
  }

  #store true data for model building/debugging
  truth <- list(y1=y1,y2=y2,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,
                N.super=N[1]+sum(N.recruit),z=z,s=s,s.cell=s.cell)

  #discard undetected or untelemetered individuals
  y.all <- apply(y1,c(1,2),sum) + apply(y2,c(1,2),sum) #N x n.year
  y.all[tel.z.states==1] <- 1 #add telemetry observation
  
  keep.idx <- which(rowSums(y.all)>0)
  y1 <- y1[keep.idx,,]
  y2 <- y2[keep.idx,,]
  tel.z.states <- tel.z.states[keep.idx,]
  eligible.states <- eligible.states[keep.idx,]
  s <- s[keep.idx,]
  s.cell <- s.cell[keep.idx]
  
  #renumber tel.ID
  tel.ID <- match(tel.ID,keep.idx)
  for(g in 1:n.year){
    if(length(tel.ID.g[[g]])>0){
      tel.ID.g[[g]] <- match(tel.ID.g[[g]],keep.idx)
    }
  }
  
  return(list(y1=y1,y2=y2,N=N,N.recruit=N.recruit,N.survive=N.survive,
              X1=X1,X2=X2,K1=K1,K2=K2,n.year=n.year,
              tel.z.states=tel.z.states,locs=locs,n.tel.inds=n.tel.inds,n.tel.years=n.tel.years,
              n.locs.ind=n.locs.ind,tel.ID=tel.ID,tel.ID.g=tel.ID.g,tel.year=tel.year,
              xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
              n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,s=s,
              D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.y1=lambda.y1,
              truth=truth))
}
