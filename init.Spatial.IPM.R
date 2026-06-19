e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.Spatial.IPM <- function(data,inits=NA,M=NA){
  n <- dim(data$y1)[1]
  if(n>M) stop("M must be larger than the number of captured/detected individuals.")
  X1 <- data$X1 
  X2 <- data$X2
  K1 <- data$K1
  K2 <- data$K2
  J1 <- unlist(lapply(X1,nrow))
  J1.max <- max(J1)
  J2 <- unlist(lapply(X2,nrow))
  J2.max <- max(J2)
  xlim <- data$xlim
  ylim <- data$ylim
  dSS <- data$dSS
  cells <- data$cells
  res <- data$res
  InSS <- data$InSS
  y1 <- array(0,dim=c(M,n.primary,J1.max))
  y1[1:n,1:n.primary,1:J1.max] <- data$y1 #all these guys must be observed
  y2 <- array(0,dim=c(M,n.primary,J2.max))
  y2[1:n,1:n.primary,1:J2.max] <- data$y2 #all these guys must be observed
  tel.z.states <- matrix(2,M,n.primary)
  tel.z.states[1:n,] <- data$tel.z.states
  tel.z.states[is.na(tel.z.states)] <- 2  #use 2 to code NA for nimble function
  
  #initialize z objects
  z.init <- matrix(0,M,n.primary)
  z.start.init <- z.stop.init <- rep(0,M)
  y12D <- apply(y1,c(1,2),sum)
  y22D <- apply(y2,c(1,2),sum)
  y2D <- y12D + y22D
  #add telemetry 
  for(i in 1:data$n.tel.inds){
    y2D[data$tel.ID[i],which(data$tel.z.states[data$tel.ID[i],]>0)] <- 1
  }
  #note, y2D is used in custom N/z update. must have scr detections and live telemetry years flagged >0
  for(i in 1:n){
    det.idx <- which(y2D[i,]>0)
    if(i %in% data$tel.ID){
      this.tel.ind <- which(data$tel.ID==i)
      det.idx <- sort(unique(c(det.idx,which(data$tel.z.states[i,]>0))))
    }
    if(length(det.idx)>0){
      z.start.init[i] <- min(det.idx)
      z.stop.init[i] <- max(det.idx)
      z.init[i,z.start.init[i]:z.stop.init[i]] <- 1
    }
  }
  z.super.init <- c(rep(1,n),rep(0,M-n))
  z.obs <- 1*(rowSums(y1)>0|rowSums(y2)>0) #indicator for "ever observed"
  
  #initialize N structures from z.init
  N.init <- colSums(z.init[z.super.init==1,])
  N.survive.init <- N.recruit.init <- rep(NA,n.primary-1)
  for(g in 2:n.primary){
    N.survive.init[g-1] <- sum(z.init[,g-1]==1&z.init[,g]==1&z.super.init==1)
    N.recruit.init[g-1] <- N.init[g]-N.survive.init[g-1]
  }
  
  #remaining SCR stuff to initialize
  #put X in ragged array
  #also make K1D, year by trap operation history, as ragged array.
  X1.nim <- array(0,dim=c(n.primary,J1.max,2))
  K1D1 <- matrix(0,n.primary,J1.max)
  X2.nim <- array(0,dim=c(n.primary,J2.max,2))
  K1D2 <- matrix(0,n.primary,J2.max)
  for(g in 1:n.primary){
    X1.nim[g,1:J1[g],1:2] <- X1[[g]]
    K1D1[g,1:J1[g]] <- rep(K1[g],J1[g])
    X2.nim[g,1:J2[g],1:2] <- X2[[g]]
    K1D2[g,1:J2[g]] <- rep(K2[g],J2[g])
  }
  
  #pull out state space with buffer around maximal trap dimensions
  s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
  idx <- which(rowSums(y1)>0|rowSums(y2)>0) #switch for those actually caught or telemetered
  idx <- sort(unique(c(idx,data$tel.ID))) #add telemetry individuals
  for(i in idx){
    trps <- matrix(0,nrow=0,ncol=2) #get locations of traps of capture across years for ind i
    for(g in 1:n.primary){
      if(sum(y1[i,g,])>0){
        trps.g <- matrix(X1.nim[g,which(y1[i,g,]>0),],ncol=2,byrow=FALSE)
        trps <- rbind(trps,trps.g)
      }
      if(sum(y2[i,g,])>0){
        trps.g <- matrix(X2.nim[g,which(y2[i,g,]>0),],ncol=2,byrow=FALSE)
        trps <- rbind(trps,trps.g)
      }
    }
    if(i%in%data$tel.ID){
      these.locs <- matrix(0,nrow=0,ncol=2) #get telemetry locs across years for tel.ID i
      this.tel.ind <- which(data$tel.ID==i)
      for(g in 1:data$n.tel.years[this.tel.ind]){
        locs.g <- data$locs[this.tel.ind,g,,]
        these.locs <- rbind(these.locs,locs.g)
      }
      trps <- rbind(trps,these.locs)
    }
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,] <- trps
    }
  }
  #in case s.init exactly on state space boundary
  eps <- 1e-5
  s.init[,1][s.init[,1]<xlim[1]] <- xlim[1] + eps
  s.init[,1][s.init[,1]>xlim[2]] <- xlim[2] - eps
  s.init[,2][s.init[,2]<ylim[1]] <- ylim[1] + eps
  s.init[,2][s.init[,2]>ylim[2]] <- ylim[2] - eps
  
  #If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
  e2dist  <-  function (x, y){
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
  }
  getCell  <-  function(s,res,cells){
    cells[trunc(s[1]/res)+1,trunc(s[2]/res)+1]
  }
  alldists <- e2dist(s.init,dSS)
  alldists[,InSS==0] <- Inf
  for(i in 1:M){
    this.cell <- cells[trunc(s.init[i,1]/res)+1,trunc(s.init[i,2]/res)+1]
    if(InSS[this.cell]==0){
      cands <- alldists[i,]
      new.cell <- which(alldists[i,]==min(alldists[i,]))
      s.init[i,] <- dSS[new.cell,]
    }
  }
  
  #check starting logProbs one session at a time
  for(g in 1:n.primary){
    #marking process
    D1 <- e2dist(s.init, X1.nim[g,1:J1[g],])
    pd1 <- p01[g]*exp(-D1*D1/(2*sigma[g]*sigma[g]))
    logProb <- array(0,dim=c(M,J1[g]))
    for(i in 1:M){
      for(j in 1:J1[g]){
        logProb[i,j] <- dbinom(y1[i,g,j],size=K1D1[g,j],prob=pd1[i,j],log=TRUE)
      }
    }
    if(!is.finite(sum(logProb)))stop(paste("Starting observation model 1 likelihood not finite, year",g))
    #sighting process
    D2 <- e2dist(s.init, X2.nim[g,1:J2[g],])
    pd2 <- p02[g]*exp(-D2*D2/(2*sigma[g]*sigma[g]))
    logProb <- array(0,dim=c(M,J2[g]))
    for(i in 1:M){
      for(j in 1:J2[g]){
        logProb[i,j] <- dbinom(y2[i,g,j],size=K1D2[g,j],prob=pd2[i,j],log=TRUE)
      }
    }
    if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood 2 not finite, year",g))
  }
  dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are

  
  return(list(s=s.init,z=z.init,N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
              N.super=n,z.start=z.start.init,z.stop=z.stop.init,z.super=z.super.init,
              z.obs=z.obs,y2D=y2D,y1=y1,y2=y2,X1=X1.nim,X2=X2.nim,K1D1=K1D1,K1D2=K1D2,
              dummy.data=dummy.data,tel.z.states=tel.z.states))
  
}