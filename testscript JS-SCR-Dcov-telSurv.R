library(nimble)
library(coda)
source("sim.JS.SCR.Dcov.telSurv.R")
source("init.Spatial.IPM.R") #initializes data for nimble
source("Nimble Model JS-SCR-Dcov-telSurv.R")
source("Nimble Functions JS-SCR-Dcov-telSurv.R") #contains custom distributions and updates
source("sSampler Dcov.R") # activity center sampler that proposes from prior when z.super=0.
source("mask.check.R")

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")

#must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)

n.year <- 7 #number of years
gamma <- rep(0.2,n.year-1) #yearly per-capita recruitment
phi <- rep(0.8,n.year-1) #yearly survival
p01 <- rep(0.1,n.year) #yearly detection probabilities at activity center, telemetry capture
p02 <- rep(0.1,n.year) #yearly detection probabilities at activity center, SCR capture
sigma <- rep(0.5,n.year) #yearly detection function scale
#Number of occasions per year per method
#to skip sampling by a method in a year, set its K=0
K1 <- c(0,2,0,0,2,0,0) #yearly sampling occasions, telemetry capture process
K2 <- c(5,0,0,5,0,0,5) #yearly sampling occasions, SCR detection process
if(length(K1)!=length(K2))stop("K1 and K2 must be same length")
if(length(K1)!=n.year)stop("K1 and K2 must be of length n.year")

#Make some traps
#Need to make both X1 and X2, the trap locations for method 1 and method 2
X1 <- vector("list",n.year) #one trapping array per year
for(g in 1:n.year){ #using same trapping array every year here
  X1[[g]] <- as.matrix(expand.grid(3:11,3:11))
}
X2 <- vector("list",n.year) #one trapping array per year
for(g in 1:n.year){ #using same trapping array every year here
  X2[[g]] <- as.matrix(expand.grid(3:11,3:11))
}

## Habitat covariate stuff##
#buffer maximal trap extent
X1.all <- X2.all <- matrix(NA,nrow=0,ncol=2)
for(g in 1:n.year){
  X1.all <- rbind(X1.all,X1[[g]])
  X2.all <- rbind(X2.all,X2[[g]])
}
X.all <- rbind(X1.all,X2.all)

buff <- 2 #state space buffer. Buffers maximal x and y dimensions of X below across years
#get x and y extent by buffering state space
xlim <- range(X.all[,1])+c(-buff,buff)
ylim <- range(X.all[,2])+c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift
for(g in 1:n.year){
  X1[[g]][,1] <- X1[[g]][,1]-x.shift
  X1[[g]][,2] <- X1[[g]][,2]-y.shift
  X2[[g]][,1] <- X2[[g]][,1]-x.shift
  X2[[g]][,2] <- X2[[g]][,2]-y.shift
}
X.all[,1] <- X.all[,1]-x.shift
X.all[,2] <- X.all[,2]-y.shift
X1.all[,1] <- X1.all[,1]-x.shift
X1.all[,2] <- X1.all[,2]-y.shift
X2.all[,1] <- X2.all[,1]-x.shift
X2.all[,2] <- X2.all[,2]-y.shift

res <- 0.25 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#simulate a D.cov, higher cov.pars for large scale cov
library(fields)
set.seed(13235) #change seed for new D.cov
grid <- list(x=x.vals,y=y.vals) 
obj <- Exp.image.cov(grid=grid,aRange=2,setup=TRUE)
D.cov <- sim.rf(obj)
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)

image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
points(X.all,pch=4)

#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(1,length(D.cov))
InSS[dSS.tmp[,1]<1&dSS.tmp[,2]<1] <- 0
InSS[dSS.tmp[,1]<1&dSS.tmp[,2]>11] <- 0
InSS[dSS.tmp[,1]>11&dSS.tmp[,2]<1] <- 0
InSS[dSS.tmp[,1]>11&dSS.tmp[,2]>11] <- 0

image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y),main="Habitat",col=cols1)

#Density covariates
D.beta0 <- -0.65
D.beta1 <- 0.5
#what is implied expected N in state space?
lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

#Telemetry settings
n.tel.locs <- 15 #number of locs per individual
#col.year.pars and col.protocol determine how collars are deployed and how long they last (not model parameters)
col.year.pars <- c(2,2,3) #parameters for truncated poisson, lambda, lower truncation, upper truncation
#data simulator requires lower bound be 1 or higher. 1 means it fails before 2nd year
#collar lifetime frequencies for col.year.pars
table(rtruncpois(10000,lambda=col.year.pars[1],lower=col.year.pars[2],upper=col.year.pars[3]))/10000
#collaring protocol: #1) never replace a collar if currently collared on capture 2) always replace
col.protocol <- 1 

#simulate some data
set.seed(133455) #change seed for new data set
data <- sim.JS.SCR.Dcov.telSurv(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,
                                gamma=gamma,n.year=n.year,phi=phi,p01=p01,p02=p02,
                                sigma=sigma,X1=X1,X2=X2,,K1=K1,K2=K2,xlim=xlim,ylim=ylim,res=res,
                                n.tel.locs=n.tel.locs,col.year.pars=col.year.pars,col.protocol)

#what is observed data
# str(data$y1) #method 1 capture history: n x n.year x J1
# str(data$y2) #method 2 detection history: n x n.year x J2
# str(data$tel.z.states) #telemetry survival observations: n x n.year
# str(data$locs) #telemetry locations: n.tel.inds x n.year x n.tel.locs x 2

#inds captured by method
colSums(apply(data$y1,c(1,2),sum)>0)
colSums(apply(data$y2,c(1,2),sum)>0)
#by all methods
colSums((apply(data$y1,c(1,2),sum)+apply(data$y2,c(1,2),sum))>0)

#spatial recaptures by method
table(apply(1*(data$y1>0),c(1,2),sum))
table(apply(1*(data$y2>0),c(1,2),sum))

colSums(apply(data$y1>0,c(1,2),sum)>0) #collars deployed per year
colSums(data$tel.z.states==1,na.rm=TRUE) #collars active per year

#visualize realized activity centers
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",col=cols1)
points(X1.all,pch=4,lwd=2) #method 1
points(X2.all,pch=4,lwd=2,col="darkred") #method 2
points(data$truth$s,pch=16)

#visualize detections and telemetry by year
#note, if X1 and X2 traps are the same X2 will plot over X1
n <- nrow(data$y1)
for(g in 1:n.year){
  image(data$x.vals,data$y.vals,matrix(data$D.cov*data$InSS,data$n.cells.x,data$n.cells.y),
        main=paste("Year",g),xlab="X",ylab="Y",col=cols1)
  if(data$K1[g]>0){
    points(data$X1[[g]],pch=4,lwd=2)
  }
  if(data$K2[g]>0){
    points(data$X2[[g]],pch=4,lwd=2,col="darkred")
  }
  points(data$truth$s[data$truth$z[,g]==1,1],data$truth$s[data$truth$z[,g]==1,2],pch=16) #activity centers
  for(i in 1:n){
    trapcaps <- which(data$y1[i,g,]>0)
    traps <-  data$X1[[g]][1:data$J1[g],][trapcaps,]
    trapcaps2 <- which(data$y2[i,g,]>0)
    traps2 <-  data$X2[[g]][1:data$J2[g],][trapcaps2,]
    traps <- rbind(traps,traps2)
    if(nrow(traps)>0){
      s <- data$s[i,]
      points(s[1],s[2],col="goldenrod",pch=16)
      if(nrow(traps)>0){
        for(j in 1:nrow(traps)){
          lines(x=c(s[1],traps[j,1]),y=c(s[2],traps[j,2]),col="goldenrod")
        }
      }
    }
  }
  if(!any(is.null(data$tel.ID.g[[g]]))){
    for(i in 1:length(data$tel.ID.g[[g]])){
      kept.idx <- data$tel.ID.g[[g]][i]  # kept individual index
      tel.idx <- which(data$tel.ID==kept.idx)  # position in tel.ID vector
      if(length(tel.idx)>0){
        s <- data$s[kept.idx,]
        # find which telemetry year index corresponds to year g
        tel.g.idx <- which(data$tel.year[tel.idx,]==g)
        if(length(tel.g.idx)>0){
          for(l in 1:data$n.locs.ind[tel.idx,tel.g.idx]){
            lines(x=c(s[1],data$locs[tel.idx,tel.g.idx,l,1]),
                  y=c(s[2],data$locs[tel.idx,tel.g.idx,l,2]),col="gray80")
          }
          points(data$locs[tel.idx,tel.g.idx,,1],
                 data$locs[tel.idx,tel.g.idx,,2],
                 pch=16,cex=0.5,col="lightblue")
          points(s[1],s[2],col="darkblue",pch=16)
        }
      }
    }
  }
}

data$truth$N.super #N.super

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
           n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
           x.vals=data$x.vals,y.vals=data$y.vals)

##Initialize##
#Hard to predict appropriate M, depends on many factors like detection prob, number of years
#level of population turnover. Maybe make sure it is at least 1.6*N.super to start
M <- 300 #data augmentation level. Check N.super posterior to make sure it never hits M
n <- nrow(data$y1)
if(n > M) stop("Must augment more than number of individuals captured")

#initialize N and z objects and activity centers
#Use reasonable inits for p0, lam0, and sigma since we check to make sure initial observation
#model likelihood is finite. initializing with 1 parameter per session, just set all to same value
inits <- list(p01=rep(0.1,n.year),p02=rep(0.1,n.year),sigma=rep(0.5,n.year))
nimbuild <- init.Spatial.IPM(data,inits,M=M)

#pull these from data (won't be in environment if not simulated directly above)
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
cellArea <- res^2
D.cov <- data$D.cov
InSS <- data$InSS
x.vals <- data$x.vals
y.vals <- data$y.vals
n.cells <- data$n.cells
n.cells.x <- data$n.cells.x
n.cells.y <- data$n.cells.y

#can plot initialized ACs with telemetry for sanity check
#I should change this to plot by session and link ACs to traps
# plot(X.all,xlim=xlim,ylim=ylim,pch=4)
# for(i in 1:data$n.tel.inds){
#   for(g in 1:data$n.tel.years[i]){
#     for(l in 1:data$n.locs.ind[i,g]){
#       lines(x=c(nimbuild$s.init[data$tel.ID[i],1],data$locs[i,g,l,1]),
#             y=c(nimbuild$s.init[data$tel.ID[i],2],data$locs[i,g,l,2]),col="gray80")
#     }
#   }
#   points(data$locs[i,,,1],data$locs[i,,,2],pch=16,cex=0.5,col="lightblue")
#   points(nimbuild$s[data$tel.ID[i],1],nimbuild$s[data$tel.ID[i],2],pch=16,col="darkblue")
# }

#deal with missing years
sampled.years.m1 <- which(K1!=0)
sampled.years.m2 <- which(K2!=0)
n.sampled.m1 <- length(sampled.years.m1)
n.sampled.m2 <- length(sampled.years.m2)

#constants for Nimble
#might want to center D.cov here. Simulated D.cov in this testscript is already effectively centered.
constants <- list(n.year=n.year,M=M,J1=J1,J2=J2,K1D1=nimbuild$K1D1,K1D2=nimbuild$K1D2,
                  xlim=xlim,ylim=ylim,D.cov=D.cov,cellArea=cellArea,n.cells=n.cells,
                  res=res,n.tel.inds=data$n.tel.inds,n.tel.years=data$n.tel.years,
                  n.locs.ind=data$n.locs.ind,tel.ID=data$tel.ID,tel.year=data$tel.year,
                  sampled.years.m1=sampled.years.m1,sampled.years.m2=sampled.years.m2,
                  n.sampled.m1=n.sampled.m1,n.sampled.m2=n.sampled.m2)

#inits for Nimble
Niminits <- list(N=nimbuild$N,N.survive=nimbuild$N.survive,N.recruit=nimbuild$N.recruit,
                 ER=nimbuild$N.recruit,N.super=nimbuild$N.super,z.super=nimbuild$z.super,
                 z=nimbuild$z,z.start=nimbuild$z.start,z.stop=nimbuild$z.stop,
                 s=nimbuild$s,D0=nimbuild$N[1]/(sum(InSS)*res^2),D.beta1=0)

#data for Nimble
dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are
Nimdata <- list(y1=nimbuild$y1,y2=nimbuild$y2,X1=nimbuild$X1,X2=nimbuild$X2,
                dummy.data=nimbuild$dummy.data,cells=cells,InSS=InSS,
                locs=data$locs,tel.z.states=nimbuild$tel.z.states)

# set parameters to monitor
parameters <- c('N','gamma.fixed','N.recruit','N.survive','N.super','lambda.y1',
                'phi.fixed','p01','p02','sigma.fixed','D0','D.beta1')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel,constants=constants,data=Nimdata,check=FALSE,inits=Niminits)
#if you add/remove parameters in model file, do so in config.nodes
config.nodes <- c('gamma.fixed','phi.fixed','p01','p02','sigma.fixed')
conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt,
                      nodes=config.nodes,useConjugacy = FALSE)

#add N/z samplers
z.super.ups <- round(M*0.25) #how many z.super update proposals per iteration? 
#25% of M seems reasonable, but optimal will depend on data set
#loop here bc potentially different numbers of traps to vectorize in each year
y1.nodes <- pd1.nodes <- c()
y2.nodes <- pd2.nodes <- c()
for(g in 1:n.year){
  y1.nodes <- c(y1.nodes,Rmodel$expandNodeNames(paste0("y1[1:",M,",",g,",1:",J1[g],"]"))) #if you change y structure, change here
  pd1.nodes <- c(pd1.nodes,Rmodel$expandNodeNames(paste0("pd1[1:",M,",",g,",1:",J1[g],"]"))) #if you change y structure, change here
  y2.nodes <- c(y2.nodes,Rmodel$expandNodeNames(paste0("y2[1:",M,",",g,",1:",J2[g],"]"))) #if you change y structure, change here
  pd2.nodes <- c(pd2.nodes,Rmodel$expandNodeNames(paste0("pd2[1:",M,",",g,",1:",J2[g],"]"))) #if you change y structure, change here
}
N.nodes <- Rmodel$expandNodeNames(paste0("N"))
N.survive.nodes <- Rmodel$expandNodeNames(paste0("N.survive[1:",n.year-1,"]"))
N.recruit.nodes <- Rmodel$expandNodeNames(paste0("N.recruit[1:",n.year-1,"]"))
ER.nodes <- Rmodel$expandNodeNames(paste0("ER[1:",n.year-1,"]"))
z.nodes <- Rmodel$expandNodeNames(paste0("z[1:",M,",1]"))
tel.z.states.nodes <- Rmodel$expandNodeNames(paste0("tel.z.states[1:",M,",1]"))
calcNodes <- c(N.nodes,ER.nodes,N.recruit.nodes,N.survive.nodes,pd1.nodes,y1.nodes,
               pd2.nodes,y2.nodes,z.nodes,tel.z.states.nodes)
conf$addSampler(target = c("z"),
                type = 'zSampler',control = list(M=M,n.year=n.year,J1=J1,J2=J2,
                                                 sampled.years.m1=sampled.years.m1,sampled.years.m2=sampled.years.m2,
                                                 n.sampled.m1=n.sampled.m1,n.sampled.m2=n.sampled.m2,
                                                 z.obs=nimbuild$z.obs,z.super.ups=z.super.ups,
                                                 y1.nodes=y1.nodes,pd1.nodes=pd1.nodes,
                                                 y2.nodes=y2.nodes,pd2.nodes=pd2.nodes,
                                                 N.nodes=N.nodes,tel.z.states.nodes=tel.z.states.nodes,
                                                 z.nodes=z.nodes,ER.nodes=ER.nodes,
                                                 N.survive.nodes=N.survive.nodes,
                                                 N.recruit.nodes=N.recruit.nodes,
                                                 y2D=nimbuild$y2D,calcNodes=calcNodes), silent = TRUE)

#activity center sampler. This sampler tunes activity centers when z.super[i]=1 and
#draws from the prior otherwise.
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSamplerDcov',control=list(i=i,res=res,n.cells.x=n.cells.x,n.cells.y=n.cells.y,
                                                 xlim=xlim,ylim=ylim),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

conf$addSampler(target = c("D0","D.beta1"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time

mvSamples <-  as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:500),]))

#reminder what the targets are
data$N
data$N.recruit
data$N.survive
data$truth$N.super #N.super


cols.deployed <- colSums(apply(data$y1>0,c(1,2),sum)>0) #collars deployed per year
cols.active <- colSums(data$tel.z.states,na.rm=TRUE) #collars active per year
methods <- ifelse(K1 > 0 & K2 > 0, "C",
                  ifelse(K1 > 0, "A",
                         ifelse(K2 > 0, "B", NA)))

library(vioplot)
vioplot(mvSamples[-c(1:500),3:(n.year+2)],ylim=c(0,200),xlim=c(-0.5,n.year+0.5))
mtext("Method(s) Used",3,at=0,line=2)
mtext(methods,3,at=1:n.year,line=2)
mtext("cols.deployed",3,at=0,line=1)
mtext(cols.deployed,3,at=1:n.year,line=1)
mtext("cols.active",3,at=0,line=0)
mtext(cols.active,3,at=1:n.year,line=0)
