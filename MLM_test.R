# Testing MLM model with correction for detection probability:
# Data from: Jackson et al. 2012

################################################
# load data
################################################
# set wd
setwd("~/Documents/Thesis Research/CA metacoms/CAmetacomsCode/MLM_Var.Select")
# use file "MLM_dataset.csv"
h <- read.csv("MLM_dataset.csv")
# remove sites with no species
h <- h[-3,]

#Normalize predictor variables:
for(i in 5:25){h[,i] <- (h[,i] - mean(h[,i]))/sd(h[,i])} 

hh <- cbind(h[,1:26],factor('ACRA')) 
names(hh)[26] <- 'PRESENCE' 

names(hh)[27] <- 'SPP' 
levels(hh$SPP) <- c("ACRA","ARTR","CATH","CHMA","DILA","GALA","GEMA","LISU","POBI","SACA","THDI","TRGR","UVGR","VIRO")

T <- dim(h)[1] 
T #[1] 54
t <- T 
t  #[1] 54
names(hh)

for(i in 27:39){
  hh[(t+1):(t+T),] <- cbind(h[,1:25], h[,i], names(h)[i])
  t <- t+T
}
t #756

hh$Elevation2 <- hh$Elevation^2
hh$LITU2 <- hh$LITU^2
hh$Ca2 <- hh$Ca^2
hh$P2 <- hh$P^2
hh$Herb2=hh$Herb^2

################################################
# end load data
################################################

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Format data for model
################################################
# All data for Bayesian model:
X <- hh[, c(5:25, 28:32)] #include all covariates
Y <- hh$PRESENCE
Species <- as.factor(as.numeric(hh$SPP))
Nspecies <- length(levels(Species))
Ncov <- ncol(X)
Nobs <- length(Y)
J <- rep(5, times=Nobs)
  
# Convert Y (pres/abs) to binomial based on simulated species detection probs.
ps <- runif(Nspecies, 0.25, 1)
pdetect <- rep(ps, each=54) # each species has 54 sites
for(i in 1:Nobs){
  if(Y[i] ==1){
    Y[i] <- rbinom(1, 5, pdetect[i])
  }
}

# Some of the X variables are highly correlated:
# Keep CA, Elevation
# Remove: OM, N, Mg, pH, K, Litu2, Ca2

X <- X[, -c(15,20,17,14,19,23)]
Ncov <- ncol(X)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# JAGS model run:
################################################
library(rjags)
library(random)
load.module("lecuyer")
load.module("glm")
# data:
jags_d <- list(Y=Y,
               X=X,
               Species=Species,
               Nspecies=Nspecies,
               Ncov=Ncov,
               Nobs=Nobs,
               J=J)

# zinit <- NULL
# zinit <- ifelse(Y > 0, 1, 0)


# set wd:
setwd("~/Documents/Thesis Research/CA metacoms/CAmetacomsCode/MLM_Var.Select")

# parameters:
params <- c("alpha", "betas", "I", "tau.beta", "p.detect", "p.include", 
            "sd.beta.post", "psi")

jinits <- function() {
  list(
    z=ifelse(Y > 0, 1, 0),
    .RNG.name=c("base::Super-Duper"),
    .RNG.seed=as.numeric(randomNumbers(n = 1, min = 1, max = 1e+06, col=1))
  )
}

# initialize model:

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

jags.parsamps <- NULL
jags.parsamps <- foreach(i=1:3, .packages=c('rjags','random')) %dopar% {
  #setwd("C:\Users\Joe\Documents\GitHub\CA_Metacoms")
  store<-1000
  nadap<-50000
  nburn<-50000
  thin<-50
  mod <- jags.model(file = "MLM_model.txt", 
                    data = jags_d, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, variable.names = params, thin=thin)
  return(out)
}

bundle <- NULL
bundle <- list(jags.parsamps[[1]][[1]],
               jags.parsamps[[2]][[1]],
               jags.parsamps[[3]][[1]])

class(bundle) <- "mcmc.list"

stopCluster(cl)

# # initialize model:
# mod <- NULL
# mod <- jags.model(file = "MLM_model.txt", 
#                   data = jags_d, n.chains = 3, n.adapt=nadap,
#                   inits = list(z=zinit))
# update(mod, n.iter=nburn)
# 
# out <- NULL
# out <- coda.samples(mod, n.iter = store*thin, variable.names = params, thin=thin)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Check Convergence:
################################################
library(ggmcmc)
library(mcmcplots)
alpha.df <- ggs(bundle, family="alpha")
beta.df <- ggs(bundle, family="betas")
I.df <- ggs(bundle, family="I")
sd.beta.df <- ggs(bundle, family="sd.beta.post")
p.detect.df <- ggs(bundle, family="p.detect")
p.include.df <- ggs(bundle, family="p.include")

ggs_Rhat(alpha.df)
ggs_Rhat(beta.df)
ggs_Rhat(tau.beta.df)
ggs_Rhat(p.detect.df)
ggs_Rhat(p.include.df)

quartz(height=4, width=11)
#x11(height=4, width=11)
caterplot(bundle, parms="betas", horizontal=F, random=50)

caterplot(bundle, parms="tau.beta", horizontal=F, val.lim=c(-1, 10))
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Check best model:
################################################

Ipost <- NULL
Ipost <- array(dim=c(store*3, Ncov)) # indicator variable array
for (i in 1:Ncov){
  string <- paste("I[", i, "]", sep="")
  Ipost[, i] <- subset(I.df, Parameter==string)$value
}

# what are the unique models that have nonzero posterior probability?
uniquemods <- unique(Ipost, MARGIN=1)
# how many do we have?
nmods <- dim(uniquemods)[1]
nmods
model.probabilities <- rep(NA, nmods)
for (i in 1:nmods){
  TFs <- apply(Ipost, 1, function(x) all(x == uniquemods[i,]))
  model.probabilities[i] <- sum(TFs) / (store*3)
}

sum(model.probabilities)
ordered.mod.probs <- model.probabilities[order(-model.probabilities)]
ordered.mods <- uniquemods[order(-model.probabilities), ]

ordered.mod.probs[1:10]
ordered.mods[1:5, ]

# Best model probability: 0.58 (followed by 0.087)
# Includes variables: c(1,10,11,14,17)
# Elevation, ACSA, LITU, Ca, Elevation2

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now refit best model:
################################################
# (I deleted most of the data, just to save on space)

# X_best for best covariates:

X_best <- X[, c(1,10,11,14,17)]
Ncov <- ncol(X_best)

# data:
jags_d_best <- list(Y=Y,
                    X=X_best,
                    Species=Species,
                    Nspecies=Nspecies,
                    Ncov=Ncov,
                    Nobs=Nobs,
                    J=J)

# parameters:
params_best <- c("alpha", "betas", "p.detect", 
                "sd.beta.post", "psi")
jinits <- function() {
  list(
    z=ifelse(Y > 0, 1, 0),
    .RNG.name=c("base::Super-Duper"),
    .RNG.seed=as.numeric(randomNumbers(n = 1, min = 1, max = 1e+06, col=1))
  )
}
# initialize model:

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

jags.parsamps <- NULL
jags.parsamps <- foreach(i=1:3, .packages=c('rjags','random')) %dopar% {
  #setwd("C:\Users\Joe\Documents\GitHub\CA_Metacoms")
  store<-1000
  nadap<-50000
  nburn<-50000
  thin<-50
  mod <- jags.model(file = "MLM_model_Best.txt", 
                    data = jags_d_best, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, 
                      variable.names = params_best, thin=thin)
  return(out)
}

bundle_best <- NULL
bundle_best <- list(jags.parsamps[[1]][[1]],
               jags.parsamps[[2]][[1]],
               jags.parsamps[[3]][[1]])

class(bundle_best) <- "mcmc.list"

stopCluster(cl)

################################################
# Check Convergence:
################################################
library(ggmcmc)
library(mcmcplots)
alpha.df.best <- ggs(bundle_best, family="alpha")
beta.df.best <- ggs(bundle_best, family="betas")
sd.beta.df.best <- ggs(bundle_best, family="sd.beta.post")
p.detect.df.best <- ggs(bundle_best, family="p.detect")
psi.df.best <- ggs(bundle_best, family="psi")

ggs_Rhat(alpha.df.best)
ggs_Rhat(beta.df.best)
ggs_Rhat(sd.beta.df.best)
ggs_Rhat(p.detect.df.best)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(bundle_best, parms="betas", horizontal=F, random=50)
caterplot(bundle_best, parms="sd.beta.post", horizontal=F)

################################################
# Check random vs. fixed:
################################################
# If 95 HDI of st.dev. of beta[j] overlaps zero, then fixed effect.
# (i.e. no significant variability in effect among species)

source(file="HDI.R")
hdi.sd <- array(0, dim=c(Ncov, 2))

for(i in 1:Ncov){
  sub <- subset(sd.beta.df.best, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd[i, ] <- hdi
}
hdi.sd
# Fixed: ACSA
# Random: Elevation, LITU, CA, 


################################################
# Extract 'linear predictor' of the model: logit(psi)
################################################

linpred.best <- NULL

for(i in 1:Nobs){
  sub <- subset(psi.df.best, Parameter==paste("psi[", i, "]", sep=""))$value
  sub <- log(sub/(1-sub))
  linpred.best[i] <- mean(sub)
}

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now fit the model that has NO RANDOM TERMS (i.e. no species-level variability in beta[j]:
################################################
# (I deleted most of the data, just to save on space)

# initialize model:

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

jags.parsamps <- NULL
jags.parsamps <- foreach(i=1:3, .packages=c('rjags','random')) %dopar% {
  #setwd("C:\Users\Joe\Documents\GitHub\CA_Metacoms")
  store<-1000
  nadap<-50000
  nburn<-50000
  thin<-50
  mod <- jags.model(file = "MLM_model_NoRandom.txt", 
                    data = jags_d_best, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, 
                      variable.names = params_best, thin=thin)
  return(out)
}

bundle_norand <- NULL
bundle_norand <- list(jags.parsamps[[1]][[1]],
                    jags.parsamps[[2]][[1]],
                    jags.parsamps[[3]][[1]])

class(bundle_norand) <- "mcmc.list"

stopCluster(cl)

################################################
# Check Convergence:
################################################
library(ggmcmc)
library(mcmcplots)
alpha.df.norand <- ggs(bundle_norand, family="alpha")
beta.df.norand <- ggs(bundle_norand, family="betas")
p.detect.df.norand <- ggs(bundle_norand, family="p.detect")
psi.df.norand <- ggs(bundle_norand, family="psi")

ggs_Rhat(alpha.df.norand)
ggs_Rhat(beta.df.norand)
ggs_Rhat(p.detect.df.norand)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(bundle_norand, parms="betas", horizontal=F)
caterplot(bundle_norand, parms="alpha", horizontal=F)

################################################
# Extract 'linear predictor' of the model: logit(psi)
################################################

linpred.norand <- NULL

for(i in 1:Nobs){
  sub <- subset(psi.df.norand, Parameter==paste("psi[", i, "]", sep=""))$value
  sub <- log(sub/(1-sub))
  linpred.norand[i] <- mean(sub)
}

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Conduct PCA to ordinate sites, map environmental effects
################################################
MLM.fitted <- array(linpred.best - linpred.norand, c(Nobs/Nspecies, Nspecies))

rownames(MLM.fitted)=c(1:54)
colnames(MLM.fitted)=c("ACRA","ARTR","CATH","CHMA","DILA","GALA","GEMA","LISU","POBI","SACA","THDI","TRGR","UVGR","VIRO")

# standardize over spp
MLM.fitted.standard <- MLM.fitted
for(j in 1:Nspecies){
  MLM.fitted.standard[,j] <- (MLM.fitted[,j]-mean(MLM.fitted[,j]))/sd(MLM.fitted[,j])
} 

ss <- cor(MLM.fitted.standard)
U <- svd(ss)
mlm.fit <- MLM.fitted.standard %*% U$v
mlm.fit <- mlm.fit[,1:2]

# environmental variables (only those with random effects)
envir.vars <- cbind(hh$Elevation, hh$LITU, hh$Ca)
mlm.envir <- NULL
for(j in 1:ncol(envir.vars)){
  mlm.envir <- cbind(mlm.envir, envir.vars[,j]*mlm.fit[,1],envir.vars[,j]*mlm.fit[,2])
}

envir.points <- t(array(colMeans(mlm.envir),c(2,dim(mlm.envir)[2]/2)))

# plot mlm
par(mfcol=c(1,1))
x11(height=5, width=5)
plot(-mlm.fit,xlab="PC1",ylab="PC2",type="n")
text(-mlm.fit,label=c(1:(Nobs/Nspecies)),cex=.5)

arrow.coordMLM <- cbind(array(0,dim(envir.points)),-envir.points)

arrows(arrow.coordMLM[,1],arrow.coordMLM[,2],arrow.coordMLM[,3],arrow.coordMLM[,4], col="black", length=0.1)

text(1.3*-envir.points,label=c("Elevation", "LITU", "Ca"),cex=.7)

