# Testing MLM model with correction for detection probability:
# Data from: Jackson et al. 2012
source("calc_waic.R")
library(mcmcplots)

################################################
# Simulate data
################################################

Nsite <- 100
Ncov <- 3
Nspecies <- 20
J <- 4

# species-specific intercepts:
alpha <- rnorm(Nspecies, 0, 1)

# covariate values
Xcov <- matrix(rnorm(Nsite*Ncov, 0, 2), 
               nrow=Nsite, ncol=Ncov)

# I'll assume 2 of the 3 covariates have effects that vary among species
Beta <- array(c(rnorm(Nspecies, 0, 2), 
                rnorm(Nspecies, -1, 1),
                rep(1, Nspecies)
                ), 
              dim=c(Nspecies, Ncov)
              )

# species-specific detection probs
p0 <- plogis(rnorm(Nspecies, 1, 0.5))
p0

#### Occupancy states ####
Yobs <- array(0, dim = c(Nspecies, Nsite)) # Simulated observations

for(n in 1:Nspecies){
  for(k in 1:Nsite){
    lpsi <- alpha[n] + Beta[n, ] %*% Xcov[k, ] # Covariate effects on occurrence 
    psi <- 1/(1+exp(-lpsi)) #anti-logit
    
    z <- rbinom(1, 1, psi) # True Occupancy
    Yobs[n, k] <- rbinom(1, J, p0[n] * z) # Observed Occupancy
  }
}

################################################
# Format data for model
################################################
# X needs to have repeated covariates for each species, long form
X <- array(0, dim=c(Nsite*Nspecies, Ncov))
t <- 1; i <- 1
TT <- Nsite
while(i <= Nspecies){
  X[t:TT, ] <- Xcov
  t <- t+Nsite
  TT <- TT + Nsite
  i <- i+1
}

# Species
Species <- rep(c(1:Nspecies), each=Nsite)

# Observations/data:
Y <- NULL
for(i in 1:Nspecies){
  Y <- c(Y, Yobs[i, ])
}

# All sites surveyed same # times:
J <- rep(J, times=Nspecies*Nsite)

# Number of total observations
Nobs <- Nspecies*Nsite
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# JAGS model run:
################################################
library(rjags)

# data:
jags_d <- list(Y=Y,
               X=X,
               Species=Species,
               Nspecies=Nspecies,
               Ncov=Ncov,
               Nobs=Nobs,
               J=J)

# parameters:
params <- c("alpha", "betas", "p.detect", "sd.beta", "z", "psi")

jinits <- function() {
  list(
    z=ifelse(Y > 0, 1, 0)
  )
}

store<-1000
nadap<-1000
nburn<-2000
thin<-7
chains <- 3
mod0 <- jags.model(file = "model_statements/MLM_model_0f.txt", 
                  data = jags_d, n.chains = chains, n.adapt=nadap,
                  inits = jinits)
update(mod0, n.iter=nburn)
out0 <- coda.samples(mod0, n.iter = store*thin, 
                     variable.names = params, thin=thin)
caterplot(out0, "betas")
caterpoints(c(Beta))
WAIC_0f <- calc_waic(posterior=out0, jags_d)

## Fit other models
mod1 <- jags.model(file = "model_statements/MLM_model_1f.txt", 
                   data = jags_d, n.chains = chains, n.adapt=nadap,
                   inits = jinits)
update(mod1, n.iter=nburn)
out1 <- coda.samples(mod1, n.iter = store*thin, variable.names = params, thin=thin)
caterplot(out1, "betas")
caterpoints(c(Beta))
WAIC_f1 <- calc_waic(posterior=out1, jags_d)

mod2 <- jags.model(file = "model_statements/MLM_model_2f.txt", 
                   data = jags_d, n.chains = chains, n.adapt=nadap,
                   inits = jinits)
update(mod2, n.iter=nburn)
out2 <- coda.samples(mod2, n.iter = store*thin, variable.names = params, thin=thin)
caterplot(out2, "betas")
caterpoints(c(Beta))
WAIC_f2 <- calc_waic(posterior=out2, jags_d)

mod3 <- jags.model(file = "model_statements/MLM_model_3f.txt", 
                   data = jags_d, n.chains = chains, n.adapt=nadap,
                   inits = jinits)
update(mod3, n.iter=nburn)
out3 <- coda.samples(mod3, n.iter = store*thin, variable.names = params, thin=thin)
caterplot(out3, "betas")
caterpoints(c(Beta))
WAIC_f3 <- calc_waic(posterior=out3, jags_d)

mod12 <- jags.model(file = "model_statements/MLM_model_12f.txt", 
                   data = jags_d, n.chains = chains, n.adapt=nadap,
                   inits = jinits)
update(mod12, n.iter=nburn)
out12 <- coda.samples(mod12, n.iter = store*thin, variable.names = params, thin=thin)
caterplot(out12, "betas")
caterpoints(c(Beta))
WAIC_f12 <- calc_waic(posterior=out12, jags_d)

mod13 <- jags.model(file = "model_statements/MLM_model_13f.txt", 
                    data = jags_d, n.chains = chains, n.adapt=nadap,
                    inits = jinits)
update(mod13, n.iter=nburn)
out13 <- coda.samples(mod13, n.iter = store*thin, variable.names = params, thin=thin)
caterplot(out13, "betas")
caterpoints(c(Beta))
WAIC_f13 <- calc_waic(posterior=out13, jags_d)

mod23 <- jags.model(file = "model_statements/MLM_model_23f.txt", 
                    data = jags_d, n.chains = chains, n.adapt=nadap,
                    inits = jinits)
update(mod23, n.iter=nburn)
out23 <- coda.samples(mod23, n.iter = store*thin, variable.names = params, thin=thin)
caterplot(out23, "betas")
caterpoints(c(Beta))
WAIC_f23 <- calc_waic(posterior=out23, jags_d)

mod123 <- jags.model(file = "model_statements/MLM_model_123f.txt", 
                    data = jags_d, n.chains = chains, n.adapt=nadap,
                    inits = jinits)
update(mod123, n.iter=nburn)
out123 <- coda.samples(mod123, n.iter = store*thin, variable.names = params, thin=thin)
caterplot(out123, "betas")
caterpoints(c(Beta))
WAIC_f123 <- calc_waic(posterior=out123, jags_d)

# compare WAIC values
vals <- c(full = WAIC_0f$WAIC, 
  fixed1 = WAIC_f1$WAIC, 
  fixed2=WAIC_f2$WAIC, 
  fixed3=WAIC_f3$WAIC, 
  fixed12 = WAIC_f12$WAIC, 
  fixed23=WAIC_f23$WAIC, 
  fixed13=WAIC_f13$WAIC,   
  allfixed=WAIC_f123$WAIC
  )

plot(factor(names(vals)), vals,  
     xlab="model", ylab="WAIC")

# correct model and full model are best supported....






################################################
# Check Convergence:
################################################
library(ggmcmc)
library(mcmcplots)
alpha.df <- ggs(out, family="alpha")
beta.df <- ggs(out, family="betas")
I.df <- ggs(out, family="I")
sd.beta.df <- ggs(out, family="sd.beta")
p.detect.df <- ggs(out, family="p.detect")

ggs_Rhat(alpha.df) + xlim(.9, 1.2) + geom_vline(xintercept=1.1, linetype="dashed")
ggs_Rhat(beta.df) + xlim(.9, 1.2) + geom_vline(xintercept=1.1, linetype="dashed")
ggs_Rhat(sd.beta.df) + xlim(.9, 1.2) + geom_vline(xintercept=1.1, linetype="dashed")
ggs_Rhat(p.detect.df) + xlim(.9, 1.2) + geom_vline(xintercept=1.1, linetype="dashed")

caterplot(out, parms="betas")
caterpoints(c(Beta))

caterplot(out, parms="sd.beta")
caterpoints(c(1, 0))



################################################
# Check random vs. fixed:
################################################
# If 95 HDI of st.dev. of beta[j] overlaps zero, then fixed effect.
# (i.e. no significant variability in effect among species)

source(file="HDI.R")
hdi.sd <- array(0, dim=c(Ncov, 2))

for(i in 1:Ncov){
  sub <- subset(ggs(out), 
                Parameter==paste("sd.beta[",i,"]", sep="")
                )$value
  hdi.sd[i, ] <- HDI(sub) #HDI of st.dev. for each covariate
}
hdi.sd

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
  nadap<-20000
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

out_norand <- NULL
out_norand <- list(jags.parsamps[[1]][[1]],
                    jags.parsamps[[2]][[1]],
                    jags.parsamps[[3]][[1]])

class(out_norand) <- "mcmc.list"

stopCluster(cl)

################################################
# Check Convergence:
################################################
library(ggmcmc)
library(mcmcplots)
alpha.df.norand <- ggs(out_norand, family="alpha")
beta.df.norand <- ggs(out_norand, family="betas")
p.detect.df.norand <- ggs(out_norand, family="p.detect")
psi.df.norand <- ggs(out_norand, family="psi")

ggs_Rhat(alpha.df.norand)
ggs_Rhat(beta.df.norand)
ggs_Rhat(p.detect.df.norand)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(out_norand, parms="betas", horizontal=F)
caterplot(out_norand, parms="alpha", horizontal=F)

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

rownames(MLM.fitted)=c(1:Nsite)
colnames(MLM.fitted)=paste0("Species", 1:Nspecies)
# standardize over spp

MLM.fitted.standard <- MLM.fitted
for(j in 1:Nspecies){
  MLM.fitted.standard[,j] <- (MLM.fitted[,j]-mean(MLM.fitted[,j]))/sd(MLM.fitted[,j])
} 

ss <- cor(MLM.fitted.standard)
U <- svd(ss)
mlm.fit <- MLM.fitted.standard %*% U$v
mlm.fit <- mlm.fit[,1:2]

# environmental variables (only those with significant random effects)
envir.vars <- Xcov[, c(1:4)]
mlm.envir <- NULL
for(j in 1:ncol(envir.vars)){
  mlm.envir <- cbind(mlm.envir, envir.vars[,j]*mlm.fit[,1],envir.vars[,j]*mlm.fit[,2])
}

envir.points <- t(array(colMeans(mlm.envir),c(2,dim(mlm.envir)[2]/2)))

# plot mlm

x11(height=5, width=5)
plot(-mlm.fit,xlab="PC1",ylab="PC2",type="n")
#text(-mlm.fit,label=c(1:(Nobs/Nspecies)),cex=.5)
points(-mlm.fit, pch=19, cex=0.5)

arrow.coordMLM <- cbind(array(0,dim(envir.points)),-envir.points)

arrows(arrow.coordMLM[,1],arrow.coordMLM[,2],arrow.coordMLM[,3],arrow.coordMLM[,4], 
       code=2, col="black", length=0.05, lwd=.8)

text(1.3*-envir.points,label=c("Cov1", "Cov2", "Cov3", "Cov4"),cex=1, font=2)

