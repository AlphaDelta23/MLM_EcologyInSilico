# Kuo and Mallick variable selection with JAGS
# MJ 10/21/13

# Example using simulated data with a Poisson r.v.
# - select important covariates
# - calculate posterior model probabilities
# - generate Bayesian model averaged predictions
library(rjags)
library(mcmcplots)
library(ggmcmc)

#---------------#
# simulate data #
#---------------#
nobs <- 150
nvars <- 20

# probability of covariate inclusion
pI <- .25 
I <- rbinom(nvars, 1, pI)

# covariate values
Xcov <- matrix(rnorm(nobs*nvars), 
               nrow=nobs, ncol=nvars)

intercept <- 1
X <- cbind(rep(1, nobs), Xcov)
slopes <- I*rnorm(nvars)
Beta <- c(intercept, slopes) 

# simulate expectations & data
lambda <- exp(X %*% Beta)
Y <- rpois(nobs, lambda)

# Also, add some missing data (for model averaging)
Y[1:3] <- NA

#--------------#
# futz w/ JAGS #
#--------------#
jd <- list(X=Xcov, Y=Y, nobs=nobs, nvars=nvars)

setwd("~/Desktop")
# specify model
cat("
    model{
# intercept
alpha ~ dnorm(0, .001)

# priors for slopes (commented out lines = optional heirarchical prior)
# sd.beta ~ dunif(0, 10)
# tau.beta <- 1 / (sd.beta * sd.beta)

for (i in 1:nvars){
#  I[i] ~ dbern(p[i])
# beta.full[i] ~ dnorm(0, tau.beta)
  I[i] ~ dbern(p)
  beta.full[i] ~ dnorm(0, .001) 
  beta[i] <- I[i] * beta.full[i]
}

p ~ dbeta(1,1)
#for(i in 1:nvars){p[i] ~ dbeta(1, 1)} # prior Pr(covariate inclusion)

# likelihood
for (i in 1:nobs){
  log(lambda[i]) <- alpha + inprod(beta[], X[i, ])
  Y[i] ~ dpois(lambda[i])
}
    }
    ", 
    file="KM.txt")

chains <- 3
iters <- 5000
nthin <- 1
store <- chains * iters / nthin

mod <- jags.model("KM.txt", data=jd, n.adapt=2000, n.chains=chains)
out <- coda.samples(mod, 
                    variable.names=c("I", "beta", "alpha", "p",
                                     "Y[1]", "Y[2]", "Y[3]"), 
                    n.iter=iters, thin=nthin)

#---------------#
# Assess output #
#---------------#
summary(out)
plot(out)

# compare estimates to known values
par(mfrow=c(1,1))
caterplot(out, "beta")
caterpoints(x=slopes)

# calculate posterior model probabilities
ggd <- ggs(out)
Ipost <- array(dim=c(store, nvars)) # indicator variable array
for (i in 1:nvars){
  string <- paste("I[", i, "]", sep="")
  Ipost[, i] <- subset(ggd, Parameter==string)$value
}

# what are the unique models that have nonzero posterior probability?
uniquemods <- unique(Ipost, MARGIN=1)
# how many do we have?
nmods <- dim(uniquemods)[1]
nmods
model.probabilities <- rep(NA, nmods)
for (i in 1:nmods){
  TFs <- apply(Ipost, 1, function(x) all(x == uniquemods[i,]))
  model.probabilities[i] <- sum(TFs) / store
}

model.probabilities

# what is most probable?
best <- which.max(model.probabilities)
mpmodel <- uniquemods[best,]
mpmodel

# How does the most probable model compare to the true model?
rbind(mpmodel, I)

# Bayesian model averaging to make predictions (already done)
ggs_density(ggd, "Y") + ggtitle("Bayesian model averaged predictions")
