Multilevel modeling to determine covariate effects on community composition
========================================================

A common goal of community ecology is to understand how and why species composition shifts across space. Common techniques to determine which environmental covariates might lead to such shifts typically rely on ordination of community data to reduce the amount of data. These techniques include redundancy analysis (RDA), canonical correspondence analysis (CCA), and nonmetric multi-dimensional scaling (NMDS), each paired with permutation tests. However, as brought to light by Jackson et al. (2012: Ecography), these ordination techniques do not discern species-level covariate effects, making it difficult to attribute community-level pattern shifts to species-level changes. Jackson et al. (2012) propose a hierarchical modeling framework as an alternative, which we extend in this post to account for the case of imperfect detection. 

Multilevel models can estimate species-level random and fixed covariate effects to determine the relative contribution of environmental covariates to changing composition across space (Jackson et al. 2012). For presence/absence data, the model looks like:
$$y_{q} \sim Bern(\psi_{q})$$
$$\psi_{q} = logit^{-1}(\alpha_{spp[q]} + b_{spp[q]}  x_{site[q]})$$
$$\alpha_{spp[q]} \sim N(\mu_{\alpha}, \sigma_{intercept}^2)$$
$$b_{spp[q]} \sim N(\mu_{b}, \sigma_{slope}^2)$$

Here $y_q$ is the presences/absence of each species at each site ($q=1, ... , NM,$ where $N$ is the number of species and $M$ the number of sites). This model can be extended to incorporate multiple covariates. 

We are concerned here with whether species respond differently to environmental gradients (e.g. elevation, temperature, precipitation). If this is the case, then we expect community composition changes along such gradients. Concretely, we are interested in whether $\sigma_{slope}^2$ for any covariate differs from zero. 

Jackson et al. (2012) provide code for a maximum likelihood implementation of their model with data from Southern Appalachian understory herbs using the *R* package *lme4*. Here we present a simple extension of Jackson and colleague's work, correcting for detection error with repeat surveys (i.e. multi-species occupancy modeling). Specifically, the model could be changed slightly to:

$$y_{q} \sim Binom(z_q p_{spp[q]}, J_q)$$
$$z_q \sim Bern(\psi_q)$$
$$\psi_{q} = logit^{-1}(\alpha_{spp[q]} + b_{spp[q]}  x_{site[q]})$$
$$\alpha_{spp[q]} \sim N(\mu_{\alpha}, \sigma_{intercept}^2)$$
$$b_{spp[q]} \sim N(\mu_{b}, \sigma_{slope}^2)$$

Now $y_q$ is the number of times each species is observed at each site over $J$ surveys. $p_{spp[q]}$ represents the species-specific probability of detection when the species is present, and $z_q$ represents the 'true' occurence of the species, a Bernoulli random variable with probability, $\psi_q$. 

To demonstrate the method, we simulate data for a 20 species community across 50 sites with 4 repeat surveys. We assume that three site level environmental covariates were measured, one of which has variable affects on occurrence probabilities, one of which has consistent effects for all species, and one of which has no effect. We also assumed that species specific detection probabilities varied, but were independent of environmental covariates.

```{r}
################################################
# Simulate data
################################################

Nsite <- 50
Ncov <- 3
Nspecies <- 20
J <- 3

# species-specific intercepts:
alpha <- rnorm(Nspecies, 0, 1)

# covariate values
Xcov <- matrix(rnorm(Nsite*Ncov), 
               nrow=Nsite, ncol=Ncov)

# I'll assume 2 of the 3 covariates have significant effects
# one of which varies among species 
Beta <- array(c(rnorm(Nspecies, 0, 1), 
                rep(1, Nspecies), 
                rep(0, Nspecies)
                ), 
              dim=c(Nspecies, Ncov))

# species-specific detection probs
p0 <- runif(Nspecies, 0.35, 1)
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
```

We fit the following model with [JAGS](http://mcmc-jags.sourceforge.net/) with vague priors. 

```{r MLM_model.txt}
model {
  ### PRIORS ###
  for(i in 1:Nspecies){
    alpha[i] ~ dnorm(logit(psi.mean), psi.tau)T(-12,12)
    lp.detect[i] ~ dnorm(logit(p.detect.mean), p.detect.tau)T(-12,12)
    p.detect[i] <-  1/(1+exp(-lp.detect[i]))
  }
  for(j in 1:Ncov){
    mean.beta[j] ~ dnorm(meta.mean.beta, meta.tau.beta)
    sd.beta[j] ~ dunif(0, 10)
    tau.beta[j] <- pow(sd.beta[j]+0.001, -2)
  }
  
  # Heirarchical Priors #
  psi.mean ~ dbeta(1,1)
  p.detect.mean ~ dbeta(1,1)

  meta.mean.beta ~ dnorm(0,0.01)
  meta.sd.beta ~ dunif(0,10)
  meta.tau.beta <- pow(meta.sd.beta, -2)

  sd.psi ~ dunif(0,10)
  psi.tau <- pow(sd.psi, -2)
  
  sd.p.detect ~ dunif(0,10)
  p.detect.tau <- pow(sd.p.detect, -2)

  sd.p.include ~ dunif(0,10)
  p.include.tau <- pow(sd.p.include, -2)

  #--------------------------------------------------------------------------
  #--------------------------------------------------------------------------

  # Likelihood model:
  for(i in 1:Nobs){
    logit(psi[i]) <- alpha[Species[i]] + inprod(betas[Species[i],], X[i, ])
    z[i] ~ dbern(psi[i])
    Y[i] ~ dbinom(z[i] * p.detect[Species[i]], J[i])
  }
  
  # Deal with covariates:
  for(j in 1:Ncov){
    for(i in 1:Nspecies){
      betas[i,j] ~ dnorm(mean.beta[j], tau.beta[j])
    }
  }
}
```

Again, our inferences are primarily based on the $\sigma_{slope}^2$ terms. 

This re-run estimated covariate 5 as having a fixed effect, as it's 95% HDI for $\sigma_{slope}$ included zero. The model does a very good job at estimating the simulated slopes:

![betas][beta.estimates]

We can show the relative contribution of covariate random effects to composition using ordination, as discussed in Jackson et al. (2012). To do this, we compare the linear predictor (i.e. $logit^{-1}(\psi_q)$) of the best model that includes significant random effects to a model that does not have any random effects. That model looks like this:

```{r}
# model {
#   ### PRIORS ###
#   for(i in 1:Nspecies){
#     alpha[i] ~ dnorm(logit(psi.mean), psi.tau)T(-12,12)
#     lp.detect[i] ~ dnorm(logit(p.detect.mean), p.detect.tau)T(-12,12)
#     p.detect[i] <-  1/(1+exp(-lp.detect[i]))
#   }
#   for(j in 1:Ncov){
#     mean.beta[j] ~ dnorm(meta.mean.beta, meta.tau.beta)
#     
#     sd.beta[j] ~ dnorm(0, .001)T(0,100) # each covariate has a uncorrelated beta stdev.
#     tau.beta[j] <- pow(sd.beta[j]+0.001, -2)
#   }
#   
#   # Heirarchical Priors #
#   psi.mean ~ dbeta(1,1)
#   p.detect.mean ~ dbeta(1,1)
#   
#   meta.mean.beta ~ dnorm(0,0.01)
#   meta.sd.beta ~ dunif(0,50)
#   meta.tau.beta <- pow(meta.sd.beta, -2)
#   
#   sd.psi ~ dunif(0,50)
#   psi.tau <- pow(sd.psi, -2)
#   
#   sd.p.detect ~ dunif(0,50)
#   p.detect.tau <- pow(sd.p.detect, -2)
#   
#   #--------------------------------------------------------------------------
#   #--------------------------------------------------------------------------
#   
#   # Likelihood model:
#   for(i in 1:Nobs){
#     logit(psi[i]) <- alpha[Species[i]] + inprod(betas[], X[i, ])
#     z[i] ~ dbern(psi[i])
#     Y[i] ~ dbinom(z[i]*p.detect[Species[i]], J[i])
#   }
#   
#   # Deal with covariates:
#   # Assume no random effects:
#   for(j in 1:Ncov){
#     betas[j] ~ dnorm(mean.beta[j], tau.beta[j])
#   }  
#   
# }
```

To extract the linear predictors from each model, I justed calculated the average $logit^{-1}(\psi_q)$ from the respective posterior of each observed site-species combination. 

Then the code to conduct a PCA on the results and visualize the ordination space looks like this (mostly borrowed from Jackson et al. (2012):

```{r}
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
envir.vars <- Xcov[, c(1,3,4)]
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

text(1.3*-envir.points,label=c("Cov1", "Cov3", "Cov4"),cex=1, font=2)

```
![pca.plot][pca]

[beta.estimates]: BetaEstimates_best.pdf "Estimated slopes:"
[pca]: Ordination.pdf "PCA ordination:"