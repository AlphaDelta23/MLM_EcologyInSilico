Multilevel modeling to determine covariate effects on community composition
========================================================

A common goal of community ecology is to understand how and why species composition shifts across space. Common techniques to determine which environmental covariates might lead to such shifts typically rely on ordination of community data to reduce the amount of data. These techniques include redundancy analysis (RDA), canonical correspondence analysis (CCA), and nonmetric multi-dimensional scaling (NMDS), each paired with permutation tests. However, as brought to light by Jackson et al. (2012: Ecography), these ordination techniques do not discern species-level covariate effects. Thus, it becomes difficult to attribute community-level pattern shifts to species-level changes.  

Jackson et al. (2012) propose the use of multilevel models that estimate species-level random and fixed covariate effects to determine the relative contribution of environmental covariates to changing composition across space. For presence/absence data, the model looks like:
$$Y_{q} \sim Bern(\psi_{q})$$
$$\psi_{q} = logit^{-1}(\alpha_{spp[q]} + b_{spp[q]}  x_{site[q]})$$
$$\alpha_{spp[q]} \sim N(\mu_{\alpha}, \sigma_{intercept}^2)$$
$$b_{spp[q]} \sim N(\mu_{b}, \sigma_{slope}^2)$$

Here $Y_q$ is the presences/absence of each species at each site ($q=1, ... , NM,$ where $N$ is the number of species and $M$ is the number of sites). This model can also be extended to incorporate any number of covariates (e.g. matrix notation $\beta_{spp[q]} X_{site[q]}$). In practice, the data set would have species in long data format and site-specific covariates repeated (i.e. row-bound) for each species. 

In this model, determining covariate effects on species composition involves testing whether $\sigma_{slope}^2$ of any covariate differs significantly from zero. In other words, the model tests whether species respond differentially to particular covariates, which means composition should change along gradients of those particular covariates. 

Jackson et al. (2012) provide code to implement their model on a data set of Southern Appalachian understory herbs using the *R* package *lme4* and maximum likelihood techniques. However, a simple extension of Jackson and colleague's work would be to correct for detection error using repeat surveys (e.g. multi-species occupancy modeling). This should reduce bias in assigning covariate effects. Specifically, the model could be changed slightly to:

$$Y_{q} \sim Binom(z_q p_{spp[q]}, J_q)$$
$$z_q \sim Bern(\psi_q)$$
$$\psi_{q} = logit^{-1}(\alpha_{spp[q]} + b_{spp[q]}  x_{site[q]})$$
$$\alpha_{spp[q]} \sim N(\mu_{\alpha}, \sigma_{intercept}^2)$$
$$b_{spp[q]} \sim N(\mu_{b}, \sigma_{slope}^2)$$

Now $Y_q$ is the number of times each species is observed at each site over $J$ surveys. $p_{spp[q]}$ represents the species-specific probability of detection when the species is present, and $z_q$ represents the 'true' occurence of the species, moderated by occurrence probability, $\psi_q$. 

Here I generate some binomially distributed observed data for 15 species surveyed across 50 with 4 re-surveys. I impose species-specific detection probabilities (ranging from ~.35 to 1). I also simulate 15 measured environmental covariates, but I assume only 5 of these have effects on occupancy probability. Furthermore 4 out of 5 of these covariates are random variables, meaning the effect of the covariate is species-specific. The fifth covariate I assume to have a fixed effect among species. 

```{r}
################################################
# Simulate data
################################################

Nsite <- 50
Ncov <- 15
Nspecies <- 15
J <- 3

# species-specific intercepts:
alpha <- rnorm(Nspecies, 0, 1)

# covariate values
Xcov <- matrix(rnorm(Nsite*Ncov), 
               nrow=Nsite, ncol=Ncov)

# I'll assume 5 of the 15 covariates have significant effects
# One out of the 5 has only a fixed effect; Four out of the 5 have random effects 
# Remember, these are on the logit scale, so very high/low values will have large effects 
# and might swamp out smaller ones in the model. 
Beta <- array(0, dim=c(Nspecies, Ncov))
Beta[, 1] <- rnorm(Nspecies, 1.75, 0.75)
Beta[, 2] <- rnorm(Nspecies, -1, 0.75)
Beta[, 3] <- rnorm(Nspecies, -.5, 0.75)
Beta[, 4] <- rnorm(Nspecies, .5, 0.75)
Beta[, 5] <- rnorm(Nspecies, 1, 0.005) #Fixed factor
Beta[, 6:Ncov] <- rnorm(Nspecies*10, 0, 0.05)

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
# END Simulate data
################################################

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

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

Then I run a Bayesian model in JAGS that includes indicator variable selection to choose the most likely suite of environmental covariate effects. This first run of the model is to determine the covariates that should be included in the subsequent analysis. The model statement looks like this:

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
#     sd.beta[j] ~ dnorm(0, .001)T(0,100)
#     tau.beta[j] <- pow(sd.beta[j]+0.001, -2)
#     
#   }
#   
#   # Heirarchical Priors #
#   psi.mean ~ dbeta(1,1)
#   p.detect.mean ~ dbeta(1,1)
#   p.include ~ dbeta(1,1)
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
#   sd.p.include ~ dunif(0,50)
#   p.include.tau <- pow(sd.p.include, -2)
#   
#   #--------------------------------------------------------------------------
#   #--------------------------------------------------------------------------
#   
#   # Likelihood model:
#   for(i in 1:Nobs){
#     logit(psi[i]) <- alpha[Species[i]] + inprod(betas[Species[i],], X[i, ])
#     z[i] ~ dbern(psi[i])
#     Y[i] ~ dbinom(z[i]*p.detect[Species[i]], J[i])
#   }
#   
#   # Deal with covariates:
#   for(j in 1:Ncov){
#     I[j] ~ dbern(p.include)
#     
#     for(i in 1:Nspecies){
#       beta.full[i,j] ~ dnorm(mean.beta[j], tau.beta[j])
#       betas[i,j] <- I[j] * beta.full[i,j]
#     }
#   }
#   
#   for(j in 1:Ncov){
#     sd.beta.post[j] <- sd(betas[, j])
#   }
#   
# }
```

After running the model in JAGS (checking convergence, etc.), the best model had a posterior probability of 0.986 and included only the first 5 covariates, matching the simulated case. 

I then re-ran the model with only the first 5 covariates to get more accurate estimates of covariate effects, this time not including indicator variable selection. The model statment also includes a posterior estimation of each covariate's $\sigma_{slope}$ so that I can test if the covariate has truly a random or fixed effect. 

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