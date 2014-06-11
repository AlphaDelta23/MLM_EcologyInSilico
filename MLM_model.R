model {
  
  # Likelihood model:
  for(i in 1:Nobs){
    logit(psi[i]) <- alpha[species[i]] + inprod(betas[species[i],], X[i, ])
    z[i] ~ dbern(psi[i])
    Y[i] ~ dbinom(z[i]*p.detect[species[i]], J[i])
  }
  
  # Deal with covariates
  for(i in 1:Nspecies){
    for(j in 1:Ncov){
      I[j] ~ dbern(p.include[i])
      beta.full[i,j] ~ dnorm(mean.beta[j], tau.beta[j])) 
      betas[i,j] <- I[j] * beta.full[i,j]
    }
  }

  ### PRIORS ###

  for(i in 1:Nspecies){
    alpha[i] ~ dnorm(logit(psi.mean), psi.tau)T(-12,12)
    logit(p.detect[i]) ~ dnorm(logit(p.detect.mean), p.detect.tau)T(-12,12)
    logit(p.include[i]) ~ dnorm(logit(p.include.mean), p.include.tau)T(-12,12)
  }
  for(j in 1:Ncov){
    mean.beta[j] ~ dnorm(meta.mean.beta, meta.tau.beta)
    sd.beta[j] ~ dunif(0, sd.beta.high)
    tau.beta[j] <- pow(sd.beta[j], -2)
  }
  
  # Heirarchical Priors #
  psi.mean ~ dbeta(1,1)
  p.detect.mean ~ dbeta(1,1)
  p.include.mean ~ dbeta(1,1)

  meta.mean.beta ~ dnorm(0,0.01)
  meta.sd.beta ~ dunif(0,100)
  meta.tau.beta <- pow(meta.sd.beta, -2)
  sd.beta.high ~ dunif(0,100)

  sd.psi ~ dunif(0,100)
  psi.tau <- pow(sd.psi, -2)
  
  sd.p.detect ~ dunif(0,100)
  p.detect.tau <- pow(sd.p.detect, -2)
  
  sd.p.include ~ dunif(0,100)
  p.include.tau <- pow(sd.p.include, -2)
   
}