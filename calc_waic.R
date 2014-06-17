# Function calculates waic for occupancy model 
# input mcmc.list with posterior sims and data.frame for model input
# output list: log posterior predictive density, p_WAIC2, WAIC, mean likelihood
calc_waic <- function(posterior, data){
  with(data,{
    chains <- length(posterior)
    store <- dim(posterior[[1]])[1]
    L <- array(dim=c(Nobs, chains, store))
    L_bar <- rep(NA, Nobs)
    var_LL <- rep(NA, Nobs)
    
    for (i in 1:Nobs){
      for (j in 1:chains){
        post_sims <- posterior[[j]]
        indx_z <- which(dimnames(post_sims)[[2]] == 
                          paste("z[", i, "]", sep=""))
        zvals <- post_sims[, indx_z]
        indx_p <- which(dimnames(post_sims)[[2]] == 
                          paste("p.detect[", Species[i], "]", sep=""))
        pvals <- post_sims[, indx_p]
        L[i, j, ] <- dbinom(rep(Y[i], store), 
                            size=J[i], 
                            prob = zvals * pvals, 
                            log=TRUE)
      }
      L_bar[i] <- mean(exp(c(L[i, , ])))
      var_LL[i] <- var(c(L[i, , ]))
    }
    
    lppd <- sum(log(L_bar))
    p_WAIC <- sum(var_LL)
    WAIC <- -2 * (lppd - p_WAIC)
    return(list(lppd=lppd, p_WAIC=p_WAIC, WAIC=WAIC, L_bar))
  })
}
