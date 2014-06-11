Multilevel modeling to determine covariate effects on community composition
========================================================

A common goal of community ecology is to understand how and why species composition shifts across space. Common techniques to discover which environmental covariates might lead to such shifts include redundancy analysis (RDA), canonical correspondence analysis (CCA), and nonmetric multi-dimensional scaling (NMDS), each paired with permutation tests. However, as brought to light by Jackson et al. (2012: Ecography), these ordination techniques do not discern species-level covariate effects. Thus, it becomes difficult to attribute community-level patterns to species-level effects.  

Jackson et al. (2012) propose the use of multilevel models that estimate species-level random and fixed covariate effects to determine the relative contribution of covariates to changing composition. For incidence data, the model looks like:
$$Y_{q} \sim Bern(\mu_{q})$$
$$\mu_{q} = logit^{-1}(\alpha_{spp[q]} + b_{spp[q]}  x_{site[q]})$$
$$\alpha_{spp[q]} \sim N(\mu_{\alpha}, \sigma_{intercept}^2)$$
$$b_{spp[q]} \sim N(\mu_{b}, \sigma_{slope}^2)$$

Here $Y_q$ is the presences/absence of each species at each site ($q=1, ... , NM,$ where $N$ is the number of species and $M$ is the number of sites). This model can also be extended to incorporate any number of covariates (e.g. matrix notation $\beta_{spp[q]} X_{site[q]}$). In practice, the data set would have species in long data format and site-specific covariates repeated (i.e. row-bound) for each species. 

In this model, determining covariate effects on species composition involves testing whether $\sigma_{slope}^2$ of any covariate differs significantly from zero. In other words, the model tests whether species respond differentially to particular covariates, which means composition should change along gradients of those particular covariates. 

Jackson et al. (2012) provide code to implement their model on a data set of Southern Appalachian understory herbs using the *R* package *lme4* and maximum likelihood techniques. However, a simple extension of Jackson and colleague's work would be to correct for detection error using repeat surveys (e.g. multi-species occupancy modeling). This should reduce bias in assigning covariate effects. Specifically, the model could be changed slightly to:

$$Y_{q} \sim Binom(z_q p_{spp[q]}, J_q)$$
$$z_q \sim Bern(\mu_q)$$
$$\mu_{q} = logit^{-1}(\alpha_{spp[q]} + b_{spp[q]}  x_{site[q]})$$
$$\alpha_{spp[q]} \sim N(\mu_{\alpha}, \sigma_{intercept}^2)$$
$$b_{spp[q]} \sim N(\mu_{b}, \sigma_{slope}^2)$$

Now $Y_q$ is the number of times each species is observed at each site over $J$ surveys. $p_{spp[q]}$ represents the species-specific probability that the species will be observed in a survey when present, and $z_q$ represents the 'true' occurence of the species, moderated by occurrence probability, $\mu_q$. 



