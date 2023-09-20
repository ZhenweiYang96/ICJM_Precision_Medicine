rm(list=ls())
library(rjags)
library(mcmcplots)
library(GLMMadaptive)
library(ggplot2)
library(tidyverse)
library(splines)
library(future)
plan(multisession, workers = 3)
load("Cleaned data/pass_cores.RData")
load("Cleaned data/pass.RData")
load("Cleaned data/pass_id.RData")
load("R script/Data analysis/ICJM1_inits/mm_inits.RData")

model.txt <- "model{
  for (i in 1:N) {
    for (j in num_psa[i]:(num_psa[i+1]-1)) {
      Y[j] ~ dt(mu[j], tau, 3)
      mu[j] <- inprod(XL_psa[j,], betaL_psa[]) + inprod(ZL_psa[j,], b[i,1:4])

    }
    
    for (j in num_cr[i]:(num_cr[i+1]-1)) {
      pos[j] ~ dbin(p[j], total[j])
      logit(p[j]) <- inprod(XL_cr[j,], betaL_cr[]) + inprod(ZL_cr[j,], b[i,5:7])
    }
    
    b[i,1:Nb] ~ dmnorm(mub[], inv_D[,])
  }
  
  inv_D[1:Nb,1:Nb] ~ dwish(4 * diagtb[,], Nb+1)
  D[1:Nb,1:Nb] <- inverse(inv_D[,])
  
  for (l in 1:Nb) {
    for (k in 1:Nb) {diagtb[l,k] <- ifelse(l == k, tb[l], 0)}
    tb[l] ~ dgamma(0.5, 0.01)
  }
 
  for (l in 1:Nbeta_psa) {betaL_psa[l] ~ dnorm(0, 0.01)}
  for (l in 1:Nbeta_cr) {betaL_cr[l] ~ dnorm(0, 0.01)}
  
  tau ~ dgamma(0.01,0.01)
}"


# Data preparation --------------------------------------------------------

# PSA part
n <- nrow(pass.id)
num_psa <- as.vector(c(1, 1+ cumsum(tapply(pass$CISNET_ID, pass$CISNET_ID, length))))
Y <- pass$PSAValue
lgt <- pass$TimeSince_Dx
age <- pass$DxAge
evt <- pass.id[,c("time.cmp1", "time.cmp2")]
XL_psa <- array(1, dim = c(nrow(pass), 5))
knot.longi <- c(min(c(as.matrix(evt)), lgt), 
                as.vector(quantile(lgt, probs = seq(0,1,length.out=4)))[2:3],
                max(c(as.matrix(evt)), lgt))

XL_psa[,2:4] <- ns(lgt, knots = knot.longi[2:3], B = knot.longi[c(1,4)])
XL_psa[,5] <- pass$DxAge
ZL_psa <- XL_psa[,1:4]

# Cores ratio part
num_cr <- as.vector(c(1, 1+ cumsum(tapply(pass.cores$CISNET_ID, pass.cores$CISNET_ID, length))))
ZL_cr <- XL_cr <- model.matrix(~ poly(TimeSince_Dx, 2, raw = T), data = pass.cores)
N <- length(unique(pass.cores$CISNET_ID))

data.list <- list(N = N, 
                  num_psa = num_psa, XL_psa = XL_psa, ZL_psa = ZL_psa,
                  Nbeta_psa = ncol(XL_psa),
                  mub = rep(0,7), Nb = ncol(ZL_psa) + ncol(ZL_cr),
                  Y = Y, 
                  num_cr = num_cr, total = pass.cores$total_cores,
                  pos = pass.cores$pos_cores, XL_cr = XL_cr, ZL_cr = ZL_cr,
                  Nbeta_cr = ncol(XL_cr))

# create inits
bin.mixed.glmm <- mixed_model(fixed = cbind(pos_cores, total_cores - pos_cores) ~ poly(TimeSince_Dx, 2, raw = T),
                              random = ~ poly(TimeSince_Dx, 2, raw = T) | CISNET_ID, data = pass.cores,
                              family = binomial(), 
                              n_phis = 1)



set.seed(2022)
rng.name <- sample(c("base::Wichmann-Hill",
                     "base::Marsaglia-Multicarry",
                     "base::Super-Duper",
                     "base::Mersenne-Twister"), 3, replace = T)
rng.num <- sample(1:100000, 3, replace = F)

inv_D_mat <- cbind(rbind(cbind(solve(mm$D),0,0,0),0,0,0))
inv_D_mat[5:7,5:7] <- solve(bin.mixed.glmm$D)

initials <- list(
  list(betaL_psa = fixef(mm),
       betaL_cr = fixef(bin.mixed.glmm),
       b = cbind(ranef(mm), ranef(bin.mixed.glmm)),
       inv_D = inv_D_mat,
       tau = 1/exp(mm$phis)^2,
       .RNG.name = rng.name[1],
       .RNG.seed = rng.num[1]),
  list(betaL_psa = fixef(mm),
       betaL_cr = fixef(bin.mixed.glmm),
       b = cbind(ranef(mm), ranef(bin.mixed.glmm)),
       inv_D = inv_D_mat,
       tau = 1/exp(mm$phis)^2,
       .RNG.name = rng.name[2],
       .RNG.seed = rng.num[2]),
  list(betaL_psa = fixef(mm),
       betaL_cr = fixef(bin.mixed.glmm),
       b = cbind(ranef(mm), ranef(bin.mixed.glmm)),
       inv_D = inv_D_mat,
       tau = 1/exp(mm$phis)^2,
       .RNG.name = rng.name[3],
       .RNG.seed = rng.num[3])
)


out <- lapply(1:3, function(i) {
  future::future({
    bin.mixed <- jags.model(file = textConnection(model.txt), data = data.list,
                            inits = initials[[i]], n.chains = 1, n.adapt = 3000,
                            quiet = T)  #
    update(bin.mixed, n.iter = 3000)
    mcmc.mvbin <- coda.samples(bin.mixed, variable.names = c("betaL_cr", "b", "betaL_psa", 
                                                             "tau", "D"),
                               n.iter = 10000, thin = 10)
    return(list(model = bin.mixed,
                mcmc = mcmc.mvbin))
  })
})
samples.binmixed <- lapply(out, future::value)
mcmc.mvmm <-  as.mcmc.list(lapply(1:3, function(x)
  as.mcmc(samples.binmixed[[x]]$mcmc)))
mcmcplot(mcmc.mvmm)
summary(mcmc.mvmm)

parm.mvmm.pool <- colMeans(do.call(rbind, mcmc.mvmm))
parm.mvmm <- list(
  betaL_psa = parm.mvmm.pool[grep("betaL_psa", names(parm.mvmm.pool))],
  betaL_cr = parm.mvmm.pool[grep("betaL_cr", names(parm.mvmm.pool))],
  D = matrix(parm.mvmm.pool[grep("D", names(parm.mvmm.pool))],7,7),
  b = matrix(parm.mvmm.pool[50:5880], 833, 7),
  tau = parm.mvmm.pool[grep("tau", names(parm.mvmm.pool))]
)

save(mcmc.mvmm, file = "R script/Data analysis/ICJM2_inits/MCMC_mvmm.RData")
save(parm.mvmm, file = "R script/Data analysis/ICJM2_inits/mvmm_inits.RData")

