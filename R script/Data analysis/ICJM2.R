library(splines)
library(rjags)
library(mcmcplots)
library(mcmcse)
library(GLMMadaptive)
library(future)
load("Cleaned data/pass.RData")
load("Cleaned data/pass_id.RData")
load("Cleaned data/pass_cores.RData")
load("R script/Data analysis/ICJM2_inits/mvmm_inits.RData")
source("R script/Functions/function.R")

plan(multisession, workers = 3)

# JAGS --------------------------------------------------------------------
## Parameters
qp = 15
kntbh_spec = "quantile"

## longitudinal 
lgt_cr <- pass.cores$TimeSince_Dx
lgt_psa <- pass$TimeSince_Dx
evt <- pass.id[,c("time.cmp1", "time.cmp2")]

# PSA
N_pt <- nrow(pass.id)
age <- pass$DxAge
num_psa <- as.vector(c(1, 1+ cumsum(tapply(pass$CISNET_ID, pass$CISNET_ID, length))))
Y <- pass$PSAValue
XL_psa <- array(1, dim = c(nrow(pass), 5))
knot.longi <- c(min(c(as.matrix(evt)), lgt_psa, lgt_cr), 
                as.vector(quantile(lgt_psa, probs = seq(0,1,length.out=4)))[2:3],
                max(c(as.matrix(evt)), lgt_psa, lgt_cr))

XL_psa[,2:4] <- ns(lgt_psa, knots = knot.longi[2:3], B = knot.longi[c(1,4)])
XL_psa[,5] <- pass$DxAge
ZL_psa <- XL_psa[,1:4]

# Cores_ratio
num_cr <- as.vector(c(1, 1+ cumsum(tapply(pass.cores$CISNET_ID, pass.cores$CISNET_ID, length))))
ZL_cr <- XL_cr <- model.matrix(~ poly(TimeSince_Dx, 2, raw = T), data = pass.cores)

## Time-to-event
evt <- pass.id[,c("time.cmp1", "time.cmp2")]
delta <- cbind(as.integer(pass.id$status.cmp == 1),
               as.integer(pass.id$status.cmp == 2),
               as.integer(pass.id$status.cmp == 0))
X <- pass.id$density
age.id <- pass.id$DxAge
f.org <- c(-0.949107912342758, -0.741531185599394, -0.405845151377397,
           0, 0.405845151377397, 0.741531185599394, 0.949107912342758,
           -0.991455371120813, -0.864864423359769, -0.586087235467691,
           -0.207784955007898, 0.207784955007898, 0.586087235467691,
           0.864864423359769, 0.991455371120813)
w <- c(0.0630920926299786, 0.140653259715526, 0.190350578064785,
       0.209482141084728, 0.190350578064785, 0.140653259715526,
       0.0630920926299786, 0.0229353220105292, 0.10479001032225,
       0.169004726639268, 0.204432940075299, 0.204432940075299,
       0.169004726639268, 0.10479001032225, 0.0229353220105292)
w <- w[order(f.org)]
f.org <- f.org[order(f.org)]

## number of knots = number of B(t,v)s  - 4L
get_knots <- function(nkn, Time, gkx) {
  knts <- quantile(Time, probs = seq(0,1, length.out = nkn + 2), names = F)
  knts <- tail(head(knts, -1), -1)
  return(sort(c(rep(range(Time, outer(Time/2, gkx + 1),0), 4),knts)))
}
knot.surv <- get_knots(8, c(as.matrix(evt)), gkx = f.org)

# left time and right time
lt <- evt[,1]; rt <- evt[,2]

## Treatment part for all three scenarios
evt.trt <- rt
evt.trt.qt <- sapply(evt.trt, function(x) {Qt(x, f.org)})
T.trt.s <- array(NA, c(N_pt, qp, 12)) # BH in the survival
for (i in 1:N_pt) {
  T.trt.s[i,,] <- bsp_dm(evt.trt.qt[,i], knot.surv)
}

T.trt.h <- bsp_dm(evt.trt, knot.surv) # BH in the hazard function: 833*12

# PSA
XL_psa.trt.s <- array(1, c(N_pt, qp, 5, 2)) # XL in the survival
# N_patient, quadrature, covartiates, slope
for (i in 1:N_pt) {
  XL_psa.trt.s[i,,2:4,1] <- ns_dm(evt.trt.qt[,i], knot.longi)
  XL_psa.trt.s[i,,2:4,2] <- ns_dm(evt.trt.qt[,i] - 1, knot.longi)
}
XL_psa.trt.s[,,5,] <- age.id
ZL_psa.trt.s <- XL_psa.trt.s[,,1:4,]  # ZL in the survival

XL_psa.trt.h <- array(1, c(N_pt, 5, 2)) # XL in the hazard
# N_patient, covartiates, slope
for (i in 1:N_pt) {
  XL_psa.trt.h[i,2:4,1] <- ns_dm(evt.trt[i], knot.longi)
  XL_psa.trt.h[i,2:4,2] <- ns_dm(evt.trt[i] - 1, knot.longi)
}
XL_psa.trt.h[,5,] <- age.id
ZL_psa.trt.h <- XL_psa.trt.h[,1:4,]

# Core ratio
XL_cr.trt.s <- array(1, c(N_pt, qp, 3)) # XL in the survival
# N_patient, quadrature, covartiates, slope
for (i in 1:N_pt) {
  XL_cr.trt.s[i,,2:3] <- poly(evt.trt.qt[,i], 2, raw = T)
}
ZL_cr.trt.s <- XL_cr.trt.s # ZL in the survival

XL_cr.trt.h <- array(1, c(N_pt, 3)) # XL in the hazard
# N_patient, covartiates, slope
for (i in 1:N_pt) {
  XL_cr.trt.h[i,2:3] <- poly(evt.trt[i], 2, raw = T)
}
ZL_cr.trt.h <- XL_cr.trt.h


## Progression for the exact survival
evt.prg <- lt
evt.prg.qt <- sapply(evt.prg, function(x) {Qt(x, f.org)})
T.prg.s <- array(NA, c(N_pt, qp, 12)) # BH in the survival
for (i in 1:N_pt) {
  T.prg.s[i,,] <- bsp_dm(evt.prg.qt[,i], knot.surv)
}

# PSA
XL_psa.prg.s <- array(1, c(N_pt, qp, 5, 2)) # XL in the survival
# N_patient, quadrature, covariates, slope
for (i in 1:N_pt) {
  XL_psa.prg.s[i,,2:4,1] <- ns_dm(evt.prg.qt[,i], knot.longi)
  XL_psa.prg.s[i,,2:4,2] <- ns_dm(evt.prg.qt[,i] - 1, knot.longi)
}
XL_psa.prg.s[,,5,] <- age.id
ZL_psa.prg.s <- XL_psa.prg.s[,,1:4,]  # ZL in the survival

# Core ratio
XL_cr.prg.s <- array(1, c(N_pt, qp, 3)) # XL in the survival
# N_patient, quadrature, covariates, slope
for (i in 1:N_pt) {
  XL_cr.prg.s[i,,2:3] <- poly(evt.prg.qt[,i], 2, raw = T)
}
ZL_cr.prg.s <- XL_cr.prg.s  # ZL in the survival

## Progression for the interval censoring part

# t.diff
evt.prginc.diff <- rt - lt

# hazard part
evt.prginc.qt <- array(NA, c(N_pt, qp)) # quadrature timepoints in the outer integral, t_{prg} to t_{trt}
T.prginc.h <- array(0, c(N_pt, qp, 12))
XL_psa.prginc.h <- array(0, c(N_pt, qp, 5, 2))
XL_cr.prginc.h <- array(0, c(N_pt, qp, 3))
for (i in 1:N_pt) {
  evt.prginc.qt[i,] <- Qt(rt = evt[i,2], lt = evt[i,1], qkq=f.org)
  T.prginc.h[i,,] <- bsp_dm(evt.prginc.qt[i,], knot.surv)
  # PSA
  XL_psa.prginc.h[i,,2:4,1] <- ns_dm(evt.prginc.qt[i,], knot.longi)
  XL_psa.prginc.h[i,,2:4,2] <- ns_dm(evt.prginc.qt[i,] - 1, knot.longi)
  XL_psa.prginc.h[i,,1,] <- 1
  XL_psa.prginc.h[i,,5,] <- age.id[i]
  
  # Core ratio
  XL_cr.prginc.h[i,,2:3] <- poly(evt.prginc.qt[i,], 2, raw = T)
  XL_cr.prginc.h[i,,1] <- 1
}
ZL_psa.prginc.h <- XL_psa.prginc.h[,,1:4,]
ZL_cr.prginc.h <- XL_cr.prginc.h

# survival part
evt.prginc.qt.qt <- array(NA, c(N_pt, qp, qp)) # first qp is the qds from the outer integral; second is from [0,qds]
T.prginc.s <- array(0, c(N_pt, qp, 12, qp)) # first qp is the qds from the outer integral; second is from [0,qds]
XL_psa.prginc.s <- array(0, c(N_pt, qp, 5, qp, 2)) # first qp is the qds from the outer integral; second is from [0,qds]
XL_cr.prginc.s <- array(0, c(N_pt, qp, 3, qp)) # first qp is the qds from the outer integral; second is from [0,qds]
for (i in 1:N_pt) {
  for (m in 1:qp) {
    evt.prginc.qt.qt[i,m,] <- Qt(rt = evt.prginc.qt[i,m], qkq=f.org)
    T.prginc.s[i,m,,] <- t(bsp_dm(evt.prginc.qt.qt[i,m,], knot.surv))
    # PSA
    XL_psa.prginc.s[i,m,2:4,,1] <- t(ns_dm(evt.prginc.qt.qt[i,m,], knot.longi))
    XL_psa.prginc.s[i,m,2:4,,2] <- t(ns_dm(evt.prginc.qt.qt[i,m,] - 1, knot.longi))
    XL_psa.prginc.s[i,,1,,] <- 1
    XL_psa.prginc.s[i,,5,,] <- age.id[i]
    # Core ratio
    XL_cr.prginc.s[i,m,2:3,] <- t(poly(evt.prginc.qt.qt[i,m,], 2, raw = T))
    XL_cr.prginc.s[i,,1,] <- 1
  }
}
ZL_psa.prginc.s <- XL_psa.prginc.s[,,1:4,,]
ZL_cr.prginc.s <- XL_cr.prginc.s


# model specification -----------------------------------------------------

model <- "model {
  for (i in 1:N_pt) {
    for (j in num_psa[i]:(num_psa[i+1]-1)) {
      Y[j] ~ dt(mu[j], tau, 3)
      mu[j] <- inprod(XL_psa[j,], betaL_psa[]) + inprod(ZL_psa[j,], b[i,1:4])

    }
    
    for (j in num_cr[i]:(num_cr[i+1]-1)) {
      pos[j] ~ dbin(p[j], total[j])
      logit(p[j]) <- inprod(XL_cr[j,], betaL_cr[]) + inprod(ZL_cr[j,], b[i,5:7])
    }
    
    #### Time-to-event part
    
    ### Treatment part
    ## Survival part
    log.s.trt.gk[i,1:qp] <- exp(1) ^ (T.trt.s[i,1:qp,] %*% gambh[,2] + X[i] * gamma[2] +
                                alpha[1,2] * (
                                      XL_psa.trt.s[i,1:qp,,1] %*% betaL_psa[] + 
                                      ZL_psa.trt.s[i,1:qp,,1] %*% b[i,1:4]
                                ) + 
                                alpha[2,2] * (
                                      XL_psa.trt.s[i,1:qp,,1] %*% betaL_psa[] + 
                                      ZL_psa.trt.s[i,1:qp,,1] %*% b[i,1:4] - 
                                      XL_psa.trt.s[i,1:qp,,2] %*% betaL_psa[] - 
                                      ZL_psa.trt.s[i,1:qp,,2] %*% b[i,1:4]
                                ) + 
                                alpha[3,2] * (
                                      XL_cr.trt.s[i,1:qp,] %*% betaL_cr[] + 
                                      ZL_cr.trt.s[i,1:qp,] %*% b[i,5:7]
                                ))
                                      
    log.s.trt[i] <- - evt.trt[i]/2 * (t(log.s.trt.gk[i,]) %*% w )
    
    ## Hazard part
    log.h.trt[i] <- inprod(T.trt.h[i, ], gambh[,2]) + X[i] * gamma[2] + 
                  alpha[1,2] * (
                    inprod(XL_psa.trt.h[i,,1], betaL_psa[]) + inprod(ZL_psa.trt.h[i,,1], b[i,1:4])) +
                  alpha[2,2] * (
                    inprod(XL_psa.trt.h[i,,1], betaL_psa[]) + inprod(ZL_psa.trt.h[i,,1], b[i,1:4]) - 
                    inprod(XL_psa.trt.h[i,,2], betaL_psa[]) - inprod(ZL_psa.trt.h[i,,2], b[i,1:4])) + 
                  alpha[3,2] * (
                    inprod(XL_cr.trt.h[i,], betaL_cr[]) + inprod(ZL_cr.trt.h[i,], b[i,5:7]))  
                    
    ### Exact progression part
    ## Survival part
    log.s.prg.gk[i,1:qp] <- exp(1) ^ (T.prg.s[i,1:qp,] %*% gambh[,1] + X[i] * gamma[1] +
                                alpha[1,1] * (
                                      XL_psa.prg.s[i,1:qp,,1] %*% betaL_psa[] + 
                                      ZL_psa.prg.s[i,1:qp,,1] %*% b[i,1:4]) + 
                                alpha[2,1] * (
                                      XL_psa.prg.s[i,1:qp,,1] %*% betaL_psa[] + 
                                      ZL_psa.prg.s[i,1:qp,,1] %*% b[i,1:4] - 
                                      XL_psa.prg.s[i,1:qp,,2] %*% betaL_psa[] - 
                                      ZL_psa.prg.s[i,1:qp,,2] %*% b[i,1:4]) + 
                                alpha[3,1] * (
                                      XL_cr.prg.s[i, 1:qp,] %*% betaL_cr[] + 
                                      ZL_cr.prg.s[i, 1:qp,] %*% b[i, 5:7])) 
    log.s.prg[i] <- - evt.prg[i]/2 * (t(log.s.prg.gk[i,]) %*% w)
    
    ### Interval-censored progression part
    ## Survival part
    for (m in 1:qp) {
      log.s.prginc.int.gk.gk[i,m,1:qp] <- exp(1) ^ (t(T.prginc.s[i,m,,1:qp]) %*% gambh[,1] + X[i] * gamma[1] + 
                                              alpha[1,1] * (t(XL_psa.prginc.s[i,m,,1:qp,1]) %*% betaL_psa[] + 
                                            t(ZL_psa.prginc.s[i,m,,1:qp,1]) %*% b[i,1:4]) + 
                                              alpha[2,1] * (
                                            t(XL_psa.prginc.s[i,m,,1:qp,1]) %*% betaL_psa[] + 
                                            t(ZL_psa.prginc.s[i,m,,1:qp,1]) %*% b[i,1:4] - 
                                            t(XL_psa.prginc.s[i,m,,1:qp,2]) %*% betaL_psa[] - 
                                            t(ZL_psa.prginc.s[i,m,,1:qp,2]) %*% b[i,1:4]) +
                                              alpha[3,1] * (t(XL_cr.prginc.s[i,m,,1:qp]) %*% betaL_cr[] + 
                                            t(ZL_cr.prginc.s[i,m,,1:qp]) %*% b[i,5:7]))
      s.prginc.int.gk[i,m] <- exp(1) ^ (- evt.prginc.qt[i,m]/2 * inprod(w, log.s.prginc.int.gk.gk[i,m,]))
    }
    
    h.prginc.int.gk[i,1:qp] <- exp(1) ^ (T.prginc.h[i,1:qp,] %*% gambh[,1] + X[i] * gamma[1] +  
                                 alpha[1,1] * (
                                   XL_psa.prginc.h[i,1:qp,,1] %*% betaL_psa[] +
                                   ZL_psa.prginc.h[i,1:qp,,1] %*% b[i,1:4]) + 
                                 alpha[2,1] * (
                                   XL_psa.prginc.h[i,1:qp,,1] %*% betaL_psa[] +
                                   ZL_psa.prginc.h[i,1:qp,,1] %*% b[i,1:4] - 
                                   XL_psa.prginc.h[i,1:qp,,2] %*% betaL_psa[] -
                                   ZL_psa.prginc.h[i,1:qp,,2] %*% b[i,1:4]) + 
                                alpha[3,1] * (
                                   XL_cr.prginc.h[i,1:qp,] %*% betaL_cr[] +
                                   ZL_cr.prginc.h[i,1:qp,] %*% b[i,5:7 ]))
                                   
    prginc.int[i] <- evt.prginc.diff[i]/2 * inprod(w, s.prginc.int.gk[i,] * h.prginc.int.gk[i,]) 
    prginc[i] <- ifelse(prginc.int[i] < 1e-16,
                        1e-16,
                        prginc.int[i])
    
    #### LIKELIHOOD
    loglike[i] <- delta[i,1] * (log(prginc[i]) + log.s.trt[i]) + 
                  delta[i,2] * (log.h.trt[i] + log.s.trt[i] + log.s.prg[i]) +
                  delta[i,3] * (log.s.trt[i] + log.s.prg[i])
    phi[i] <- 5000 - loglike[i]
    zeros[i] ~ dpois(phi[i])
    
    b[i,1:Nb] ~ dmnorm(mub[], inv_D[,])
  }
  
  #### PRIORS
  tau ~ dgamma(0.01,0.01)
  inv_D[1:Nb,1:Nb] ~ dwish(4 * diagtb[,], Nb+1)
  for (l in 1:Nb) {
    for (k in 1:Nb) {diagtb[l,k] <- ifelse(l == k, tb[l], 0)}
    tb[l] ~ dgamma(0.5, 0.01)
  }
  
  D[1:Nb,1:Nb] <- inverse(inv_D[,])
  
  for (l in 1:Nbeta_psa) {betaL_psa[l] ~ dnorm(0, 0.01)}
  for (l in 1:Nbeta_cr) {betaL_cr[l] ~ dnorm(0, 0.01)}
  for (k in 1:Nrisk) {
    gambh[1:Ngambh,k] ~ dmnorm(mu.gambh[], tau.smooth[k] * tau.gambh[,])
    tau.smooth[k] ~ dgamma(5, 0.5)
    alpha[1,k] ~ dnorm(0, 0.01) # tau_alpha[1,k] * psi_alpha[1,k]
    alpha[2,k] ~ dnorm(0, 0.01) #tau_alpha[2,k] * psi_alpha[2,k]
    alpha[3,k] ~ dnorm(0, 0.01) #tau_alpha[3,k] * psi_alpha[3,k]
    
    #tau_alpha[1,k] ~ dgamma(0.1, 0.1)
    #tau_alpha[2,k] ~ dgamma(0.1, 0.1)
    #tau_alpha[3,k] ~ dgamma(0.1, 0.1)
    
    #psi_alpha[1,k] ~ dgamma(1, 0.01)
    #psi_alpha[2,k] ~ dgamma(1, 0.01)
    #psi_alpha[3,k] ~ dgamma(1, 0.01)
    
    gamma[k] ~ dnorm(0, 0.01)
  }
  
}"

DD <- diag(12)

data <- list(N_pt = N_pt, qp = qp,
             num_psa = num_psa, Y =Y, XL_psa = XL_psa, ZL_psa = ZL_psa,
             num_cr = num_cr, total = pass.cores$total_cores,
             pos = pass.cores$pos_cores, XL_cr = XL_cr, ZL_cr = ZL_cr,
             Nrisk = ncol(delta)-1, 
             T.trt.s = T.trt.s, X = X, XL_psa.trt.s = XL_psa.trt.s, ZL_psa.trt.s = ZL_psa.trt.s,
             XL_cr.trt.s = XL_cr.trt.s, ZL_cr.trt.s = ZL_cr.trt.s,
             evt.trt = evt.trt, w = w, 
             T.trt.h = T.trt.h, XL_psa.trt.h = XL_psa.trt.h, ZL_psa.trt.h = ZL_psa.trt.h,
             XL_cr.trt.h = XL_cr.trt.h, ZL_cr.trt.h = ZL_cr.trt.h,
             T.prg.s = T.prg.s, XL_psa.prg.s = XL_psa.prg.s, ZL_psa.prg.s = ZL_psa.prg.s, 
             XL_cr.prg.s = XL_cr.prg.s, ZL_cr.prg.s = ZL_cr.prg.s, 
             evt.prg= evt.prg,
             T.prginc.s = T.prginc.s, XL_psa.prginc.s = XL_psa.prginc.s, ZL_psa.prginc.s = ZL_psa.prginc.s,
             XL_cr.prginc.s = XL_cr.prginc.s, ZL_cr.prginc.s = ZL_cr.prginc.s,
             evt.prginc.qt = evt.prginc.qt,
             T.prginc.h = T.prginc.h, XL_psa.prginc.h = XL_psa.prginc.h, ZL_psa.prginc.h = ZL_psa.prginc.h,
             XL_cr.prginc.h = XL_cr.prginc.h, ZL_cr.prginc.h = ZL_cr.prginc.h,
             evt.prginc.diff = evt.prginc.diff, 
             delta = delta, zeros = rep(0, N_pt),
             mub = rep(0,dim(ZL_psa)[2] + dim(ZL_cr)[2]), 
             Nb = dim(ZL_psa)[2] + dim(ZL_cr)[2], 
             Nbeta_psa = dim(XL_psa)[2], Nbeta_cr = dim(XL_cr)[2], 
             Ngambh = 12,
             mu.gambh = rep(0,12), tau.gambh = crossprod(diff(DD, diff = 2)) + 1e-6* DD)

set.seed(100)
rng.name <- sample(c("base::Wichmann-Hill",
                     "base::Marsaglia-Multicarry",
                     "base::Super-Duper",
                     "base::Mersenne-Twister"), 3, replace = T)
rng.num <- sample(1:100000, 3, replace = F)

initials <- list(
  list(alpha = matrix(rnorm(6), 3, 2),
       #tau_alpha = matrix(runif(8), 4, 2),
       #psi_alpha = matrix(runif(8), 4, 2),
       betaL_psa= parm.mvmm$betaL_psa,
       betaL_cr= parm.mvmm$betaL_cr,
       gambh = matrix(rnorm(24),12,2),
       gamma = rnorm(2),
       b = parm.mvmm$b,
       inv_D = solve(parm.mvmm$D),
       tau = parm.mvmm$tau,
       .RNG.name = rng.name[1],
       .RNG.seed = rng.num[1]),
  list(alpha = matrix(rnorm(6), 3, 2),
       #tau_alpha = matrix(runif(8), 4, 2),
       #psi_alpha = matrix(runif(8), 4, 2),
       betaL_psa= parm.mvmm$betaL_psa,
       betaL_cr= parm.mvmm$betaL_cr,
       gambh = matrix(rnorm(24),12,2),
       gamma = rnorm(2),
       b = parm.mvmm$b,
       inv_D = solve(parm.mvmm$D),
       tau = parm.mvmm$tau,
       .RNG.name = rng.name[2],
       .RNG.seed = rng.num[2]),
  list(alpha = matrix(rnorm(6), 3, 2),
       #tau_alpha = matrix(runif(8), 4, 2),
       #psi_alpha = matrix(runif(8), 4, 2),
       betaL_psa= parm.mvmm$betaL_psa,
       betaL_cr= parm.mvmm$betaL_cr,
       gambh = matrix(rnorm(24),12,2),
       gamma = rnorm(2),
       b = parm.mvmm$b,
       inv_D = solve(parm.mvmm$D),
       tau = parm.mvmm$tau,
       .RNG.name = rng.name[3],
       .RNG.seed = rng.num[3])
)

time.start <- Sys.time()
out <- lapply(1:3, function(i) {
  future::future({
    icjm.psa.cr <- jags.model(file = textConnection(model), data = data,
                              inits = initials[[i]], n.chains = 1, n.adapt = 5000,
                              quiet = T)  #
    update(icjm.psa.cr, n.iter = 5000)
    mcmc <- coda.samples(icjm.psa.cr, variable.names = c("alpha", "gambh", "D",
                                                         "gamma", "tau", 
                                                         "betaL_psa", "betaL_cr"),
                         n.iter = 15000, thin = 50)
    return(list(model = icjm.psa.cr,
                mcmc = mcmc))
  })
})
samples.icjm.psa.cr <- lapply(out, future::value)
time.end <- Sys.time()

time.end - time.start

mcmc.mat <- do.call(rbind, lapply(samples.icjm.psa.cr, 
                                  function(x) {as.matrix(x$mcmc)}))
mcmc_all <- as.mcmc.list(lapply(1:3, function(x) 
  as.mcmc(samples.icjm.psa.cr[[x]]$mcmc)))

ICJM2 <- list(mcmc = lapply(samples.icjm.psa.cr, function(x) {x$mcmc}),
              model_info = list(var_names = list(id = "CISNET_ID",
                                                 longitime = "TimeSince_Dx",
                                                 survtime = c("time.cmp1", "time.cmp2"),
                                                 longi.bsVar = "DxAge",
                                                 surv.bsVar = "density"),
                                knots = list(knot.longi = knot.longi,
                                             knot.surv = knot.surv)),
              inits = initials,
              summary = list(coef = apply(mcmc.mat, 2, mean),
                             coef.sd = apply(mcmc.mat, 2, sd),
                             gelman_rubin = my.gelman.diag(mcmc_all),
                             se = mcse.mat(mcmc.mat)[,2]))
#jags_model = lapply(samples.icjm.psa.cr, function(x) {x$model}))
mcmcplot(ICJM2$mcmc)

save(ICJM2, file = "Output/Data analysis/ICJM2.RData")
ICJM2$summary$coef
