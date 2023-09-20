iccsjm <- function(df = 3, data, summary = T, n.adapt = 3000,
                   n.burnin = 3000, n.iter = 10000, 
                   qp = 15, kntbh_spec = "quantile", 
                   reduced_mcmc = FALSE, itersaved_perchain = 200,
                   seed) {
  
  # Load packages -----------------------------------------------------------
  library(splines)
  library(GLMMadaptive)
  library(rjags)
  library(mcmcse)
  library(future)
  source("R script/Functions/function.R")
  
  
  # Check argument ----------------------------------------------------------
  
  if (!(qp %in% c(7, 15))) {
    stop("The number of quadrature points should either be 7 or 15!")
  }
  
  if (!(kntbh_spec %in% c("quantile", "equal"))) {
    stop("Knots in the baseline hazard function should either be 'quantile' or 'equal'!")
  }
  source("R script/Functions/function.R")
  
  # Start with the mixed model ----------------------------------------------
  bound <- range(data$TimeSince_Dx, data$time.cmp1, data$time.cmp2)
  data[,(ncol(data)+1):(ncol(data)+3)] <- ns(data$TimeSince_Dx, 3, B = bound)
  colnames(data)[(ncol(data)-2):ncol(data)] <- c("TimeSince_Dx.1",
                                                 "TimeSince_Dx.2",
                                                 "TimeSince_Dx.3")
  mm <- mixed_model(fixed = PSAValue ~ TimeSince_Dx.1 + TimeSince_Dx.2 + 
                      TimeSince_Dx.3 + DxAge,
                    data = data, 
                    random = ~ TimeSince_Dx.1 + TimeSince_Dx.2 + 
                      TimeSince_Dx.3 | CISNET_ID,
                    family = students.t(df = df), n_phis = 1,
                    initial_values = list("betas" = gaussian()))
  #print("mixed model done!")
  
  # Data preparation --------------------------------------------------------
  ## Quadrature knots
  f.org <- JMbayes:::gaussKronrod(qp)$sk
  w <- JMbayes:::gaussKronrod(qp)$wk
  w <- w[order(f.org)]
  f.org <- f.org[order(f.org)]
  
  ## data frames
  data.id <- data[!duplicated(data$CISNET_ID),]
  N_pt <- nrow(data.id)
  lgt <- data$TimeSince_Dx
  evt <- data.id[,c("time.cmp1", "time.cmp2")]
  age <- data$DxAge
  num <- as.vector(c(1, 1+ cumsum(tapply(data$CISNET_ID, data$CISNET_ID, length))))
  Y <- data$PSAValue
  XL <- array(1, dim = c(nrow(data), 5))
  knot.longi <- c(min(c(as.matrix(evt)), lgt), 
                  as.vector(quantile(lgt, probs = seq(0,1,length.out=4)))[2:3],
                  max(c(as.matrix(evt)), lgt))
  
  XL[,2:4] <- ns(lgt, knots = knot.longi[2:3], B = knot.longi[c(1,4)])
  XL[,5] <- data$DxAge
  ZL <- XL[,1:4]
  
  ## event_time
  delta <- cbind(as.integer(data.id$status.cmp == 1),
                 as.integer(data.id$status.cmp == 2),
                 as.integer(data.id$status.cmp == 0))
  X <- data.id$density
  age.id <- data.id$DxAge
  knot.surv <- get_knots(8, c(as.matrix(evt)), gkx = f.org, bh = kntbh_spec)
  
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
  
  XL.trt.s <- array(1, c(N_pt, qp, 5, 2)) # XL in the survival
  # N_patient, quadrature, covartiates, slope
  for (i in 1:N_pt) {
    XL.trt.s[i,,2:4,1] <- ns_dm(evt.trt.qt[,i], knot.longi)
    XL.trt.s[i,,2:4,2] <- ns_dm(evt.trt.qt[,i] - 1, knot.longi)
  }
  XL.trt.s[,,5,] <- age.id
  ZL.trt.s <- XL.trt.s[,,1:4,]  # ZL in the survival
  
  XL.trt.h <- array(1, c(N_pt, 5, 2)) # XL in the hazard
  # N_patient, covartiates, slope
  for (i in 1:N_pt) {
    XL.trt.h[i,2:4,1] <- ns_dm(evt.trt[i], knot.longi)
    XL.trt.h[i,2:4,2] <- ns_dm(evt.trt[i] - 1, knot.longi)
  }
  XL.trt.h[,5,] <- age.id
  ZL.trt.h <- XL.trt.h[,1:4,]
  
  ## Progression for the exact survival
  evt.prg <- lt
  evt.prg.qt <- sapply(evt.prg, function(x) {Qt(x, f.org)})
  T.prg.s <- array(NA, c(N_pt, qp, 12)) # BH in the survival
  for (i in 1:N_pt) {
    T.prg.s[i,,] <- bsp_dm(evt.prg.qt[,i], knot.surv)
  }
  
  XL.prg.s <- array(1, c(N_pt, qp, 5, 2)) # XL in the survival
  # N_patient, quadrature, covariates, slope
  for (i in 1:N_pt) {
    XL.prg.s[i,,2:4,1] <- ns_dm(evt.prg.qt[,i], knot.longi)
    XL.prg.s[i,,2:4,2] <- ns_dm(evt.prg.qt[,i] - 1, knot.longi)
  }
  XL.prg.s[,,5,] <- age.id
  ZL.prg.s <- XL.prg.s[,,1:4,]  # ZL in the survival
  
  ## Progression for the interval censoring part
  
  # t.diff
  evt.prginc.diff <- rt - lt
  
  # hazard part
  evt.prginc.qt <- array(NA, c(N_pt, qp)) # quadrature timepoints in the outer integral, t_{prg} to t_{trt}
  T.prginc.h <- array(0, c(N_pt, qp, 12))
  XL.prginc.h <- array(0, c(N_pt, qp, 5, 2))
  for (i in 1:N_pt) {
    evt.prginc.qt[i,] <- Qt(rt = evt[i,2], lt = evt[i,1], qkq=f.org)
    T.prginc.h[i,,] <- bsp_dm(evt.prginc.qt[i,], knot.surv)
    XL.prginc.h[i,,2:4,1] <- ns_dm(evt.prginc.qt[i,], knot.longi)
    XL.prginc.h[i,,2:4,2] <- ns_dm(evt.prginc.qt[i,] - 1, knot.longi)
    XL.prginc.h[i,,1,] <- 1
    XL.prginc.h[i,,5,] <- age.id[i]
  }
  ZL.prginc.h <- XL.prginc.h[,,1:4,]
  
  # survival part
  evt.prginc.qt.qt <- array(NA, c(N_pt, qp, qp)) # first qp is the qds from the outer integral; second is from [0,qds]
  T.prginc.s <- array(0, c(N_pt, qp, 12, qp)) # first qp is the qds from the outer integral; second is from [0,qds]
  XL.prginc.s <- array(0, c(N_pt, qp, 5, qp, 2)) # first qp is the qds from the outer integral; second is from [0,qds]
  for (i in 1:N_pt) {
    for (m in 1:qp) {
      evt.prginc.qt.qt[i,m,] <- Qt(rt = evt.prginc.qt[i,m], qkq=f.org)
      T.prginc.s[i,m,,] <- t(bsp_dm(evt.prginc.qt.qt[i,m,], knot.surv))
      XL.prginc.s[i,m,2:4,,1] <- t(ns_dm(evt.prginc.qt.qt[i,m,], knot.longi))
      XL.prginc.s[i,m,2:4,,2] <- t(ns_dm(evt.prginc.qt.qt[i,m,] - 1, knot.longi))
      XL.prginc.s[i,,1,,] <- 1
      XL.prginc.s[i,,5,,] <- age.id[i]
    }
  }
  ZL.prginc.s <- XL.prginc.s[,,1:4,,]
  
  #print("data preparation done!")
  
  # Start JAGS --------------------------------------------------------------
  model <- "model {
  for (i in 1:N_pt) {
    for (j in num[i]:(num[i+1]-1)) {
      Y[j] ~ dt(mu[j], tau, 3)
      mu[j] <- inprod(XL[j,], betaL[]) + inprod(ZL[j,], b[i,])
    }
    
    #### Time-to-event part
    
    ### Treatment part
    ## Survival part
    log.s.trt.gk[i,1:qp] <- exp(1) ^ (T.trt.s[i,1:qp,] %*% gambh[,2] + X[i] * gamma[2] +
                                alpha[1,2] * (
                                      XL.trt.s[i,1:qp,,1] %*% betaL[] + 
                                      ZL.trt.s[i,1:qp,,1] %*% b[i,]) + 
                                alpha[2,2] * (
                                      XL.trt.s[i,1:qp,,1] %*% betaL[] + 
                                      ZL.trt.s[i,1:qp,,1] %*% b[i,] - 
                                      XL.trt.s[i,1:qp,,2] %*% betaL[] - 
                                      ZL.trt.s[i,1:qp,,2] %*% b[i,]))
                                      
    log.s.trt[i] <- - evt.trt[i]/2 * (t(log.s.trt.gk[i,]) %*% w )
    
    ## Hazard part
    log.h.trt[i] <- inprod(T.trt.h[i, ], gambh[,2]) + X[i] * gamma[2] + 
                  alpha[1,2] * (
                    inprod(XL.trt.h[i,,1], betaL[]) + inprod(ZL.trt.h[i,,1], b[i,])) +
                  alpha[2,2] * (
                    inprod(XL.trt.h[i,,1], betaL[]) + inprod(ZL.trt.h[i,,1], b[i,]) - 
                    inprod(XL.trt.h[i,,2], betaL[]) - inprod(ZL.trt.h[i,,2], b[i,]))
                    
    ### Exact progression part
    ## Survival part
    log.s.prg.gk[i,1:qp] <- exp(1) ^ (T.prg.s[i,1:qp,] %*% gambh[,1] + X[i] * gamma[1] +
                                alpha[1,1] * (
                                      XL.prg.s[i,1:qp,,1] %*% betaL[] + 
                                      ZL.prg.s[i,1:qp,,1] %*% b[i,]) + 
                                alpha[2,1] * (
                                      XL.prg.s[i,1:qp,,1] %*% betaL[] + 
                                      ZL.prg.s[i,1:qp,,1] %*% b[i,] - 
                                      XL.prg.s[i,1:qp,,2] %*% betaL[] - 
                                      ZL.prg.s[i,1:qp,,2] %*% b[i,]))
    log.s.prg[i] <- - evt.prg[i]/2 * (t(log.s.prg.gk[i,]) %*% w)
    
    ### Interval-censored progression part
    ## Survival part
    for (m in 1:qp) {
      log.s.prginc.int.gk.gk[i,m,1:qp] <- exp(1) ^ (t(T.prginc.s[i,m,,1:qp]) %*% gambh[,1] + X[i] * gamma[1] + 
                                              alpha[1,1] * (t(XL.prginc.s[i,m,,1:qp,1]) %*% betaL[] + 
                                            t(ZL.prginc.s[i,m,,1:qp,1]) %*% b[i,]) + 
                                              alpha[2,1] * (
                                            t(XL.prginc.s[i,m,,1:qp,1]) %*% betaL[] + 
                                            t(ZL.prginc.s[i,m,,1:qp,1]) %*% b[i,] - 
                                            t(XL.prginc.s[i,m,,1:qp,2]) %*% betaL[] - 
                                            t(ZL.prginc.s[i,m,,1:qp,2]) %*% b[i,]))
      s.prginc.int.gk[i,m] <- exp(1) ^ (- evt.prginc.qt[i,m]/2 * inprod(w, log.s.prginc.int.gk.gk[i,m,]))
    }
    
    h.prginc.int.gk[i,1:qp] <- exp(1) ^ (T.prginc.h[i,1:qp,] %*% gambh[,1] + X[i] * gamma[1] +  
                                 alpha[1,1] * (
                                   XL.prginc.h[i,1:qp,,1] %*% betaL[] +
                                   ZL.prginc.h[i,1:qp,,1] %*% b[i,]) + 
                                 alpha[2,1] * (
                                   XL.prginc.h[i,1:qp,,1] %*% betaL[] +
                                   ZL.prginc.h[i,1:qp,,1] %*% b[i,] - 
                                   XL.prginc.h[i,1:qp,,2] %*% betaL[] -
                                   ZL.prginc.h[i,1:qp,,2] %*% b[i,]))
                                   
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
  
  for (l in 1:Nbeta) {betaL[l] ~ dnorm(0, 0.01)}
  for (k in 1:Nrisk) {
    gambh[1:Ngambh,k] ~ dmnorm(mu.gambh[], tau.smooth[k] * tau.gambh[,])
    tau.smooth[k] ~ dgamma(5, 0.5)
    alpha[1,k] ~ dnorm(0, 0.01)
    alpha[2,k] ~ dnorm(0, 0.01)
    gamma[k] ~ dnorm(0, 0.01)
  }
  
}"
  
  DD <- diag(12)
  
  data <- list(N_pt = N_pt, qp = qp, num = num, Y = Y, XL = XL, ZL = ZL,
               T.trt.s = T.trt.s, X = X, XL.trt.s = XL.trt.s, ZL.trt.s = ZL.trt.s,
               evt.trt = evt.trt, w = w, 
               T.trt.h = T.trt.h, XL.trt.h = XL.trt.h, ZL.trt.h = ZL.trt.h,
               T.prg.s = T.prg.s, XL.prg.s = XL.prg.s, ZL.prg.s = ZL.prg.s, 
               evt.prg= evt.prg,
               T.prginc.s = T.prginc.s, XL.prginc.s = XL.prginc.s, ZL.prginc.s = ZL.prginc.s,
               evt.prginc.qt = evt.prginc.qt,
               T.prginc.h = T.prginc.h, XL.prginc.h = XL.prginc.h, ZL.prginc.h = ZL.prginc.h,
               evt.prginc.diff = evt.prginc.diff, 
               delta = delta, zeros = rep(0, N_pt),
               mub = rep(0,dim(ZL)[2]), Nb = dim(ZL)[2], Nbeta = dim(XL)[2], 
               Nrisk = ncol(delta) - 1, Ngambh = 12,
               mu.gambh = rep(0,12), tau.gambh = crossprod(diff(DD, diff = 2)) + 1e-6* DD)
  
  set.seed(seed)
  rng.name <- sample(c("base::Wichmann-Hill",
                       "base::Marsaglia-Multicarry",
                       "base::Super-Duper",
                       "base::Mersenne-Twister"), 3, replace = T)
  rng.num <- sample(1:100000, 3, replace = F)
  
  initials <- list(
    list(alpha = matrix(rnorm(4), 2, 2),
         betaL= fixef(mm),
         gambh = matrix(rnorm(24),12,2),
         gamma = rnorm(2),
         b = ranef(mm),
         inv_D = solve(mm$D),
         tau = 1/exp(mm$phis)^2,
         .RNG.name = rng.name[1],
         .RNG.seed = rng.num[1]),
    list(alpha = matrix(rnorm(4), 2, 2),
         betaL= fixef(mm),
         gambh = matrix(rnorm(24),12,2),
         gamma = rnorm(2),
         b = ranef(mm),
         inv_D = solve(mm$D),
         tau = 1/exp(mm$phis)^2,
         .RNG.name = rng.name[2],
         .RNG.seed = rng.num[2]),
    list(alpha = matrix(rnorm(4), 2, 2),
         betaL= fixef(mm),
         gambh = matrix(rnorm(24),12,2),
         gamma = rnorm(2),
         b = ranef(mm),
         inv_D = solve(mm$D),
         tau = 1/exp(mm$phis)^2,
         .RNG.name = rng.name[3],
         .RNG.seed = rng.num[3])
  )
  
  
  out <- lapply(1:3, function(i) {
    future::future({
      icjm <- jags.model(file = textConnection(model), data = data,
                         inits = initials[[i]], n.chains = 1, n.adapt = n.adapt,
                         quiet = T)  #
      update(icjm, n.iter = n.burnin)
      mcmc <- coda.samples(icjm, variable.names = c("alpha", "gambh", "D",
                                                    "gamma", "tau", "betaL"),
                           n.iter = n.iter, thin = 10)
      return(list(model = icjm,
                  mcmc = mcmc))
    })
  })
  
  samples.icjm <- lapply(out, future::value)
  #print("JAGS model done!")
  iter.perchain <- nrow(as.matrix(samples.icjm[[1]]$mcmc))
  if (reduced_mcmc) {
    mcmcoutput <- lapply(samples.icjm, function(x) {
      as.mcmc(as.matrix(x)[(iter.perchain-itersaved_perchain+1):iter.perchain,])
    })
  } else {
    mcmcoutput <- lapply(samples.icjm, function(x) {x$mcmc})
  }
  
  if (summary == F) {
    return(list(mcmc = mcmcoutput,
                model_info = list(var_names = list(id = "CISNET_ID",
                                                   longitime = "TimeSince_Dx",
                                                   survtime = c("time.cmp1", "time.cmp2"),
                                                   longi.bsVar = "DxAge",
                                                   surv.bsVar = "density"),
                                  knots = list(knot.longi = knot.longi,
                                               knot.surv = knot.surv)),
                #jags_model = lapply(samples.icjm, function(x) {x$model}),
                inits = initials))
  } else if (summary == T) {
    mcmc.mat <- do.call(rbind, lapply(samples.icjm, 
                                      function(x) {as.matrix(x$mcmc)}))
    mcmc_all <- as.mcmc.list(lapply(1:3, function(x) 
      as.mcmc(samples.icjm[[x]]$mcmc)))
    return(list(mcmc = mcmcoutput,
                model_info = list(var_names = list(id = "CISNET_ID",
                                                   longitime = "TimeSince_Dx",
                                                   survtime = c("time.cmp1", "time.cmp2"),
                                                   longi.bsVar = "DxAge",
                                                   surv.bsVar = "density"),
                                  knots = list(knot.longi = knot.longi,
                                               knot.surv = knot.surv)),
                #jags_model = lapply(samples.icjm, function(x) {x$model}),
                inits = initials,
                summary = list(coef = apply(mcmc.mat, 2, mean),
                               coef.sd = apply(mcmc.mat, 2, sd),
                               gelman_rubin = my.gelman.diag(mcmc_all),
                               se = mcse.mat(mcmc.mat)[,2])))
  } 
  
}
