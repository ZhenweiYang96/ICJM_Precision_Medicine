# use JAGS results --------------------------------------------------------
rm(list=ls())
library(MASS)
library(JMbayes)
library(rjags)
library(Matrix)
library(splines)
library(truncnorm)
library(future)
library(doParallel)

cl <- detectCores() - 1
plan(multisession, workers = cl) 
registerDoParallel(cl)  

# load full model ---------------------------------------------------------
load("Output/Data analysis/ICJM1.RData")
mcmc <- do.call(rbind, lapply(1:3, function(x) {
  as.matrix(ICJM1$mcmc[[x]])
}))
# get parameters together
param <- apply(mcmc, 2, mean)

# covaraince matrix for random effects
D_c3 <- matrix(param[grep("D", names(param))], 4, 4)

# Prepare seeds -----------------------------------------------------------
set.seed(2022)
random.seed <- sample(c(5001:1000000), 1000, replace = FALSE)
n.dataset <- 200
save(random.seed, file = "R script/Simulation/seed/seed_data.RData")


# Simulate for 200 test sets -----------------------------------------------

t1 <- Sys.time()
out <- lapply(1:n.dataset, function(datanum) {
  future({
    load("Cleaned data/pass.RData")
  load("Cleaned data/pass_id.RData")
  rm(list=setdiff(ls(), c("D_c3", "random.seed", "datanum", "n.dataset", "param",
                          "pass", "pass.id", "ICJM1","t1")))
  
  # quadrature
  f.org <- JMbayes:::gaussKronrod()$sk
  w <- JMbayes:::gaussKronrod()$wk
  knot.longi <- ICJM1$model_info$knots$knot.longi
  knot.surv <- ICJM1$model_info$knots$knot.surv
  
  # seed
  random_seed_inuse <- random.seed[datanum]
  set.seed(random_seed_inuse)
  
  # start preparation
  n <- 1000
  t.max <- max(pass$TimeSince_Dx, pass$time.cmp1, pass$time.cmp2) # 12.48
  
  betas <- param[grep("betaL", names(param))]
  
  sigma.y <- sqrt(1/param[grep("tau", names(param))])
  
  gammas <- param[grep("gamma", names(param))]
  alpha <- t(matrix(param[grep("alpha", names(param))],2,2)) # row = risk
  
  gambh <- matrix(param[grep("gambh", names(param))],12,2)
  mean.Cens <-  mean(pass.id$time.cmp2[pass.id$status.cmp == 0]) # 5.114849
  
  V <- D_c3
  
  V <- nearPD(V)$mat
  D <- V
  
  ################################################

  tm <- t.max %/% 0.25
  times <- c(replicate(n, 
                       c(0, seq(0.25, tm * 0.25, 0.25) + rtruncnorm(tm, -0.125, 0.125, 0, 0.036), 
                         t.max)))
  density <- rnorm(n, mean(pass.id$density), sd(pass.id$density))
  age <- rnorm(n, mean(pass.id$DxAge), sd(pass.id$DxAge))
  
  K <- length(times)/n
  DF <- data.frame(TimeSince_Dx = times, density = rep(density, each = K), DxAge = rep(age, each = K))
  X <- model.matrix(~ ns(TimeSince_Dx, knots = knot.longi[2:3], B = knot.longi[c(1,4)]) + DxAge  , data = DF)
  Z <- model.matrix(~ ns(TimeSince_Dx, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), data = DF)
  
  
  ################################################
  
  b <- mvrnorm(n, rep(0, nrow(D)), D)
  
  id <- rep(1:n, each = K)
  eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, 1:4])) # linear predictor
  
  y <- rgt(n * K, mu =  eta.y, sigma = sigma.y, df =  3)

  eta.t <- cbind(as.vector(density * gammas[1]),
                 as.vector(density * gammas[2]))
  
  invS <- function (t, u, i) {
    h <- function (s, k) {
      age.i <- age[i]
      XX <- cbind(1, ns(s, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), age.i)
      XX.back <- cbind(1, ns(s - 1, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), age.i)
      ZZ <- cbind(1, ns(s, knots = knot.longi[2:3], B = knot.longi[c(1,4)]))
      ZZ.back <- cbind(1, ns(s - 1, knots = knot.longi[2:3], B = knot.longi[c(1,4)]))
      f <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), 1:4]))
      f.back <- as.vector(XX.back %*% betas + rowSums(ZZ.back * b[rep(i, nrow(ZZ)), 1:4]))
      bh <- splineDesign(knots = knot.surv, s, ord = 4L, outer.ok = T)
      
      exp(bh %*% gambh[,k]  + eta.t[i,k] + f * alpha[k,1] + (f - f.back) * alpha[k,2])
      
    }
    integrate(h, lower = 0, upper = t, k = 1)$value + integrate(h, lower = 0, upper = t, k = 2)$value + log(u)
  }
  u <- runif(n)
  trueTimes <- matrix(NA, n)
  
  for (i in 1:n) {
    Up <- 50
    tries <- 5
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
      tries <- tries - 1
      Up <- Up + 200
      Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    }
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
  }
  trueTimes[is.na(trueTimes)] <- 1e06
  
  # simulate which event is
  # calculate the instantaneous hazard
  h <- function (s, k ,i) {
    age.i <- age[i]
    XX <- cbind(1, ns(s, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), age.i)
    XX.back <- cbind(1, ns(s - 1, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), age.i)
    ZZ <- cbind(1, ns(s, knots = knot.longi[2:3], B = knot.longi[c(1,4)]))
    ZZ.back <- cbind(1, ns(s - 1, knots = knot.longi[2:3], B = knot.longi[c(1,4)]))
    f <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), 1:4]))
    f.back <- as.vector(XX.back %*% betas + rowSums(ZZ.back * b[rep(i, nrow(ZZ)), 1:4]))
    bh <- splineDesign(knots = knot.surv, s, ord = 4L, outer.ok = T)
    
    return(ifelse(exp(bh %*% gambh[,k]  + eta.t[i,k] + f * alpha[k,1] + (f - f.back) * alpha[k,2]) < 1e-16,
                  1e-16,
                  exp(bh %*% gambh[,k]  + eta.t[i,k] + f * alpha[k,1] + (f - f.back) * alpha[k,2])))
    
  }
  
  event.2 <- sapply(1:n, function(i) {
    rbinom(1, 1, h(trueTimes[i], 2, i)/(h(trueTimes[i], 2, i) + h(trueTimes[i], 1, i))) + 1  # in rbinom 1 = event 2
  })
  
  # simulate censoring times from an exponential distribution,
  # and calculate the observed event times, i.e., min(true event times, censoring times)
  Ctimes <- runif(n, 0, 2 * mean.Cens)
  Time <- pmin(Ctimes, trueTimes)
  
  event <- sapply(1:n, function(i) {ifelse(trueTimes[i] < Ctimes[i], event.2[i], 0)}) # event indicator
  #table(event)
  
  fixed_visits <- matrix(times, ncol = n)
  fv_idx <- c(1,3,5, seq(9,K, 8))
  if (K %in% fv_idx) {
    fixed_visits <- fixed_visits[fv_idx,]
  } else {
    fixed_visits <- fixed_visits[c(fv_idx, K),]
  }
  
  Time1 <- sapply(1:n, function(i) {
    max(fixed_visits[fixed_visits[,i] <= Time[i],i])
  })
  
  Time2 <- sapply(1:n, function(i) {
    ifelse(event[i] == 1, min(fixed_visits[fixed_visits[,i] >= Time[i], i]), Time[i])
  })
  
  # sum(Time2 < Time1)
  # sum(Time1 <0 & Time2 <0)
  
  times.mat <- matrix(times, ncol = n)
  ind <- sapply(1:n, function(i) {
    if(event[i] != 1) {
      times.mat[,i] <= rep(Time[i], K)
    } else {
      if (i <= 300) {
        times.mat[,i] <= rep(Time2[i], K)
      } else {
        rep(T, K)
      }
    }
  })
  
  y.obs <- y[ind]

  id.obs <- id[ind]
  id.obs <- match(id.obs, unique(id.obs))

  dat <- DF[ind, ]
  dat$CISNET_ID <- id.obs
  dat$PSAValue <- y.obs

  dat$time.cmp1 <- Time1[id.obs]
  dat$time.cmp2 <- Time2[id.obs]
  dat$time <- Time[id.obs]
  dat$status.cmp <- event[id.obs]
  dat$b1 <- b[id.obs,1]
  dat$b2 <- b[id.obs,2]
  dat$b3 <- b[id.obs,3]
  dat$b4 <- b[id.obs,4]

  dat <- dat[c("CISNET_ID", "TimeSince_Dx", "PSAValue", "time", "time.cmp1",
               "time.cmp2","status.cmp", "density","DxAge",
               "b1", "b2", "b3", "b4")]

  dat.id <- dat[!duplicated(dat$CISNET_ID), ]
  
  # Generate 1000 subjects, the first 300 in training set,
  # split the data
  train.dat <- dat[dat$CISNET_ID %in% 1:300,]
  train.dat.id <- dat.id[dat.id$CISNET_ID %in% 1:300,]
  save(train.dat, file = paste0("Output/Simulation/Simulated datasets/Training sets/trainingset_", datanum,".RData"))
  save(train.dat.id, file = paste0("Output/Simulation/Simulated datasets/Training sets/trainingset_id_", datanum,".RData"))
  
  
  store <- dat[dat$CISNET_ID %in% 301:1000 & dat$status.cmp == 1, ][seq(41,41 + 51 * 150, 51),]
  test.idx <- c(store[store$time < store$TimeSince_Dx,]$CISNET_ID[1:100],
                dat.id[dat.id$CISNET_ID %in% 301:1000 &
                         dat.id$status.cmp == 2,]$CISNET_ID[1:30],
                dat.id[dat.id$CISNET_ID %in% 301:1000 &
                         dat.id$status.cmp == 0,]$CISNET_ID[1:70])
  
  test.dat <- dat[dat$CISNET_ID %in% test.idx,]
  test.dat.id <- dat.id[dat.id$CISNET_ID %in% test.idx,]
  save(test.dat, file = paste0("Output/Simulation/Simulated datasets/Test sets/testset_", datanum,".RData"))
  save(test.dat.id, file = paste0("Output/Simulation/Simulated datasets/Test sets/testset_id_", datanum,".RData"))
  
  # Complete data
  dat.cmpl <- DF
  dat.cmpl$CISNET_ID <- id
  dat.cmpl$PSAValue <- y
  
  dat.cmpl$time.cmp1 <- Time1[id]
  dat.cmpl$time.cmp2 <- Time2[id]
  dat.cmpl$time <- Time[id]
  dat.cmpl$status.cmp <- event[id]
  dat.cmpl$b1 <- b[id,1]
  dat.cmpl$b2 <- b[id,2]
  dat.cmpl$b3 <- b[id,3]
  dat.cmpl$b4 <- b[id,4]
  dat.cmpl <- dat.cmpl[c("CISNET_ID", "TimeSince_Dx", "PSAValue", "time", "time.cmp1", 
               "time.cmp2","status.cmp", "density","DxAge",
               "b1", "b2", "b3", "b4")]
  
  test.dat.cmpl <- dat.cmpl[dat.cmpl$CISNET_ID %in% test.idx,]
  save(test.dat.cmpl, file = paste0("Output/Simulation/Simulated datasets/Test sets/testset_complete_", datanum,".RData"))
  print(datanum)
  })
})
res <- lapply(out, future::value)
t2 <- Sys.time()
t2-t1


