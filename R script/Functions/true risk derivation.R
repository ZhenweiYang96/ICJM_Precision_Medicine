get_risk <- function(object, Survtimes, newdata, idVar = "CISNET_ID") {
  library(splines)
  test <- newdata
  test.id <- test[!duplicated(test[[idVar]]),]
  
  # Storage of the risks
  step <- Survtimes
  risk <- data.frame(t.hor = step,
                     true_risk = rep(NA, length(step)))
  risk.store <- numeric(length(step))
  
  # Extract the true parameters
  mcmc <- do.call("rbind", lapply(object$mcmc,as.matrix))
  parameters <- colMeans(mcmc)
  beta <- parameters[grep("beta", names(parameters))]
  b <- as.numeric(newdata[1,grep("b", colnames(newdata))])
  gam_bh <- matrix(parameters[grep("gambh", names(parameters))],ncol = 2)
  gamma <- matrix(parameters[grep("gamma", names(parameters))],ncol = 2)
  alpha <- matrix(parameters[grep("alpha", names(parameters))],ncol = 2) # col = risk
  knot.longi <- object$model_info$knots$knot.longi
  knot.surv <- object$model_info$knots$knot.surv
  survtimeVar <- object$model_info$var_names$survtime
  longitimeVar <- object$model_info$var_names$longitime
  longi.bsVar <- object$model_info$var_names$longi.bsVar
  surv.bsVar <- object$model_info$var_names$surv.bsVar
  X <- test.id[[surv.bsVar]]
  
  u <- max(Survtimes)
  t <- min(Survtimes)
  
  # Prepare for the calculation
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
  
  
  # Start calculation
  overall.surv <- numeric(length(step))
  
  # denominator - survival part
  t.start <- step[1]
  XL.hs.now <- matrix(1,15, length(beta))
  ZL.hs.now <- matrix(1,15, length(b))
  XL.hs.now[,2:length(b)] <- ZL.hs.now[,2:length(b)] <- t(sapply(1:15, function(n) {
    ns(f.org[n] * t.start/2 + t.start/2, 
       knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,length(b))])
  }))
  XL.hs.now[,(length(b)+1):length(beta)] <- test.id[[longi.bsVar]]
  
  XL.hs.back <- matrix(1, 15, length(beta))
  ZL.hs.back <- matrix(1, 15, length(b))
  XL.hs.back[,2:length(b)] <- ZL.hs.back[,2:length(b)] <- t(sapply(1:15, function(n) {
    ns(f.org[n] * t.start/2 + t.start/2 - 1, 
       knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,length(b))])
  }))
  XL.hs.back[,(length(b)+1):length(beta)] <- test.id[[longi.bsVar]]
  
  T.hs <- splineDesign(knot.surv, (f.org + 1) * t/2, ord = 4L)
  
  hazard.survival <- exp(T.hs[1:15,] %*% gam_bh[,1:2] + rep(gamma[,1:2] * X, each = 15) +
                                         (XL.hs.now[1:15,] %*% beta + ZL.hs.now[1:15,] %*% b) %*% alpha[1,1:2]  + 
                                         (XL.hs.now[1:15,] %*% beta + ZL.hs.now[1:15,] %*% b - 
                                            XL.hs.back[1:15,] %*% beta - ZL.hs.back[1:15,] %*% b) %*% alpha[2,1:2])
  surv <-  exp(- t.start/2 * (rep(1,2) %*% (t(hazard.survival[1:15,1:2]) %*% w)))
  surv[surv < 1e-16] <- 1e-16
  
  for (i in 1:length(step)) {
    
    u <- step[i]
    t <- ifelse(i == 1, step[i], step[i-1])
    t.step <- f.org * (u-t)/2 + (u+t)/2
    T.h <- splineDesign(knot.surv, x = f.org * (u-t)/2 + (u+t)/2, ord = 4L)
    
    XL.h.now <- matrix(1, 15, length(beta))
    ZL.h.now <- matrix(1, 15, length(b))
    XL.h.now[,2:length(b)] <- ZL.h.now[,2:length(b)] <- ns(f.org * (u-t)/2 + (u+t)/2, 
                                                           knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,length(b))])
    XL.h.now[,(length(b) + 1):length(beta)] <- test.id[[longi.bsVar]]
    
    XL.h.back <- matrix(1, 15, length(beta))
    ZL.h.back <- matrix(1, 15, length(b))
    XL.h.back[,2:length(b)] <- ZL.h.back[,2:length(b)] <- ns(f.org * (u-t)/2 + (u+t)/2 - 1, 
                                                                         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,length(b))])
    XL.h.back[,(length(b) + 1):length(beta)] <- test.id[[longi.bsVar]]
    
    T.s <- array(1, dim = c(15, nrow(gam_bh), 15)) # second 15 is the second integral
    XL.s.now <- array(1, dim = c(15, length(beta), 15))
    ZL.s.now <- array(1, dim = c(15, length(b), 15))
    XL.s.now[,(length(b) + 1):length(beta),] <- test.id[[longi.bsVar]]
    
    for (j in 1:15) {
      XL.s.now[j,2:length(b),] <- ZL.s.now[j,2:length(b),] <- t(ns(f.org * t.step[j]/2 + t.step[j]/2,
                                                                               knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,length(b))]))
    }
    
    XL.s.back <- array(1, dim = c(15, length(beta), 15))
    ZL.s.back <- array(1, dim = c(15, length(b), 15))
    XL.s.back[,(length(b) + 1):length(beta),] <- test.id[[longi.bsVar]]
    
    for (j in 1:15) {
      XL.s.back[j,2:length(b),] <- ZL.s.back[j,2:length(b),] <- t(ns(f.org * t.step[j]/2 + t.step[j]/2 - 1,
                                                                                 knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,length(b))]))
    }
    
    T.s.prep <- splineDesign(knot.surv, x = outer(f.org + 1, t.step/2), ord = 4L) # 1:15 - step[1]:f.org[1:15]
    gk <- split(1:225, 1:15)  # f.org group
    for (j in 1:15) {
      T.s[,,j] <- T.s.prep[gk[[j]],]
    }
    
    # In addition (overall survival till different time points)
    XL.os.now <- matrix(1,15, length(beta))
    ZL.os.now <- matrix(1,15, length(b))
    XL.os.now[,2:length(b)] <- ZL.os.now[,2:length(b)] <- t(sapply(1:15, function(n) {
      ns(f.org[n] * u/2 + u/2, 
         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,length(b))])
    }))
    XL.os.now[,(length(b) + 1):length(beta)] <- test.id[[longi.bsVar]]
    
    XL.os.back <- matrix(1,15, length(beta))
    ZL.os.back <- matrix(1,15, length(b))
    XL.os.back[,2:length(b)] <- ZL.os.back[,2:length(b)] <- t(sapply(1:15, function(n) {
      ns(f.org[n] * u/2 + u/2 - 1, 
         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,length(b))])
    }))
    XL.os.back[,(length(b) + 1):length(beta)] <- test.id[[longi.bsVar]]
    
    
    T.os <- splineDesign(knot.surv, (f.org + 1) * u/2, ord = 4L)
    
    
    survival <- numeric(15)
    hazard <- numeric(15)
    hazard.overall.survival <- array(NA, dim = c(15, 2))
    
    for (m in 1:15) {
      survival[m] <- exp(
        - t.step[m]/2 * sum(w %*% exp(t(T.s[m,,1:15]) %*% gam_bh[,1:2]  + rep(gamma[,1:2] * X, each = 15) + 
                                        (t(XL.s.now[m,,1:15]) %*% beta + t(ZL.s.now[m,,1:15]) %*% b) %*% alpha[1,1:2] + 
                                        (t(XL.s.now[m,,1:15]) %*% beta + t(ZL.s.now[m,,1:15]) %*% b - 
                                            t(XL.s.back[m,,1:15]) %*% beta - t(ZL.s.back[m,,1:15]) %*% b) %*% alpha[2,1:2]))
      )
      
    }
    hazard.overall.survival <- exp(T.os[,] %*% gam_bh[,1:2] + rep(gamma[,1:2] * X, each = 15) + 
                                     (XL.os.now[,] %*% beta + ZL.os.now[,] %*% b) %*% alpha[1,1:2] + 
                                            (XL.os.now[,] %*% beta + ZL.os.now[,] %*% b - 
                                               XL.os.back[,] %*% beta - ZL.os.back[,] %*% b) %*% alpha[2,1:2])
    hazard[1:15] <- exp(T.h[1:15,] %*% gam_bh[,1] + rep(gamma[,1] * X, each =15) + 
                            alpha[1,1] * (XL.h.now[1:15,] %*% beta + ZL.h.now[1:15,] %*% b) + 
                              alpha[2,1] * (XL.h.now[1:15,] %*% beta + ZL.h.now[1:15,] %*% b - 
                                                XL.h.back[1:15,] %*% beta - ZL.h.back[1:15,] %*% b))
    risk.store[i] <- (u-t)/2 * ((hazard * survival) %*% w) /
      c(surv)
    
    overall.surv[i] <- exp(- u/2 * (hazard.overall.survival[,1] %*% w + 
                                             hazard.overall.survival[,2] %*% w))
    overall.surv[overall.surv[i] < 1e-16] <- 1e-16
    risk[i,2] <- sum(risk.store[1:i])
  }
  
  return(list(risk = risk,
              surv = surv,
              parameters = list(beta = beta,
                                b = b,
                                gam_bh = gam_bh,
                                gamma = gamma,
                                alpha = alpha),
              overall.surv = overall.surv))
}

