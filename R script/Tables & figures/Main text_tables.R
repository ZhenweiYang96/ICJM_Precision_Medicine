library(rjags)

# Table 1 -----------------------------------------------------------------
load("Output/Data analysis/ICJM1.RData")
load("Output/Data analysis/ICJM2.RData")

### ICJM 1
mcmc1 <- do.call(rbind, lapply(ICJM1$mcmc, as.matrix))
# Gamma
round(exp(ICJM1$summary$coef[grepl("gamma",names(ICJM1$summary$coef))]),2)
round(apply(exp(mcmc1[,grepl("gamma",names(ICJM1$summary$coef))]), 
            2, function(x) {quantile(x, 0.025)}),2)
round(apply(exp(mcmc1[,grepl("gamma",names(ICJM1$summary$coef))]), 
            2, function(x) {quantile(x, 0.975)}),2)

# Alpha [row=forms, col=event type]
round(exp(ICJM1$summary$coef[grepl("alpha",names(ICJM1$summary$coef))]),2)
round(apply(exp(mcmc1[,grepl("alpha",names(ICJM1$summary$coef))]), 
            2, function(x) {quantile(x, 0.025)}),2)
round(apply(exp(mcmc1[,grepl("alpha",names(ICJM1$summary$coef))]), 
            2, function(x) {quantile(x, 0.975)}),2)


### ICJM 2
mcmc2 <- do.call(rbind, lapply(ICJM2$mcmc, as.matrix))
# Gamma
round(exp(ICJM2$summary$coef[grepl("gamma",names(ICJM2$summary$coef))]),2)
round(apply(exp(mcmc2[,grepl("gamma",names(ICJM2$summary$coef))]), 
            2, function(x) {quantile(x, 0.025)}),2)
round(apply(exp(mcmc2[,grepl("gamma",names(ICJM2$summary$coef))]), 
            2, function(x) {quantile(x, 0.975)}),2)

# Alpha [row=forms, col=event type]
round(exp(ICJM2$summary$coef[grepl("alpha",names(ICJM2$summary$coef))]),2)
round(apply(exp(mcmc2[,grepl("alpha",names(ICJM2$summary$coef))]), 
            2, function(x) {quantile(x, 0.025)}),2)
round(apply(exp(mcmc2[,grepl("alpha",names(ICJM2$summary$coef))]), 
            2, function(x) {quantile(x, 0.975)}),2)
