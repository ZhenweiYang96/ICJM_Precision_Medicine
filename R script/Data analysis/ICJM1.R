library(future)
source("R script/Functions/ICCSJoint model function JAGS.R")
load("Cleaned data/pass.RData")

plan(multisession, workers = 3)

t0 <- Sys.time()
ICJM1 <- iccsjm(data = pass, n.adapt = 5000, 
                n.burnin = 5000, n.iter = 10000, 
                qp = 15, kntbh_spec = "quantile",
                seed = 2022)
t1 <- Sys.time()
t1 - t0

save(ICJM1,
     file = "work 2021 second half/Joint model final/ICJM1.RData")
ICJM1$summary$coef
mcmcplots::mcmcplot(ICJM1$mcmc)
