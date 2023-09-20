library(future)
source("R script/Functions/ICCSJoint model function JAGS.R")
load("Cleaned data/pass.RData")

plan(multisession, workers = 3)

iccsjm_obs <- iccsjm(data = pass, n.adapt = 3000, 
                     n.burnin = 3000, n.iter = 10000, seed = 2022)
save(iccsjm_obs,
     file = "Output/Observed data results/ICJM1.RData")