# load parallel packages -------------------------------------
library(future)
source("R script/Functions/ICCSJoint model function JAGS.R")
load("R script/Simulation/seed/seed_data.RData")

plan(
  list(
    tweak(multisession, workers = 5),
    tweak(multisession, workers = 3)
  )
)

t1 <- Sys.time()
out <- lapply(1:100, function(datanum) {
  future({
    
    # load data file
    load(paste0("Output/Simulation/Simulated datasets/Training sets/trainingset_",
                datanum, ".RData"))
    
    # start modelling
    iccsjm.model <- iccsjm(data = train.dat, n.adapt = 3000, 
                     n.burnin = 3000, n.iter = 10000, seed = random.seed[datanum])
    save(iccsjm.model,
         file = paste0("Output/Simulation/Model results/Joint model_",
                       datanum, ".RData"))
    # 
    plan(
      list(
        tweak(multisession, workers = 5),
        tweak(multisession, workers = 3)
      )
    )
    return(paste0(datanum,"done!"))
  })
})
res <- lapply(out, future::value)
t2 <- Sys.time()

t2 - t1

