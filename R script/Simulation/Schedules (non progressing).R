library(future)
source("R script/Functions/screening schedule planning function.R")
source("R script/Functions/dynpred_Risk prediction function.R")
source("R script/Functions/function.R")
source("R script/Functions/Final personalized schedule evaluation function.R")

set.seed(4321)
schedule_seed_nonpro <- matrix(sample(1:1000000, 1000 * 100), 1000, 100)
save(schedule_seed_nonpro, file = "R script/Simulation/seed/schedule generating seed (non progressing patients).RData")

plan(multisession, workers = 15)

# Looping around all the test set [non-progression]-----------------------------------------
for (i in 1:200) {
  load(paste0("Output/Simulation/Simulated datasets/Test sets/testset_", i, ".RData"))
  load(paste0("Output/Simulation/Model results/Joint model_", i, ".RData"))
  test.dat <- subset(test.dat, status.cmp != 1)
  id <- unique(test.dat$CISNET_ID)
  out <- lapply(1:100, function(j) {
    future({
      # Windows cannot load packages in nodes
      library(rjags)
      
      testset <- test.dat[test.dat$CISNET_ID == id[j] & test.dat$TimeSince_Dx <= 10.125,]
      
      fv <- seq(0.5,10, 0.5)
      len <- length(testset$TimeSince_Dx)
      if (len < 3) {
        visits <- fv
      } else {
        if (length(seq(3, len,2)) < length(fv)) {
          visits <- c(testset$TimeSince_Dx[seq(3,len, 2)],
                      fv[(length(testset$TimeSince_Dx[seq(3,len, 2)]) + 1):length(fv)])
        } else {
          visits <- testset$TimeSince_Dx[seq(3,len, 2)]
        }
      }
      

      flexible_schedule <- personalizedSchedule.iccsjm_evaluation(object = iccsjm.model, newdata = testset,
                                                                  last_test_time = 0,
                                                                  fixed_grid_visits = visits, 
                                                                  detection_delay_limit = 1.5,
                                                                  risk_range = c(seq(0.003, 0.3, 0.003), seq(0.31,1,0.01)),
                                                                  seed = schedule_seed_nonpro[i,j])
      save(flexible_schedule, file = paste0("Output/Simulation/Personalized schedule result [non progressing]/detail/schedule_dataset_",
                                            i,
                                            "_subject_", 
                                            j, 
                                            "_seed_", 
                                            schedule_seed_nonpro[i,j], 
                                            ".RData"))
      num_biopsy_flexible <- sum(flexible_schedule < testset$time[1]) + 1
      
      fixed_schedule_1 <- c(1,seq(2,10,2))
      fixed_schedule_2 <- seq(1,10,1)
      num_biopsy_fixed_1 <- sum(fixed_schedule_1 < testset$time[1]) + 1
      num_biopsy_fixed_2 <- sum(fixed_schedule_2 < testset$time[1]) + 1
      
      return(list(num_biopsy_flexible = num_biopsy_flexible,
                  num_biopsy_fixed_1 = num_biopsy_fixed_1,
                  num_biopsy_fixed_2 = num_biopsy_fixed_2))
    })
  })
  result <- lapply(out, future::value)
  save(result, file = paste0("Output/Simulation/Personalized schedule result [non progressing]/indicators_", i, ".RData"))
}

