library(future)
plan(multisession, workers = 15)

set.seed(1234)
schedule_seed_pro <- matrix(sample(1:1000000, 1000 * 100), 1000, 100)
save(schedule_seed_pro, file = "R script/Simulation/seed/schedule generating seed (progressing patients).RData")

source("R script/Functions/screening schedule planning function.R")
source("R script/Functions/[notrt]dynpred_Risk prediction function.R")
source("R script/Functions/function.R")
source("R script/Functions/Final personalized schedule evaluation function.R")

# Looping around all the test set -----------------------------------------
for (i in 1:200) {
  load(paste0("Output/Simulation/Simulated datasets/Test sets/testset_", i, ".RData"))
  load(paste0("Output/Simulation/Model results/Joint model_", i, ".RData"))
  test.dat <- subset(test.dat, status.cmp == 1)
  id <- unique(test.dat$CISNET_ID)
  out <- lapply(1:100, function(j) {
    future({
      # Windows cannot load packages in nodes
      library(rjags)
      
      testset <- test.dat[test.dat$CISNET_ID == id[j] & test.dat$TimeSince_Dx <= 10.125,]
      len <- length(testset$TimeSince_Dx)
      visits <- testset$TimeSince_Dx[seq(3,len, 2)]

      flexible_schedule <- personalizedSchedule.iccsjm_evaluation(object = iccsjm.model, newdata = testset,
                                                                  last_test_time = 0,
                                                                  fixed_grid_visits = visits, 
                                                                  detection_delay_limit = 1.5,
                                                                  risk_range = c(seq(0.003, 0.3, 0.003), seq(0.31,1,0.01)),
                                                                  seed = schedule_seed_pro[i,j])
      save(flexible_schedule, file = paste0("Output/Simulation/Personalized schedule result [progressing]/detail/schedule_dataset_",
                                            i,
                                            "_subject_", 
                                            j, 
                                            "_seed_", 
                                            schedule_seed_pro[i,j], 
                                            ".RData"))

      num_biopsy_flexible <- sum(flexible_schedule < testset$time[1]) + 2
      treatment_delay_flexible <- flexible_schedule[num_biopsy_flexible - 1] - testset$time[1]
      
      fixed_schedule_1 <- c(1,seq(2,10,2))
      fixed_schedule_2 <- seq(1,10,1)
      num_biopsy_fixed_1 <- sum(fixed_schedule_1 < testset$time[1]) + 2
      num_biopsy_fixed_2 <- sum(fixed_schedule_2 < testset$time[1]) + 2
      treatment_delay_fixed_1 <- fixed_schedule_1[num_biopsy_fixed_1 - 1] - testset$time[1]
      treatment_delay_fixed_2 <- fixed_schedule_2[num_biopsy_fixed_2 - 1] - testset$time[1]
      
      return(list(num_biopsy_flexible = num_biopsy_flexible,
                  num_biopsy_fixed_1 = num_biopsy_fixed_1,
                  num_biopsy_fixed_2 = num_biopsy_fixed_2,
                  treatment_delay_flexible = treatment_delay_flexible,
                  treatment_delay_fixed_1 = treatment_delay_fixed_1,
                  treatment_delay_fixed_2 = treatment_delay_fixed_2))
    })
  })
  result <- lapply(out, future::value)
  save(result, file = paste0("Output/Simulation/Personalized schedule result [progressing]/indicators_", i, ".RData"))
}

