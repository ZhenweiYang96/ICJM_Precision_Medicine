library(future)
plan(multisession, workers = 15)

load("Output/Data analysis/ICJM1.RData")
source("R script/Functions/dynpred_Risk prediction function.R")
source("R script/Functions/function.R")
source("R script/Functions/true risk derivation.R")

set.seed(4567)
seed_modeltest <- array(sample(501:1000000000, 1000 * 200 * 6), dim = c(1000,200,6))
save(seed_modeltest, file = "R script/Simulation/seed/seed_modeltest.RData")

# 11hrs - 100 
results_risk <- list()
for (i in 1:200) {
  load(paste0("Output/Simulation/Simulated datasets/Test sets/testset_complete_",i,".RData"))
  load(paste0("Output/Simulation/Model results/Joint model_",i,".RData"))
  ids <- unique(test.dat.cmpl$CISNET_ID)
  out <- lapply(1:200, function(j) {
    future({
      library(rjags)
      risk_pred <- data.frame(start_time = c(0,1,2,3,4,6),
                              true_risk = 0,
                              est_risk = 0)
      st.vec <- c(0,1,2,3,4,6)
      for (m in seq_len(length(st.vec))) {
        nd = test.dat.cmpl[test.dat.cmpl$CISNET_ID == ids[j] & test.dat.cmpl$TimeSince_Dx <= st.vec[m],]
        true_risk <- get_risk(ICJM1, Survtimes = c(st.vec[m], st.vec[m]+2), newdata = nd)
        est_risk <- csdypred(object = iccsjm.model, iter = 150, n.adapt = 200, 
                             Survtimes = seq(st.vec[m], st.vec[m]+2, length.out=  100),
                             newdata = nd, seed = seed_modeltest[i,j,m])
        risk_pred[m,2] <- true_risk$risk[2,2]
        risk_pred[m,3] <- est_risk$risk[nrow(est_risk$risk),2]
      }
      save(risk_pred, file = paste0("Output/Simulation/Model test/detail/risk comparison_dataset",
                                    i,
                                    "_subject",
                                    j,
                                    "_seed",
                                    seed_modeltest[i,j,m],
                                    ".RData"))
      risk_pred
    })
  })
  result <- lapply(out, future::value)
  save(result, file = paste0("Output/Simulation/Model test/risk comparison_dataset_",i,".RData"))
  results_risk[[i]] <- result
  print(paste(i, "done"))
}
save(results_risk, file = "Output/Simulation/Model test/risk comparison_full.RData")


# Plots -------------------------------------------------------------------
# for (i in 1:200) {
#   for (j in 1:200) {
#     results_risk[[i]][[j]]$MSE = (results_risk[[i]][[j]]$true_risk - 
#                                             results_risk[[i]][[j]]$est_risk)^2
#     results_risk[[i]][[j]]$absE = abs(results_risk[[i]][[j]]$true_risk - 
#                                     results_risk[[i]][[j]]$est_risk)
#     results_risk[[i]][[j]]$E = results_risk[[i]][[j]]$true_risk - 
#                                     results_risk[[i]][[j]]$est_risk
#   }
# }
# 
# result.pool <- list()
# for(i in 1:200){
#   result.pool[[i]] <- do.call(rbind,results_risk[[i]])
# }
# comparison <- do.call(rbind, result.pool)
# comparison$dataset <- rep(1:200, each = 200 * 6)
# comparison$id <- rep(rep(1:200, each = 6), 200)
# comparison$type <- "quantile knots"
# save(comparison, file = "Output/Simulation/Model test/comparison.RData")

# library(ggplot2)
# library(tidyverse)
# boxplot.e <- comparison %>% 
#   #gather(year, MSE, `0`:`6`) %>% 
#   ggplot(aes(x = as.factor(start_time), y = E)) + 
#   stat_boxplot(geom = "errorbar", width = 0.5) +  
#   geom_boxplot() +
#   theme_bw() +
#   xlab("Start year of prediction") +
#   ylab("Error (True - Estimated risk)") + 
#   ggtitle("Precision of 2-year risk prediction with different start years")
# boxplot.e
# ggsave(boxplot.e, file = "work 2021 second half/Simulation/Real simulation/plots/boxplot_Error.png")
# 
# boxplot.absE <- comparison %>% 
#   #gather(year, MSE, `0`:`6`) %>% 
#   ggplot(aes(x = as.factor(start_time), y = absE)) + 
#   stat_boxplot(geom = "errorbar", width = 0.5) +  
#   geom_boxplot() +
#   theme_bw() +
#   xlab("Start year of prediction") +
#   ylab("Absolute error") + 
#   ggtitle("Precision of 2-year risk prediction with different start years")
# boxplot.absE
# ggsave(boxplot.absE, file = "work 2021 second half/Simulation/Real simulation/plots/boxplot_Aboslute_Error.png")
# 
# boxplot.sE <- comparison %>% 
#   #gather(year, MSE, `0`:`6`) %>% 
#   ggplot(aes(x = as.factor(start_time), y = MSE)) + 
#   stat_boxplot(geom = "errorbar", width = 0.5) +  
#   geom_boxplot() +
#   theme_bw() +
#   xlab("Start year of prediction") +
#   ylab("Squared error") + 
#   ggtitle("Precision of 2-year risk prediction with different start years")
# boxplot.sE
# ggsave(boxplot.sE, file = "work 2021 second half/Simulation/Real simulation/plots/boxplot_square_Error.png")
# 
# boxplot.mse <- comparison %>% 
#   group_by(dataset, start_time) %>% 
#   summarise(mean = mean(MSE)) %>% 
#   ggplot(aes(x = as.factor(start_time), y = mean)) + 
#   stat_boxplot(geom = "errorbar", width = 0.5) +  
#   geom_boxplot() +
#   theme_bw() +
#   xlab("Start year of prediction") +
#   ylab("Mean Square Error") + 
#   ggtitle("Precision of 2-year risk prediction with different start years")
# boxplot.mse
# ggsave(boxplot.mse, file = "work 2021 second half/Simulation/Real simulation/plots/boxplot_MSE.png")


