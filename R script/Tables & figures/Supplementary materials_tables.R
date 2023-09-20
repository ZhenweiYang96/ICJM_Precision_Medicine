

# Table S2 ----------------------------------------------------------------

store <- data.frame(i = 1:200, 
                    censoring = 0,
                    progression = 0, 
                    treamtent = 0)
for (i in 1:200) {
  load(paste0("Output/Simulation/Simulated datasets/Training sets/trainingset_id_", 
              i,".RData"))
  store[i,2:4] <- prop.table(table(train.dat.id$status.cmp))
}

load("Cleaned data/pass_id.RData")
prop.table(table(pass.id$status.cmp))

tab.s2 <- data.frame(events = c("Cancer progression", "Treatment", "Censoring"),
                     `Simulated data` = round(colMeans(store)[c(3,4,2)]*100, 2),
                     `Observed data` = round(as.vector(prop.table(table(pass.id$status.cmp)))[c(2,3,1)]* 100, 2))
tab.s2


# Table S3 ----------------------------------------------------------------

load("Output/Data analysis/ICJM1.RData")
load("Output/Data analysis/ICJM2.RData")
library(rjags)
mcmc.icjm1 <- do.call(rbind, lapply(1:3, function(i) {as.matrix(ICJM1$mcmc[[i]])}))
mcmc.icjm1 <- mcmc.icjm1[,grep("betaL|gamma|alpha", colnames(mcmc.icjm1))]
colnames(mcmc.icjm1)
summary(mcmc.icjm1)
mcmc.icjm1 <- mcmc.icjm1[,c("betaL[1]", "betaL[2]", "betaL[3]",
                            "betaL[4]", "betaL[5]", "gamma[1]",
                            "alpha[1,1]", "alpha[2,1]", "gamma[2]",
                            "alpha[1,2]", "alpha[2,2]")]
tab.icjm1 <- data.frame(parameters = c("intercept_PSA", "Time 1_PSA", "Time 2_PSA",
                                       "Time 3_PSA", "Age", "log(PSA density)_pro",
                                       "log2(PSA + 1) value_pro", "log2(PSA + 1) yearly change_pro",
                                       "log(PSA density)_trt",
                                       "log2(PSA + 1) value_trt", "log2(PSA + 1) yearly change_trt"),
                        estimates = round(colMeans(mcmc.icjm1), 2),
                        CI = sapply(seq_len(ncol(mcmc.icjm1)), function(i) {
                          paste0("[", round(quantile(mcmc.icjm1[,i], probs = 0.025),2), 
                                ", ", round(quantile(mcmc.icjm1[,i], probs = 0.975),2), "]")
                        }))
tab.icjm1


mcmc.icjm2 <- do.call(rbind, lapply(1:3, function(i) {as.matrix(ICJM2$mcmc[[i]])}))
mcmc.icjm2 <- mcmc.icjm2[,grep("betaL|gamma|alpha", colnames(mcmc.icjm2))]
colnames(mcmc.icjm2)
summary(mcmc.icjm2)
mcmc.icjm2 <- mcmc.icjm2[,c("betaL_psa[1]", "betaL_psa[2]", "betaL_psa[3]",
                            "betaL_psa[4]", "betaL_psa[5]", 
                            "betaL_cr[1]", "betaL_cr[2]", "betaL_cr[3]",
                            "gamma[1]", "alpha[1,1]", "alpha[2,1]", "alpha[3,1]",
                            "gamma[2]", "alpha[1,2]", "alpha[2,2]", "alpha[3,2]")]
tab.icjm2 <- data.frame(parameters = c("intercept_PSA", "Time 1_PSA", "Time 2_PSA",
                                       "Time 3_PSA", "Age", 
                                       "intercept_cr", "Time 1_cr", "Time 2_cr",
                                       "log(PSA density)_pro",
                                       "log2(PSA + 1) value_pro", "log2(PSA + 1) yearly change_pro",
                                       "logit(cr)_pro",
                                       "log(PSA density)_trt",
                                       "log2(PSA + 1) value_trt", "log2(PSA + 1) yearly change_trt",
                                       "logit(cr)_trt"),
                        estimates = round(colMeans(mcmc.icjm2), 2),
                        CI = sapply(seq_len(ncol(mcmc.icjm2)), function(i) {
                          paste0("[", round(quantile(mcmc.icjm2[,i], probs = 0.025),2), 
                                 ", ", round(quantile(mcmc.icjm2[,i], probs = 0.975),2), "]")
                        }))
tab.icjm2
