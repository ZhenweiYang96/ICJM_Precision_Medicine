library(ggplot2)
library(tidyverse)
library(splines)
library(rjags)
load("Cleaned data/pass_id.RData")
load("Cleaned data/pass.RData")
load("Cleaned data/pass_cores.RData")
load("R script/Data analysis/ICJM2_inits/MCMC_mvmm.RData")

# Figure S1 ---------------------------------------------------------------
set.seed(2021)
ids <- sort(sample(pass.id$CISNET_ID, 20))

s1.1 <- pass %>% 
  filter(CISNET_ID %in% ids) %>% 
  ggplot() + 
  geom_point(aes(x = TimeSince_Dx, y = PSAValue), size = 0.5) + 
  geom_line(aes(x = TimeSince_Dx, y = PSAValue)) + 
  facet_wrap(~ CISNET_ID, nrow = 4, ncol = 5) + 
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", size = 1, fill= NA),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=16)) + 
  ylab(expression(log[2](PSA + 1))) + 
  xlab("Time (years)")
s1.1
ggsave(s1.1, file = "Output/Tables & figures/Supplementary/FigureS1_1_plot_PSA.png")

s1.2 <- pass.cores %>% 
  filter(CISNET_ID %in% ids) %>% 
  ggplot() + 
  geom_point(aes(x = TimeSince_Dx, y = cores_ratio), size = 0.5) + 
  geom_line(aes(x = TimeSince_Dx, y = cores_ratio)) + 
  facet_wrap(~ CISNET_ID, nrow = 4, ncol = 5) + 
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", size = 1, fill= NA),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=16)) + 
  scale_x_continuous(breaks = seq(0,9,3))+
  ylab("Core ratio (%)") + 
  xlab("Time (years)")
s1.2
ggsave(s1.2, file = "Output/Tables & figures/Supplementary/FigureS1_2_plot_coresratio.png")


# Figure S2 Model evaluation error ----------------------------------------
load("Output/Simulation/Model test/risk comparison_full.RData")
for (i in 1:200) {
  for (j in 1:200) {
    results_risk[[i]][[j]]$MSE = (results_risk[[i]][[j]]$true_risk - 
                                    results_risk[[i]][[j]]$est_risk)^2
    results_risk[[i]][[j]]$absE = abs(results_risk[[i]][[j]]$true_risk - 
                                        results_risk[[i]][[j]]$est_risk)
    results_risk[[i]][[j]]$E = results_risk[[i]][[j]]$true_risk - 
      results_risk[[i]][[j]]$est_risk
  }
}
result.pool <- list()
for(i in 1:200){
  result.pool[[i]] <- do.call(rbind,results_risk[[i]])
}
comparison <- do.call(rbind, result.pool)
boxplot.e <- comparison %>% 
  #gather(year, MSE, `0`:`6`) %>% 
  ggplot(aes(x = as.factor(start_time), y = E)) + 
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot() +
  theme_bw() +
  xlab("Start year of prediction") +
  ylab("Error (True - Estimated risk)") + 
  theme(text = element_text(size=20))
boxplot.e
ggsave(boxplot.e, file = "Output/Tables & figures/Supplementary/FigureS2_Evaluation_error_plot.png")


# Figure S3 Baseline hazard plot ------------------------------------------
source("R script/Functions/function.R")
load("Output/Data analysis/ICJM1.RData")
# Extract the estimated baseline hazard functions
est.bhzard <- data.frame(dataset = rep(1:200, each = 500),
                         time = rep(seq(0, 10, length.out = 500), 200),
                         risk = 0)
for (i in 1:200) {
  load(paste0("Output/Simulation/Model results/Joint model_", i, ".RData"))
  gam_bh <- matrix(iccsjm.model$summary$coef[grep("gambh", names(iccsjm.model$summary$coef))], 12, 2)
  dm_gambh <- splineDesign(seq(0, 10, length.out = 500),
                           knots = iccsjm.model$model_info$knots$knot.surv, outer.ok = TRUE)
  est.bhzard[est.bhzard$dataset == i,3] <- dm_gambh %*% gam_bh[,1]
}
save(est.bhzard, file="Output/Tables & figures/Supplementary/estbh.RData")

# Extract the true baseline hazard function
true.bhzard <- data.frame(dataset = rep(0, 500),
                          time = seq(0,10, length.out = 500),
                          risk = 0)
gam_bh <- matrix(ICJM1$summary$coef[grep("gambh", names(ICJM1$summary$coef))], 12, 2)
dm_gambh <- splineDesign(seq(0, 10, length.out = 500),
                         knots = ICJM1$model_info$knots$knot.surv, outer.ok = TRUE)
true.bhzard[,3] <- dm_gambh %*% gam_bh[,1]
save(true.bhzard, file="Output/Tables & figures/Supplementary/truebh.RData")

load("Output/Tables & figures/Supplementary/truebh.RData")
load("Output/Tables & figures/Supplementary/estbh.RData")

plot <- ggplot() + 
  geom_line(aes(x = est.bhzard$time, y = est.bhzard$risk, group = est.bhzard$dataset), alpha = 0.5, linewidth = 1) + 
  geom_line(aes(x = true.bhzard$time, y = true.bhzard$risk), color = "red", linewidth = 1.5) + 
  ylab("log(baseline hazard)") + 
  xlab("Time (years)") + 
  theme_bw()+ 
  theme(text = element_text(size=20))
plot

ggsave(plot, file = "Output/Tables & figures/Supplementary/FigureS3_baseline_hazard.png")


# Figure S4 fitness of the mixed model -------------------------------------
load("R script/Data analysis/ICJM2_inits/MCMC_mvmm.RData")
load("Output/Data analysis/ICJM2.RData")
n_sample <- dim(mcmc.mvmm[[1]])[1] * 3
knot.longi <- ICJM2$model_info$knots$knot.longi
ids.index <- sapply(ids, function(x) {which(pass.id$CISNET_ID == x)})
#ids.index[1]=72
para <- do.call(rbind, mcmc.mvmm)
b <- array(NA, c(833, 7, n_sample))
for (i in 1:n_sample) {
  b[,,i] <- matrix(para[i,50:5880], 833, 7)
}
para.sub <- array(NA, c(20, 15, n_sample))
dimnames(para.sub)[[1]] <- sapply(ids, function(i){as.character(i)})
dimnames(para.sub)[[2]] <- c(sapply(1:5, function(i) {paste0("betaL_psa_",i)}),
                             sapply(1:4, function(i) {paste0("b_psa_",i)}),
                             sapply(1:3, function(i) {paste0("betaL_cr_",i)}),
                             sapply(1:3, function(i) {paste0("b_cr_",i)}))
for (i in 1:20) {
  for (j in 1:5) {
    para.sub[i,j,] <- para[,grep("betaL_psa",colnames(para))[j]]
  }
  for (j in 6:9) {
    para.sub[i,j,] <- b[ids.index[i], j-5,]
  }
  for (j in 10:12) {
    para.sub[i,j,] <- para[,grep("betaL_cr",colnames(para))[j-9]]
  }
  for (j in 13:15) {
    para.sub[i,j,] <- b[ids.index[i], j-8,]
  }
}
para.sub

mvmm.fit.psa <- data.frame(
  CISNET_ID = rep(ids, each = 100),
  time = c(sapply(ids, function(x) {seq(0, max(pass$TimeSince_Dx[pass$CISNET_ID==x]) + 1, 
                                        length.out=100)})),
  age = rep(sapply(ids, function(x) {pass.id$DxAge[pass.id$CISNET_ID == x]}), each = 100),
  psafit = 0,
  psafit_lower = 0,
  psafit_upper = 0
)

mvmm.fit.cr <- data.frame(
  CISNET_ID = rep(ids, each = 100),
  time = c(sapply(ids, function(x) {seq(0, max(pass.cores$TimeSince_Dx[pass.cores$CISNET_ID==x]) + 1, 
                                        length.out=100)})),
  crfit = 0,
  crfit_lower = 0,
  crfit_upper = 0
)

for (i in 1:length(ids)) {
  Mat_psa <- model.matrix(~ ns(time, knots = knot.longi[2:3],
                               B = knot.longi[c(1,4)]) + age, data = mvmm.fit.psa[mvmm.fit.psa$CISNET_ID==ids[i],])
  Mat_cr <- model.matrix(~ poly(time, 2, raw = TRUE), data = mvmm.fit.cr[mvmm.fit.cr$CISNET_ID==ids[i],])
  lp_psa <- cbind(Mat_psa, Mat_psa[,1:4]) %*% para.sub[i,1:9,]
  lp_cr <- cbind(Mat_cr, Mat_cr) %*% para.sub[i,10:15,]
  
  
  mvmm.fit.psa[mvmm.fit.psa$CISNET_ID==ids[i], "psafit"] <- rowMeans(lp_psa)
  mvmm.fit.psa[mvmm.fit.psa$CISNET_ID==ids[i], "psafit_lower"] <- apply(lp_psa, 1, function(x) {quantile(x, 0.025)})
  mvmm.fit.psa[mvmm.fit.psa$CISNET_ID==ids[i], "psafit_upper"] <- apply(lp_psa, 1, function(x) {quantile(x, 0.975)})
  mvmm.fit.cr[mvmm.fit.cr$CISNET_ID==ids[i], "crfit"] <- rowMeans(lp_cr)
  mvmm.fit.cr[mvmm.fit.cr$CISNET_ID==ids[i], "crfit_lower"] <- apply(lp_cr, 1, function(x) {quantile(x, 0.025)})
  mvmm.fit.cr[mvmm.fit.cr$CISNET_ID==ids[i], "crfit_upper"] <- apply(lp_cr, 1, function(x) {quantile(x, 0.975)})
}
mvmm.fit.cr$crfit <- exp(mvmm.fit.cr$crfit)/(1+exp(mvmm.fit.cr$crfit))
mvmm.fit.cr$crfit_lower <- exp(mvmm.fit.cr$crfit_lower)/(1+exp(mvmm.fit.cr$crfit_lower))
mvmm.fit.cr$crfit_upper <- exp(mvmm.fit.cr$crfit_upper)/(1+exp(mvmm.fit.cr$crfit_upper))
pass.sub <- pass %>% 
  filter(CISNET_ID %in% ids)
draw_psa <- data.frame(
  ID = c(pass.sub$CISNET_ID, mvmm.fit.psa$CISNET_ID),
  time = c(pass.sub$TimeSince_Dx, mvmm.fit.psa$time),
  psa = c(pass.sub$PSAValue, mvmm.fit.psa$psafit),
  psa_lower = c(pass.sub$PSAValue, mvmm.fit.psa$psafit_lower),
  psa_upper = c(pass.sub$PSAValue, mvmm.fit.psa$psafit_upper),
  type = c(rep("observed",nrow(pass.sub)),
           rep("fit",nrow(mvmm.fit.psa)))
)
s5.1 <- mvmm.fit.psa %>% 
  ggplot() + 
  geom_point(data = pass.sub, aes(x = TimeSince_Dx, y = PSAValue), size = 0.5) + 
  geom_line(data = pass.sub, aes(x = TimeSince_Dx, y = PSAValue)) +
  geom_line(aes(x = time, y = psafit), color="red") + 
  geom_ribbon(aes(x=time, ymin=psafit_lower, ymax = psafit_upper), fill="red", alpha = 0.3) +
  facet_wrap(~ CISNET_ID, nrow = 4, ncol = 5) + 
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", size = 1, fill= NA),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=16)) +
  ylab(expression(log[2](PSA + 1))) + 
  scale_x_continuous(breaks = seq(0,7.5,2.5)) +
  xlab("Time (years)")
s5.1
ggsave("Output/Tables & figures/Supplementary/FigureS4_1_fittness_of_PSA.png",
       s5.1)


pass_core.sub <- pass.cores %>% 
  filter(CISNET_ID %in% ids)
s5.2 <- mvmm.fit.cr %>% 
  ggplot() + 
  geom_point(data = pass_core.sub, aes(x = TimeSince_Dx, y = cores_ratio), size = 0.5) + 
  geom_line(data = pass_core.sub, aes(x = TimeSince_Dx, y = cores_ratio)) +
  geom_line(aes(x = time, y = crfit*100), color="red") + 
  geom_ribbon(aes(x=time, ymin=crfit_lower*100, ymax = crfit_upper*100), fill="red", alpha = 0.3) +
  facet_wrap(~ CISNET_ID, nrow = 4, ncol = 5) + 
  ylim(c(0,100)) + 
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", size = 1, fill= NA),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=16)) + 
  ylab("Core ratio (%)") + 
  xlab("Time (years)")
s5.2
ggsave("Output/Tables & figures/Supplementary/FigureS4_2_fittness_of_core_ratio.png",
       s5.2)


# Figure S5a Expected trajectory of PSA ---------------------------------------------------------------------
mcmc <- do.call("rbind",lapply(1:3, function(x) {as.matrix(ICJM2$mcmc[[x]])}))
pars <- apply(mcmc, 2, mean)
pars.lower <- apply(mcmc, 2, function(x) {quantile(x, prob = 0.025)})
pars.upper <- apply(mcmc, 2, function(x) {quantile(x, prob = 0.975)})

knot.longi <- ICJM2$model_info$knots$knot.longi
T <- ns(seq(0,10, length.out =  100), knots = knot.longi[2:3], B = knot.longi[c(1,4)])
X.t <- cbind(1, T, 0)
Y <- 2^(X.t %*% pars[grep("beta", names(pars), value = T)][4:8]) - 1
data <- data.frame(time = seq(0,10, length.out =  100),
                   PSA = Y,
                   PSA.lower = 2^(X.t %*% pars.lower[grep("beta", names(pars), value = T)][4:8]) - 1,
                   PSA.upper = 2^(X.t %*% pars.upper[grep("beta", names(pars), value = T)][4:8]) - 1)

plot <- ggplot(data) +
  geom_line(aes(x = time,
                y = PSA), size = 1.1)+
  geom_ribbon(aes(ymax = PSA.upper, ymin = PSA.lower, x = time), fill = "gray70", alpha = 0.5) +
  xlab("Time (year)") +
  ylim(c(2.5,10)) +
  ylab("Expected PSA level (ng/ml)") + 
  theme_bw()+
  theme(text = element_text(size=20))
plot
ggsave(plot, filename = "Output/Tables & figures/Supplementary/FigureS5a_Effect_plot_PSA.png")


# Figure S5b Expected trajectory of core ratio ------------------------------
T.cr <- poly(seq(0,10, length.out =  100), degree = 2, raw = TRUE)
X.t.cr <- cbind(1, T.cr)
Y.cr <- X.t.cr %*% pars[grep("beta", names(pars), value = T)][1:3]
data.cr <- data.frame(time = seq(0,10, length.out =  100),
                      cr = exp(Y.cr)/(1 + exp(Y.cr)) * 100,
                      cr.lower = exp(X.t.cr %*% pars.lower[grep("beta", names(pars), value = T)][1:3])/
                        (1 + exp(X.t.cr %*% pars.lower[grep("beta", names(pars), value = T)][1:3])) * 100,
                      cr.upper = exp(X.t.cr %*% pars.upper[grep("beta", names(pars), value = T)][1:3])/
                        (1 + exp(X.t.cr %*% pars.upper[grep("beta", names(pars), value = T)][1:3])) * 100)

plot <- ggplot(data.cr) +
  geom_line(aes(x = time,
                y = cr), size = 1.1)+
  geom_ribbon(aes(ymax = cr.upper, ymin = cr.lower, x = time), fill = "gray70", alpha = 0.5) +
  xlab("Time (year)") +
  ylim(c(0,100)) +
  ylab("Expected core ratio (%)") + 
  theme_bw()+
  theme(text = element_text(size=20))
plot
ggsave(plot, filename = "Output/Tables & figures/Supplementary/FigureS5b_Effect_plot_core_ratio.png")
  
# Figure 6 HR for core ratio ----------------------------------------------
mean(data.cr$cr) # 15

cr_vec <- seq(0.075 ,0.3, 0.005)
which(cr_vec == 0.15)
hazard.value.cr <- exp(pars["alpha[3,1]"] * log(cr_vec/(1-cr_vec)))
hazard.value.cr.full <- exp(mcmc[,"alpha[3,1]"] %*% matrix(log(cr_vec/(1-cr_vec)), nrow = 1))
hazard.value.ratio.cr.full <- apply(hazard.value.cr.full, 1, function(x) {x/x[16]})
data.cr <- data.frame(cr = cr_vec * 100,
                      hazard.ratio = hazard.value.cr/hazard.value.cr[16],
                      hazard.ratio.upper = apply(hazard.value.ratio.cr.full, 1,
                                                 function(x) quantile(x, probs = 0.975)),
                      hazard.ratio.lower = apply(hazard.value.ratio.cr.full, 1,
                                                 function(x) quantile(x, probs = 0.025)))
plot.value.cr <- ggplot(data.cr) +
  geom_line(aes(x = cr,
                y = hazard.ratio), size = 1.1)+
  geom_ribbon(aes(ymax = hazard.ratio.upper, ymin = hazard.ratio.lower, x = cr), fill = "gray70", alpha = 0.5) +
  #geom_point(aes(x = PSA,
  #               y = hazard.ratio), size = 3) +
  xlab("Core ratio (%)") +
  ylab("Hazard ratio") +
  theme_bw() +
  theme(text = element_text(size=25)) +
  ylim(c(0, 4)) +
  geom_hline(aes(yintercept = 1),  linetype = "dashed")
plot.value.cr
ggsave(plot.value.cr, filename = "Output/Tables & figures/Supplementary/FigureS6_effect_plot_core_ratio_value.png")
round(max(data.cr$hazard.ratio), 2)
round(min(data.cr$hazard.ratio), 2)
  

