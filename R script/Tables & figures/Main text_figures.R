library(ggplot2)
library(tidyverse)
library(splines)
library(latex2exp)
library(cowplot)
library(rjags)
load("Cleaned data/pass_id.RData")
load("Cleaned data/pass.RData")
load("Cleaned data/pass_cores.RData")


# Figure 2  --------------------------------------------------------------------
load("Output/Data analysis/ICJM1.RData")
source("R script/Functions/dynpred_Risk prediction function.R")
source("R script/Functions/function.R")
source("R script/Functions/screening schedule planning function.R")

test_163 <- pass[pass$CISNET_ID == 163,]
testset <- test_163[test_163$TimeSince_Dx <1.5,]

result <- personalizedSchedule.iccsjm(object = ICJM1, newdata = testset, 
                                      last_test_time = 0.5, fixed_grid_visits = seq(1.5,9,0.5), 
                                      iter = 200)

loss.df <- data.frame(threshold = numeric(101),
                      expected_num_tests = 0,
                      expected_detection_delay = 0,
                      euclidean_distance = 0)
for (i in 1:101) {
  loss.df$threshold[i] = result$all_schedules[[i]]$cumulative_risk_threshold
  loss.df$expected_num_tests[i] = result$all_schedules[[i]]$expected_num_tests
  loss.df$expected_detection_delay[i] = result$all_schedules[[i]]$expected_detection_delay
  loss.df$euclidean_distance[i] = result$all_schedules[[i]]$euclidean_distance
}

loss.optimal <- loss.df[which.min(loss.df$euclidean_distance),]
loss.df <- loss.df[-which.min(loss.df$euclidean_distance),]
loss.df <- loss.df[!duplicated(loss.df$expected_num_tests),]

plot.loss <- ggplot() +
  geom_segment(aes(x = 1, xend = loss.df$expected_num_tests, y = 0, yend = loss.df$expected_detection_delay), linetype = "dashed", color ="gray", alpha = 0.3) + 
  geom_segment(aes(x = 1, xend = loss.optimal$expected_num_tests, y = 0, yend = loss.optimal$expected_detection_delay), color = "red", linetype = "dashed") +
  geom_point(aes(x = loss.df[loss.df$expected_detection_delay<1.5,]$expected_num_tests, 
                 y = loss.df[loss.df$expected_detection_delay<1.5,]$expected_detection_delay), size = 8) + 
  geom_point(aes(x = loss.df[loss.df$expected_detection_delay>=1.5,]$expected_num_tests, 
                 y = loss.df[loss.df$expected_detection_delay>=1.5,]$expected_detection_delay), size = 8, alpha = 0.4, fill = "lightgray") +
  geom_point(aes(loss.optimal$expected_num_tests, loss.optimal$expected_detection_delay), size = 10, shape = 23, fill="red") + 
  geom_hline(aes(yintercept = 1.5), linetype = "dashed", color = "orange") + 
  geom_point(aes(x=  1, y = 0), shape = 24, fill = "darkred", size = 6, color ="white") +
  geom_label(aes(x = 2.1, y = 0.1, label ="Ideal schedule"), size = 6, fontface = "bold") + # , family = "CMU Serif"
  geom_label(aes(x = 5.1, y = 1.5, label = "Maximum clinically acceptable delay"), color = "orange", size = 7, fontface = "bold") + # , family = "CMU Serif"
  geom_label(aes(x = 7.3, y = 0.7, label = "\u03D5 = 0%"), size = 6, fill ="black", color= "white") + #\u03BA # , family = "Arial"
  geom_label(aes(x = 1.9, y = 4.55, label = "\u03D5 = 100%"), size = 6, fill = "gray", color= "white") + # , family = "Arial"
  geom_label(aes(x = 1.7, y = 0.75, label = "\u03D5* = 14%"), size = 6, fill = "red", color ="white") + # , family = "Arial"
  scale_x_continuous(breaks = seq(1,9,1)) + 
  ylab("Expected detection delay (years)") + 
  xlab("Expected number of biopsies") + 
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", size = 1, fill= NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=18), # , family = "LM Roman 10"
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))
plot.loss
ggsave(plot.loss, file = "Output/Tables & figures/Main text/Figure2_loss function.png")

# Figure 3a Effect plot for PSA value ---------------------------------------------------------------------
load("Output/Data analysis/ICJM2.RData")
mcmc <- do.call("rbind",lapply(1:3, function(x) {as.matrix(ICJM2$mcmc[[x]])}))
pars <- apply(mcmc, 2, mean)
pars.lower <- apply(mcmc, 2, function(x) {quantile(x, prob = 0.025)})
pars.upper <- apply(mcmc, 2, function(x) {quantile(x, prob = 0.975)})

knot.longi <- ICJM2$model_info$knots$knot.longi
T <- ns(seq(0,10, length.out =  100), knots = knot.longi[2:3], B = knot.longi[c(1,4)])
X.t <- cbind(1, T, 0)
Y <- 2^(X.t %*% pars[grep("betaL_psa", names(pars), value = T)]) - 1
data <- data.frame(time = seq(0,10, length.out =  100),
                   PSA = Y,
                   PSA.lower = 2^(X.t %*% pars.lower[grep("beta", names(pars), value = T)][4:8]) - 1,
                   PSA.upper = 2^(X.t %*% pars.upper[grep("beta", names(pars), value = T)][4:8]) - 1)
mean(data$PSA) # mean value is 5ng/ml

psa.vec <- seq(2.5 ,10, 0.1)
which(psa.vec == 5)
hazard.value <- exp(pars["alpha[1,1]"] * log(psa.vec+1,base = 2))
hazard.value.full <- exp(mcmc[,"alpha[1,1]"] %*% matrix(log(psa.vec+1,base = 2), nrow = 1))
hazard.value.ratio.full <- apply(hazard.value.full, 1, function(x) {x/x[26]})
data.psa <- data.frame(PSA = psa.vec,
                       hazard.ratio = hazard.value/hazard.value[26],
                       hazard.ratio.upper = apply(hazard.value.ratio.full, 1, 
                                                  function(x) quantile(x, probs = 0.975)),
                       hazard.ratio.lower = apply(hazard.value.ratio.full, 1, 
                                                  function(x) quantile(x, probs = 0.025)))
plot.value <- ggplot(data.psa) +
  geom_line(aes(x = PSA,
                y = hazard.ratio), size = 1.1)+
  geom_ribbon(aes(ymax = hazard.ratio.upper, ymin = hazard.ratio.lower, x = PSA), fill = "gray70", alpha = 0.5) +
  #geom_point(aes(x = PSA,
  #               y = hazard.ratio), size = 3) +
  xlab("PSA level (ng/ml)") +
  ylab("Hazard ratio") + 
  theme_bw() +
  theme(text = element_text(size=23)) +  # , family = "LM Roman 10"
  geom_hline(aes(yintercept = 1), linetype = "dashed")
plot.value
ggsave(plot.value, filename = "Output/Tables & figures/Main text/Figure3a_effect plot - PSA value.png", 
       width = 22, height = 18, unit = "cm")
round(max(data.psa$hazard.ratio), 2)
round(c(max(data.psa$hazard.ratio.lower), max(data.psa$hazard.ratio.upper)), 2)
round(min(data.psa$hazard.ratio), 2)
round(c(min(data.psa$hazard.ratio.lower), min(data.psa$hazard.ratio.upper)), 2)

# Figure 3b Effect plot for PSA yearly change ---------------------------------------------------------------
psa <- data[,2]
idx = 1:90
mean(psa[idx + 10] - psa[idx]) # mean change is 0.33
psa.diff.vec <- seq(0.15, 0.6, 0.01)

which(psa.diff.vec == 0.3)
hazard.diff <- exp(pars["alpha[2,1]"] * 
                     (log(5 + 1, base  = 2) - log(5 - psa.diff.vec + 1, base = 2)))
hazard.diff.full <- exp(mcmc[,"alpha[2,1]"] %*% matrix((log(5 + 1, base  = 2) - log(5 - psa.diff.vec + 1, base = 2)), 
                                                       nrow = 1))
hazard.diff.ratio.full <- apply(hazard.diff.full, 1, function(x) {x/x[16]})

data.psadiff <- data.frame(PSA = psa.diff.vec,
                           hazard.ratio = hazard.diff/hazard.diff[16],
                           hazard.ratio.upper = apply(hazard.diff.ratio.full, 1, 
                                                      function(x) quantile(x, probs = 0.975)),
                           hazard.ratio.lower = apply(hazard.diff.ratio.full, 1, 
                                                      function(x) quantile(x, probs = 0.025)))
plot.diff <- ggplot(data.psadiff) +
  geom_line(aes(x = PSA,
                y = hazard.ratio), size = 1.1)+
  geom_ribbon(aes(ymax = hazard.ratio.upper, ymin = hazard.ratio.lower, x = PSA), fill = "gray70", alpha = 0.5) +
  xlab("Change in PSA level over the previous year (ng/ml)") +
  ylim(c(0.7,1.5)) + 
  ylab("Hazard ratio") + 
  theme_bw() +
  theme(text = element_text(size=23)) +  # , family = "LM Roman 10"
  geom_hline(aes(yintercept = 1), linetype = "dashed")
plot.diff
ggsave(plot.diff, filename = "Output/Tables & figures/Main text/Figure3b_Effect plot - PSA diff.png",
       width = 22, height = 18, unit = "cm")
round(max(data.psadiff$hazard.ratio), 2)
round(c(max(data.psadiff$hazard.ratio.lower), max(data.psadiff$hazard.ratio.upper)), 2)
round(min(data.psadiff$hazard.ratio), 2)
round(c(min(data.psadiff$hazard.ratio.lower), min(data.psadiff$hazard.ratio.upper)), 2)

# Figure 4 HR for core ratio ----------------------------------------------
T.cr <- poly(seq(0,10, length.out =  100), degree = 2, raw = TRUE)
X.t.cr <- cbind(1, T.cr)
Y.cr <- X.t.cr %*% pars[grep("beta", names(pars), value = T)][1:3]
data.cr <- data.frame(time = seq(0,10, length.out =  100),
                      cr = exp(Y.cr)/(1 + exp(Y.cr)) * 100,
                      cr.lower = exp(X.t.cr %*% pars.lower[grep("beta", names(pars), value = T)][1:3])/
                        (1 + exp(X.t.cr %*% pars.lower[grep("beta", names(pars), value = T)][1:3])) * 100,
                      cr.upper = exp(X.t.cr %*% pars.upper[grep("beta", names(pars), value = T)][1:3])/
                        (1 + exp(X.t.cr %*% pars.upper[grep("beta", names(pars), value = T)][1:3])) * 100)
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
ggsave(plot.value.cr, filename = "Output/Tables & figures/Main text/Figure4_effect plot - core ratio value.png")
round(max(data.cr$hazard.ratio), 2)
round(c(max(data.cr$hazard.ratio.lower), max(data.cr$hazard.ratio.upper)), 2)
round(min(data.cr$hazard.ratio), 2)
round(c(min(data.cr$hazard.ratio.lower), min(data.cr$hazard.ratio.upper)), 2)

# Figure 5a number of biopsies ------------------------------------------
comparison <- data.frame(dataset = rep(1:200, each = 100), 
                         id = rep(1:100, 200),
                         num_biopsy_flexible = 0,
                         num_biopsy_fixed_1 = 0,
                         num_biopsy_fixed_2 = 0,
                         detection_delay_flexible = 0,
                         detection_delay_fixed_1 = 0,
                         detection_delay_fixed_2 = 0,
                         patients = "progressed")

comparison.nonpro <- data.frame(dataset = rep(1:200, each = 100), 
                                id = rep(1:100, 200),
                                num_biopsy_flexible = 0,
                                num_biopsy_fixed_1 = 0,
                                num_biopsy_fixed_2 = 0,
                                patients = "non-progressed")

for (i in 1:200) {
  load(paste0("Output/Simulation/Personalized schedule result [progressing]/indicators_",i,".RData"))
  comparison$num_biopsy_flexible[comparison$dataset == i] <- sapply(1:100, function(j) {result[[j]]$num_biopsy_flexible})
  comparison$num_biopsy_fixed_1[comparison$dataset == i] <- sapply(1:100, function(j) {result[[j]]$num_biopsy_fixed_1})
  comparison$num_biopsy_fixed_2[comparison$dataset == i] <- sapply(1:100, function(j) {result[[j]]$num_biopsy_fixed_2})
  comparison$detection_delay_flexible[comparison$dataset == i] <- sapply(1:100, function(j) {result[[j]]$treatment_delay_flexible})
  comparison$detection_delay_fixed_1[comparison$dataset == i] <- sapply(1:100, function(j) {result[[j]]$treatment_delay_fixed_1})
  comparison$detection_delay_fixed_2[comparison$dataset == i] <- sapply(1:100, function(j) {result[[j]]$treatment_delay_fixed_2})
  load(paste0("Output/Simulation/Personalized schedule result [non progressing]/indicators_",i,".RData"))
  comparison.nonpro$num_biopsy_flexible[comparison$dataset == i] <- sapply(1:100, function(j) {result[[j]]$num_biopsy_flexible})
  comparison.nonpro$num_biopsy_fixed_1[comparison$dataset == i] <- sapply(1:100, function(j) {result[[j]]$num_biopsy_fixed_1})
  comparison.nonpro$num_biopsy_fixed_2[comparison$dataset == i] <- sapply(1:100, function(j) {result[[j]]$num_biopsy_fixed_2})
} 

boxplot.num.all <- comparison %>% 
  dplyr::select(-detection_delay_flexible:-detection_delay_fixed_2) %>% 
  gather(Schedule_type, num_biopsies, num_biopsy_flexible:`num_biopsy_fixed_2`) %>% 
  rbind(gather(comparison.nonpro, Schedule_type, num_biopsies, num_biopsy_flexible:`num_biopsy_fixed_2`)) %>% 
  ggplot(aes(x = Schedule_type, y = num_biopsies, color = as.factor(patients))) +
  stat_boxplot(geom = "errorbar") +  
  geom_boxplot() +
  xlab("Schedule type") +
  ylab("Number of biopsies conducted") + 
  theme_bw() +
  theme(legend.position = 'top',
        text = element_text(size=  18), # , family = "LM Roman 10", face = "bold"
        panel.border = element_rect(color = "black", size = 1, fill= NA)) + 
  scale_color_grey(name = "Type of patients", start = 0.1, end = 0.6) +
  scale_x_discrete(labels=c("Fixed\n(PASS)","Fixed\n(Annual)","Personalized")) + 
  ylim(c(0, 11))
boxplot.num.all
ggsave(boxplot.num.all, file = "Output/Tables & figures/Main text/Figure5a_num_biopsies.png")

# PROGRESSION - number of biopsies - number
median(comparison$num_biopsy_flexible)
median(comparison$num_biopsy_fixed_1)
median(comparison$num_biopsy_fixed_2)
round(mean(comparison$num_biopsy_fixed_1 - comparison$num_biopsy_flexible), 2)
round(mean(comparison$num_biopsy_fixed_2 - comparison$num_biopsy_flexible), 2)
round(mean((comparison$num_biopsy_fixed_1 - comparison$num_biopsy_flexible)/comparison$num_biopsy_fixed_1), 2)
round(mean((comparison$num_biopsy_fixed_2 - comparison$num_biopsy_flexible)/comparison$num_biopsy_fixed_2), 2)

# NON-PROGRESSION - number of biopsies - number
median(comparison.nonpro$num_biopsy_flexible)
median(comparison.nonpro$num_biopsy_fixed_1)
median(comparison.nonpro$num_biopsy_fixed_2)
round(mean(comparison.nonpro$num_biopsy_fixed_1 - comparison.nonpro$num_biopsy_flexible), 2)
round(mean(comparison.nonpro$num_biopsy_fixed_2 - comparison.nonpro$num_biopsy_flexible), 2)
round(mean((comparison.nonpro$num_biopsy_fixed_1 - comparison.nonpro$num_biopsy_flexible)/comparison.nonpro$num_biopsy_fixed_1), 2)
round(mean((comparison.nonpro$num_biopsy_fixed_2 - comparison.nonpro$num_biopsy_flexible)/comparison.nonpro$num_biopsy_fixed_2), 2)


# Figure 5b progressing - detection delay ---------------------------------
boxplot.delay.all <- comparison %>% 
  filter(!is.na(detection_delay_fixed_1)) %>% 
  gather(Schedule_type, delay, detection_delay_flexible:`detection_delay_fixed_2`) %>% 
  ggplot(aes(x = Schedule_type, y = delay)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot() +
  xlab("Schedule type") +
  ylab("Detection delay (year)") + 
  theme_bw() +
  theme(text = element_text(size=  18),
        panel.border = element_rect(color = "black", size = 1, fill= NA)) + 
  scale_x_discrete(labels=c("Fixed\n(PASS)","Fixed\n(Annual)","Personalized"))
boxplot.delay.all
ggsave(boxplot.delay.all, file = "Output/Tables & figures/Main text/Figure5b_detection_delay (progressing).png")

# PROGRESSION - detection delay - number
median(comparison$detection_delay_flexible, na.rm = TRUE)
median(comparison$detection_delay_fixed_1[!is.na(comparison$detection_delay_fixed_1)])
median(comparison$detection_delay_fixed_2[!is.na(comparison$detection_delay_fixed_2)])



#### Plot 3a and 4b together

fig5 <- plot_grid(boxplot.num.all, boxplot.delay.all, align = "h", axis = "bt", rel_widths = c(1, 1))
ggsave(fig5, file = "Output/Tables & figures/Main text/Figure5_ddnb.png",
       width = 30, height = 15, units = "cm")


