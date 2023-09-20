load("Original data/Biopsy.RData") # 27 missing
load("Original data/baselinedata.RData") # 17 missing
load("Cleaned data/pass_id.RData") 
library(tidyverse)

ids <- pass.id$CISNET_ID
pass.id$TimeSince_Dx = 0
colnames(pass.id)[16] <- "cores_ratio"
biopsy.nw <- rbind(biopsy[,c("CISNET_ID", "TimeSince_Dx", "cores_ratio")],
                   pass.id[,c("CISNET_ID", "TimeSince_Dx", "cores_ratio")]) %>% 
  filter(CISNET_ID %in% ids) %>% 
  arrange(CISNET_ID, TimeSince_Dx)


pass.cores <- merge(biopsy.nw, pass.id %>% select(- TimeSince_Dx, -PSAValue, 
                                                  -cores_ratio)) %>% 
  filter(TimeSince_Dx <= time.cmp2,
         !is.na(cores_ratio)) %>% 
  arrange(CISNET_ID, TimeSince_Dx)

reference.df <- matrix(NA, 100,100)
for (i in 1:100) {
  reference.df[i,1:i] <- round((1:i)/i * 100,2)
} 
pos_cores <- numeric(2724)
total_cores <- numeric(2724)

# all number can be found
`%!in%` <- Negate(`%in%`)
sum(pass.cores$cores_ratio %!in% c(reference.df, 0))

# find the one with fewest cores as possible
for (i in 1:2724) {
  cr <- pass.cores$cores_ratio[i]
  if (cr %in% reference.df[12,]) {
    total_cores[i] <- 12
    pos_cores[i] <- which(reference.df[12,] == cr)
  } else if (cr %in% reference.df) {
    #total_cores_possible <- do.call(c, sapply(1:100, function(j) {which(reference.df[,j] == cr)}))
    total_cores[i] <- 100
    pos_cores[i] <- round(cr)
  } else if (cr == 0) {
    total_cores[i] <- 12
    pos_cores[i] <- 0
  } else {
    total_cores[i] <- NA
    pos_cores[i] <- NA
  }
}

pass.cores$total_cores <- total_cores
pass.cores$pos_cores <- pos_cores

save(pass.cores, file = "Cleaned data/pass_cores.RData")
