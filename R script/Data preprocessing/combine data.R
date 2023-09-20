############################
## Notice: for now, the time variable for people without event should be their last negative biopsy (standard)
############################


# data
load("Original data/baselinedata.RData")
load("Original data/PSA.RData")
load("Original data/Biopsy.RData")

# packages
library(tidyverse)


# event -------------------------------------------------------------------
biopsy$sumscore <- biopsy$primary_gleason + biopsy$secondary_gleason
res <- NULL
progress_time <- NULL
for (i in 1:nrow(biopsy)) {
  score <- biopsy$sumscore[biopsy$CISNET_ID == biopsy$CISNET_ID[i]]
  time <- biopsy$TimeSince_Dx[biopsy$CISNET_ID == biopsy$CISNET_ID[i]]
  if (length(which(score[complete.cases(score)] >= 7)) > 1) {
    res <- c(res, "progressed-middle")
    progress_time <- c(progress_time, time[which(score >=7)[1]])
  } else if (all(length(which(score >= 7)) == 1, which.max(score) != length(score))) {
    res <- c(res, "progressed-middle")
    progress_time <- c(progress_time, time[which(score >=7)[1]])
  } else if (any(all(is.na(score)), score[which.max(score)] == 6)) {
    res <- c(res, "not progressed")
    progress_time <- c(progress_time, max(time))
  } else if (all(length(which(score >= 7)) == 1, which.max(score) == length(score))) {
    res <- c(res, "progressed-last")
    progress_time <- c(progress_time, time[length(time)])
  }
}
biopsy$progress_dtl<- res
biopsy$progress_time <- progress_time

biopsy <- biopsy %>% mutate(progressed = case_when(progress_dtl %in% c("progressed-middle", "progressed-last") ~ 1, progress_dtl == "not progressed" ~ 0))

# add progressed and progressed_dtl to baseline data
pass <- merge(bldata, biopsy[!duplicated(biopsy$CISNET_ID),c(1,9,10,11)])

pass <- pass %>% 
  mutate(status.cmp = case_when(
    progressed == 0 & event_trt == 0 ~ 0,
    progressed == 0 & event_trt == 1 ~ 2,
    T ~ 1
  )) %>% 
  mutate(status = case_when(status.cmp >= 1 ~ 1,
                            T ~ 0))
# end event


# start - end time set aside for now **-------------------------------------------------------

start <- NULL; end <- NULL
for(i in 1:850) {
  subdt <- biopsy[biopsy$CISNET_ID == i,]
  if (subdt$progress_dtl[1] == "progressed-middle") {
    start[i] = ifelse(which(subdt$sumscore >= 7)[1] !=1, subdt$TimeSince_Dx[which(subdt$sumscore >= 7)[1] -1], 0.001) 
    end[i] =  subdt$TimeSince_Dx[which(subdt$sumscore >= 7)[1]]
  } else if (subdt$progress_dtl[1] == "not progressed") {
    start[i] = NA#max(subdt$TimeSince_Dx)
    end[i] = NA
  } else if (subdt$progress_dtl[1] == "progressed-last") {
    start[i] = ifelse(nrow(subdt) !=1, subdt$TimeSince_Dx[which.max(subdt$TimeSince_Dx) - 1], 0)
    end[i] = max(subdt$TimeSince_Dx)
  }
}
pass$start_time <- start
pass$end_time <- end

lastsecond <- NULL
for(i in 1:850) {
  subdt <- biopsy[biopsy$CISNET_ID == i,]
  if (subdt$progressed[1] == 1) {
    lastsecond[i] = ifelse(which(subdt$sumscore >= 7)[1] !=1, subdt$TimeSince_Dx[which(subdt$sumscore >= 7)[1] -1], 0.001) 
  } else if (subdt$progressed[1] == 0) {
    lastsecond[i] = NA
  }
}
pass$lastsecond_progress_time <- lastsecond

# time to event -----------------------------------------------------------

pass <- pass %>% 
  mutate(time = case_when(status.cmp == 1 ~ progress_time,
                          status.cmp == 2 ~ TimeOnAS,
                          status.cmp == 0 ~ progress_time),
         time.cmp1 = case_when(status.cmp == 1 ~ start_time,
                               status.cmp == 2 ~ progress_time,
                               status.cmp == 0 ~ progress_time),
         time.cmp2 = case_when(status.cmp == 1 ~ end_time,
                               status.cmp == 2 ~ TimeOnAS,
                               status.cmp == 0 ~ TimeOnAS))


# density -----------------------------------------------------------------

pass <- pass %>% 
  mutate(density = log(dx_psa/prostate_size))

# Merge with PSA ----------------------------------------------------------

psa <- psa[psa$PSAValue >=0,]
psa$PSAValue <- log(psa$PSAValue+1, base = 2)

# merge
pass <- merge(psa, pass)

# exclude those after the event
pass <- pass[pass$TimeSince_Dx < pass$time,]

# pass id -----------------------------------------------------------------

pass.id <- pass[!duplicated(pass$CISNET_ID),]

# Center Age --------------------------------------------------------------
pass$DxAge <- pass$DxAge - median(pass.id$DxAge)
pass.id$DxAge <- pass.id$DxAge - median(pass.id$DxAge)

# Save the data -----------------------------------------------------------

save(pass, file = "Cleaned data/pass.RData")
save(pass.id, file = "Cleaned data/pass_id.RData")

