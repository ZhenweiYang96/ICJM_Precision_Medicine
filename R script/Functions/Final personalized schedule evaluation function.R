personalizedSchedule.iccsjm_evaluation <- function(object, newdata, last_test_time = 0,
                                                   fixed_grid_visits, gap = 0, detection_delay_limit=Inf, total_tests_limit=Inf,
                                                   mcmcchainlen = 500,
                                                   seed = 100L, iter = 100, n.adapt = 200, cache_size=100, 
                                                   risk_range = seq(0, 1, 0.01)) {
  # get information 
  time_event <- newdata$time.cmp2[1]
  
  biopsy_schedules <- numeric(0)
  num_biopsy <- 0
  time_last_biopsy <- last_test_time
  time_current <- fixed_grid_visits[1]
  clinical_visits_future <- fixed_grid_visits[fixed_grid_visits > time_last_biopsy]
  next_biopsy <- Inf
  for (i in 1:(length(fixed_grid_visits)-1)) {
    # print(i)
    # print(fixed_grid_visits[i])
    # t1 <- Sys.time()
    testset <- newdata[newdata$TimeSince_Dx <= time_current,]
    proposed_schedules <- personalizedSchedule.iccsjm(object = object, newdata = testset, 
                                                        last_test_time = time_last_biopsy,
                                                        fixed_grid_visits = clinical_visits_future,
                                                        gap = gap,
                                                        detection_delay_limit = detection_delay_limit, 
                                                        iter = iter, cache_size = cache_size, 
                                                        risk_range = risk_range, seed = seed,
                                                        mcmcchainlen = mcmcchainlen) 
    # t2 <- Sys.time()
    # print(t2 - t1)
    #print(proposed_schedules$optimal_schedule$planned_test_schedule)
    #print(proposed_schedules$next_biopsy)
    next_biopsy <- proposed_schedules$next_biopsy
    if (next_biopsy == time_current) {
      num_biopsy <- num_biopsy + 1
      biopsy_schedules[num_biopsy] <- next_biopsy
      time_last_biopsy <- time_current
      clinical_visits_future <- fixed_grid_visits[fixed_grid_visits > time_current]
      generated_schedules <- proposed_schedules$optimal_schedule$planned_test_schedule
      next_biopsy <- generated_schedules[generated_schedules > next_biopsy][1]
      if (time_current > time_event) {
        biopsy_schedules <- c(biopsy_schedules, 
                              generated_schedules[generated_schedules > time_current])
        break
      }
    }
    # if (time_current > time_event) {
    #   biopsy_schedules <- c(biopsy_schedules, 
    #                         generated_schedules[generated_schedules > time_current])
    #   break
    # }
    clinical_visits_future <- fixed_grid_visits[fixed_grid_visits > time_current]
    time_current <- fixed_grid_visits[i+1]
    # print(proposed_schedules$optimal_schedule$planned_test_schedule)
  }
  
  if (length(biopsy_schedules) == 0) {
    biopsy_schedules <- max(fixed_grid_visits)
  }
  if (max(biopsy_schedules) != max(fixed_grid_visits)) {
    biopsy_schedules <- c(biopsy_schedules, max(fixed_grid_visits))
  }
  return(biopsy_schedules)
  
}
