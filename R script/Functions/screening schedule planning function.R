personalizedSchedule.iccsjm <- function (object, newdata, idVar = "CISNET_ID", last_test_time = NULL,
                                          fixed_grid_visits, gap = 0, detection_delay_limit=Inf, total_tests_limit=Inf,
                                          mcmcchainlen = 500,
                                          seed = 100L, iter = 500L, n.adapt = 200, cache_size=100, 
                                          risk_based_schedules_only=T, 
                                          risk_range = seq(0, 1, 0.01), ...) {
  if(length(unique(newdata[[idVar]])) > 1){
    stop("Personalized schedule can be created for only one patient at a time.\n")
  }
  
  timeVar <- object$model_info$var_names$longitime
  
  if(is.null(last_test_time)){
    last_test_time <- max(newdata[[timeVar]], na.rm = T)
  }
  
  cur_visit_time = fixed_grid_visits[1]#max(newdata[[timeVar]], na.rm = T)
  fixed_grid_visits = c(cur_visit_time, fixed_grid_visits[fixed_grid_visits > cur_visit_time])
  horizon = max(fixed_grid_visits)
  
  #Step 1: Create a CACHE of predicted survival probabilities to speed up operation
  RISK_CACHE_TIMES <- seq(from = last_test_time, to = horizon, length.out = cache_size)
  RISK_CACHE_TIMES <- unique(sort(c(fixed_grid_visits, RISK_CACHE_TIMES), decreasing = F))
  risk_pred <- csdypred(object = object, Survtimes = RISK_CACHE_TIMES, iter = iter,
                        newdata = newdata, n.adapt = n.adapt, seed = seed, 
                        mcmcchainlen = mcmcchainlen)
  
  OVERALL_SURV_CACHE_FULL <- t(risk_pred$overall.surv)
  OVERALL_SURV_CACHE_FULL_rescaled <- apply(OVERALL_SURV_CACHE_FULL,2, function(x) {x/x[1]})
  RISK_CACHE_FULL <- do.call('rbind', lapply(risk_pred$full.results, "[[", 1))
  RISK_CACHE_FULL_rescaled <- apply(RISK_CACHE_FULL,2, function(x) {x/x[length(x)]})
  
  
  ##############################################
  #Definition of the function for fixed schedule
  ##############################################
  fixedSchedule <- function(risk_threshold){
    proposed_test_times <- c()
    previous_test_time <- last_test_time
    risk_schedule_temp <- RISK_CACHE_FULL
    #risk.store <- list() # risk results at time of scheduled test
    #num.biopsy <- 0
    #biopsy.do <- rep(0, length(fixed_grid_visits)) # indicator
    
    for(i in 1:length(fixed_grid_visits)){
      if(fixed_grid_visits[i] - previous_test_time>=gap){
        risk_cache_index = which.min(abs(RISK_CACHE_TIMES-fixed_grid_visits[i])) # the time points vectors are updated
        #RISK_CACHE_TIMES[risk_cache_index]
        risk_row_visit = risk_schedule_temp[risk_cache_index,]
        if(mean(risk_row_visit, na.rm = T) >= risk_threshold){
          #num.biopsy <- num.biopsy  + 1
          #biopsy.do[i] <- 1
          #print(risk.start.time)
          proposed_test_times <- c(proposed_test_times, fixed_grid_visits[i])
          previous_test_time = fixed_grid_visits[i]
          #risk.store[[num.biopsy]] <- risk_schedule_temp[risk_schedule_temp[,1] <= fixed_grid_visits[i],]
          risk_schedule_temp <- apply(RISK_CACHE_FULL, 2, FUN = function(x){(x - x[risk_cache_index])}) / 
            matrix(rep(OVERALL_SURV_CACHE_FULL[risk_cache_index,], each = nrow(RISK_CACHE_FULL)), 
                        nrow(RISK_CACHE_FULL),
                        iter)
          
        }
      } 
    }
    
    
    
    if(length(proposed_test_times)==0){
      proposed_test_times = horizon
    }
    
    # return(list(risk_at_biopsy = risk.store,
    #             schedule = data.frame(visit_time = fixed_grid_visits,
    #                                   biopsy_scheduled = biopsy.do),
    #             proposed_test_times = proposed_test_times))
    return(proposed_test_times)
    
  }
  # ptt <- res$proposed_test_times
  #Function to create all possible test schedules 
  allTestSchedules = function(cur_visit_time, min_test_gap, horizon){
    times = list()
    count = 0
    
    fixed_grid_visits = c(-1, fixed_grid_visits)
    
    recursive = function(last_test_time, current_visit_time, cur_path){
      if(current_visit_time < cur_visit_time){
        count <<- count + 1
        times[[count]] <<- as.numeric(strsplit(cur_path, split = "-")[[1]])
      }else{
        if(last_test_time - current_visit_time>=min_test_gap){
          #case 1 do a test
          recursive(current_visit_time,
                    tail(fixed_grid_visits[fixed_grid_visits<current_visit_time],1),
                    paste0(current_visit_time,"-",cur_path))
        }
        
        #case 2 don't do a test
        recursive(last_test_time,
                  tail(fixed_grid_visits[fixed_grid_visits<current_visit_time],1),
                  cur_path)
        
      }
    }
    recursive(horizon, tail(fixed_grid_visits[fixed_grid_visits<horizon],1), 
              as.character(horizon)) # in our case, count = 2^13
    
    #I will sort them by first test
    first_test_time = sapply(times, min)
    total_biopsies = sapply(times, length)
    times = times[order(first_test_time, total_biopsies, decreasing = F)]
    return(times)
    
    # times = do.call("c", lapply(seq_along(fixed_grid_visits[-length(fixed_grid_visits)]), 
    #                     function(i) combn(fixed_grid_visits[-length(fixed_grid_visits)], i, FUN = list)))
    # 
    # times = lapply(times, function(x) {c(x, horizon)})
  }
  
  wk = JMbayes2:::gaussKronrod()$wk
  sk = JMbayes2:::gaussKronrod()$sk
  getGaussianQuadWeightsPoints = function(lower_upper_limit){
    p1 = 0.5 * sum(lower_upper_limit)
    p2 = 0.5 * diff(lower_upper_limit)
    
    return(list(weights = p2*wk, points = p2 * sk + p1))
  } 
  
  
  ###############################################
  # Consequences
  ###############################################
  getConsequences = function(proposed_test_times){
    planned_test_schedule = proposed_test_times
    #Make the last test at horizon
    # if the last two tests are too close, just do at the last time point
    if(horizon - tail(planned_test_schedule,1) <= gap){
      planned_test_schedule = planned_test_schedule[-length(planned_test_schedule)]
    }
    planned_test_schedule = c(planned_test_schedule, horizon)
    
    #Now we create intervals in which to integrate the conditional surv prob
    if(length(planned_test_schedule)==1){
      test_intervals = list(c(last_test_time, planned_test_schedule))
    }else{
      test_intervals = c(last_test_time,
                         rep(planned_test_schedule[-length(planned_test_schedule)],each=2),
                         planned_test_schedule[length(planned_test_schedule)])
      test_intervals = split(test_intervals, rep(1:(length(test_intervals)/2), each=2))
    }
    
    interval_res = vector("list", length(test_intervals))
    exp_num_tests = rep(0, iter)
    for(j in 1:length(test_intervals)){
      wt_points=getGaussianQuadWeightsPoints(test_intervals[[j]])
      lower_limit = test_intervals[[j]][1]
      upper_limit = test_intervals[[j]][2]
      
      lower_limit_nearest_index = which.min(abs(lower_limit - RISK_CACHE_TIMES))
      upper_limit_nearest_index = which.min(abs(upper_limit - RISK_CACHE_TIMES))
      interval_res[[j]]$cum_risk_interval = RISK_CACHE_FULL_rescaled[upper_limit_nearest_index,] - 
        RISK_CACHE_FULL_rescaled[lower_limit_nearest_index,]
      
      cond_expected_fail_time = sapply(1:length(wt_points$points), function(i){
        overall_surv_at_points = RISK_CACHE_FULL_rescaled[upper_limit_nearest_index,] - 
          RISK_CACHE_FULL_rescaled[which.min(abs(wt_points$points[i]-RISK_CACHE_TIMES)),]
          
          # OVERALL_SURV_CACHE_FULL_rescaled[which.min(abs(wt_points$points[i]-RISK_CACHE_TIMES)),] - 
          # OVERALL_SURV_CACHE_FULL_rescaled[upper_limit_nearest_index,]
        scaled_cum_surv_at_points = overall_surv_at_points/interval_res[[j]]$cum_risk_interval
        
        return(wt_points$weights[i] * scaled_cum_surv_at_points)
      })
      
      interval_res[[j]]$delay = upper_limit - (lower_limit + apply(cond_expected_fail_time, 1, sum)) 
      
      #### impose a restriction (need inspection later!!!!!!!!!!!!!)
      #interval_res[[j]]$delay = sapply(interval_res[[j]]$delay, function(x) {ifelse(x < 0, 0,x)})
      
      exp_num_tests = exp_num_tests + j * interval_res[[j]]$cum_risk_interval
    }
    # need investigate (origin uses survival probability)
    #exp_num_tests = exp_num_tests + (j+1) * OVERALL_SURV_CACHE_FULL_rescaled[upper_limit_nearest_index,] #RISK_CACHE_FULL_rescaled
    exp_num_tests = mean(exp_num_tests, na.rm = T)
    
    expected_detection_delay = mean(apply(sapply(interval_res, FUN = function(x){
      x$delay * x$cum_risk_interval
    }),1, sum, na.rm=T),na.rm=T)
    
    #proposed test times always has a final test at horizon
    return(list(expected_detection_delay = expected_detection_delay,
                expected_num_tests = exp_num_tests,
                planned_test_schedule = planned_test_schedule))
  }
  
  
  # Now we create schedules
  if(risk_based_schedules_only==T){
    risk_thresholds = risk_range
    schedule_plans = lapply(risk_thresholds, fixedSchedule)
    all_schedules <- vector("list", length=length(schedule_plans))
    names(all_schedules) <- c(paste("Risk:", risk_range))
  }else{
    schedule_plans = allTestSchedules(cur_visit_time, gap, horizon)
    all_schedules <- vector("list", length=length(schedule_plans))
  }
  
  for(i in 1:length(all_schedules)){
    all_schedules[[i]] = getConsequences(schedule_plans[[i]])
    all_schedules[[i]]$euclidean_distance = sqrt(all_schedules[[i]]$expected_detection_delay^2 + (all_schedules[[i]]$expected_num_tests-1)^2)
    
    if(risk_based_schedules_only==T){
      all_schedules[[i]]$cumulative_risk_threshold = risk_thresholds[i]
    }
  }
  
  expected_delays = sapply(all_schedules, "[[", "expected_detection_delay")
  expected_num_tests = sapply(all_schedules, "[[", "expected_num_tests")
  dist = sapply(all_schedules, "[[", "euclidean_distance")
  
  ############
  # this is the code for applying constraints before final optimization
  ###########
  flag = 0
  old_total_test_limit = total_tests_limit
  while(sum(expected_num_tests<=total_tests_limit)==0){
    flag = 1
    total_tests_limit = total_tests_limit + 1
  }
  
  if(flag==1){
    print(paste0("Could not find any schedule with a total_tests_limit=", old_total_test_limit))
    print(paste0("Using new total_tests_limit=",total_tests_limit))
  }
  
  flag = 0
  old_detection_delay_limit = detection_delay_limit
  delta_delay = diff(range(expected_delays))/100
  while(sum(expected_delays<=detection_delay_limit) == 0){
    flag = 1
    detection_delay_limit = detection_delay_limit + delta_delay
  }
  
  if(flag==1){
    print(paste0("Could not find any schedule with a detection_delay_limit=",old_detection_delay_limit))
    print(paste0("Using new detection_delay_limit=",detection_delay_limit))
  }
  
  filter = expected_delays<=detection_delay_limit & expected_num_tests<=total_tests_limit
  filtered_schedules = all_schedules[filter]
  filtered_dist = dist[filter]
  
  #There are many, but we select only one for now
  optimal_threshold_index = which.min(filtered_dist)[1]
  ret = filtered_schedules[[optimal_threshold_index]]
  
  full_data = list("optimal_schedule"=ret, "all_schedules"=all_schedules,
                   "next_biopsy" = ret$planned_test_schedule[1])
  
  class(full_data) <- "personalizedSchedule.iccsjm"
  return(full_data)
}



