TestpersonalizedScreening <- function(object, newdata, idVar = "CISNET_ID", 
                                      last_test_time = NULL, fixed_grid_visits, 
                                      seed = 100L, iter = 200L, n.adapt = 200, cache_size=100,
                                      risk_threshold) { 
  
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
  risk_pred <- csdypred(object = object, Survtimes = RISK_CACHE_TIMES, iter = iter, seed = seed, 
                        newdata = newdata, n.adapt = n.adapt)
  
  OVERALL_SURV_CACHE_FULL <- t(risk_pred$overall.surv)
  OVERALL_SURV_CACHE_FULL_rescaled <- apply(OVERALL_SURV_CACHE_FULL,2, function(x) {x/x[1]})
  RISK_CACHE_FULL <- do.call('rbind', lapply(risk_pred$full.results, "[[", 1))
  RISK_CACHE_FULL_rescaled <- apply(RISK_CACHE_FULL,2, function(x) {x/x[length(x)]})
  
  # store list of calculated risk 
  proposed_test_times <- c()
  previous_test_time <- last_test_time
  risk_schedule_temp <- RISK_CACHE_FULL
  risk.store <- list() # risk results at time of scheduled test
  num.biopsy <- 0
  biopsy.do <- rep(0, length(fixed_grid_visits)) # indicator
  
  for(i in 1:length(fixed_grid_visits)){
    risk_cache_index = which.min(abs(RISK_CACHE_TIMES-fixed_grid_visits[i])) # the time points vectors are updated
    #RISK_CACHE_TIMES[risk_cache_index]
    risk_row_visit = risk_schedule_temp[risk_cache_index,]
    if(mean(risk_row_visit, na.rm = T) >= risk_threshold){
      num.biopsy <- num.biopsy  + 1
      biopsy.do[i] <- 1
      #print(risk.start.time)
      proposed_test_times <- c(proposed_test_times, fixed_grid_visits[i])
      previous_test_time = fixed_grid_visits[i]
      
      # store the risk curve
      start_index <- which(risk_schedule_temp[,1] == 0)
      risk_schedule_cache <- risk_schedule_temp[start_index:risk_cache_index,]
      risk_pred_store <- data.frame(t.hor = risk_pred$risk[start_index:risk_cache_index,1],
                                    mean = apply(risk_schedule_cache,1,mean),
                                    upper = apply(risk_schedule_cache,1,function(x) {quantile(x, probs = 0.975)}),
                                    lower = apply(risk_schedule_cache,1,function(x) {quantile(x, probs = 0.025)}))
      
      risk.store[[num.biopsy]] <- risk_pred_store
      
      # update the temp 
      risk_schedule_temp <- apply(RISK_CACHE_FULL, 2, FUN = function(x){(x - x[risk_cache_index])}) / 
        matrix(rep(OVERALL_SURV_CACHE_FULL[risk_cache_index,], each = nrow(RISK_CACHE_FULL)), 
               nrow(RISK_CACHE_FULL),
               iter)
      
    }
  }
  
  if(length(proposed_test_times)==0){
    proposed_test_times = horizon
  }
  
  if (max(proposed_test_times) != horizon) {
    num.biopsy <- num.biopsy  + 1
    proposed_test_times <- c(proposed_test_times, fixed_grid_visits[i])
    previous_test_time = fixed_grid_visits[i]
    
    # store the risk curve
    start_index <- which(risk_schedule_temp[,1] == 0)
    risk_schedule_cache <- risk_schedule_temp[start_index:risk_cache_index,]
    risk_pred_store <- data.frame(t.hor = risk_pred$risk[start_index:risk_cache_index,1],
                                  mean = apply(risk_schedule_cache,1,mean),
                                  upper = apply(risk_schedule_cache,1,function(x) {quantile(x, probs = 0.975)}),
                                  lower = apply(risk_schedule_cache,1,function(x) {quantile(x, probs = 0.025)}))
    
    risk.store[[num.biopsy]] <- risk_pred_store
  }  
  
  return(list(risk_at_biopsy = risk.store,
              schedule = data.frame(visit_time = fixed_grid_visits,
                                    biopsy_scheduled = biopsy.do),
              PSA_predicted = risk_pred$longi,
              proposed_test_times = proposed_test_times))
}
