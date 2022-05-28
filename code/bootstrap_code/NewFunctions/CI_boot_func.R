# require input 'data_obs' is of the form c("date", "observations)

# need functions: 'Get_Re_est_func' and 'Para_model_func'

# deleted things about semi-parametric and other methods e.g., union...

##########

CI_boot_func <- function(data_obs, alpha=0.95, n_boot, block_size=10, 
                         starting_threshold=1,
                         IncubationParams, OnsetToCountParams, 
                         smooth_para_boot=21, smooth_para_deConv=21){
  
  start_date_index <- which((data_obs$observations>=starting_threshold)==TRUE)[1]
  data_obs <- data_obs[start_date_index:nrow(data_obs),]
  
  ####
  data_obs_log <- data_obs
  data_obs_log$log_observations <- log(data_obs_log$observations + 1)
  
  smoothed_incidence_data <- data_obs_log %>%
    complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(log_observations = 0)) %>%
    mutate(log_loess = getLOESSCases(dates = date, count_data = log_observations, days_incl=smooth_para_boot),
           log_diff = log_observations - log_loess)
  log_diff <- smoothed_incidence_data$log_diff
  
  # plot(log_diff)
  
  ### bootstrap sample and Re_boot
  count_boot <- 0
  all_boot_R_mean <- data.frame()
  
  while (count_boot < n_boot) {
    
    log_diff_boot <- Block_boot_overlap_func(log_diff, block_size)
    # plot(log_diff_boot)
    
    log_smoothed_data <- smoothed_incidence_data$log_loess
    
    # ts_boot <- exp(log_diff_boot + log_smoothed_data)  
    ts_boot <- exp(log_diff_boot + log_smoothed_data) - 1
    
    ts_boot[ts_boot<0] <- 0
    ts_boot <- round(ts_boot)
    
    data_obs_boot_temp <- smoothed_incidence_data[,c(1,2)]
    data_obs_boot_temp$observations <- ts_boot
    
    # par(mfrow=c(2,1))
    # plot(data_obs)
    # plot(ts_boot)
    
    ##### Calculate Re posterior based on the bootstrap sample
    allRe_posterior_boot <- Get_allRe_func(data_obs_boot_temp, 
                                           IncubationParams,OnsetToCountParams,
                                           smooth_para_deConv,
                                           timevarying = FALSE, 
                                           delay = 0)
                                      
    R_mean_boot <- allRe_posterior_boot %>% filter(variable == "R_mean") %>% dplyr::select(date, value)
    
    count_boot <- count_boot + 1
    all_boot_R_mean = bind_rows(all_boot_R_mean, R_mean_boot)
    
  } # end while
  
  ######################################################
  low_quan <- (1-alpha)/2
  high_quan <- 1-(1-alpha)/2
  
  Summarize_boot_R_mean <- all_boot_R_mean %>% group_by(date) %>% 
    summarize( sd = sd(value), 
               var = var(value), 
               mean = mean(value), .groups = "drop")
  
  ######################################################
  ###### Constrcut CI based on various ways
  
  # R_est based on the original data: used as middle point for some CI
  allRe_posterior_ori <- Get_allRe_func(data_obs, IncubationParams, OnsetToCountParams,
                                        smooth_para_deConv,
                                        timevarying = FALSE, 
                                        delay = 0)
  
  R_est_ori <- allRe_posterior_ori %>%
    filter(variable == "R_mean") %>%
    dplyr::select(date, value)
  
  # below might be not needed if the date would be definitely the same, think and ask.
  R_est_ori_common <- R_est_ori %>% filter(date >= max( min(R_est_ori$date), min( Summarize_boot_R_mean$date)) &
                                             date <= min( max(R_est_ori$date), max( Summarize_boot_R_mean$date)))
  
  Summarize_boot_R_mean_common <- Summarize_boot_R_mean %>% 
    filter(date >= max( min(R_est_ori$date), min( Summarize_boot_R_mean$date)) &
             date <= min( max(R_est_ori$date), max( Summarize_boot_R_mean$date)))
  
  #####
  R_est_combine <-  R_est_ori_common
  names(R_est_combine) <- c("date", "R_est_ori")
  R_est_combine$R_mean_mean <- Summarize_boot_R_mean_common$mean
  R_est_combine$R_mean_var <- Summarize_boot_R_mean_common$var        # estimated variance via bootstrap

  ### Normal CI with 'R_est_ori' as middle point and sd based on 'R_mean'
  CI_lower <- R_est_combine$R_est_ori - qnorm(high_quan)*Summarize_boot_R_mean_common$sd
  CI_upper <- R_est_combine$R_est_ori + qnorm(high_quan)*Summarize_boot_R_mean_common$sd
  CI_upper[CI_upper<0] <- CI_lower[CI_lower<0] <- 0
  R_est_combine$R_est_ori_l <- CI_lower
  R_est_combine$R_est_ori_u <- CI_upper
  
  ### Normal CI with 'R_mean_mean' as middle point and sd based on 'R_mean'
  CI_lower <- R_est_combine$R_mean_mean - qnorm(high_quan)*Summarize_boot_R_mean_common$sd
  CI_upper <- R_est_combine$R_mean_mean + qnorm(high_quan)*Summarize_boot_R_mean_common$sd
  CI_upper[CI_upper<0] <- CI_lower[CI_lower<0] <- 0
  R_est_combine$R_mean_mean_l <- CI_lower
  R_est_combine$R_mean_mean_u <- CI_upper
  
  return(R_est_combine)
}

#############################################################
### overlapping block bootstrap
Block_boot_overlap_func <- function(ts, block_size){
  
  # get the weekdays for each position of ts
  weekdays_index <- (1:length(ts)) %% 7
  weekdays_index[which(weekdays_index==0)] <- 7
  
  ts_boot <-c()
  last_day_index <- 7
  
  ###### get the ts_boot: make sure glue wrt the correct days
  while(length(ts_boot) < length(ts)){
    start_index <- sample(1:(length(ts)-block_size+1), 1)
    sampled_index <- start_index:(start_index+block_size-1)
    sampled_weekdays <- weekdays_index[sampled_index]
    
    # make sure the day related to the first sample is after the previous ts_boot
    first_day_index <- which(sampled_weekdays==last_day_index)[1] + 1
    ts_boot_index <- sampled_index[first_day_index:block_size]
    
    last_day_index <- tail(weekdays_index[ts_boot_index],1)
    
    ts_boot <- c(ts_boot, ts[ts_boot_index])
  }
  
  # take the same length as previous ts
  ts_boot <- ts_boot[1:length(ts)]
  
  return(ts_boot)
}