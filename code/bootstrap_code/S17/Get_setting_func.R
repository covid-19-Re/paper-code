### Get R_est, noise model and init_infection_vec based on real data

Get_setting_func <- function(dataType,
                             dateType,
                             country_val,
                             Region=NULL,
                             R_true_type,
                             noise_model_type,
                             smooth_para_resi, 
                             smooth_para_deConv,
                             n_boot){
  
  if(is.null(Region)) {Region=country_val}
  
  ### Load real data
  LoadPath <- paste0('newdata/', country_val, '-Data.rds')
  Data <- readRDS(LoadPath)
  
  data_obs <- Data %>%
    filter(data_type == dataType,
           date_type == dateType, 
           region == Region) %>%
    dplyr::select(date, value)
  
  data_obs <- data_obs %>% complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0))

  # plot(data_obs)
  
  ####### Get R_est, fitted noise model and initial vector from original data in order to simulate date
  IncubationParams <- getGammaParams(meanParam = 5.3, sdParam = 3.2)
  if(dataType == 'Confirmed cases') {OnsetToCountParams = getGammaParams(5.5, 3.8)}
  if(dataType == 'Hospitalized patients') {OnsetToCountParams = OnsetToCountParams = getGammaParams(5.1, 4.2)}
  if(dataType == 'Deaths') {OnsetToCountParams = OnsetToCountParams = getGammaParams(15, 6.9)}
  
  ############################ get R_est
  data_obs_observations <- data.frame("date"=data_obs$date, "observations"=data_obs$value)
  R_est_CI_combine <- CI_boot_func(data_obs_observations, alpha=0.95, n_boot, block_size=10, starting_threshold=20,
                                   IncubationParams, OnsetToCountParams,
                                   smooth_para_boot=smooth_para_resi, 
                                   smooth_para_deConv=smooth_para_deConv)
  
  # if(R_true_type=="step"){
  #   R_est <- data.frame("date"=R_est_CI_combine$date, "value"=R_est_CI_combine$R_mean_mean)
  #   
  #   R_step <- rep(c(1.5, 0.5, 1.5, 1, 0.7, 1, 1.3), each=length(R_est$value)/8)
  #   R_step <- c(R_step, rep(1, length(R_est$value)-length(R_step)))
  #   R_est$value <- R_step
  # }
  
  if(R_true_type=="smooth"){
    R_est <- data.frame("date"=R_est_CI_combine$date, "value"=R_est_CI_combine$R_mean_mean)
    
    dates <- R_est$date
    count_data <- R_est$value
    
    n_points <- length(unique(dates))
    sel_span <- 21 / n_points
    n_pad <- round(length(count_data) * sel_span * 0.5)
    c_data <- data.frame(value = count_data, date_num = as.numeric(dates))
    c_data.lo <- loess(value ~ date_num, data = c_data, span = sel_span, degree = 1)
    smoothed <- predict(c_data.lo)
    smoothed[smoothed < 0] <- 0
    
    # plot(R_est$value)
    # plot(smoothed)
    
    R_est$value <- smoothed
  }
  
  #plot(R_est, type='l', ylim=c(0,4))
   
  ######################################## get residual
  # get residual on log-scale
  data_obs_log <- data_obs
  data_obs_log$value <- log(data_obs_log$value + 1)
  
  names(data_obs_log) <- c("date", "log_value")
  
  smoothed_incidence_data <- data_obs_log %>%
    complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(log_value = 0)) %>%
    mutate(log_loess = getLOESSCases(dates = date, count_data = log_value, days_incl=smooth_para_resi),
           log_diff = log_value - log_loess)
  
  # take the same length as R_est (need noiseless observations has same length as simulate(noise_model))
  smoothed_incidence_data_cut <- smoothed_incidence_data %>% filter(date >= min(R_est$date))
  
  ### Plots
  # Smooth_df <- as.data.frame(smoothed_incidence_data)
  # gg_diff <- ggplot() + geom_point(data=Smooth_df, aes(x=date, y=log_diff))
  # gg_smooth <- ggplot() + geom_point(data=Smooth_df, aes(x=date, y=log_value)) + geom_line(data=Smooth_df, aes(x=date, y=log_loess), color="red")
  # gg_smooth_diff <- ggplot() + geom_point(data=Smooth_df, aes(x=log_loess, y=log_diff))
  # gg_smooth
  # gg_diff
  # gg_smooth_diff
  
  # get noise model on log_diff
  log_diff_ts <- ts(smoothed_incidence_data_cut$log_diff, frequency=7)
  
  if(noise_model_type=="CHE_confirmed_21_21"){
    fitted_noise_model <- Arima(log_diff_ts, order = c(3, 0, 4), seasonal = c(1, 0, 2),
                                include.mean = FALSE)
  }
  
  ### Plots: fit the model such that the acf and pacf of the residuals are withing the confidence bind as much as possible.
  # gg_acf <- ggAcf(fitted_noise_model$residuals)
  # gg_pacf <- ggPacf(fitted_noise_model$residuals)
  # gg_acf
  # gg_pacf

  ####
  estimatedInfections_originalData <- estimateInfectionTS(data_obs_observations, 
                                                          IncubationParams, OnsetToCountParams,
                                                          smooth_param = T,
                                                          smooth_para_deConv,
                                                          timevarying = FALSE, 
                                                          n_boot)   
  ### INIT INFECTION:
  init_infection_vec <- estimatedInfections_originalData %>%
    filter(data_type == 'infection_Simulated',
           region == 'Simulated') %>%
    group_by(date) %>%
    summarise(mean_val = round(mean(value)),
              .groups = 'drop') %>%
    filter(date <= min(R_est$date)) %>%
    pull(mean_val)
  
  return(list("R_est"=R_est, "fitted_noise_model"=fitted_noise_model, "init_infection_vec"=init_infection_vec, "data_obs"=data_obs))
}