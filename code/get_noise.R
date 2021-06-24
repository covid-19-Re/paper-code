### Get noise model based on real data

library(forecast)

get_noise <- function(country_val = 'CHE', 
                      dateType = 'report', 
                      dataType = 'Confirmed cases', 
                      smooth_days = 21){
  #country_val = 'CHN'
  
  ### Load real data
  LoadPath <- paste0('../data/', country_val, '-Data.rds')
  Data <- readRDS(LoadPath)
  
  data_obs <- Data %>%
    filter(data_type == dataType,
           date_type == dateType, 
           region == country_val) %>%
    dplyr::select(date, value)
  
  data_obs <- data_obs %>% complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0))
  
  
  ######################################## get residual
  # get residual on log-scale
  data_obs_log <- data_obs
  data_obs_log$value <- log(data_obs_log$value + 1)
  
  names(data_obs_log) <- c("date", "log_value")
  
  smoothed_incidence_data <- data_obs_log %>%
    complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(log_value = 0)) %>%
    mutate(log_loess = getLOESSCases(dates = date, count_data = log_value, days_incl=smooth_days),
           log_diff = log_value - log_loess)
  
  # get noise model on log_diff
  log_diff_ts <- ts(smoothed_incidence_data$log_diff, frequency=7)
  
  ###### Test section #####
  #fitted_noise_model <- auto.arima(log_diff_ts)
  # fitted_noise_model <- Arima(log_diff_ts, order = c(1, 0, 1), seasonal = c(0, 0, 0),
  #                             include.mean = FALSE)
  # fitted_noise_model
  # ggAcf(fitted_noise_model$residuals)
  # ggPacf(fitted_noise_model$residuals)
  
  ###### Set values #####
  if(country_val=="CHE"){
    fitted_noise_model <- Arima(log_diff_ts, order = c(2, 0, 1), seasonal = c(0, 1, 1),
                                include.mean = FALSE)
  } else if (country_val =="USA"){
    fitted_noise_model <- Arima(log_diff_ts, order = c(4, 0, 0), seasonal = c(0, 0, 0),
                                include.mean = FALSE)
  } else if (country_val == 'NZL'){
    fitted_noise_model <- Arima(log_diff_ts, order = c(4, 0, 1), seasonal = c(1, 0, 0),
                                include.mean = FALSE)
  } else if (country_val == 'FRA'){
    fitted_noise_model <- Arima(log_diff_ts, order = c(0, 0, 6), seasonal = c(0, 1, 1),
          include.mean = FALSE)
  } else if (country_val == 'CHN'){
    fitted_noise_model <- Arima(log_diff_ts, order = c(1, 0, 1), seasonal = c(0, 0, 0),
                                include.mean = FALSE)
  }
  
  return(fitted_noise_model)
}
