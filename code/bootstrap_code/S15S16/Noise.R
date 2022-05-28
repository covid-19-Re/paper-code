source('LoadPackagesAndFunctions.R')

################## Load real data
Country_vec <- c("CHE", "CHN", "FRA", "NZL", "USA")
gg_acf_pacf_res_list <- list() 

for (i in 1:length(Country_vec)) {
  country_val <- Country_vec[i]
  Region <- Country_vec[i]
  dataType = 'Confirmed cases'
  dateType <- 'report'
  
  LoadPath <- paste0('newdata/', country_val, '-Data.rds')
  Data <- readRDS(LoadPath)
  
  data_obs <- Data %>%
    filter(data_type == dataType,
           date_type == dateType, 
           region == Region) %>%
    dplyr::select(date, value)
  
  data_obs <- data_obs %>% complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0))
  
  # get residual on log-scale
  data_obs_log <- data_obs
  data_obs_log$value <- log(data_obs_log$value + 1)
  
  names(data_obs_log) <- c("date", "log_value")
  
  smoothed_incidence_data <- data_obs_log %>%
    complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(log_value = 0)) %>%
    mutate(log_loess = getLOESSCases(dates = date, count_data = log_value, days_incl=21),
           log_diff = log_value - log_loess)
  
  # get noise model on log_diff
  original_residuals <- ts(smoothed_incidence_data$log_diff, frequency=7)
  
  ### Plots: Look at the acf and pacf plots: indicating that the residuals are auto-correlated; and the fitted noise model
  ### is better
  gg_acf_res <- ggAcf(original_residuals)
  gg_pacf_res <- ggPacf(original_residuals)
  gg_acf_pacf_res <- ggarrange(gg_acf_res, gg_pacf_res, ncol = 2, nrow = 1)
  gg_acf_pacf_res_list[[i]] <- annotate_figure(gg_acf_pacf_res, 
                                               top = text_grob(paste0("ACF and PACF plots of ", Country_vec[i])))
  
  if(country_val=="CHE"){
    gg_acf_res_swiss <- ggAcf(original_residuals)
    gg_pacf_res_swiss <- ggPacf(original_residuals)
    
    #arima_fit <- Arima(original_residuals, order = c(3, 0, 4), seasonal = c(1, 0, 2), include.mean = FALSE)
    arima_fit <- Arima(original_residuals, order = c(2, 0, 1), seasonal = c(0, 1, 1), include.mean = FALSE)
    summary(arima_fit)
    gg_acf_fit <- ggAcf(arima_fit$residuals)
    gg_pacf_fit <- ggPacf(arima_fit$residuals)
  }
}

gg_acf_pacf_res <- ggarrange(gg_acf_pacf_res_list[[1]], gg_acf_pacf_res_list[[2]],
                             gg_acf_pacf_res_list[[3]], gg_acf_pacf_res_list[[4]],
                             gg_acf_pacf_res_list[[5]],
                             ncol = 1, nrow = 5)
gg_acf_pacf_Swiss <- ggarrange(gg_acf_res_swiss, gg_pacf_res_swiss, gg_acf_fit, gg_pacf_fit,
                     ncol = 2, nrow = 2)

ggsave("S15S16/acf_pacf_res.pdf", plot=gg_acf_pacf_res, width = 8, height = 14)
ggsave("S15S16/gg_acf_pacf_Swiss.pdf", plot=gg_acf_pacf_Swiss, width = 10, height = 10)
