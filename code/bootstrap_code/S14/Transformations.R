library(ggpubr)
library(ggplot2)
source('LoadPackagesAndFunctions.R')

Country_vec <- c("CHE", "CHN", "FRA", "NZL", "USA")

Trans_met_vec <- c("Square-root-transformation", "Log-transformation")

Reply_plot_country_list <- list()
set.seed(2021)
for (country_index in 1:length(Country_vec)) {
  
  print(country_index)
  
  ################## Load real data
  country_val <- Country_vec[country_index]
  Region <- Country_vec[country_index]
  dataType = 'Confirmed cases'
  #dataType = 'Hospitalized patients'
  #dataType = 'Deaths'
  dateType <- 'report'
  
  LoadPath <- paste0('newdata/', country_val, '-Data.rds')
  Data <- readRDS(LoadPath)
  
  data_obs <- Data %>%
    filter(data_type == dataType,
           date_type == dateType, 
           region == Region) %>%
    dplyr::select(date, value)
  
  data_obs <- data_obs %>% complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0))
  
  Reply_plot_list <- list()
  for (trans_Resi_index in 1:length(Trans_met_vec)) {
    
    trans_resi_method <- Trans_met_vec[trans_Resi_index]
    
    ########## Transform data
    data_obs_trans <- data_obs
    
    if(trans_resi_method == "Square-root-transformation"){
      data_obs_trans$trans_value <- sqrt(data_obs_trans$value)
      
      smoothed_incidence_data <- data_obs_trans %>%
        complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>%
        mutate(trans_loess = getLOESSCases(dates = date, count_data = trans_value, days_incl=21),
               trans_diff = trans_value - trans_loess)  
      
      trans_diff <- smoothed_incidence_data$trans_diff
      trans_loess <- smoothed_incidence_data$trans_loess
    }
    
    if(trans_resi_method == "Log-transformation"){
      data_obs_trans$trans_value <- log(data_obs_trans$value + 1)
      
      smoothed_incidence_data <- data_obs_trans %>%
        complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>%
        mutate(trans_loess = getLOESSCases(dates = date, count_data = trans_value, days_incl=21),
               trans_diff = trans_value - trans_loess)
      
      trans_diff <- smoothed_incidence_data$trans_diff
      trans_loess <- smoothed_incidence_data$trans_loess
    }
    
    ########## Plots
    Smooth_df <- as.data.frame(smoothed_incidence_data)
    gg_smooth <- ggplot() + geom_point(data=Smooth_df, aes(x=date, y=trans_value)) + 
      geom_line(data=Smooth_df, aes(x=date, y=trans_loess), color="red") + 
      labs(x = 'date', y = 'transformed observation') 
    gg_smooth_diff <- ggplot() + geom_point(data=Smooth_df, aes(x=trans_loess, y=trans_diff)) + 
      labs(x = 'smoothed value', y = 'residual') 
    
    Resi_plot <- ggarrange(gg_smooth, gg_smooth_diff, ncol = 2, nrow = 1)
    
    Reply_plot_list[[trans_Resi_index]] <- annotate_figure(Resi_plot, top = text_grob(trans_resi_method))
  }
  
  #####
  Comb_reply_plot <-  ggarrange(Reply_plot_list[[1]], Reply_plot_list[[2]],
                               ncol = 2, nrow = 1)
  Reply_plot_country_list[[country_index]] <- annotate_figure(Comb_reply_plot, top = text_grob(country_val))
}

all_reply_plot <-  ggarrange(Reply_plot_country_list[[1]], Reply_plot_country_list[[2]], Reply_plot_country_list[[3]],
                             Reply_plot_country_list[[4]], Reply_plot_country_list[[5]],
                             ncol = 1, nrow = 5)

all_reply_plot

ggsave(paste0("S14/Transformations-all.pdf"), plot=all_reply_plot, width = 10, height = 15)
