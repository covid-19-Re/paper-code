source('LoadPackagesAndFunctions.R')
source('S17/Get_setting_func.R')

#####
country_val <- 'CHE'
Region <- "CHE"
dataType = 'Confirmed cases'
dateType <- 'report'

### Get R_est, noise model and init_infection_vec based on real data
set.seed(111)

noise_model_type = "CHE_confirmed_21_21"
R_true_type="smooth"

setting_list <- Get_setting_func(dataType, dateType, country_val, Region, R_true_type, noise_model_type, 
                                 smooth_para_resi=21, smooth_para_deConv=21, n_boot=100)
R_est <- setting_list$R_est
fitted_noise_model <- setting_list$fitted_noise_model
init_infection_vec <- setting_list$init_infection_vec
data_obs <- setting_list$data_obs

noise=list('fitted_noise_model'=fitted_noise_model)

#####
RealData_plot <- ggplot() + geom_point(data=data_obs, aes(x=date, y=value)) +
  labs(x = 'date', y = 'observation from real data') +
  ylim(c(0, max(data_obs[,2])*3))

#############################################################
####### simulated data
IncubationParams <- getGammaParams(meanParam = 5.3, sdParam = 3.2)
if(dataType == 'Confirmed cases') {OnsetToCountParams = getGammaParams(5.5, 3.8)}

list_noiseless <- list_noise <- plot_noiseless <- plot_noise <- list()
for (i in 1:4) {
  print(i)
  simulation <- Para_model_func(R_hat = R_est$value,
                                R_dates = R_est$date,
                                Obs_dates = data_obs$date,
                                IncubationParams, OnsetToCountParams,
                                noise,
                                init_infection = init_infection_vec)
  list_noise[[i]] <- simulation[,c(5,3)]
  list_noiseless[[i]] <- simulation[,c(5,4)]
  
  plot_noiseless[[i]] <- ggplot() + geom_point(data=list_noiseless[[i]], aes(x=date, y=Observations_noiseless)) +
    labs(x = 'date', y = 'simulated noiseless observation') +
    ylim(c(0, max(data_obs[,2])*3))
  
  plot_noise[[i]] <- ggplot() + geom_point(data=list_noise[[i]], aes(x=date, y=observations)) +
    labs(x = 'date', y = 'simulated observation') +
    ylim(c(0, max(data_obs[,2])*3))
}

library(ggpubr)

comb_noiseless <- ggarrange(RealData_plot, plot_noiseless[[1]], plot_noiseless[[2]], 
                            plot_noiseless[[3]],  plot_noiseless[[4]], 
                         ncol = 5, nrow = 1)
comb_noiseless <- annotate_figure(comb_noiseless, top = text_grob("Simulated noiseless observations"))
comb_noise <- ggarrange(RealData_plot, plot_noise[[1]], plot_noise[[2]], 
                        plot_noise[[3]], plot_noise[[4]], 
                         ncol = 5, nrow = 1)
comb_noise <- annotate_figure(comb_noise, top = text_grob("Simulated observation"))

Comb_simu_plot <-  ggarrange(comb_noiseless, comb_noise,
                             ncol = 2, nrow = 1)
Comb_simu_plot

ggsave("S17/SimulatedData.png", plot=Comb_simu_plot, width = 15, height = 8)

##############################
# Redo figure
noise_df = bind_rows(list_noise, .id = 'rep') %>%
  rename(value = observations)
noiseless_df = bind_rows(list_noiseless, .id = 'rep') %>%
  rename(value = Observations_noiseless)

raw_df <- data_obs %>%
  mutate(rep = '0')

full_df <- bind_rows(raw_df, raw_df, noise_df, noiseless_df, .id = 'data_type') %>%
  mutate(data_type = case_when(
    data_type %in% c('1', '3') ~ 'With noise',
    data_type %in% c('2', '4') ~ 'Noiseless'),
    col = ifelse(rep == 0, 'data', 'simulation'))

ggplot(full_df) + 
  geom_point(aes(x=date, y=value, colour = col), show.legend = F) +
  labs(x = 'Date', y = 'Observations') +
  scale_color_manual(values = c('red', 'black')) +
  facet_grid(rows = vars(rep), cols = vars(data_type)) + 
  theme_minimal() +
  theme(strip.text.y = element_blank(),
        text = element_text(size = 20),
        strip.text.x = element_text(size = 25))

ggsave("S17/SimulatedData.pdf", height = 15, width = 12)

