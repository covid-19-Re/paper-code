###########################################################
## paperSimulations.R
## author: J.S. Huisman
## Last update: March 2021
###########################################################

# Note: it takes quite long to calculate each of 
# the simulation grids, so these commands have been commented out. 
# Instead the files with simulation results are included in the 
# git repository. 

# Also the plotting takes some time because boxplots across
# all simulations are drawn for every timepoint. (expect ~1 min per plot)
# To speed up the plotting, one can choose to replace the boxplots
# by ribbons showing the mean instead.

###########################################################

library(tidyverse)
library(lubridate)

library(ggplot2)
library(patchwork)
library(viridis)

theme_set(theme_minimal() +
            theme(
              strip.text = element_text(size=25),
              axis.text= element_text(size=20),
              axis.title =  element_text(size=25),
              legend.text= element_text(size=20),
              legend.title= element_text(size=25)
            ))

###########################################################
simDir = '../simulations'

source('generateSimulations.R')
source('get_noise.R')

plot_path = '../figures'
###########################################################
## Necessary Functions #####
getCountParams <- function(obs_type, misdelay = 0){
  switch(obs_type,
         incubation = getGammaParams(5.3, 3.2),
         zero = list(shape = 0, scale = 0),
         death = getGammaParams(15.0 + misdelay, 6.9),
         hospitalisation = getGammaParams(5.14 + misdelay, 4.2),
         confirmed = getGammaParams(5.5 + misdelay, 3.8) )
}

###########################################################
# per scenario at least 100 simulations; 
# ensemble of 100 bootstrap values in the estimation
n_sim = 100

###########################################################
## In Parallel ####
library("parallel")

# we parallelise the computation of the simulation ensemble

## Parallel  Functions ####
sim_Re <- function(i){
  
  simulation <- simulateTS(param_list$shift_times, param_list$R_levels,
                           getCountParams("incubation"),
                           getCountParams(param_list$data_type),
                           init_infection = param_list$init_I,
                           noise = param_list$noise,
                           timevarying = param_list$timevar_sim)
  simulation['replicate'] = i
  
  return(simulation)
}

estimate_Re <- function(i, simulation){

  simulation_i <- simulation %>% filter(replicate == i)
  
  if(param_list$epiestim){
    param_list$n_boot <- 0
  }

  estimatedInfections <- estimateInfectionTS(simulation_i,
                                             getCountParams("incubation"),
                                             getCountParams(as.character(param_list$data_type), 
                                                            misdelay = param_list$delay),
                                             smooth_param = param_list$smooth,
                                             fixed_shift = param_list$fixed_shift,
                                             n_boot = param_list$n_boot,
                                             timevarying = param_list$timevar_est)

  estimatedInfections <- estimatedInfections %>%
    mutate(local_infection = T)
  
  estimatedRe <- estimateReTS(estimatedInfections)

  if(param_list$epiestim | param_list$n_boot == 0){
    summarisedRe <- as_tibble(estimatedRe) %>%
      pivot_wider(names_from = "variable", values_from = "value") %>%
      rename(median_R_mean = R_mean,
             median_R_highHPD = R_highHPD,
             median_R_lowHPD = R_lowHPD)
  } else {
    summarisedRe <- cleanCountryReEstimate(estimatedRe, method = 'bootstrap',
                                           rename_types = F) %>%
      dplyr::select(-c(country, region, source, data_type, estimate_type)) %>%
      mutate(replicate = i)
  }
  
  return(summarisedRe)
}

parallel_simulation <- function(param_list,
                                ensemble_size = 100, number_of_cores = 8){

  num.cores <- number_of_cores - 1

  # Initialise the cluster
  cl <- makeCluster(num.cores, type="FORK")
  # Set a seed for the random number generator, which then seeds all new jobs
  # to ensure no random sequences are used on multiple cores
  clusterSetRNGStream(cl, 150)
  clusterExport(cl=cl, varlist = c("param_list"),
                envir=environment())

  # Re simulation
  sim_result <- parLapply(cl, seq_len(ensemble_size), sim_Re)
  sim_result <- bind_rows(sim_result)
  
  # Re estimation
  Re_result <- parLapply(cl, seq_len(ensemble_size), estimate_Re, sim_result)
  Re_result <- bind_rows(Re_result)
  
  # Finally, close the cluster we created
  stopCluster(cl)
  # clean up memory because otherwise cluster connection problems can occur
  gc()
  return(list(sim = sim_result, Re = Re_result))
}

###########################################################
## Generally useful function ########
compute_exp_metrics <- function(exp_result){

  rmse_df <- exp_result %>%
    mutate(sre = ((median_R_mean - Re)/Re)^2 ) %>%
    mutate(se = (median_R_mean - Re)^2 ) %>%
    group_by(date, exp_name) %>%
    summarize(
      rmse = sqrt(sum(se)/n()),
      rmsre = sqrt(sum(sre)/n()),
      .groups = "drop"
    ) 
  
  emp_cov <- exp_result %>%
    group_by(date, exp_name) %>%
    summarize(
      empirical_cov = sum(in_CI)/n(),
      .groups = "drop"
    ) 
  
  above_1_df <- exp_result %>%
    group_by(Re, exp_name) %>%
    summarize(
      above_1 = sum(above_1)/n(),
      below_1 = sum(below_1)/n(),
      .groups = "drop"
    ) 
  return(list(rmse = rmse_df, cov = emp_cov, above_1 = above_1_df))
}

######################################################################################################################
# Run 1 scenario 1000 times (to test; not included in paper) ####
param_list <- list(
shift_times = c(0, 20, 50, 60, 80, 90, 100, 160, 180, 190, 210, 250, 270, 290, 310, 340, 365, 374),
R_levels = c(2.5, 0.5, 0.7, 1.1, 1.3, 0.9, 1.1, 0.8, 1.0),
init_I = 100, noise = list("fitted_noise_model" = get_noise(country_val = "CHE")),
data_type = "confirmed", n_boot = 100,
smooth = TRUE,
epiestim = FALSE,
delay = 0,
fixed_shift = FALSE,
timevar_sim = FALSE,
timevar_est = FALSE)

tictoc::tic()
par_Re_est <- parallel_simulation(param_list, ensemble_size = 100, number_of_cores = 8)
tictoc::toc()
# 3095.391 sec = 52 min; for 1000 ensemble
# 662.532 sec  = 11 min; for 100 ensemble

## Very rough plotting example
ggplot(par_Re_est$sim) +
  geom_line(aes(x = date, y = infections, colour = replicate))

ggplot(par_Re_est$Re) +
  geom_line(aes(x = date, y = median_R_mean, colour = replicate))

######################################################################################################################
## Fig 1, S4: one scenario, various noise
config_df <- expand.grid(init = c(100),
                         noise_country = c('CHE', 'FRA', 'USA', 'CHN', 'NZL')
                         )

## Commented out because it takes long to run! ##
# for (i in 1:nrow(config_df)){
#   param_list <- list(
#     shift_times = c(0, 20, 50, 60, 80, 90, 100, 160, 180, 190, 210, 250, 270, 290, 310, 340, 365, 374),
#     R_levels = c(2.5, 0.5, 0.7, 1.1, 1.3, 0.9, 1.1, 0.8, 1.0),
#     init_I = config_df[i, 'init'], 
#     noise = list(fitted_noise_model = get_noise(country_val = config_df[i, 'noise_country']) ),
#     data_type = "confirmed", n_boot = 100,
#     smooth = TRUE,
#     delay = 0,
#     fixed_shift = FALSE,
#     timevar_sim = FALSE,
#     timevar_est = FALSE)
#   
#   par_Re_est <- parallel_simulation(param_list, ensemble_size = 100, number_of_cores = 8)
#   
#   dir.create(file.path(simDir, 'Fig1_noise'))
#   config_result <- par_Re_est$Re %>%
#     left_join(par_Re_est$sim, by = c("replicate", "date")) %>%
#     mutate(in_CI = (Re > median_R_lowHPD) & (Re < median_R_highHPD),
#            above_1 = (median_R_lowHPD > 1),
#            below_1 = (median_R_highHPD < 1)) %>%
#     bind_cols(config_df[i, ]) %>%
#     mutate(exp_name = paste0(init, '_', noise_country))
#   
#   write_csv(config_result, paste0(simDir, '/Fig1_noise/all_',
#                                   config_df[i, 'noise_country'], '.csv'))
#   
# }

## Reading in #####

exp_result <- data.frame()
for (i in 1:nrow(config_df)){
  config_result <- read_csv(paste0(simDir, '/Fig1_noise/all_', 
                         config_df[i, 'noise_country'],'.csv'))
  
  exp_result <- bind_rows(exp_result, config_result)
}


## Plotting Supplementary Fig S4; Same Scenario with different Noise ####

plot_Re_est <- exp_result %>%
  pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
               names_to = 'name', values_to = 'value') %>%
  mutate(grouping = paste0(date, name),
         plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

ggplot() +
  geom_boxplot(data = plot_Re_est, 
               aes(x = date, y = value, group = grouping, colour = plot_colour), 
               position = 'identity', outlier.shape = NA )+
  geom_line(data = plot_Re_est %>% filter(replicate == 1), aes(x = date, y = Re), size = 1, colour = 'black') +
  facet_wrap(vars(noise_country), ncol = 1)+
  coord_cartesian(ylim = c(0,5)) +
  labs(x = 'Date', y = 'Re') +
  scale_colour_manual(values = c(viridis(2)[c(2,1)]), 
                      labels = c('Point estimate', 'CI'),
                      breaks = c('Point estimate', 'CI'),
                      name = 'Variable') +
  theme(legend.position = 'bottom')

ggsave(paste0(plot_path, '/', 'Supp_Re_noise.png'), height = 15, width = 12)

## Plotting Fig 1 ####

plot_Re_est <- exp_result %>%
  filter(exp_name == '100_CHE') %>%
  pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
               names_to = 'name', values_to = 'value') %>%
  mutate(grouping = paste0(date, name),
         plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

exp_metrics <- compute_exp_metrics(exp_result)

Re_plot <- ggplot()+
  geom_boxplot(data = plot_Re_est,
               aes(x = date, y = value, group = grouping, colour = plot_colour), 
               position = 'identity', outlier.shape = NA ) +
  geom_line(data = plot_Re_est %>% filter(replicate == 1), aes(x = date, y = Re), size = 1, colour = 'black') +
  labs(x = 'Date', y = 'Re') +
  scale_colour_manual(values = c(viridis(2)[c(2,1)]), 
                      labels = c('Point estimate', 'CI'),
                      breaks = c('Point estimate', 'CI'),
                      name = 'Variable') +
  theme(legend.position = 'bottom')

Cov_plot <- ggplot(exp_metrics$cov %>%
                     filter(exp_name == '100_CHE') ) +
  geom_point(aes(x = date, y = empirical_cov)) +
  geom_hline(yintercept = 0.95, colour = 'red', linetype = 'dashed', alpha = 0.5) +
  labs(x = 'Date', y = 'Empirical Coverage') 

rmse_plot <- ggplot(exp_metrics$rmse %>%
                      filter(exp_name == '100_CHE') ) +
  geom_point(aes(x = date, y = rmse)) +
  labs(x = 'Date', y = 'RMSE') 

##
cross_df <- exp_metrics$above_1 %>%
  filter(exp_name == '100_CHE') %>%
  pivot_longer(cols = c('above_1', 'below_1'))

cross_plot <- ggplot(cross_df ) +
  geom_point(aes(x = Re, y = value, colour = name)) +
  scale_colour_manual(values = c(viridis(2)), 
                      labels = c('Signif. above 1', 'Signif. below 1'),
                      breaks = c('above_1', 'below_1'),
                      name = 'Variable') +
  labs(x = 'Re', y = 'Fraction of simulations') +
  theme(legend.position = 'bottom')

##### this uses patchwork
  Re_plot  + 
  Cov_plot + 
  rmse_plot + 
  cross_plot +
  plot_layout(ncol=1) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25),
        plot.tag.position = 'topright')

ggsave(paste0(plot_path, '/', 'Fig1.png'), height = 15, width = 12)



###########################################################
### Fig1: various scenarios, one noise #####
config_df <- expand.grid(R1 = c(2.5, 3.0, 3.5),
                         init = c(10, 100, 1000),
                         t3 = c(7, 14, 28),
                         noise_country = 'CHE')

## Commented out because it takes long to run! ##
# for (i in 1:nrow(config_df)){
#   t3 <- config_df[i, 't3']
#   
#   param_list <- list(
#     shift_times = c(0, 20, 20+t3, 60, 60+t3, 90, 90+t3, 160, 160+t3, 190, 
#                     190+t3, 250, 250+t3, 290, 290+t3, 340, 340+t3, 374),
#     R_levels = c(config_df[i, 'R1'], 0.5, 0.7, 1.1, 1.3, 0.9, 1.1, 0.8, 1.0),
#     init_I = config_df[i, 'init'],
#     noise = list(fitted_noise_model = get_noise(country_val = config_df[i, 'noise_country']) ),
#     data_type = "confirmed", n_boot = 100,
#     smooth = TRUE,
#     delay = 0,
#     fixed_shift = FALSE,
#     timevar_sim = FALSE,
#     timevar_est = FALSE)
# 
#   par_Re_est <- parallel_simulation(param_list, ensemble_size = 100, number_of_cores = 8)
# 
#   dir.create(file.path(simDir, 'Fig1_noise'))
#   config_result <- par_Re_est$Re %>%
#     left_join(par_Re_est$sim, by = c("replicate", "date")) %>%
#     mutate(in_CI = (Re > median_R_lowHPD) & (Re < median_R_highHPD),
#            above_1 = (median_R_lowHPD > 1),
#            below_1 = (median_R_highHPD < 1)) %>%
#     bind_cols(config_df[i, ]) %>%
#     mutate(exp_name = paste0(R1, '_', init, '_', t3))
# 
#   write_csv(config_result, paste0(simDir, '/Fig1_noise/all_',
#                                   config_df[i, 'R1'],'_',
#                                   config_df[i, 'init'], '_',
#                                   config_df[i, 't3'], '.csv'))
# 
# }

## Reading in #####

exp_result <- data.frame()
for (i in 1:nrow(config_df)){
  config_result <- read_csv(paste0(simDir, '/Fig1_noise/all_', 
                                   config_df[i, 'R1'],'_',
                                   config_df[i, 'init'], '_', 
                                   config_df[i, 't3'], '.csv'))

  exp_result <- bind_rows(exp_result, config_result)
}

exp_result <- exp_result %>%
  mutate(above_1 = (median_R_lowHPD > 1),
         below_1 = (median_R_highHPD < 1) )

#####################################
## Fig S2: Many Scenarios - one Noise ####

plot_Re_est <- exp_result %>%
  pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
               names_to = 'name', values_to = 'value') %>%
  mutate(grouping = paste0(date, name),
         plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

exp_metrics <- compute_exp_metrics(exp_result)


Re_plot <- ggplot() +
  geom_boxplot(data = plot_Re_est %>% filter(init == 100, R1 == 3), 
               aes(x = date, y = value, group = grouping, colour = plot_colour), 
               position = 'identity', outlier.shape = NA )+
  geom_line(data = plot_Re_est %>% filter(replicate == 1,init == 100, R1 == 3), 
            aes(x = date, y = Re), size = 1, colour = 'black') +
  coord_cartesian(ylim =c(0,5)) +
  facet_grid(cols = vars(t3), rows = vars(R1))+
  labs(x = 'Date', y = 'Re') +
  scale_colour_manual(values = c(viridis(2)[c(2,1)]), 
                      labels = c('Point estimate', 'CI'),
                      breaks = c('Point estimate', 'CI'),
                      name = 'Variable') +
  theme(legend.position = 'bottom')

Cov_plot <- ggplot(exp_metrics$cov %>%
                     separate(exp_name, into = c("R1", "init", "t3"), sep = "_")%>%
                     filter(init == 100, R1 == 3) %>%
                     mutate(t3 = factor(t3, levels = c(7, 14, 28)))) +
  geom_point(aes(x = date, y = empirical_cov)) +
  facet_grid(cols = vars(t3), rows = vars(R1))+
  geom_hline(yintercept = 0.95, colour = 'red', linetype = 'dashed', alpha = 0.5) +
  labs(x = 'Date', y = 'Empirical Coverage') 

RMSE_plot <- ggplot(exp_metrics$rmse %>%
                     separate(exp_name, into = c("R1", "init", "t3"), sep = "_")%>%
                     filter(init == 100, R1 == 3) %>%
                     mutate(t3 = factor(t3, levels = c(7, 14, 28)))) +
  geom_point(aes(x = date, y = rmse)) +
  facet_grid(cols = vars(t3), rows = vars(R1))+
  labs(x = 'Date', y = 'RMSE') 

cross_df <- exp_metrics$above_1 %>%
  separate(exp_name, into = c("R1", "init", "t3"), sep = "_")%>%
  filter(init == 100, R1 == 3) %>%
  pivot_longer(cols = c('above_1', 'below_1')) %>%
  mutate(t3 = factor(t3, levels = c(7, 14, 28)))

cross_plot <- ggplot(cross_df ) +
  geom_point(aes(x = Re, y = value, colour = name)) +
  #geom_line(aes(x = Re, y = value, colour = name), linetype = 'dashed') +
  facet_grid(cols = vars(t3), rows = vars(R1))+
  scale_colour_manual(values = c(viridis(2)), 
                      labels = c('Signif. above 1', 'Signif. below 1'),
                      breaks = c('above_1', 'below_1'),
                      name = 'Variable') +
  labs(x = 'Re', y = 'Fraction of simulations') +
  theme(legend.position = 'bottom')


#####

combi_plot <- Re_plot  + Cov_plot  + RMSE_plot + cross_plot +
  plot_layout(ncol=1) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25),
        strip.text.y = element_blank())
combi_plot

ggsave(plot = combi_plot, paste0(plot_path, '/', 'Supp_Re_scenarios.png'), height = 20, width = 17)


###########################################################
##### Fig S3: One scenario, various noise, but no smoothing ######
config_df <- expand.grid(init = c(100),
                         noise_country = c('CHE', 'FRA', 'USA', 'CHN', 'NZL'),
                         smooth = FALSE)

## Commented out because it takes long to run! ##
# for (i in 1:nrow(config_df)){
#   param_list <- list(
#     shift_times = c(0, 20, 50, 60, 80, 90, 100, 160, 180, 190, 210, 250, 270, 290, 310, 340, 365, 374),
#     R_levels = c(2.5, 0.5, 0.7, 1.1, 1.3, 0.9, 1.1, 0.8, 1.0),
#     init_I = config_df[i, 'init'],
#     noise = list(fitted_noise_model = get_noise(country_val = config_df[i, 'noise_country']) ),
#     data_type = "confirmed", n_boot = 100,
#     smooth = config_df[i, 'smooth'],
#     delay = 0,
#     fixed_shift = FALSE,
#     timevar_sim = FALSE,
#     timevar_est = FALSE)
# 
#   par_Re_est <- parallel_simulation(param_list, ensemble_size = 100, number_of_cores = 8)
# 
#   dir.create(file.path(simDir, 'Fig1_noise_no_smooth'))
#   config_result <- par_Re_est$Re %>%
#     left_join(par_Re_est$sim, by = c("replicate", "date")) %>%
#     mutate(in_CI = (Re > median_R_lowHPD) & (Re < median_R_highHPD),
#            above_1 = (median_R_mean > 1) ) %>%
#     bind_cols(config_df[i, ]) %>%
#     mutate(exp_name = noise_country)
# 
#   write_csv(config_result, paste0(simDir, '/Fig1_noise_no_smooth/all_',
#                                   config_df[i, 'noise_country'],'.csv'))
# 
# }

## Reading in #####

exp_result <- data.frame()
for (i in 1:nrow(config_df)){
  config_result <- read_csv(paste0(simDir, '/Fig1_noise_no_smooth/all_', 
                                   config_df[i, 'noise_country'], '.csv'))
  
  exp_result <- bind_rows(exp_result, config_result)
}

exp_result <- exp_result %>%
  mutate(above_1 = (median_R_lowHPD > 1),
         below_1 = (median_R_highHPD < 1) )

## Plotting #####
plot_Re_est <- exp_result %>%
  pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
               names_to = 'name', values_to = 'value') %>%
  mutate(grouping = paste0(date, name),
         plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

ggplot() +
  geom_boxplot(data = plot_Re_est, 
               aes(x = date, y = value, group = grouping, colour = plot_colour), 
               position = 'identity', outlier.shape = NA )+
  geom_line(data = plot_Re_est %>% filter(replicate == 1), 
            aes(x = date, y = Re), size = 1, colour = 'black') +
  facet_wrap(vars(noise_country), ncol = 1 )+
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = 'Date', y = 'Re') +
  scale_colour_manual(values = c(viridis(2)[c(2,1)]), 
                      labels = c('Point estimate', 'CI'),
                      breaks = c('Point estimate', 'CI'),
                      name = 'Variable') +
  theme(legend.position = 'bottom')

ggsave(paste0(plot_path, '/', 'Supp_Re_noise_non_smooth.png'), height = 14, width = 12)

######################################################################################################################
##### Fig. S5: Dependence on population size ######
config_df <- expand.grid(R1 = c(0.8, 1, 1.5),
                         init = c(10, 100, 1000, 5000, 10000),
                         epiestim = FALSE)

## Commented out because it takes long to run! ##
# for (i in 1:nrow(config_df)){
#   param_list <- list(
#     shift_times = c(0, 100),
#     R_levels = c(config_df[i, 'R1']),
#     init_I = config_df[i, 'init'],
#     noise = list(fitted_noise_model = get_noise(country_val = "CHE") ),
#     data_type = "confirmed", n_boot = 100,
#     smooth = TRUE,
#     delay = 0,
#     fixed_shift = FALSE,
#     timevar_sim = FALSE,
#     timevar_est = FALSE,
#     epiestim = config_df[i, 'epiestim'])
#   
#   par_Re_est <- parallel_simulation(param_list, ensemble_size = 100, number_of_cores = 7)
#   
#   dir.create(file.path(simDir, 'Pop_size'))
#   config_result <- par_Re_est$Re %>%
#     left_join(par_Re_est$sim, by = c("replicate", "date")) %>%
#     mutate(in_CI = (Re > median_R_lowHPD) & (Re < median_R_highHPD),
#            above_1 = (median_R_lowHPD > 1),
#            below_1 = (median_R_highHPD < 1) ) %>%
#     bind_cols(config_df[i, ]) %>%
#     mutate(exp_name = paste0(R1, '_', init))
#   
#   if (config_df[i, 'epiestim']){
#     write_csv(config_result, paste0(simDir, '/Pop_size/all_',
#                                     config_df[i, 'R1'],'_',
#                                     config_df[i, 'init'], '_',
#                                     config_df[i, 'epiestim'], '.csv'))
#   } else{
#     write_csv(config_result, paste0(simDir, '/Pop_size/all_',
#                                     config_df[i, 'R1'],'_',
#                                     config_df[i, 'init'],'.csv'))
#   }
#   
# }

## Reading in #####

exp_result <- data.frame()
for (i in 1:nrow(config_df)){

  if (config_df[i, 'epiestim']){
    config_result <- read_csv(paste0(simDir, '/Pop_size/all_',
                                    config_df[i, 'R1'],'_',
                                    config_df[i, 'init'], '_',
                                    config_df[i, 'epiestim'], '.csv'))
  } else{
    config_result <- read_csv(paste0(simDir, '/Pop_size/all_',
                                    config_df[i, 'R1'],'_',
                                    config_df[i, 'init'],'.csv'))
  }
  
  exp_result <- bind_rows(exp_result, config_result)
}

exp_result <- exp_result %>%
  mutate(in_CI = (R1 >= median_R_lowHPD) & (R1 <= median_R_highHPD))

## Plotting #####
plot_Re_est <- exp_result %>%
  pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
               names_to = 'name', values_to = 'value') %>%
  mutate(grouping = paste0(date, name),
         plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI') )

Re_plot <- ggplot() +
  geom_boxplot(data = plot_Re_est , 
               aes(x = date, y = value, group = grouping, colour = plot_colour), 
               position = 'identity', outlier.shape = NA )+
  geom_line(data = plot_Re_est %>% filter(replicate %in% c(0, 1)), 
            aes(x = date, y = R1), size = 1, colour = 'black') +
  facet_grid(rows = vars(R1), cols = vars(init), scale = 'free_y' )+
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = 'Date', y = 'Re') +
  scale_colour_manual(values = c(viridis(2)[c(2,1)]), 
                      breaks = c('Point estimate', 'CI'),
                      labels = c('Point estimate', 'CI'),
                      name = 'Variable') +
  theme(legend.position = 'bottom',
        panel.spacing = unit(2, 'lines'))
Re_plot

ggsave(paste0(plot_path, '/', 'Supp_Re_pop_size.png'), height = 10, width = 18)
#ggsave(paste0(plot_path, '/', 'Supp_Re_pop_size_epiestim.png'), height = 10, width = 18)

exp_metrics <- compute_exp_metrics(exp_result)

Cov_plot <- ggplot(exp_metrics$cov %>%
                  separate(exp_name, into = c('R1', 'init'), sep = '_') %>%
                  mutate(init = factor(init, levels = c(10, 100, 1000, 5000, 10000)) ) ) +
  geom_point(aes(x = date, y = empirical_cov)) +
  facet_grid(rows = vars(R1), cols = vars(init), scale = 'free_y' )+
  geom_hline(yintercept = 0.95, colour = 'red', linetype = 'dashed', alpha = 0.5) +
  labs(x = 'Date', y = 'Empirical Coverage') +
  theme(panel.spacing = unit(2, 'lines'))
Cov_plot

ggsave(paste0(plot_path, '/', 'Supp_Re_pop_size_cov.png'), height = 10, width = 18)
#ggsave(paste0(plot_path, '/', 'Supp_Re_pop_size_epiestim_cov.png'), height = 10, width = 18)

######################################################################################################################
## Fig. S6: Fixed delay Grid #####
config_df <- expand.grid(fixed = c(TRUE, FALSE))

## Commented out because it takes long to run! ##
# for (i in 1:nrow(config_df)){
#   param_list <- list(
#     shift_times = c(0, 20, 50, 60, 80, 90, 100, 160, 180, 190, 210, 250, 270, 290, 310, 340, 365, 374),
#     R_levels = c(2.5, 0.5, 0.7, 1.1, 1.3, 0.9, 1.1, 0.8, 1.0),
#     init_I = 100,
#     noise = list(fitted_noise_model = get_noise(country_val = "CHE") ),
#     data_type = "confirmed", n_boot = 100,
#     smooth = TRUE,
#     delay = 0,
#     fixed_shift = config_df[i, 'fixed'],
#     timevar_sim = FALSE,
#     timevar_est = FALSE)
# 
#   par_Re_est <- parallel_simulation(param_list, ensemble_size = 100, number_of_cores = 8)
# 
#   dir.create(file.path(simDir, 'Fixed'))
#   write_csv(par_Re_est$sim, paste0(simDir, '/Fixed/sim_',
#                                    config_df[i, 'fixed'], '.csv'))
#   write_csv(par_Re_est$Re, paste0(simDir, '/Fixed/est_',
#                                   config_df[i, 'fixed'],'.csv'))
# 
# }

## Reading in #####

exp_result <- data.frame()
for (i in 1:nrow(config_df)){
  sim <- read_csv(paste0(simDir, '/Fixed/sim_', 
                         config_df[i, 'fixed'],'.csv'))
  Re <- read_csv(paste0(simDir, '/Fixed/est_', 
                        config_df[i, 'fixed'],'.csv'))
  
  config_result <- Re %>%
    left_join(sim, by = c("replicate", "date")) %>%
    mutate(in_CI = (Re > median_R_lowHPD) & (Re < median_R_highHPD),
           above_1 = (median_R_lowHPD > 1),
           below_1 = (median_R_highHPD < 1)) %>%
    mutate(fixed = config_df[i, ]) 
  
  exp_result <- bind_rows(exp_result, config_result)
}

## Plotting #####

plot_Re_est <- exp_result %>%
  pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
               names_to = 'name', values_to = 'value') %>%
  mutate(grouping = paste0(date, name),
         plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

exp_metrics <- compute_exp_metrics(exp_result %>% mutate(exp_name = fixed))


Re_plot <- ggplot() +
  geom_boxplot(data = plot_Re_est, 
               aes(x = date, y = value, group = grouping, colour = plot_colour), 
               position = 'identity', outlier.shape = NA )+
  geom_line(data = plot_Re_est %>% filter(replicate == 1), 
            aes(x = date, y = Re), size = 1, colour = 'black') +
  facet_grid(cols = vars(fixed),
             labeller = labeller(fixed = setNames(c("Deconvolution", "Fixed Shift"), c(FALSE, TRUE)) ) )+
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = 'Date', y = 'Re') +
  scale_colour_manual(values = c(viridis(2)[c(2,1)]), 
                      labels = c('Point estimate', 'CI'),
                      breaks = c('Point estimate', 'CI'),
                      name = 'Variable') +
  theme(legend.position = 'bottom')

Cov_plot <- ggplot(exp_metrics$cov) +
  geom_point(aes(x = date, y = empirical_cov)) +
  facet_grid(cols = vars(exp_name),
             labeller = labeller(exp_name = setNames(c("Deconvolution", "Fixed Shift"), c(FALSE, TRUE)) ) )+
  geom_hline(yintercept = 0.95, colour = 'red', linetype = 'dashed', alpha = 0.5) +
  labs(x = 'Date', y = 'Empirical Coverage') +
  theme(strip.text.x = element_blank())

rmse_plot <- ggplot(exp_metrics$rmse ) +
  geom_point(aes(x = date, y = rmse)) +
  facet_grid(cols = vars(exp_name),
             labeller = labeller(exp_name = setNames(c("Deconvolution", "Fixed Shift"), c(FALSE, TRUE)) ) )+
  labs(x = 'Date', y = 'RMSE') +
  theme(strip.text.x = element_blank())

Re_plot  + Cov_plot  + rmse_plot +
  plot_layout(ncol=1) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))

ggsave(paste0(plot_path, '/', 'Supp_Re_fixed.png'), height = 15, width = 12)


######################################################################################################################
## Fig. S7: Delay misspecification Grid ####
config_df <- expand.grid(misdelay = c(-2, -1, 0, 1, 2),
                         delay_opt = c('confirmed', 'death'))
config_df[,'delay_opt'] <- as.character(config_df[,'delay_opt'])

## Commented out because it takes long to run! ##
# for (i in 1:nrow(config_df)){
#   param_list <- list(
#     shift_times = c(0, 20, 50, 60, 80, 90, 100, 160, 180, 190, 210, 250, 270, 290, 310, 340, 365, 374),
#     R_levels = c(2.5, 0.5, 0.7, 1.1, 1.3, 0.9, 1.1, 0.8, 1.0),
#     init_I = 100,
#     noise = list(fitted_noise_model = get_noise(country_val = "CHE") ),
#     data_type = config_df[i, "delay_opt"], n_boot = 100,
#     smooth = TRUE,
#     delay = config_df[i, "misdelay"],
#     fixed_shift = FALSE,
#     timevar_sim = FALSE,
#     timevar_est = FALSE)
#   
#   par_Re_est <- parallel_simulation(param_list, ensemble_size = 100, number_of_cores = 8)
#   
#   dir.create(file.path(simDir, 'Misdelay'))
#   write_csv(par_Re_est$sim, paste0(simDir, '/Misdelay/sim_',
#                                    config_df[i, 'delay_opt'],'_',
#                                    config_df[i, 'misdelay'], '.csv'))
#   write_csv(par_Re_est$Re, paste0(simDir, '/Misdelay/est_',
#                                   config_df[i, 'delay_opt'],'_',
#                                   config_df[i, 'misdelay'], '.csv'))
#   
# }

## Reading in #####

exp_result <- data.frame()
for (i in 1:nrow(config_df)){
  sim <- read_csv(paste0(simDir, '/Misdelay/sim_', 
                         config_df[i, 'delay_opt'],'_',
                         config_df[i, 'misdelay'],'.csv'))
  Re <- read_csv(paste0(simDir, '/Misdelay/est_', 
                        config_df[i, 'delay_opt'],'_',
                        config_df[i, 'misdelay'],'.csv'))
  
  config_result <- Re %>%
    left_join(sim, by = c("replicate", "date")) %>%
    mutate(in_CI = (Re > median_R_lowHPD) & (Re < median_R_highHPD),
           above_1 = (median_R_lowHPD > 1),
           below_1 = (median_R_highHPD < 1)) %>%
    bind_cols(config_df[i, ]) %>%
    mutate(exp_name = paste0(delay_opt, '_', misdelay))
  
  exp_result <- bind_rows(exp_result, config_result)
}

## Plotting #####
exp_result <- exp_result %>%
  filter(delay_opt == "confirmed")

plot_Re_est <- exp_result %>%
  pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
               names_to = 'name', values_to = 'value') %>%
  mutate(grouping = paste0(date, name),
         plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

exp_metrics <- compute_exp_metrics(exp_result)

Re_plot <- ggplot() +
  geom_boxplot(data = plot_Re_est, 
               aes(x = date, y = value, group = grouping, colour = plot_colour), 
               position = 'identity', outlier.shape = NA )+
  geom_line(data = plot_Re_est %>% filter(replicate == 1), 
            aes(x = date, y = Re), size = 1, colour = 'black') +
  facet_grid(cols = vars(misdelay), rows = vars(delay_opt),
             labeller = labeller(delay_opt = setNames(c("Confirmed Cases", "Deaths"),
                                                      c("confirmed", "death")) ))+
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = 'Date', y = 'Re') +
  scale_colour_manual(values = c(viridis(2)[c(2,1)]), 
                      labels = c('Point estimate', 'CI'),
                      breaks = c('Point estimate', 'CI'),
                      name = 'Variable') +
  theme(legend.position = 'bottom',
        strip.text.y = element_blank())


Cov_plot <- ggplot(exp_metrics$cov %>% 
                     separate(exp_name, into = c("delay_opt", "misdelay"), sep = "_")) +
  geom_point(aes(x = date, y = empirical_cov)) +
  facet_grid(cols = vars(misdelay), rows = vars(delay_opt),
             labeller = labeller(delay_opt = setNames(c("Confirmed Cases", "Deaths"),
                                                      c("confirmed", "death")) ))+
  geom_hline(yintercept = 0.95, colour = 'red', linetype = 'dashed', alpha = 0.5) +
  labs(x = 'Date', y = 'Empirical Coverage') +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank())

rmse_plot <- ggplot(exp_metrics$rmse %>% 
                      separate(exp_name, into = c("delay_opt", "misdelay"), sep = "_")) +
  geom_point(aes(x = date, y = rmse)) +
  facet_grid(cols = vars(misdelay), rows = vars(delay_opt),
             labeller = labeller(delay_opt = setNames(c("Confirmed Cases", "Deaths"),
                                                      c("confirmed", "death")) ))+
  labs(x = 'Date', y = 'RMSE') +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank())

Re_plot  + Cov_plot  + rmse_plot +
  plot_layout(ncol=1) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25),
        plot.tag.position = 'topright',
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(plot_path, '/', 'Supp_Re_misdelay_confirmed.png'), height = 15, width = 17)

######################################################################################################################
## Fig. S8: Time-varying delay Grid #####
config_df <- expand.grid(timevar_est = c(TRUE, FALSE),
                             delay_opt = c('confirmed', 'death'))
config_df[,'delay_opt'] <- as.character(config_df[,'delay_opt'])

## Commented out because it takes long to run! ##
# for (i in 1:nrow(config_df)){
#   param_list <- list(
#     shift_times = c(0, 20, 50, 60, 80, 90, 100, 160, 180, 190, 210, 250, 270, 290, 310, 340, 365, 374),
#     R_levels = c(2.5, 0.5, 0.7, 1.1, 1.3, 0.9, 1.1, 0.8, 1.0),
#     init_I = 100,
#     noise = list(fitted_noise_model = get_noise(country_val = "CHE") ),
#     data_type = config_df[i, "delay_opt"], 
#     n_boot = 100,
#     smooth = TRUE,
#     delay = 0,
#     fixed_shift = FALSE,
#     timevar_sim = TRUE,
#     timevar_est = config_df[i, "timevar_est"])
#   
#   par_Re_est <- parallel_simulation(param_list, ensemble_size = 100, number_of_cores = 8)
#   
#   dir.create(file.path(simDir, 'Timevar'))
#   write_csv(par_Re_est$sim, paste0(simDir, '/Timevar/sim_', 
#                                    config_df[i, 'delay_opt'],'_',
#                                    config_df[i, 'timevar_est'], '.csv'))
#   write_csv(par_Re_est$Re, paste0(simDir, '/Timevar/est_', 
#                                   config_df[i, 'delay_opt'],'_',
#                                   config_df[i, 'timevar_est'], '.csv'))
#   
# }


## Reading in #####

exp_result <- data.frame()
for (i in 1:nrow(config_df)){
  sim <- read_csv(paste0(simDir, '/Timevar/sim_', 
                         config_df[i, 'delay_opt'],'_',
                         config_df[i, 'timevar_est'],'.csv'))
  Re <- read_csv(paste0(simDir, '/Timevar/est_', 
                        config_df[i, 'delay_opt'],'_',
                        config_df[i, 'timevar_est'],'.csv'))
  
  config_result <- Re %>%
    left_join(sim, by = c("replicate", "date")) %>%
    mutate(in_CI = (Re > median_R_lowHPD) & (Re < median_R_highHPD),
           above_1 = (median_R_lowHPD > 1),
           below_1 = (median_R_highHPD < 1)) %>%
    bind_cols(config_df[i, ]) %>%
    mutate(exp_name = paste0(delay_opt, '_', timevar_est))
  
  exp_result <- bind_rows(exp_result, config_result)
}

## Plotting #####
exp_result <- exp_result %>%
  #filter(delay_opt == 'confirmed')
  filter(delay_opt == 'death')

plot_Re_est <- exp_result %>%
  pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
               names_to = 'name', values_to = 'value') %>%
  mutate(grouping = paste0(date, name),
         plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

exp_metrics <- compute_exp_metrics(exp_result)

Re_plot <- ggplot() +
  geom_boxplot(data = plot_Re_est, 
               aes(x = date, y = value, group = grouping, colour = plot_colour), 
               position = 'identity', outlier.shape = NA )+
  geom_line(data = plot_Re_est %>% filter(replicate == 1), 
            aes(x = date, y = Re), size = 1, colour = 'black') +
  facet_grid(cols = vars(timevar_est), rows = vars(delay_opt),
             labeller = labeller(delay_opt = setNames(c("Confirmed Cases", "Deaths"),
                                                      c("confirmed", "death")),
                                 timevar_est = setNames(c("Estimated with time variation",
                                                          "Estimated without time variation"),
                                                        c(TRUE, FALSE)) ))+
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = 'Date', y = 'Re') +
  scale_colour_manual(values = c(viridis(2)[c(2,1)]), 
                      labels = c('Point estimate', 'CI'),
                      breaks = c('Point estimate', 'CI'),
                      name = 'Variable') +
  theme(legend.position = 'bottom',
        strip.text.y = element_blank())


Cov_plot <- ggplot(exp_metrics$cov %>% 
                     separate(exp_name, into = c("delay_opt", "timevar_est"), sep = "_")) +
  geom_point(aes(x = date, y = empirical_cov)) +
  facet_grid(cols = vars(timevar_est), rows = vars(delay_opt),
             labeller = labeller(delay_opt = setNames(c("Confirmed Cases", "Deaths"),
                                                      c("confirmed", "death")),
                                 timevar_est = setNames(c("Estimated with time variation",
                                                          "Estimated without time variation"),
                                                        c(TRUE, FALSE)) ))+
  geom_hline(yintercept = 0.95, colour = 'red', linetype = 'dashed', alpha = 0.5) +
  labs(x = 'Date', y = 'Empirical Coverage') +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank())

rmse_plot <- ggplot(exp_metrics$rmse %>% 
                      separate(exp_name, into = c("delay_opt", "timevar_est"), sep = "_")) +
  geom_point(aes(x = date, y = rmse)) +
  facet_grid(cols = vars(timevar_est), rows = vars(delay_opt),
             labeller = labeller(delay_opt = setNames(c("Confirmed Cases", "Deaths"),
                                                      c("confirmed", "death")),
                                 timevar_est = setNames(c("Estimated with time variation",
                                                          "Estimated without time variation"),
                                                        c(TRUE, FALSE)) ))+
  labs(x = 'Date', y = 'RMSE') +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank())

Re_plot  + Cov_plot  + rmse_plot +
  plot_layout(ncol=1) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25),
        plot.tag.position = 'topright')

#ggsave(paste0(plot_path, '/', 'Supp_Re_timevar.png'), height = 14, width = 12)
ggsave(paste0(plot_path, '/', 'Supp_Re_timevar_death.png'), height = 14, width = 12)



