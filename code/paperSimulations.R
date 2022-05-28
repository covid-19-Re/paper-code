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
  
  if(!is.null(param_list$epiestim)){
    if(param_list$epiestim){
    param_list$n_boot <- 0
    }
  }

  estimatedInfections <- estimateInfectionTS(simulation_i,
                                             getCountParams("incubation"),
                                             getCountParams(as.character(param_list$data_type), 
                                                            misdelay = param_list$delay),
                                             smooth_param = param_list$smooth,
                                             fixed_shift = param_list$fixed_shift,
                                             n_boot = param_list$n_boot,
                                             timevarying = param_list$timevar_est,
                                             days_incl = ifelse(is.null(param_list$days_incl), 
                                                                21, param_list$days_incl),
                                             block_size = ifelse(is.null(param_list$block_size), 
                                                                10, param_list$block_size))

  estimatedInfections <- estimatedInfections %>%
    mutate(local_infection = T)
  
  estimatedRe <- estimateReTS(estimatedInfections)

  if(param_list$n_boot == 0){
    summarisedRe <- as_tibble(estimatedRe) %>%
      pivot_wider(names_from = "variable", values_from = "value") %>%
      rename(median_R_mean = R_mean,
             median_R_highHPD = R_highHPD,
             median_R_lowHPD = R_lowHPD)
  } else if(!is.null(param_list$all_cis)){
    summarisedRe <- extractAllCIs(estimatedRe) %>%
      dplyr::select(-c(country, region, source, data_type, estimate_type)) %>%
      mutate(replicate = i)
  }
  else {
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

### All CIs
extractAllCIs <- function(countryEstimatesRaw,
                                   alpha=0.95){
  
  cleanEstimate <- as_tibble(countryEstimatesRaw)
  high_quan <- 1 - (1 - alpha) / 2
  low_quan <- (1 - alpha) / 2

  # EpiEstim-based
  # legacy_ReEstimates <- cleanEstimate %>%
  #   pivot_wider(names_from = "variable", values_from = "value") %>%
  #   dplyr::group_by(date, country, region, data_type, source, estimate_type) %>%
  #   dplyr::summarize(
  #     median_R_mean = median(R_mean),
  #     median_R_highHPD = median(R_highHPD),
  #     median_R_lowHPD = median(R_lowHPD),
  #     .groups = "keep"
  #   ) %>%
  #   dplyr::select(country, region, source, data_type, estimate_type, date,
  #                 median_R_mean, median_R_highHPD, median_R_lowHPD) %>%
  #   arrange(country, region, source, data_type, estimate_type, date) %>%
  #   ungroup()
  
  ##################################
    orig_ReEstimate <- cleanEstimate %>%
      filter(replicate == 0 ) %>%
      pivot_wider(names_from = "variable", values_from = "value") %>%
      rename(median_R_mean = R_mean,
             median_R_highHPD = R_highHPD,
             median_R_lowHPD = R_lowHPD) %>%
    dplyr::select(country, region, source, data_type, estimate_type, date,
                  median_R_mean, median_R_highHPD, median_R_lowHPD)
    # this is called median to be compatible with legacy code
    
    MM_ReEstimates <- cleanEstimate %>%
      filter(replicate != 0 ) %>%
      pivot_wider(names_from = "variable", values_from = "value") %>%
      dplyr::group_by(date, country, region, data_type, source, estimate_type) %>%
      dplyr::summarize(
        sd_mean = sd(R_mean), #across all bootstrap replicates
        .groups = "drop"
      ) %>%
      right_join(orig_ReEstimate, by = c('date', 'country', 'region',
                                         'data_type', 'source', 'estimate_type')) %>%
      dplyr::mutate(median_R_highHPD = median_R_mean + qnorm(high_quan)*sd_mean,
                    median_R_lowHPD = median_R_mean - qnorm(high_quan)*sd_mean
      ) %>%
      mutate(median_R_highHPD = ifelse(median_R_highHPD <0, 0, median_R_highHPD),
             median_R_lowHPD = ifelse(median_R_lowHPD <0, 0, median_R_lowHPD)
      )
    # we add the estimate-type extension later, because simpleUnion and
    # wideHPDs still derive from this df
    
    simple_Union <- MM_ReEstimates %>%
      left_join(orig_ReEstimate, by = c('date', 'country', 'region',
                                        'data_type', 'source', 'estimate_type')) %>%
      rowwise() %>%
      mutate(median_R_mean = median_R_mean.x,
             median_R_highHPD = max(median_R_highHPD.x, median_R_highHPD.y),
             median_R_lowHPD = min(median_R_lowHPD.x, median_R_lowHPD.y)) %>%
      dplyr::select(country, region, source, data_type, estimate_type, date,
                    median_R_mean, median_R_highHPD, median_R_lowHPD)
    
    ##################################
    # quantile and reverse quantile
    
    quantile <- cleanEstimate %>%
      filter(replicate != 0 ) %>%
      pivot_wider(names_from = "variable", values_from = "value") %>%
      dplyr::group_by(date, country, region, data_type, source, estimate_type) %>%
      dplyr::summarize(
        median_R_highHPD = quantile(R_mean, probs=high_quan, na.rm=T), #across all bootstrap replicates
        median_R_lowHPD = quantile(R_mean, probs=low_quan, na.rm=T),
        .groups = "drop"
      ) %>%
      right_join(orig_ReEstimate, by = c('date', 'country', 'region',
                                         'data_type', 'source', 'estimate_type')) %>%
      rename(median_R_highHPD = median_R_highHPD.x, median_R_lowHPD = median_R_lowHPD.x) %>%
      dplyr::select(country, region, source, data_type, estimate_type, date,
                    median_R_mean, median_R_highHPD, median_R_lowHPD)
    
    rev_quantile <- cleanEstimate %>%
      filter(replicate != 0 ) %>%
      pivot_wider(names_from = "variable", values_from = "value") %>%
      right_join(orig_ReEstimate, by = c('date', 'country', 'region',
                                         'data_type', 'source', 'estimate_type')) %>%
      # everything prefaced by median stems from orig data
      mutate(R_rev = R_mean - median_R_mean) %>%
      dplyr::group_by(date, country, region, data_type, source, estimate_type) %>%
      dplyr::summarize(
        median_R_highHPD = quantile(R_rev, probs=low_quan, na.rm=T), #across all bootstrap replicates
        median_R_lowHPD = quantile(R_rev, probs=high_quan, na.rm=T),
        #median_R_mean = median_R_mean,
        .groups = "drop"
      ) %>% 
      right_join(orig_ReEstimate, by = c('date', 'country', 'region',
                                         'data_type', 'source', 'estimate_type')) %>%
      mutate(median_R_highHPD = median_R_mean - median_R_highHPD.x,
             median_R_lowHPD = median_R_mean - median_R_lowHPD.x) %>%
      dplyr::select(country, region, source, data_type, estimate_type, date,
                    median_R_mean, median_R_highHPD, median_R_lowHPD)
    
    ##################################
    # combine all Re estimates
    unsortedReEstimates <- bind_rows(#EpiEstim_boot = legacy_ReEstimates, 
                                     EpiEstim_single = orig_ReEstimate,
                                     estimateR = simple_Union, 
                                     quant = quantile,
                                     rev_quant = rev_quantile, .id = 'CI_method') 
    
    ReEstimates <- unsortedReEstimates  %>%
        arrange(country, region, source, data_type, estimate_type, date) %>%
        ungroup()
    
  return(ReEstimates)
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
#     all_cis = TRUE, ## New 3.5.2022! 
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

# plot_Re_est <- exp_result %>%
#   pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
#                names_to = 'name', values_to = 'value') %>%
#   mutate(grouping = paste0(date, name),
#          plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

plot_Re_est <- exp_result %>%
  group_by(CI_method, date, noise_country, exp_name) %>%
  summarise(R_mean = mean(median_R_mean),
            sd_R_mean = sd(median_R_mean),
            R_highHPD = mean(median_R_highHPD),
            sd_highHPD = sd(median_R_highHPD),
            R_lowHPD = mean(median_R_lowHPD),
            sd_lowHPD = sd(median_R_lowHPD),
            Re = Re,
            .groups = 'drop') 

ggplot(plot_Re_est) +
  geom_ribbon(aes(x = date, ymin = R_highHPD - sd_highHPD, ymax = R_highHPD + sd_highHPD),
                  fill = viridis(2)[1])+
  geom_ribbon(aes(x = date, ymin = R_lowHPD - sd_lowHPD, ymax = R_lowHPD + sd_lowHPD),
              fill = viridis(2)[1])+
  geom_ribbon(aes(x = date, ymin = R_mean - sd_R_mean, ymax = R_mean + sd_R_mean),
              fill = viridis(2)[2])+
  geom_line(aes(x = date, y = Re), size = 1, colour = 'black') +
  #facet_wrap(vars(noise_country), ncol = 1)+
  facet_grid(cols = vars(CI_method), rows = vars(noise_country),
             labeller = labeller(CI_method = setNames(c("EpiEstim", "estimateR", "Quantiles",
                                                        "Reverse Quantiles"), 
                                  c("EpiEstim_single", "estimateR", "quant", "rev_quant")) ))+
  coord_cartesian(ylim = c(0,5)) +
  labs(x = 'Date', y = 'Re', colour = 'Variable') +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.spacing = unit(2, 'lines'))

ggsave(paste0(plot_path, '/', 'Supp_Re_noise.pdf'), height = 12, width = 14)

## Additional metrics; specifically coverage
up_exp_result <- exp_result %>%
  mutate(exp_name = paste0(exp_name, '_', CI_method))
exp_metrics <- compute_exp_metrics(up_exp_result)

cov_df <- exp_metrics$cov %>%
  separate(exp_name, into = c('init', 'noise_country', 'CI_method'), sep = c(4, 8)) %>%
  mutate(noise_country = sub('_', '', noise_country))

ggplot(cov_df) +
  geom_point(aes(x = date, y = empirical_cov)) +
  facet_grid(cols = vars(CI_method), rows = vars(noise_country),
             labeller = labeller(CI_method = setNames(c("EpiEstim", "estimateR", "Quantiles",
                                                        "Reverse Quantiles"), 
                                                      c("EpiEstim_single", "estimateR", "quant", "rev_quant")) ))+
  geom_hline(yintercept = 0.95, colour = 'red', linetype = 'dashed', alpha = 0.5) +
  labs(x = 'Date', y = 'Empirical Coverage') +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.spacing = unit(2, 'lines'))

ggsave(paste0(plot_path, '/', 'Supp_Re_noise_cov.pdf'), height = 12, width = 14)


## Plotting Fig 1 ####

plot_Re_est <- exp_result %>%
  filter(exp_name == '100_CHE') %>%
  pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
               names_to = 'name', values_to = 'value') %>%
  mutate(grouping = paste0(date, name),
         plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

# plot_Re_est <- exp_result %>%
#   group_by(CI_method, date, noise_country, exp_name) %>%
#   summarise(R_mean = mean(median_R_mean),
#             R_highHPD = mean(median_R_highHPD),
#             sd_highHPD = sd(median_R_highHPD),
#             R_lowHPD = mean(median_R_lowHPD),
#             sd_lowHPD = sd(median_R_lowHPD),
#             .groups = 'drop') 

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

ggsave(paste0(plot_path, '/', 'Figure1.pdf'), height = 15, width = 12)



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

# plot_Re_est <- exp_result %>%
#   group_by(CI_method, date, noise_country, exp_name) %>%
#   summarise(R_mean = mean(median_R_mean),
#             R_highHPD = mean(median_R_highHPD),
#             sd_highHPD = sd(median_R_highHPD),
#             R_lowHPD = mean(median_R_lowHPD),
#             sd_lowHPD = sd(median_R_lowHPD),
#             .groups = 'drop') 

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

ggsave(plot = combi_plot, paste0(plot_path, '/', 'Supp_Re_scenarios.pdf'), height = 20, width = 17)


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

# plot_Re_est <- exp_result %>%
#   group_by(CI_method, date, noise_country, exp_name) %>%
#   summarise(R_mean = mean(median_R_mean),
#             R_highHPD = mean(median_R_highHPD),
#             sd_highHPD = sd(median_R_highHPD),
#             R_lowHPD = mean(median_R_lowHPD),
#             sd_lowHPD = sd(median_R_lowHPD),
#             .groups = 'drop') 

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

ggsave(paste0(plot_path, '/', 'Supp_Re_noise_non_smooth.pdf'), height = 14, width = 12)

######################################################################################################################
##### Fig. S5: Dependence on population size ######
config_df <- expand.grid(R1 = c(0.8, 1, 1.5),
                         init = c(10, 100, 1000, 5000, 10000),
                         epiestim = FALSE)
                         #epiestim = TRUE)

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

# plot_Re_est <- exp_result %>%
#   group_by(CI_method, date, noise_country, exp_name) %>%
#   summarise(R_mean = mean(median_R_mean),
#             R_highHPD = mean(median_R_highHPD),
#             sd_highHPD = sd(median_R_highHPD),
#             R_lowHPD = mean(median_R_lowHPD),
#             sd_lowHPD = sd(median_R_lowHPD),
#             .groups = 'drop') 

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

ggsave(paste0(plot_path, '/', 'Supp_Re_pop_size.pdf'), height = 10, width = 18)
#ggsave(paste0(plot_path, '/', 'Supp_Re_pop_size_epiestim.pdf'), height = 10, width = 18)

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

ggsave(paste0(plot_path, '/', 'Supp_Re_pop_size_cov.pdf'), height = 10, width = 18)
#ggsave(paste0(plot_path, '/', 'Supp_Re_pop_size_epiestim_cov.pdf'), height = 10, width = 18)


Re_plot  + 
  Cov_plot + 
  plot_layout(ncol=1) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25),
        plot.tag.position = 'topleft')

ggsave(paste0(plot_path, '/', 'Supp_Re_pop_size_combi.pdf'), height = 22, width = 18)
#ggsave(paste0(plot_path, '/', 'Supp_Re_pop_size_epiestim_combi.pdf'), height = 22, width = 18)


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

# plot_Re_est <- exp_result %>%
#   group_by(CI_method, date, noise_country, exp_name) %>%
#   summarise(R_mean = mean(median_R_mean),
#             R_highHPD = mean(median_R_highHPD),
#             sd_highHPD = sd(median_R_highHPD),
#             R_lowHPD = mean(median_R_lowHPD),
#             sd_lowHPD = sd(median_R_lowHPD),
#             .groups = 'drop') 

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

ggsave(paste0(plot_path, '/', 'Supp_Re_fixed.pdf'), height = 15, width = 12)

# stats
exp_metrics$rmse %>%
  group_by(exp_name) %>%
  summarise(mean = mean(rmse),
            cum = sum(rmse))
#c("Deconvolution", "Fixed Shift")
#c(FALSE, TRUE))
exp_metrics$cov %>%
  group_by(exp_name) %>%
  summarise(mean = mean(empirical_cov))

######################################################################################################################
## Fig. S7: Delay misspecification Grid ####
config_df <- expand.grid(misdelay = c(-5,-2, -1, 0,1,2, 5, 10),
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
  #filter(delay_opt == "death") # does not provide a lot of additional information

# plot_Re_est <- exp_result %>%
#   pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
#                names_to = 'name', values_to = 'value') %>%
#   mutate(grouping = paste0(date, name),
#          plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

plot_Re_est <- exp_result %>%
  group_by(date, misdelay, exp_name) %>%
  summarise(R_mean = mean(median_R_mean),
            R_highHPD = mean(median_R_highHPD),
            sd_highHPD = sd(median_R_highHPD),
            R_lowHPD = mean(median_R_lowHPD),
            sd_lowHPD = sd(median_R_lowHPD),
            .groups = 'drop')

exp_metrics <- compute_exp_metrics(exp_result)

Re_plot <- ggplot(plot_Re_est) +
  geom_ribbon(aes(x = date, ymin = R_highHPD - sd_highHPD, 
                  ymax = R_highHPD + sd_highHPD), fill = viridis(2)[1])+
  geom_ribbon(aes(x = date, ymin = R_lowHPD - sd_lowHPD, 
                  ymax = R_lowHPD + sd_lowHPD), fill = viridis(2)[1])+
  geom_line(aes(x = date, y = R_mean), size = 1, colour = viridis(2)[2]) +
  facet_wrap(vars(misdelay), ncol= 4) +
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = 'Date', y = 'Re', colour = 'Variable', fill = 'Variable') +
  theme(legend.position = 'bottom',
        strip.text.y = element_blank())
Re_plot
ggsave(paste0(plot_path, '/', 'Supp_Re_misdelay_A.pdf'), height = 9, width = 16)

Cov_plot <- ggplot(exp_metrics$cov %>% 
                     separate(exp_name, into = c("delay_opt", "misdelay"), sep = "_")%>%
                     mutate(misdelay = factor(misdelay, levels = c('-5','-2', '-1', '0', '1', '2', '5', '10'))) ) +
  geom_point(aes(x = date, y = empirical_cov)) +
  # facet_grid(cols = vars(misdelay), rows = vars(delay_opt),
  #            labeller = labeller(delay_opt = setNames(c("Confirmed Cases", "Deaths"),
  #                                                     c("confirmed", "death")) ))+
  facet_wrap(vars(misdelay), ncol= 4) +
  geom_hline(yintercept = 0.95, colour = 'red', linetype = 'dashed', alpha = 0.5) +
  labs(x = 'Date', y = 'Empirical Coverage') +
  theme(#strip.text.x = element_blank(),
        strip.text.y = element_blank())
Cov_plot
ggsave(paste0(plot_path, '/', 'Supp_Re_misdelay_B.pdf'), height = 9, width = 16)

rmse_plot <- ggplot(exp_metrics$rmse %>% 
                      separate(exp_name, into = c("delay_opt", "misdelay"), sep = "_") %>%
                      mutate(misdelay = factor(misdelay, levels = c('-5','-2', '-1', '0', '1', '2', '5', '10'))) ) +
  geom_point(aes(x = date, y = rmse)) +
  # facet_grid(cols = vars(misdelay), rows = vars(delay_opt),
  #            labeller = labeller(delay_opt = setNames(c("Confirmed Cases", "Deaths"),
  #                                                     c("confirmed", "death")) ))+
  facet_wrap(vars(misdelay), ncol= 4) +
  labs(x = 'Date', y = 'RMSE') +
  theme(#strip.text.x = element_blank(),
        strip.text.y = element_blank())
rmse_plot
ggsave(paste0(plot_path, '/', 'Supp_Re_misdelay_C.pdf'), height = 9, width = 16)

Re_plot  + Cov_plot  + rmse_plot +
  plot_layout(ncol=1) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25),
        plot.tag.position = 'topright',
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(plot_path, '/', 'Supp_Re_misdelay_confirmed.pdf'), height = 20, width = 16)
#ggsave(paste0(plot_path, '/', 'Supp_Re_misdelay_death.png'), height = 15, width = 17)

##### Change layout #######
# Re_plot_conf <- Re_plot
# Cov_plot_conf <- Cov_plot
# rmse_plot_conf <- rmse_plot
# 
# rmse_plot_conf  + rmse_plot +
#   plot_layout(ncol=1) + 
#   plot_annotation(tag_levels = 'A') & 
#   theme(plot.tag = element_text(size = 25),
#         plot.tag.position = 'topright',
#         axis.text.x = element_text(angle = 45, hjust = 1))

#ggsave(paste0(plot_path, '/', 'Supp_Re_misdelay_A.png'), height = 15, width = 17)

# stats
exp_metrics$rmse %>%
  group_by(exp_name) %>%
  summarise(mean = mean(rmse),
            cum = sum(rmse))

exp_metrics$cov %>%
  group_by(exp_name) %>%
  summarise(mean = mean(empirical_cov))

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

# plot_Re_est <- exp_result %>%
#   group_by(CI_method, date, noise_country, exp_name) %>%
#   summarise(R_mean = mean(median_R_mean),
#             R_highHPD = mean(median_R_highHPD),
#             sd_highHPD = sd(median_R_highHPD),
#             R_lowHPD = mean(median_R_lowHPD),
#             sd_lowHPD = sd(median_R_lowHPD),
#             .groups = 'drop') 

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

#ggsave(paste0(plot_path, '/', 'Supp_Re_timevar.pdf'), height = 14, width = 12)
ggsave(paste0(plot_path, '/', 'Supp_Re_timevar_death.pdf'), height = 14, width = 12)

# stats
exp_metrics$rmse %>%
  group_by(exp_name) %>%
  summarise(mean = mean(rmse),
            cum = sum(rmse))

exp_metrics$cov %>%
  group_by(exp_name) %>%
  summarise(mean = mean(empirical_cov))

######################################################################################################################
# 7, 14, or 21 days of smoothing

config_df <- expand.grid(smooth_day = c(7, 14, 21))
  
## Commented out because it takes long to run! ##
# for (i in 1:nrow(config_df)){
#   param_list <- list(
#     shift_times = c(0, 20, 50, 60, 80, 90, 100, 160, 180, 190, 210, 250, 270, 290, 310, 340, 365, 374),
#     R_levels = c(2.5, 0.5, 0.7, 1.1, 1.3, 0.9, 1.1, 0.8, 1.0),
#     init_I = 100,
#     noise = list(fitted_noise_model = get_noise(country_val = "CHE") ),
#     data_type = "confirmed",
#     n_boot = 100,
#     smooth = TRUE,
#     days_incl = config_df[i, 'smooth_day'],
#     delay = 0,
#     fixed_shift = FALSE,
#     timevar_sim = FALSE,
#     timevar_est = FALSE)
# 
#   par_Re_est <- parallel_simulation(param_list, ensemble_size = 1, number_of_cores = 8)
# 
#   dir.create(file.path(simDir, 'smooth_day'))
#   config_result <- par_Re_est$Re %>%
#     left_join(par_Re_est$sim, by = c("replicate", "date")) %>%
#     mutate(in_CI = (Re > median_R_lowHPD) & (Re < median_R_highHPD),
#            above_1 = (median_R_mean > 1) ) %>%
#     bind_cols(smooth_day = config_df[i, ]) %>%
#     mutate(exp_name = smooth_day)
# 
#   write_csv(config_result, paste0(simDir, '/smooth_day/all_',
#                                   config_df[i, 'smooth_day'],'.csv'))
# 
# }

## Reading in #####

exp_result <- data.frame()
for (i in 1:nrow(config_df)){
  config_result <- read_csv(paste0(simDir, '/smooth_day/all_', 
                         config_df[i, 'smooth_day'],'.csv')) %>%
    mutate(below_1 = (median_R_highHPD < 1))
  
  exp_result <- bind_rows(exp_result, config_result)
}

## Plotting #####
plot_Re_est <- exp_result %>%
  pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
               names_to = 'name', values_to = 'value') %>%
  mutate(grouping = paste0(date, name),
         plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

# plot_Re_est <- exp_result %>%
#   group_by(CI_method, date, noise_country, exp_name) %>%
#   summarise(R_mean = mean(median_R_mean),
#             R_highHPD = mean(median_R_highHPD),
#             sd_highHPD = sd(median_R_highHPD),
#             R_lowHPD = mean(median_R_lowHPD),
#             sd_lowHPD = sd(median_R_lowHPD),
#             .groups = 'drop') 

exp_metrics <- compute_exp_metrics(exp_result)

Re_plot <- ggplot(exp_result) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD, ymax = median_R_highHPD), fill = viridis(2)[1])+
  geom_line(data = exp_result %>% filter(replicate == 1), 
            aes(x = date, y = median_R_mean), size = 1, colour = viridis(2)[2]) +
  facet_wrap(vars(smooth_day), ncol =1)+
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = 'Date', y = 'Re') +
  theme(legend.position = 'bottom',
        strip.text.y = element_blank())

Re_plot
#ggsave(paste0(plot_path, '/', 'Supp_smooth_days.png'), height = 6, width = 16)
ggsave(paste0(plot_path, '/', 'Supp_smooth_days.pdf'), height = 12, width = 10)


######################################################################################################################
# Different block sizes

config_df <- expand.grid(block_size = seq(7,14))

## Commented out because it takes long to run! ##
# for (i in 1:nrow(config_df)){
#   param_list <- list(
#     shift_times = c(0, 20, 50, 60, 80, 90, 100, 160, 180, 190, 210, 250, 270, 290, 310, 340, 365, 374),
#     R_levels = c(2.5, 0.5, 0.7, 1.1, 1.3, 0.9, 1.1, 0.8, 1.0),
#     init_I = 100,
#     noise = list(fitted_noise_model = get_noise(country_val = "CHE") ),
#     data_type = "confirmed",
#     n_boot = 100,
#     smooth = TRUE,
#     block_size = config_df[i, 'block_size'],
#     delay = 0,
#     fixed_shift = FALSE,
#     timevar_sim = FALSE,
#     timevar_est = FALSE)
#   
#   par_Re_est <- parallel_simulation(param_list, ensemble_size = 1, number_of_cores = 8)
#   
#   dir.create(file.path(simDir, 'block_size'))
#   config_result <- par_Re_est$Re %>%
#     left_join(par_Re_est$sim, by = c("replicate", "date")) %>%
#     mutate(in_CI = (Re > median_R_lowHPD) & (Re < median_R_highHPD),
#            above_1 = (median_R_mean > 1) ) %>%
#     bind_cols(block_size = config_df[i, ]) %>%
#     mutate(exp_name = block_size)
#   
#   write_csv(config_result, paste0(simDir, '/block_size/all_',
#                                   config_df[i, 'block_size'],'.csv'))
#   
# }

## Reading in #####

exp_result <- data.frame()
for (i in 1:nrow(config_df)){
  config_result <- read_csv(paste0(simDir, '/block_size/all_', 
                                   config_df[i, 'block_size'],'.csv')) %>%
    mutate(below_1 = (median_R_highHPD < 1))
  
  exp_result <- bind_rows(exp_result, config_result)
}

## Plotting #####
plot_Re_est <- exp_result %>%
  pivot_longer(cols = c('median_R_mean', 'median_R_highHPD', 'median_R_lowHPD'),
               names_to = 'name', values_to = 'value') %>%
  mutate(grouping = paste0(date, name),
         plot_colour = ifelse(name == 'median_R_mean', 'Point estimate', 'CI'))

# plot_Re_est <- exp_result %>%
#   group_by(CI_method, date, noise_country, exp_name) %>%
#   summarise(R_mean = mean(median_R_mean),
#             R_highHPD = mean(median_R_highHPD),
#             sd_highHPD = sd(median_R_highHPD),
#             R_lowHPD = mean(median_R_lowHPD),
#             sd_lowHPD = sd(median_R_lowHPD),
#             .groups = 'drop') 

exp_metrics <- compute_exp_metrics(exp_result)

# Re_plot <- ggplot(exp_result) +
#   geom_ribbon(aes(x = date, ymin = median_R_lowHPD, ymax = median_R_highHPD), fill = viridis(2)[1])+
#   geom_line(data = exp_result %>% filter(replicate == 1), 
#             aes(x = date, y = median_R_mean), size = 1, colour = viridis(2)[2]) +
#   facet_wrap(vars(block_size))+
#   coord_cartesian(ylim = c(0, 5)) +
#   labs(x = 'Date', y = 'Re') +
#   theme(legend.position = 'bottom',
#         strip.text.y = element_blank())

Re_plot <- ggplot(exp_result %>% mutate(block_size = as.factor(block_size))) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD, ymax = median_R_highHPD, 
                  group = block_size, colour = block_size, fill = block_size), alpha = 0.3) +
  scale_colour_manual(values = viridis(8)) +
  scale_fill_manual(values = viridis(8)) +
  coord_cartesian(ylim = c(0, 3.5)) +
  labs(x = 'Date', y = 'Re', colour = 'Block Size', fill = 'Block Size') +
  theme(legend.position = 'bottom',
        strip.text.y = element_blank())

Re_plot
ggsave(paste0(plot_path, '/', 'Supp_block_size.pdf'), height = 10, width = 14)

