###########################################################
## paperSimulationGrid.R
## author: J.S. Huisman
###########################################################

# Note: it takes quite long (~2 hours) to calculate each of 
# the simulation grids (and produces a lot of csv files),
# so these commands have been commented out. 
# Instead I supplied the files with final
# values ("plot_results"), shown in the paper. 
# This code could probably be parallelised easily, but that
# did not have priority at the time.

library(tidyverse)
library(cowplot)
library(viridis)

simDir = '../simulations'
plot_path = '../figures'

source('generateSimulations.R')
source('compareSimulations.R')

###########################################################
getCountParams <- function(obs_type, misdelay = 0){
  switch(obs_type,
         incubation = getGammaParams(5.3, 3.2),
         zero = list(shape = 0, scale = 0),
         death = getGammaParams(15.0 + misdelay, 6.9),
         hospitalisation = getGammaParams(6.1 + misdelay, 4.7),
         confirmed = getGammaParams(4.5 + misdelay, 4.9) )
}

###########################################################
simulateGrid <- function(cond_grid, smooth_R = FALSE){
  for (row in 1:nrow(cond_grid)){
    
    condDir = cond_grid[row, 'simulationDir']
    
    if(!dir.exists(condDir)){
      dir.create(condDir, recursive = TRUE)
    }
    
    condFile = paste0(condDir, '/', cond_grid[row, 'filename'])
    
    if(!file.exists(condFile)){
      
      simulation <- simulateTS(as.numeric(cond_grid[row, paste0('t', 1:6)]), 
                               as.numeric(cond_grid[row, paste0('R', 1:3)]),
                               getCountParams("incubation"), 
                               getCountParams(as.character(cond_grid[row, 'delay_opt'])),
                               init_infection = cond_grid[row, 'init'],
                               noise = list('weekly' = cond_grid[row, 'weekly'],
                                            'gaussian' = cond_grid[row, 'gaussian']),
                               smooth_R = smooth_R,
                               timevarying = cond_grid[row, 'timevar_sim'])
      
      write_csv(simulation, condFile)
      
    }
  }
  
}

estimateGrid <- function(cond_grid){
  
  for (row in 1:nrow(cond_grid)){
    
    condDir = cond_grid[row, 'estimationDir']
    
    if(!dir.exists(condDir)){
      dir.create(condDir, recursive = TRUE)
    }
    
    condFile = cond_grid[row, 'filename']
    simFile = paste0(cond_grid[row, 'simulationDir'],
                     '/', condFile)

    
    if(!file.exists(paste0(condDir, cond_grid[row, 'infection_file']))){
    
      simulation <- read_csv(simFile)
      OnsetParams <- getCountParams(as.character(cond_grid[row, 'delay_opt']),
                                    misdelay = cond_grid[row, 'misdelay'])
      
      estimatedInfections <- try(estimateInfectionTS(simulation, getCountParams("incubation"), 
                                                     OnsetParams,
                                                     smooth_param = (cond_grid[row, 'smooth'] == 'smooth'),
                                                     fixed_shift = cond_grid[row, 'fixed'],
                                                     timevarying = cond_grid[row, 'timevar_est']))
      if('try-error' %in% class(estimatedInfections)){
        next
      }
      write_csv(estimatedInfections, paste0(condDir, cond_grid[row, 'infection_file']))
      
      estimatedRe <- estimateReTS(estimatedInfections, delay = 0)
      write_csv(estimatedRe, paste0(condDir, cond_grid[row, 're_file']))
    }
  }
  
}

plotSimulations <- function(valid_cond_grid){
  for (row_id in 1:nrow(valid_cond_grid)){
    
    simulation <- read_csv(paste0(cond_grid[row_id, 'simulationDir'], '/', 
                                  cond_grid[row_id, 'filename']))
    
    estimatedInfections <- read_csv(paste0(cond_grid[row_id, 'estimationDir'],
                                           cond_grid[row_id, 'infection_file']))
    
    estimatedRe <- read_csv(paste0(cond_grid[row_id, 'estimationDir'], 
                                   cond_grid[row_id, 're_file']))
    
    longSim <- simulation %>%
      pivot_longer(cols = c(Re, infections, observations),
                   names_to = 'type') %>%
      mutate(type = factor(type, levels = c('Re', 'infections', 'observations')),
             source_type = 'simulation')
    
    longInfections <- estimatedInfections %>%
      dplyr::select(c(date, replicate, value)) %>% 
      group_by(date) %>%
      summarize(
        median_val = median(value),
        high_quant = quantile(value, probs=0.975, na.rm=T),
        low_quant = quantile(value, probs=0.025, na.rm=T),
        .groups = "keep"
      ) %>%
      ungroup() %>%
      mutate(type = factor('infections', levels = c('Re', 'infections', 'observations')),
             source_type = 'estimation')
    
    longRe <- cleanReTSestimate(estimatedRe) %>%
      mutate(type = factor('Re', levels = c('Re', 'infections', 'observations')),
             source_type = 'estimation')
    
    p <- ggplot() +
      geom_line(data = longSim, aes(x = date, y = value, colour = type)) +
      geom_ribbon(data = longRe, aes(x = date, ymin = median_R_lowHPD,
                                     ymax = median_R_highHPD), alpha = 0.7) +
      geom_line(data = longRe, aes(x = date, y = median_R_mean), size=1.1) +
      geom_ribbon(data = longInfections, aes(x = date, ymin = low_quant,
                                             ymax = high_quant), alpha = 0.7) +
      geom_line(data = longInfections, aes(x = date, y = median_val), size=1.1) +
      facet_grid(rows = vars(type), scale = 'free') +
      theme_minimal()
    
    plotPath <- paste0(cond_grid[row_id, 'simulationDir'], '/', 
                       valid_cond_grid[row_id, 'smooth'], '_',
                       gsub('.csv' , '.png', valid_cond_grid[row_id,'filename']))
    ggsave(plotPath, plot = p, width = 8, height = 8)
    
  }
  
}

###########################################################
defaults = list(t1 = 0, t2 = 10, t3 = 20, t4 = 10,
                t5 = 30, t6 = 50,
                R1 = 3.5, R2 = 0.5, R3 = 1.2,
                smooth = 'smooth', delay_opt = 'confirmed',
                init = 10, 
                timevar_est = FALSE, timevar_sim = FALSE,
                fixed = FALSE, misdelay = 0,
                weekly = 1, gaussian = 0)

complete_cond_grid <- function(raw_cond_grid, defaults){
  
  for (var in names(defaults)){
    if(!var %in% colnames(raw_cond_grid) ){
      if ( str_detect(var, "^t[1-9]") & !str_detect(var, "t[1-3]$") ){
        raw_cond_grid[var] = raw_cond_grid[,'t3'] + defaults[[var]]
      } else {
        raw_cond_grid[var] = defaults[var]
      }
      
    }
  }
  
  cond_grid <- raw_cond_grid %>%
    mutate(filename = paste0('t_', t2, '_', t3,
                             '_R_', R1, '_init_', init,
                             '_week_', weekly,
                             '_norm_', gaussian,
                             '_timevar_', timevar_sim,
                             '_', timevar_est,
                             '_fixed_', fixed,
                             '_misdelay_', misdelay,
                             '.csv'),
           estimationDir = paste0(simDir, '/estimated/', smooth,
                                  '/', delay_opt),
           simulationDir = paste0(simDir, '/simulated/', delay_opt),
           infection_file = paste0('/infection_', filename),
           re_file = paste0('/Re_', filename)
    )
  return(cond_grid)
}

get_valid_cond_grid <- function(cond_grid){
  valid_cond_grid <- cond_grid %>%
    mutate(infection_exists = file.exists(paste0(estimationDir, infection_file)),
           re_exists = file.exists(paste0(estimationDir, re_file)),
           infection_size = file.size(paste0(estimationDir, infection_file)),
           re_size = file.size(paste0(estimationDir, re_file)))
  
  valid_cond_grid <- valid_cond_grid %>% 
    filter(re_exists & (re_size > 0))
  
  return(valid_cond_grid)
}

###########################################################
getMaxIncidence <- function(x, valid_cond_grid){
  simulation <- read_csv(paste0(simDir, '/simulated/', 
                                valid_cond_grid[x, 'delay_opt'], '/', 
                                valid_cond_grid[x,'filename']))
  
  max_inc <- max(simulation$infections)
  return(max_inc)
}

getPlotResults <- function(valid_cond_grid){
  ReError = data.frame()
  for (row_id in 1:nrow(valid_cond_grid)){
    
    simulation <- try(read_csv(paste0(valid_cond_grid[row_id, 'simulationDir'], '/', 
                                  valid_cond_grid[row_id, 'filename'])))
    if('try-error' %in% class(simulation)){next}
    
    estimatedInfections <- try(read_csv(paste0(valid_cond_grid[row_id, 'estimationDir'],
                                           valid_cond_grid[row_id, 'infection_file'])))
    if('try-error' %in% class(estimatedInfections)){next}
    
    estimatedRe <- try(read_csv(paste0(valid_cond_grid[row_id, 'estimationDir'], 
                                   valid_cond_grid[row_id, 're_file'])))
    if('try-error' %in% class(estimatedRe)){next}
    
    new_error <- getReError(simulation, estimatedRe)
    ReError <- bind_rows(ReError, new_error)
    
  }
  
  maxInc <- sapply(1:nrow(valid_cond_grid), FUN = getMaxIncidence, valid_cond_grid)
  SlopeError <- getMeanSlopeError(valid_cond_grid)
  
  val_results <- bind_cols(valid_cond_grid, ReError, SlopeError, maxInc = maxInc)
  plot_results <- val_results %>%
    mutate(t3diff = t3 - t2) %>%
    mutate(across(c("ReRMSE", "slopeError", "RootDiff", "EmpCoverage",
                    "OneTrans"), as.numeric))
  return(plot_results)
}

###########################################################
# For the plotting
strip_plot <- function(input_plot){ input_plot + theme(legend.position = "none",
                                                       axis.text.x = element_blank(), 
                                                       axis.ticks.x = element_blank(),
                                                       axis.title.x = element_blank(),
                                                       strip.text.x = element_blank(),
                                                       strip.text.y = element_blank()) }

get_plotlist <- function(long_plot_results, ylimits, col_var = NULL) {
  stat_plots <- list()
  for (stat_i in c('ReRMSE', 'slopeError', 'RootDiff', 
                   'OneTrans', 'EmpCoverage')){
    
    subset_plot_results <- long_plot_results %>%
      filter(measure == stat_i) %>%
      mutate(t3 = factor(t3))
    
    scaleFUN <- function(x) sprintf("%.2f", x)
    
    p <- ggplot(subset_plot_results) +
      geom_point(aes(x = log10(maxInc), y = measure_val, colour = t3, 
                     shape = as.factor(R1)), size = 2, show.legend = TRUE ) +
      geom_line(aes(x = log10(maxInc), y = measure_val, 
                    colour = t3, group = t3), linetype = 'dotted',  show.legend = FALSE ) +
      scale_y_continuous(labels=scaleFUN, limits = ylimits[[stat_i]]) +
      labs(x = 'Peak infection incidence (log)', y = 'Measure value',
           colour = 'Time of 1st \ndescent (T3)', shape = 'First Re \nlevel (R1)') +
      scale_colour_viridis(discrete = T) +
      theme_bw() +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        axis.text.y= element_text(size=14),
        axis.text.x= element_text(size=14),
        axis.title.x =  element_text(size=17),
        axis.title.y =  element_blank(),
        legend.title = element_text(size=17),
        legend.text = element_text(size=15),
        strip.placement = "outside"
      )
    
    if (is.null(col_var)){
      p <- p +
        facet_grid( rows = vars(measure), 
                    scale = 'free', switch = 'y',
                    labeller = as_labeller(method_labels))
    } else {
      p <- p +
        facet_grid(cols = vars(get(col_var)), rows = vars(measure), 
                   scale = 'free', switch = 'y')
    }
    
    stat_plots <- rlist::list.append(stat_plots, p)
  }
  return(stat_plots)
}

get_plotgrid <- function(stat_plots){
  
  plot_legend <- get_legend(stat_plots[[1]])
  all_plots <- plot_grid( NULL, stat_plots[[1]] + theme(legend.position = "none",
                                                        axis.text.x = element_blank(), 
                                                        axis.ticks.x = element_blank(),
                                                        axis.title.x = element_blank(),
                                                        strip.text.y = element_blank()) ,
                          NULL, strip_plot(stat_plots[[2]]) ,
                          NULL, strip_plot(stat_plots[[3]]) ,
                          NULL, strip_plot(stat_plots[[4]]) ,
                          NULL, stat_plots[[5]] + theme(legend.position = "none",
                                                        strip.text.x = element_blank(),
                                                        strip.text.y = element_blank()) ,
                          align = "v",
                          labels = c('', 'RMSE', 
                                     '', 'Relative error of the slope', 
                                     '', 'R=1 crossing (days)',
                                     '', 'Correct R>1 / R<1',
                                     '', 'Coverage' ), label_size = 15, 
                          vjust = c(0, +0.1, 
                                    0, +0.1, 
                                    0, +0.1, 
                                    0, +0.1, 
                                    0, +0.1),
                          hjust = c(0, -1.24, 
                                    0, -0.29, 
                                    0, -0.37, 
                                    0, -0.41, 
                                    0, -0.75), 
                          rel_heights = c(0.1, 1, 
                                          0.075, 1,
                                          0.075, 1,
                                          0.075, 1,
                                          0.075, 1),
                          nrow = 10, ncol = 1)
  
  final_plot <- plot_grid(all_plots, plot_legend, ncol = 2, rel_widths = c(1, .2), align = "h", axis = "t")
  
  return(final_plot)
}
###########################################################
### FIG 1 + S2
###########################################################
#### Set up simulation grid #####

raw_cond_grid <- expand.grid(R1 = c(2.5, 3.0, 3.5),
                         init = c(10, 100, 1000),
                          t3 = c(21, 26, 31, 36),
                          t2 = 20,
                          weekly = c(0.25, 0.75, 1.0), 
                          gaussian = c(0.0, 0.1, 0.2) )

#### Run simulation grid #####
cond_grid <- complete_cond_grid(raw_cond_grid, defaults)

#simulateGrid(cond_grid, smooth_R = TRUE)
#estimateGrid(cond_grid)

##### Plot simulation traces ######

#valid_cond_grid <- get_valid_cond_grid(cond_grid)
#write_csv(valid_cond_grid, paste0(simDir, '/valid_cond_grid.csv'))

#valid_cond_grid = read_csv(paste0(simDir, '/valid_cond_grid.csv'))

#plotSimulations(valid_cond_grid)

##### calculate performance metrics #####

#plot_results <- getPlotResults(valid_cond_grid)
#write_csv(plot_results, paste0(simDir, '/plot_results.csv'))

plot_results <- read_csv(paste0(simDir, '/plot_results.csv'))

long_plot_results <- plot_results %>%
  dplyr::select(-filename, -estimationDir, -simulationDir,
                -infection_file, -re_file) %>%
  pivot_longer(cols = c(ReRMSE, RootDiff, 
                        EmpCoverage, slopeError, OneTrans),
               names_to = 'measure',
               values_to = 'measure_val') %>%
  mutate(slope = (R2-R1)/(t3-t2),
         noise = paste0('week: ', weekly, '\n norm:', gaussian)) 

method_labels <- setNames(c('RMSE', 'R=1 difference', 
                            'Coverage', 'Slope error', 
                            'Correct R>1 / R<1'),
                          c('ReRMSE', 'RootDiff', 
                            'EmpCoverage', 'slopeError', 'OneTrans') )


##### plot Fig. 1 #####

nonoise_plot_results <- long_plot_results %>%
  filter(weekly == 1,
         gaussian == 0)

ylimits <- data.frame(list('ReRMSE' = c(0,1), 'RootDiff' = c(-2, 4), 
                           'EmpCoverage' = c(0,1), 'slopeError' = c(0, 1), 
                           'OneTrans' = c(0.8, 1)))

stat_plots <- get_plotlist(nonoise_plot_results, ylimits)
final_plot <- get_plotgrid(stat_plots)

plotPath <- paste0(plot_path, '/Fig1B.pdf')
ggsave(plot = final_plot, plotPath, width = 9, height = 12)

##### S2 (Fig 1. with noise) ####

ylimits <- data.frame(list('ReRMSE' = c(0,1), 'RootDiff' = c(-2, 4), 
                           'EmpCoverage' = c(0,1), 'slopeError' = c(0, 1), 
                           'OneTrans' = c(0.7, 1)))



stat_plots <- get_plotlist(long_plot_results, ylimits, col_var = 'noise')
final_plot <- get_plotgrid(stat_plots)

plotPath <- paste0(plot_path, '/Supp_noise.pdf')
ggsave(plot = final_plot, plotPath, width = 15, height = 12)

###########################################################
### Figure S1 (non-smooth version of Fig. 1) #####

raw_cond_grid <- expand.grid(R1 = c(2.5, 3.0, 3.5),
                             init = c(10, 100, 1000),
                             t3 = c(21, 26, 31, 36),
                             t2 = 20,
                             weekly = c(0.25, 0.75, 1.0), 
                             gaussian = c(0.0, 0.1, 0.2),
                             smooth = 'non_smooth')

cond_grid <- complete_cond_grid(raw_cond_grid, defaults)
#estimateGrid(cond_grid)

#non_smooth_valid_cond_grid <- get_valid_cond_grid(cond_grid)
#write_csv(non_smooth_valid_cond_grid, paste0(simDir, '/non_smooth_valid_cond_grid.csv'))

#plot_results <- getPlotResults(non_smooth_valid_cond_grid)
#write_csv(plot_results, paste0(simDir, '/non_smooth_plot_results.csv'))
plot_results <- read_csv(paste0(simDir, '/non_smooth_plot_results.csv'))

long_plot_results <- plot_results %>%
  pivot_longer(cols = c(ReRMSE, RootDiff, 
                        EmpCoverage, slopeError, OneTrans),
               names_to = 'measure',
               values_to = 'measure_val') %>%
  mutate(slope = (R2-R1)/(t3-t2),
         noise = paste0('week: ', weekly, '\n norm:', gaussian)) 

##############
ylimits <- data.frame(list('ReRMSE' = c(0,3), 'RootDiff' = c(-25, 10), 
                           'EmpCoverage' = c(0,0.5), 'slopeError' = c(-3, 3), 
                           'OneTrans' = c(0.5, 1)))

stat_plots <- get_plotlist(long_plot_results, ylimits, col_var = 'noise')
final_plot <- get_plotgrid(stat_plots)

plotPath <- paste0(plot_path, '/Supp_noise_non_smooth.pdf')
ggsave(plot = final_plot, plotPath, width = 15, height = 12)


###########################################################
### Investigating Delay Misspecification - Supplement S3
#####
raw_cond_grid <- expand.grid(R1 = c(2.5, 3.0, 3.5),
                             init = c(10, 100, 1000),
                             t3 = c(21, 26, 31, 36),
                             t2 = 20,
                             weekly = 0.75, 
                             gaussian = 0.1,
                             misdelay = c(-2, -1, 0, 1, 2),
                             delay_opt = c('confirmed', 'death'))

cond_grid <- complete_cond_grid(raw_cond_grid, defaults)
#simulateGrid(cond_grid, smooth_R = TRUE)
#estimateGrid(cond_grid)

#misdelay_valid_cond_grid <- get_valid_cond_grid(cond_grid)
#write_csv(misdelay_valid_cond_grid, paste0(simDir, '/misdelay_valid_cond_grid.csv'))
#misdelay_valid_cond_grid <- read_csv(paste0(simDir, '/misdelay_valid_cond_grid.csv'))

#plot_results <- getPlotResults(misdelay_valid_cond_grid)
#write_csv(plot_results, paste0(simDir, '/misdelay_plot_results.csv'))

plot_results <- read_csv(paste0(simDir, '/misdelay_plot_results.csv'))

#plotSimulations(misdelay_valid_cond_grid)

##############
ylimits <- data.frame(list('ReRMSE' = c(0,1), 'RootDiff' = c(-2, 5), 
                           'EmpCoverage' = c(0,1), 'slopeError' = c(0, 1), 
                           'OneTrans' = c(0.8, 1)))

for (delay_opt_i in c('confirmed', 'death')){
  long_plot_results <- plot_results %>%
    pivot_longer(cols = c(ReRMSE, RootDiff, 
                          EmpCoverage, slopeError, OneTrans),
                 names_to = 'measure',
                 values_to = 'measure_val') %>%
    mutate(slope = (R2-R1)/(t3-t2),
           noise = paste0('week: ', weekly, '\n norm:', gaussian)) %>%
    filter(delay_opt == delay_opt_i)
  
  stat_plots <- get_plotlist(long_plot_results, ylimits, col_var = 'misdelay')
  

  (final_plot <- get_plotgrid(stat_plots))
  
  plotPath <- paste0(plot_path, '/Supp_delay_', delay_opt_i, '.pdf')
  ggsave(plotPath, width = 15, height = 12)
}
###########################################################
### Fixed Delay - Supplement S5
######
raw_cond_grid <- expand.grid(R1 = c(2.5, 3.0, 3.5),
                             init = c(10, 100, 1000),
                             t3 = c(21, 26, 31, 36),
                             t2 = 20,
                             weekly = 0.75, 
                             gaussian = 0.1,
                             misdelay = 0,
                             fixed = c(TRUE, FALSE),
                             delay_opt = 'confirmed')

cond_grid <- complete_cond_grid(raw_cond_grid, defaults)

# The fixed delay is used in the estimation step,
# but does not change the simulations -> thus the filename
# should refer to the location of the fixed_FALSE simulation
# files. This line does not change the location of the infection and Re files
cond_grid <- cond_grid %>%
  mutate(filename = sub('fixed_TRUE', 'fixed_FALSE', filename))

#estimateGrid(cond_grid)

#fixed_valid_cond_grid <- get_valid_cond_grid(cond_grid)
#write_csv(fixed_valid_cond_grid, paste0(simDir, '/fixed_valid_cond_grid.csv'))
#fixed_valid_cond_grid <- read_csv(paste0(simDir, '/fixed_valid_cond_grid.csv'))

#plot_results <- getPlotResults(fixed_valid_cond_grid)
#write_csv(plot_results, paste0(simDir, '/fixed_plot_results.csv'))
plot_results <- read_csv(paste0(simDir, '/fixed_plot_results.csv'))

#plotSimulations(fixed_valid_cond_grid)

long_plot_results <- plot_results %>%
  pivot_longer(cols = c(ReRMSE, RootDiff, 
                        EmpCoverage, slopeError, OneTrans),
               names_to = 'measure',
               values_to = 'measure_val') %>%
  mutate(slope = (R2-R1)/(t3-t2),
         noise = paste0('week: ', weekly, '\n norm:', gaussian))

##############
ylimits <- data.frame(list('ReRMSE' = c(0,1), 'RootDiff' = c(0, 6), 
                           'EmpCoverage' = c(0,1), 'slopeError' = c(0, 1), 
                           'OneTrans' = c(0.7, 1)))


stat_plots <- get_plotlist(long_plot_results, ylimits, col_var = 'fixed')

  
(final_plot <- get_plotgrid(stat_plots))
  
plotPath <- paste0(plot_path, '/Supp_fixed.pdf')
ggsave(plotPath, width = 13, height = 12)

###########################################################
### Timevarying Delay - Supplement S4
#######
raw_cond_grid <- expand.grid(R1 = c(2.5, 3.0, 3.5),
                             init = c(10, 100, 1000),
                             t3 = c(21, 26, 31, 36),
                             t2 = 20,
                             weekly = 0.75, 
                             gaussian = 0.1,
                             misdelay = 0,
                             timevar_sim = TRUE,
                             timevar_est = c(TRUE, FALSE),
                             delay_opt = c('confirmed', 'death'))

cond_grid <- complete_cond_grid(raw_cond_grid, defaults)

#simulateGrid(cond_grid)
#estimateGrid(cond_grid)

#timevar_valid_cond_grid <- get_valid_cond_grid(cond_grid)
#write_csv(timevar_valid_cond_grid, paste0(simDir, '/timevar_valid_cond_grid.csv'))
#timevar_valid_cond_grid <- read_csv(paste0(simDir, '/timevar_valid_cond_grid.csv'))

#plot_results <- getPlotResults(timevar_valid_cond_grid)
#write_csv(plot_results, paste0(simDir, '/timevar_plot_results.csv'))
plot_results <- read_csv(paste0(simDir, '/timevar_plot_results.csv'))

#plotSimulations(timevar_valid_cond_grid)

##############################################

ylimits <- data.frame(list('ReRMSE' = c(0,1), 'RootDiff' = c(-3, 3), 
                           'EmpCoverage' = c(0,1), 'slopeError' = c(-1, 1), 
                           'OneTrans' = c(0.7, 1)))

for (delay_opt_i in c('confirmed', 'death')){
  long_plot_results <- plot_results %>%
    pivot_longer(cols = c(ReRMSE, RootDiff, 
                          EmpCoverage, slopeError, OneTrans),
                 names_to = 'measure',
                 values_to = 'measure_val') %>%
    mutate(slope = (R2-R1)/(t3-t2),
           noise = paste0('week: ', weekly, '\n norm:', gaussian)) %>%
    filter(delay_opt == delay_opt_i)
  
  stat_plots <- get_plotlist(long_plot_results, ylimits, col_var = 'timevar_est')
  
  (final_plot <- get_plotgrid(stat_plots))
  
  plotPath <- paste0(plot_path, '/Supp_timevar_', delay_opt_i, '.pdf')
  ggsave(plotPath, width = 15, height = 12)
  
}
