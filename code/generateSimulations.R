###########################################################
## generateSimulations.R
## author: J.S. Huisman; some functions by Jinzhou Li
###########################################################

#plot_path = '/Users/jana/Documents/PhD/covid/Re_project/figures'
appDir = "/Users/jana/Documents/PhD/covid/covid-19-re-shiny-app"

library("fitdistrplus")
library(tidyverse)
library(EpiEstim)

source(paste0(appDir, '/app/otherScripts/2_utils_getInfectionIncidence.R'))
source(paste0(appDir, '/app/otherScripts/3_utils_doReEstimates.R'))

# The most important function is simulateTS - used to simulate
# infections and observations from a given Re trajectory/scenario

# To estimate infections from observations (using deconvolution)
# and Re from infections; we use the functions from our general method
# these are read in from the scripts '2_...' and '3_...' above;
# and require the wrappers estimateInfectionTS and estimateReTS
# to work with the simulated data

###########################################################
###### Generate Re TS ######
# input:
# - t1-t6
# - Re1-3
# output:
# - Re timeseries
getReTS <- function(shift_times, R_levels){
  t_norm = shift_times - shift_times[1]
  n_t = length(t_norm)
  
  ReTS <- vector(mode = 'numeric', length = t_norm[n_t])
  
  R_count = 1
  for (t_i in 1:(n_t-1)){
    if (t_i %% 2 != 0){
      # odd
      ReTS[(t_norm[t_i]+1):t_norm[t_i+1]] <- R_levels[R_count]
      R_count = R_count + 1
    } else{
      # even
      ReTS[t_norm[t_i]:(t_norm[t_i+1]+1)] <- seq(from = R_levels[R_count-1], to = R_levels[R_count], 
                                           length.out = (t_norm[t_i+1]-t_norm[t_i]+2)) 
    }
  }
  
  return(ReTS)
}

# shift_times = c(0, 30, 40, 50, 60, 80)
# R_levels = c(4, 0.5, 1.2)
# ReTS <- getReTS(shift_times, R_levels)
# plot(ReTS)


getSmoothReTS <- function(shift_times, R_levels, days_incl = 14) {
  
  ReTS <- getReTS(shift_times, R_levels)

  n_points <- length(ReTS)
  sel_span <- days_incl / n_points
  
  df <- data.frame(date = 1:n_points, R = ReTS)
  c_data.lo <- loess(R ~ date, data = df, span = sel_span, degree = 1)
  
  smoothed_R <- predict(c_data.lo)

  return(smoothed_R)
}


###### Generate Infection TS ######

getDiscreteSIperDay <- function(k, shapeG=2.73, scaleG=1.39) {
  ### Expression from Cori et al. 2013, Web appendix 11
  #serial interval SI is such that SI -1 is Gamma distributed
  wk <- k * pgamma(k, shape=shapeG, scale=scaleG) +
    (k-2) * pgamma(k-2, shape=shapeG, scale=scaleG) +
    (-2) * (k-1) * pgamma(k-1, shape=shapeG, scale=scaleG) +
    shapeG * scaleG * (2 * pgamma(k-1, shape=shapeG+1, scale=scaleG) - 
                           pgamma(k-2, shape=shapeG+1, scale=scaleG) - 
                           pgamma(k, shape=shapeG+1, scale=scaleG)) 
  return(wk)
}

getDiscreteSerialInterval <- function(shapeG=2.73, scaleG=1.39){
  longDSI <- sapply(0:1000, getDiscreteSIperDay, shapeG, scaleG)
  cutDSI <- min(which(longDSI[2:1000]==0))
  DiscreteSerialInterval <- longDSI[1:cutDSI]
  
  return(DiscreteSerialInterval)
}

getInfectionsDayT <- function(RT, InfectionsBeforeT, DiscreteSerialInterval, migration = T){
  # Assumption: DSI starts from 0; here we use from day 1
  memory = min(length(InfectionsBeforeT), length(DiscreteSerialInterval)-1)
  LambdaDayT = RT * sum(rev(InfectionsBeforeT)[1:memory]*
                              DiscreteSerialInterval[2:(memory+1)])
  
  InfectionsDayT = rpois(1, LambdaDayT)
  if(is.na(InfectionsDayT)){
    InfectionsDayT = 0
    #print(paste0('error for ', LambdaDayT))
    # happens occassionally for very small negative values (numerical error)
  }
  
  # if (migration & (runif(1, min = 0, max = 1) < 0.01)){
  #   InfectionsDayT = InfectionsDayT + 1
  # }
  
  return(InfectionsDayT)
}

# input:
# - timeseries of Re
# - serial interval distribution
# output:
# - timeseries of infections
getInfectionTS <- function(ReTS, init_infection = 1, ...){
  # can supply shapeG and scaleG to pass to DSI function
  n_ts = length(ReTS)
  n_init = length(init_infection)
  DiscreteSerialInterval = getDiscreteSerialInterval(...)
  
  infectionTS = vector(mode="numeric", length = n_ts + n_init - 1)
  infectionTS[1:n_init] = init_infection
  
  for(i in 2:n_ts){
    infectionTS[n_init-1+i] = getInfectionsDayT(ReTS[i], 
                                       infectionTS[1:(n_init+i-2)], 
                                       DiscreteSerialInterval)
  }
  
  infectionTS = infectionTS[n_init:(n_init+n_ts-1)]
  return(infectionTS)
}

# infectionTS <- getInfectionTS(ReTS)
# plot(infectionTS)

###### Generate Observation TS ######

getGammaParams <- function(meanParam, sdParam){
  shapeParam <- meanParam^2 / (sdParam^2)
  scaleParam <- (sdParam^2) / meanParam
  return(list(shape = shapeParam, scale = scaleParam))
}

getInvGammaParams <- function(shapeParam, scaleParam){
  meanParam <- scaleParam * shapeParam
  sdParam <- sqrt(scaleParam^2 * shapeParam)
  return(list(mean = meanParam, sd = sdParam))
}

drawDoubleGamma <- function(n_draws, IncubationParams, OnsetToCountParams){
  draws <- round(rgamma(n_draws, shape = OnsetToCountParams$shape,
           scale = OnsetToCountParams$scale) +
     rgamma(n_draws, shape = IncubationParams$shape,
             scale = IncubationParams$scale))
  return(draws)
}

#library(actuar)
# makeGammaECDF <- function(GammaParams1, GammaParams2, numberOfSamples = 1E6) {
#   draws <- rgamma(numberOfSamples, shape = GammaParams1$shape, scale = GammaParams1$scale) +
#     rgamma(numberOfSamples, shape = GammaParams2$shape, scale = GammaParams2$scale)
#   return(Vectorize(ecdf(draws)))
# }

# drawDoubleGamma <- function(n_draws, IncubationParams, OnsetToCountParams){
#   Fhat <- makeGammaECDF(IncubationParams, OnsetToCountParams)
#   
#   # the first entry corresponds to delay = 0
#   discrete_double_gamma <- discretize(Fhat, method = "upper",
#                                       from = -0.5, to = 101, step = 1)
#   
#   draws = sample(x = 0:100, size = n_draws, replace = TRUE, prob = discrete_double_gamma)
#   
#   return(draws)
# }



addTSnoise <- function(observationTS, noise){
  origObsTS <- observationTS
  allObservationDates = 1:length(origObsTS)
  
  if ('weekly' %in% names(noise)){
    # we just care about the 7-day pattern, not which day exactly is a weekend
    
    saturdays = which((allObservationDates %% 7) == 0)
    sat_reduct <- rnorm(length(saturdays), noise$weekly, noise$weekly/10)
    sat_reduct[sat_reduct>1] = 1
    sat_reduct[sat_reduct<0] = 0
    
    observationTS[saturdays] = sat_reduct * origObsTS[saturdays]
    observationTS[saturdays+2] = origObsTS[saturdays+2] + (1 - sat_reduct) * origObsTS[saturdays]
    
    sundays = which((allObservationDates %% 7) == 1)
    sun_reduct <- rnorm(length(sundays), noise$weekly, noise$weekly/10)
    sun_reduct[sun_reduct>1] = 1
    sun_reduct[sun_reduct<0] = 0
    
    observationTS[sundays] = sun_reduct * origObsTS[sundays]
    observationTS[sundays+2] = origObsTS[sundays+2] + (1 - sun_reduct) * origObsTS[sundays]
  }
  
  if ('gaussian' %in% names(noise)){
    mult_noise <- rnorm(length(observationTS), mean = 1, sd = noise$gaussian)
    observationTS = mult_noise * observationTS
  }
  
  if ('fitted_noise_model' %in% names(noise)){
    mult_noise <- simulate(noise$fitted_noise_model)     # noise model is on log-scale
    
    if(length(mult_noise) < length(observationTS)) {print("length(mult_noise) < observationTS in the fitted time series model!!!")}
    
    # y  =  mu * residual
    observationTS = observationTS * exp(mult_noise[1:length(observationTS)])
  }
  
  if ('iid_noise_sd' %in% names(noise)){
    mult_noise <- rnorm(length(observationTS), mean = 0, sd = noise$iid_noise_sd)
    
    # y  =  mu * residual
    observationTS = observationTS * exp(mult_noise)     # so the error is iid log-normal
  }
  
  if ('noiseless' %in% names(noise)){   # used for SemiPara Boot
    observationTS = observationTS
  }
  
  observationTS = round(observationTS)
  observationTS[observationTS < 0] = 0
  
  return(observationTS)
}


# timevaryingDelayDist <- function(incidenceOnDay, OnsetToCountParams,
#                                  first = 1.2, second = 1.5){
#   
#   invParams = getInvGammaParams(shapeParam = OnsetToCountParams$shape, 
#                     scaleParam = OnsetToCountParams$scale)
#   
#   # if (incidenceOnDay > 1000){
#   #   #mean is 1.2 times greater
#   #   NewOnsetToCountParams = list(shape = (first^2)*OnsetToCountParams$shape, 
#   #                                scale = OnsetToCountParams$scale/first)
#   # } else if (incidenceOnDay > 10000){
#   #   NewOnsetToCountParams = list(shape = (second^2)*OnsetToCountParams$shape, 
#   #                                scale = OnsetToCountParams$scale/second)
#   # } else{
#   #   NewOnsetToCountParams <- OnsetToCountParams
#   # }
#   
#   # 2 days longer processing time at 10000 cases
#   newMean = invParams$mean + 2e-04*incidenceOnDay
#   
#   NewOnsetToCountParams <- getGammaParams(newMean, invParams$sd)
#   
#   return(NewOnsetToCountParams)
# }

timevaryingDelayDist <- function(onset_date, OnsetToCountParams){
  
  invParams = getInvGammaParams(shapeParam = OnsetToCountParams$shape, 
                                scaleParam = OnsetToCountParams$scale)
  
  # 1/20 less per day; no lower than 2 days
  newMean = max(invParams$mean - onset_date/20.0, 2.0)
  
  NewOnsetToCountParams <- getGammaParams(newMean, invParams$sd)
  
  return(NewOnsetToCountParams)
}

# input
# - timeseries of infections
# - delay distribution
# output
# - time series of observations
getObservationTS <- function(infectionTS, IncubationParams, OnsetToCountParams,
                             noise = list('weekly' = 0.3, 'gaussian' = 0.1),
                             timevarying = FALSE, truncate = 'empirical'){
  
  n_ts = length(infectionTS)

  n_infect_tot <- sum(infectionTS)
  ObservationDates <- c()
  #infection_date = 1
  for (infection_date in 1:n_ts) {
     
    if (infectionTS[infection_date] > 0){
      if (timevarying){
        sampledOnsets <- round(rgamma(infectionTS[infection_date], 
                                      shape = IncubationParams$shape, 
                                      scale = IncubationParams$scale))
        
        drawnOnsetDates <- infection_date + sampledOnsets
          
        sampledDelays <-  sapply(drawnOnsetDates, function(onset_date){
          delay_dist = timevaryingDelayDist(onset_date, OnsetToCountParams)
          return( round(rgamma(1, 
                       shape = delay_dist$shape, 
                       scale = delay_dist$scale)) )
          })
        
        drawnCountDates <- drawnOnsetDates + sampledDelays
        
      } else {
        sampledDelays <- drawDoubleGamma(infectionTS[infection_date],
                                         IncubationParams,
                                         OnsetToCountParams)
        
        drawnCountDates <- infection_date + sampledDelays
      }
      ObservationDates <- c(ObservationDates, drawnCountDates)
    }
  }
  
  allObservationDates <- seq(1, max(ObservationDates, n_ts))
  # -1 because we add all dates
  observationTS <- unname( table( c(ObservationDates, allObservationDates) ) ) -1
  observationTS <- as.numeric(observationTS)
  
  # observationTS <- c()
  # for(i in 1:length(allObservationDates)) {
  #   observationTS[i] <- sum(ObservationDates == allObservationDates[i])
  # }

  #origObsTS <- observationTS
  
  #observationTS <- addTSnoise(observationTS, noise)
  observationTS_cut <- observationTS[1:length(infectionTS)]
  observationTS <- addTSnoise(observationTS_cut, noise)
  
  #Cut-off the timeseries a bit earlier than the last observation
  # if (truncate == 'empirical'){
  #   delay_sample <- drawDoubleGamma(10000,
  #                                   IncubationParams,
  #                                   OnsetToCountParams)
  #   extension = min(n_ts + median(delay_sample), length(origObsTS))
  #   observationTS = observationTS[1:extension]
  # }

  
  return(observationTS)
}


getEmpCountDelayDist <- function(infectionTS, IncubationParams, OnsetToCountParams,
                                 subsample = 0.4){
  n_ts = length(infectionTS)
  
  #delays_data_path <- "/Users/jana/Documents/PhD/covid/covid-19-re-shiny-app/app/data/CH/FOPH_data_delays.csv"
  delay_df <- data.frame()
  for (infection_date in 1:n_ts) {
    
    # not sure this is the right distribution to use
    # just want to make sure it is occasionally above 0.5 such that
    # days with incidence of 1 stand a chance of being sampled
    subsample_frac = min(abs(rnorm(1, mean = subsample, sd = 1)), 1)
    
    if (round(infectionTS[infection_date]*subsample_frac) > 0){
      sampledOnsets <- round(rgamma(infectionTS[infection_date], 
                                    shape = IncubationParams$shape, 
                                    scale = IncubationParams$scale))
      
      drawnOnsetDates <- infection_date + sampledOnsets
      
      sampledDelays <-  sapply(drawnOnsetDates, function(onset_date){
        delay_dist = timevaryingDelayDist(onset_date, OnsetToCountParams)
        return( round(rgamma(1, 
                             shape = delay_dist$shape, 
                             scale = delay_dist$scale)) )
      })
      
      drawnCountDates <- drawnOnsetDates + sampledDelays
      
      
      new_delay_df <- data.frame(data_type = 'Simulated',
                                 onset_date = as_date("2020-02-01") + drawnOnsetDates,
                                 count_date = as_date("2020-02-01") + drawnCountDates,
                                 delay = sampledDelays,
                                 source = 'ETH',
                                 region = 'Simulated',
                                 country = 'Simulated',
                                 local_infection = TRUE)
      
      delay_df = bind_rows(delay_df, new_delay_df)
    }
  }
  return(delay_df)
}


# IncubationParams <- getGammaParams(meanParam = 5.3, sdParam = 3.2)
# # onset to death: mean =15.0 sd=6.9 (Linton et al. best gamma distr fit)
# OnsetToCountParams = getGammaParams(15.0, 6.9)
# 
# ObservationTS <- getObservationTS(infectionTS, 
#                                   IncubationParams, OnsetToCountParams)
# 
# plot(ObservationTS)

######################################################################################################################
###### Simulation Master Script ######
simulateTS <- function(shift_times, R_levels,
                       IncubationParams, OnsetToCountParams, noise = list(),
                       timevarying = FALSE, smooth_R = FALSE,
                       ...){
  # With the dots one can specify shapeG, scaleG
  
  if (smooth_R){
    ReTS <- getSmoothReTS(shift_times, R_levels)
  } else{
    ReTS <- getReTS(shift_times, R_levels)
  }
  
  infectionTS <- getInfectionTS(ReTS, ...)
  ObservationTS <- getObservationTS(infectionTS, 
                                    IncubationParams, OnsetToCountParams, noise = noise,
                                    timevarying = timevarying)
  
  ## the observation timeseries will be longer by a few days; we cut this off
  n_re <- length(ReTS)
  #n_obs <- length(ObservationTS)
  result <- data.frame(list(Re = ReTS, 
                            infections = infectionTS, 
                            observations = ObservationTS[1:n_re],
                            date = as_date("2020-02-01") + 1:n_re) )
  return(result)
}

# simulation <- simulateTS(t, R,
#           IncubationParams, OnsetToCountParams, init_infection = 1000)
# 
# longSim <- simulation %>%
#   pivot_longer(cols = c(Re, infections, observations),
#                names_to = 'type') %>%
#   mutate(type = factor(type, levels = c('Re', 'infections', 'observations')))
# 
# ggplot(longSim) +
#   geom_line(aes(x = date, y = value, colour = type)) +
#   facet_grid(rows = vars(type), scale = 'free')

######################################################################################################################
###### Validate Deconvolution ######

getSimIncidence <- function(simulation){
  
  observation_df <- simulation %>%
    dplyr::select(date, value = observations) %>%
    mutate(data_type = 'Simulated',
           source = 'ETH',
           variable = 'incidence',
           region = 'Simulated',
           country = 'Simulated',
           date_type = 'report',
           local_infection = TRUE)
  
  return(observation_df)
}




get_infection_incidence_by_shift <- function(
  data_subset,
  constant_delay_distribution,
  smooth_incidence = T,
  empirical_delays  = tibble(),
  n_bootstrap = 5) {
  
  is_empirical = (nrow(empirical_delays) > 0)
  
  # use mode of delay_distribution. -1 because indices are offset by one as the delay can be 0.
  if(is_empirical) {
    delay_dist_matrix <- get_matrix_empirical_waiting_time_distr(
      empirical_delays, seq(min(data_subset$date)-24, max(data_subset$date), by = "days") )
    
    #lapply(1:ncol(delay_dist_matrix), function(x){which.max(delay_dist_matrix[, x]) - 1})
    
    # this should be improved! now it is not different for every shift
    n_shift = which.max(delay_dist_matrix[, 1])  - 1
  } else {
    n_shift = which.max(constant_delay_distribution)  - 1
  }
  
  results <- list(tibble())
  for (bootstrap_replicate_i in 0:n_bootstrap) {
    
    if (bootstrap_replicate_i == 0) {
      time_series <- data_subset
    } else {
      time_series <- get_bootstrap_replicate(data_subset)
    }
    
    if (smooth_incidence == T) {
      smoothed_incidence_data <- time_series %>%
        complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>% 
        mutate(value = getLOESSCases(dates = date, count_data = value))
      
      raw_total_incidence <- sum(time_series$value, na.rm = TRUE)
      smoothed_total_incidence <- sum(smoothed_incidence_data$value, na.rm = T)
      
      if (smoothed_total_incidence > 0) {
        smoothed_incidence_data <- smoothed_incidence_data %>%
          mutate(value = value * raw_total_incidence / smoothed_total_incidence)
      }
      
    } else {
      smoothed_incidence_data <- time_series  %>%
        complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0))
    }
    
    shifted_infections <- smoothed_incidence_data %>%
      mutate(date = date - n_shift)
    
    ## prepare metadata for result tibble
    data_type_subset <- unique(time_series$data_type)[1]
    data_type_name <- paste0("infection_", data_type_subset)
    
    ## dataframe containing results
    shifted_infections <- tibble(
      date = shifted_infections$date,
      region = unique(time_series$region)[1],
      country = unique(time_series$country)[1],
      source = unique(time_series$source)[1],
      data_type = data_type_name,
      replicate = bootstrap_replicate_i,
      value = shifted_infections$value,
      variable = "incidence",
      local_infection = T
    )
    
    results <- c(results, list(shifted_infections))
  }
  
  return(bind_rows(results))
}

###########################################################
# test deconvolution against "real" infection curve
estimateInfectionTS <- function(simulation, IncubationParams, OnsetToCountParams,
                                smooth_param = FALSE, fixed_shift = FALSE,
                                timevarying = FALSE, n_boot = 50){
  infection_df <- getSimIncidence(simulation)
  
  constant_delay_distributions <- list("Simulated" = get_vector_constant_waiting_time_distr(
    IncubationParams$shape, IncubationParams$scale,
    OnsetToCountParams$shape, OnsetToCountParams$scale),
    "Symptoms" = get_vector_constant_waiting_time_distr(
      IncubationParams$shape, IncubationParams$scale,
      0, 0))
  
  # load empirical delays
  if (timevarying){
    delays_onset_to_count <- getEmpCountDelayDist(infection_df$value, IncubationParams, OnsetToCountParams)
  } else {
    delays_onset_to_count <- tibble()  
  }
  
  if (fixed_shift){
    estimatedInfections <- get_infection_incidence_by_shift(
      infection_df,
      constant_delay_distribution = constant_delay_distributions[['Simulated']],
      smooth_incidence = smooth_param,
      empirical_delays  = delays_onset_to_count,
      n_bootstrap = n_boot)
    
  } else{
    estimatedInfections <- get_infection_incidence_by_deconvolution(
      infection_df,
      is_local_cases = T,
      constant_delay_distribution = constant_delay_distributions[['Simulated']],
      constant_delay_distribution_incubation = constant_delay_distributions[["Symptoms"]],
      max_iterations = 100,
      smooth_incidence = smooth_param,
      empirical_delays = delays_onset_to_count,
      n_bootstrap = n_boot,
      verbose = FALSE)
    
  }
  
  return(estimatedInfections)
}

#estimatedInfections <- estimateInfectionTS(simulation, IncubationParams, OnsetToCountParams)
######################################################################################################################
###### Validate Re Estimation ######


# test Re estimates against "real" Re estimates

estimateReTS <- function(estimatedInfections, delay = 0){
  all_delays <- list("infection_Simulated" = c(Cori = delay))
  
  truncations <- list(left = c(Cori = 5), right = c(Cori = 0))
  
  
  rawEstimatedRe <- doAllReEstimations(
    estimatedInfections,
    slidingWindow = 3,
    #slidingWindow = 1,
    methods = "Cori",
    variationTypes = c("slidingWindow"),
    all_delays = all_delays,
    truncations = truncations,
    interval_ends = max(estimatedInfections$date),
    additional_interval_ends = data.frame(region = "Simulated",
                                          data_type = "Simulated",
                                          date = max(estimatedInfections$date)))
  
  estimatedRe <- as_tibble(rawEstimatedRe)

  return(estimatedRe)
}

cleanReTSestimate <- function(rawEstimatedRe){
  estimatedRe <- as_tibble(rawEstimatedRe) %>%
    pivot_wider(names_from = "variable", values_from = "value") %>%
    group_by(date, country, region, data_type, source, estimate_type) %>%
    summarize(
      median_R_mean = median(R_mean),
      median_R_highHPD = median(R_highHPD),
      median_R_lowHPD = median(R_lowHPD),
      .groups = "keep"
    ) %>%
    dplyr::select(country, region, source, data_type, estimate_type, date,
           median_R_mean, median_R_highHPD, median_R_lowHPD) %>%
    arrange(country, region, source, data_type, estimate_type, date) %>%
    ungroup()
  return(estimatedRe)
}

#estimatedRe <- estimateReTS(estimatedInfections, delay = 0)

###########################################################

# getSimulationResults <- function(shift_times, R_levels, 
#                                  IncubationParams, OnsetToCountParams,
#                                  smooth_param = TRUE, delay = 0){
#   # Simulate and Estimate
#   simulation <- simulateTS(shift_times, R_levels,
#                            IncubationParams, OnsetToCountParams)
#   estimatedInfections <- estimateInfectionTS(simulation, IncubationParams, 
#                                              OnsetToCountParams,
#                                              smooth_param) 
#   estimatedRe <- estimateReTS(estimatedInfections, delay)
#   
#   return(list(simulation = simulation,
#               estimatedInfections = estimatedInfections,
#               estimatedRe = estimatedRe))
# }
######################################################################################################################
###### Compare Simulation and Estimation results ######

# getReRMSE <- function(ReCompare){
#   ReCompare <- ReCompare %>%
#     #mutate(SE = (Re - median_R_mean)^2)
#     mutate(SE = ( median_R_mean - Re)^2)
#   
#   rmse <- sqrt(sum(ReCompare$SE)/length(ReCompare$SE))
#   norm_rmse <- rmse / mean(ReCompare$Re)
#   #return(rmse)
#   return(norm_rmse)
# }
# 
# getRelError <- function(ReCompare){
#   ReCompare <- ReCompare %>%
#     filter(Re != 0 ) %>%
#     mutate(rel_err = (Re - median_R_mean)/Re)
#   
#   rel_err <- mean(ReCompare$rel_err)
#   return(rel_err)
# }
# 
# # getRootDiff <- function(ReCompare, timewindow = NULL){
# #   if(is.null(timewindow)){
# #     timewindow = min(which(diff(ReCompare$Re) != 0)):length(ReCompare$date)
# #   }
# #   
# #   sim_root <- ReCompare$date[timewindow][min(which(ReCompare$Re[timewindow] < 1))]
# #   est_root <- ReCompare$date[timewindow][min(which(ReCompare$median_R_mean[timewindow] < 1))]
# #   
# #   root_diff <- as.numeric(est_root - sim_root)
# #   return(root_diff)
# # }
# 
# getRootDiff <- function(simulation, estimatedRe, all = FALSE){
# 
#   sim_change_points = diff(sign(simulation$Re -1))  
#   sim_roots = simulation$date[which(sim_change_points != 0 )]
# 
#   if (! ('median_R_mean' %in% colnames(estimatedRe))){
#     estimatedRe <- cleanReTSestimate(estimatedRe)
#   }
#   
#   est_change_points = diff(sign(estimatedRe$median_R_mean -1))  
#   est_roots = estimatedRe$date[which(est_change_points != 0 )]
#   
#   root_diffs <- as.numeric(est_roots - sim_roots)
#   
#   if (all){
#     root_diff = sum(abs(root_diffs))
#   } else{
#     root_diff = root_diffs[1]
#   }
#   return( root_diff )
# }
# 
# getEmpCoverage <- function(ReCompare){
#   covCompare <- ReCompare %>%
#     mutate(coverage = (Re >= median_R_lowHPD) & (Re <= median_R_highHPD) ,
#            CI_width = median_R_highHPD - median_R_lowHPD)
# 
#   frac_coverage <- sum(covCompare$coverage)/nrow(covCompare)
#   median_CI_width <- median(covCompare$CI_width)
#   
#   return(list(frac_coverage = frac_coverage, CI_width = median_CI_width))
# }
# 
# library(dtwclust)
# library(SimilarityMeasures)
# 
# getReTSdist <- function(ReCompare, Fdist = FALSE){
#   Re_est_ts <- ts(ReCompare %>% pull(median_R_mean), start = min(ReCompare$date), end = max(ReCompare$date))
#   Re_sim_ts <- ts(ReCompare %>% pull(Re), start = min(ReCompare$date), end = max(ReCompare$date))
#   
#   DTWdist <- dtw(Re_sim_ts, Re_est_ts, distance.only = TRUE)$normalizedDistance
#   
#   if(!Fdist){
#     return(list(DTWdist = DTWdist))
#   } else {
#     FrechetDist <- Frechet(matrix(ReCompare %>% pull(median_R_mean)), 
#                            matrix(ReCompare %>% pull(Re)))
#     
#     return(list(DTWdist = DTWdist, Fdist = FrechetDist))
#   }
#   
# }
# 
# getReError <- function(simulation, estimatedRe, Fdist = FALSE){
#   
#   if (! ('median_R_mean' %in% colnames(estimatedRe))){
#     estimatedRe <- cleanReTSestimate(estimatedRe)
#   }
#   
#   #Compare
#   ReCompare <- simulation %>%
#     full_join(estimatedRe, by = c('date')) %>%
#     filter(!is.na(median_R_mean))
#   
#   ReRMSE <- getReRMSE(ReCompare)
#   RelError <- getRelError(ReCompare)
#   RootDiff <- getRootDiff(simulation, estimatedRe, all = FALSE)
#   EmpCoverage <- getEmpCoverage(ReCompare)
#   ReTSdist <- getReTSdist(ReCompare, Fdist)
#   OneTrans <- getOneTransitions(simulation, estimatedRe)
#   
#   if (!Fdist){
#     return(list(ReRMSE = ReRMSE, RelError = RelError, RootDiff = RootDiff, 
#                 EmpCoverage = EmpCoverage$frac_coverage,
#                 CIWidth = EmpCoverage$CI_width,
#                 OneTrans = OneTrans,
#                 DTWdist = ReTSdist$DTWdist))
#   } else {
#     return(list(ReRMSE = ReRMSE, RelError = RelError, RootDiff = RootDiff, 
#                 EmpCoverage = EmpCoverage$frac_coverage,
#                 CIWidth = EmpCoverage$CI_width,
#                 OneTrans = OneTrans,
#                 DTWdist = ReTSdist$DTWdist,
#                 Fdist = ReTSdist$Fdist))
#   }
#  
# }
#  
# 
# getSlopeError <- function(simulation, estimatedRe, valid_cond){
#   date_t3 = simulation$date[valid_cond$t3]
#   date_t2 = simulation$date[valid_cond$t2]
#   date_first = max(min(estimatedRe$date), date_t2)
#   date_int = as.numeric(date_t3 - date_first)
#   
#   slope_sim = (simulation$Re[simulation$date == date_t3] -
#                  simulation$Re[simulation$date == date_first])/date_int
#   
#   meanEstRe <- estimatedRe %>%
#     filter(variable == 'R_mean') %>%
#     dplyr::select(date, replicate, value) %>%
#     group_by(replicate) %>%
#     summarise(slope_est = (value[date == date_t3] -
#                              value[date == date_first])/date_int,
#               .groups = "drop")
#   
#   abs_slope_error = slope_sim - meanEstRe$slope_est
#   rel_slope_error = (abs_slope_error/slope_sim)
#   #return(abs_slope_error)
#   return(rel_slope_error)
# }
# 
# 
# getMeanSlopeError <- function(valid_cond_grid){
#   SlopeError = data.frame()
#   for (row_id in 1:nrow(valid_cond_grid)){
#     
#     simulation <- read_csv(paste0(valid_cond_grid[row_id, 'simulationDir'], '/', 
#                                   valid_cond_grid[row_id, 'filename']))
#     
#     estimatedInfections <- read_csv(paste0(valid_cond_grid[row_id, 'estimationDir'],
#                                            valid_cond_grid[row_id, 'infection_file']))
#     
#     estimatedRe <- read_csv(paste0(valid_cond_grid[row_id, 'estimationDir'], 
#                                    valid_cond_grid[row_id, 're_file']))
#     
#     new_error <- getSlopeError(simulation, estimatedRe, valid_cond_grid[row_id, ])
#     SlopeError <- bind_rows(SlopeError, list(slopeError = mean(new_error)) )
#   }
#   return(SlopeError)
# }
# 
# 
# getOneTransitions <- function(simulation, estimatedRe){
#   simulation <- simulation %>%
#     mutate(OneTrans = ifelse(Re >= 1, 1, 0))
#   
#   if (! ('median_R_mean' %in% colnames(estimatedRe))){
#     estimatedRe <- cleanReTSestimate(estimatedRe)
#   }
#   
#   one_transitions = estimatedRe %>%
#     mutate(EstOneTrans = case_when(median_R_lowHPD >= 1 ~ 1,
#                                     median_R_highHPD < 1 ~ 0)) %>%
#     full_join(simulation, by = c('date')) %>%
#     dplyr::select(date, Re, OneTrans, EstOneTrans) %>%
#     filter(!is.na(EstOneTrans)) %>%
#     mutate(compare = OneTrans == EstOneTrans)
#   
#   perc_transitions = sum(one_transitions$compare)/nrow(one_transitions)
#     
#   return(perc_transitions)
# }
# 
# 
# ###### Investigate Bootstrap Residuals ######
# plotBootstrapResid <- function(simulation, estimatedRe){
#   sim_subset <- simulation %>%
#     dplyr::select(-c(infections, observations))
#   
#   est_subset <- estimatedRe %>%
#     filter(variable == 'R_mean') %>%
#     dplyr::select(date, replicate, Re_est = value) %>%
#     left_join(sim_subset, by = 'date') %>%
#     mutate(resid = sqrt(nrow(.))*(Re_est - Re) )
# 
#   first_boot_set <- estimatedRe %>%
#     filter(variable == 'R_mean',
#            replicate == 0) %>%
#     dplyr::select(date, first_boot = value)
#   
#   boot_set <- estimatedRe %>%
#     filter(variable == 'R_mean',
#            replicate != 0) %>%
#     dplyr::select(date, replicate, value) %>%
#     left_join(first_boot_set, by = 'date') %>%
#     mutate(resid = sqrt(nrow(.))*(value - first_boot) )
#   
#   p <- ggplot() +
#     geom_histogram(data = boot_set,
#                    aes(resid, stat(density)), binwidth = 1, alpha = 0.4, fill = 'black') +
#   geom_histogram(data = est_subset,
#                     aes(resid, stat(density)), binwidth = 1, alpha = 0.4, fill = 'lightblue')
#   return(p)
# }
# 
# checkBootstrapNormal <- function(simulation, estimatedRe){
#   sim_subset <- simulation %>%
#     dplyr::select(-c(infections, observations))
#   
#   est_subset <- estimatedRe %>%
#     filter(variable == 'R_mean') %>%
#     dplyr::select(date, replicate, Re_est = value) %>%
#     left_join(sim_subset, by = 'date') %>%
#     mutate(resid = sqrt(nrow(.))*(Re_est - Re) )
# 
#   
#   first_boot_set <- estimatedRe %>%
#     filter(variable == 'R_mean',
#            replicate == 0) %>%
#     dplyr::select(date, first_boot = value)
#   
#   boot_set <- estimatedRe %>%
#     filter(variable == 'R_mean',
#            replicate != 0) %>%
#     dplyr::select(date, replicate, value) %>%
#     left_join(first_boot_set, by = 'date') %>%
#     mutate(resid = sqrt(nrow(.))*(value - first_boot) )
#   
#   
#   relevant_dates <- est_subset %>%
#     group_by(date) %>%
#     count() %>%
#     filter(n > 10) %>%
#     pull(date)
#   
#   get_pval <- function(date_id){shapiro.test(est_subset %>% filter(date == date_id) %>% pull(resid))$p.value}
#   
#   p_vals <- sapply(relevant_dates, get_pval )
#   # p-value below 0.05 indicates not normal  
#   frac_normal = sum(p_vals > 0.05)/length(p_vals)
#   dates_not_normal = relevant_dates[which(p_vals <= 0.05)]
#   
#   if (length(est_subset %>% pull(resid)) < 3){
#     overall_normal = 0
#   } else if (length(est_subset %>% pull(resid)) > 5000) {
#     resid = est_subset %>% pull(resid)
#     overall_normal = shapiro.test(resid[1:5000])$p.value > 0.05 
#   } else {
#     overall_normal = shapiro.test(est_subset %>% pull(resid))$p.value > 0.05 
#   }
#   
#   return(list(frac_normal = frac_normal, # dates_not_normal = dates_not_normal,
#               overall_normal = overall_normal))
# }
# 
# ###### Infection Incidence Estimation Error ######
# 
# getInfectionCoverage <- function(simulation, estimatedInfections){
#   
#   InfectCompare <- estimatedInfections %>%
#     dplyr::select(-c(country, region, source, data_type, variable)) %>%
#     group_by(date) %>%
#     summarize(
#       median_val = median(value),
#       high_quant = quantile(value, probs=0.975, na.rm=T),
#       low_quant = quantile(value, probs=0.025, na.rm=T),
#       .groups = "keep"
#     ) %>%
#     ungroup() %>%
#     full_join(simulation, by = c('date')) %>%
#     filter(!is.na(infections) & !is.na(median_val)) %>%
#     dplyr::select(-c(Re, observations)) %>%
#     mutate(coverage = (infections >= low_quant) & (infections <= high_quant) ,
#            CI_width = high_quant - low_quant )
#   
#   frac_coverage <- sum(InfectCompare$coverage)/nrow(InfectCompare)
#   median_CI_width <- median(InfectCompare$CI_width)
#   
#   return(list(frac_coverage = frac_coverage, CI_width = median_CI_width))
# }