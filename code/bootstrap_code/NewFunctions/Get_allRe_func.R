### A function returns 4 quantities related to the poeterior of Re:
### R_mean, R_highHPD, R_lowHPD, R_sample

### combine functions 'estimateInfectionTS' and 'estimateReTS' from Jana and Jeremie's code
### require input 'data_obs' contains c("date", "observations")

Get_allRe_func <- function(data_obs, IncubationParams, OnsetToCountParams,
                           smooth_para_deConv,
                           timevarying = FALSE, 
                           delay = 0){
                            
  
  # 1. estimatedInfections
  infection_df <- getSimIncidence(data_obs)
  
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
  
  estimatedInfections <- get_infection_incidence_by_deconvolution(
    infection_df,
    is_local_cases = T,
    constant_delay_distribution = constant_delay_distributions[['Simulated']],
    constant_delay_distribution_incubation = constant_delay_distributions[["Symptoms"]],
    max_iterations = 100,
    smooth_incidence = T,
    smooth_para_deConv = smooth_para_deConv,
    empirical_delays = delays_onset_to_count,
    n_bootstrap = 0,
    verbose = FALSE)
  
  # 2. estimateReTS
  
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