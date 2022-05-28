### Parametric model: input Re and noise model, return observations (the same as the 'simBoot' function in Jana&Jeremie's code)

Para_model_func <- function(R_hat, R_dates, Obs_dates, IncubationParams, OnsetToCountParams,
                    init_infection, noise){
  
  ## the observation timeseries will be longer by a few days
  n_re <- length(R_hat)
  n_obs <- n_re + as.numeric(max(Obs_dates) - max(R_dates))     # only add at tail, not at head
  
  last_R_hat <- R_hat[length(R_hat)]
  R_hat <- c(R_hat, rep(last_R_hat, n_obs-n_re))
  
  infectionTS <- getInfectionTS(R_hat, init_infection)
  ObservationTS <- getObservationTS(infectionTS, 
                                    IncubationParams, OnsetToCountParams,
                                    noise,
                                    timevarying = FALSE)
  ObservationTS_noiseless <- getObservationTS(infectionTS, 
                                    IncubationParams, OnsetToCountParams,
                                    noise="noiseless",
                                    timevarying = FALSE)
  
  result <- data.frame(list(Re = R_hat, 
                            infections = infectionTS, 
                            observations = ObservationTS[1:n_obs],
                            Observations_noiseless = ObservationTS_noiseless[1:n_obs]))
  if (n_obs > n_re){
    result['date'] = c(R_dates, 
                       max(R_dates) + 1:(n_obs-n_re))
  } else{
    result['date'] = R_dates
  }
  
  return(result)
}