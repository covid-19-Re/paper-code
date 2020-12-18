###########################################################
## compareSimulations.R
## author: J.S. Huisman
###########################################################

# this library is optional, as the DTW distance
# wasn't used for the paper
#library(dtwclust)

###### Compare Simulation and Estimation results ######

getReRMSE <- function(ReCompare){
  ReCompare <- ReCompare %>%
    mutate(SE = ( median_R_mean - Re)^2)
  
  rmse <- sqrt(sum(ReCompare$SE)/length(ReCompare$SE))
  norm_rmse <- rmse / mean(ReCompare$Re)
  return(norm_rmse)
}

getRootDiff <- function(simulation, estimatedRe, all = FALSE){
  
  sim_change_points = diff(sign(simulation$Re -1))  
  sim_roots = simulation$date[which(sim_change_points != 0 )]
  
  if (! ('median_R_mean' %in% colnames(estimatedRe))){
    estimatedRe <- cleanReTSestimate(estimatedRe)
  }
  
  est_change_points = diff(sign(estimatedRe$median_R_mean -1))  
  est_roots = estimatedRe$date[which(est_change_points != 0 )]
  
  root_diffs <- as.numeric(est_roots - sim_roots)
  
  if (all){
    root_diff = sum(abs(root_diffs))
  } else{
    root_diff = root_diffs[1]
  }
  return( root_diff )
}

getEmpCoverage <- function(ReCompare){
  covCompare <- ReCompare %>%
    mutate(coverage = (Re >= median_R_lowHPD) & (Re <= median_R_highHPD) ,
           CI_width = median_R_highHPD - median_R_lowHPD)
  
  frac_coverage <- sum(covCompare$coverage)/nrow(covCompare)
  median_CI_width <- median(covCompare$CI_width)
  
  return(list(frac_coverage = frac_coverage, CI_width = median_CI_width))
}

getOneTransitions <- function(simulation, estimatedRe){
  simulation <- simulation %>%
    mutate(OneTrans = ifelse(Re >= 1, 1, 0))
  
  if (! ('median_R_mean' %in% colnames(estimatedRe))){
    estimatedRe <- cleanReTSestimate(estimatedRe)
  }
  
  one_transitions = estimatedRe %>%
    mutate(EstOneTrans = case_when(median_R_lowHPD >= 1 ~ 1,
                                   median_R_highHPD < 1 ~ 0)) %>%
    full_join(simulation, by = c('date')) %>%
    dplyr::select(date, Re, OneTrans, EstOneTrans) %>%
    filter(!is.na(EstOneTrans)) %>%
    mutate(compare = OneTrans == EstOneTrans)
  
  perc_transitions = sum(one_transitions$compare)/nrow(one_transitions)
  
  return(perc_transitions)
}

getReError <- function(simulation, estimatedRe){
  
  if (! ('median_R_mean' %in% colnames(estimatedRe))){
    estimatedRe <- cleanReTSestimate(estimatedRe)
  }
  
  #Compare
  ReCompare <- simulation %>%
    full_join(estimatedRe, by = c('date')) %>%
    filter(!is.na(median_R_mean))
  
  ReRMSE <- getReRMSE(ReCompare)
  RootDiff <- getRootDiff(simulation, estimatedRe, all = FALSE)
  EmpCoverage <- getEmpCoverage(ReCompare)
  OneTrans <- getOneTransitions(simulation, estimatedRe)
  
  return(list(ReRMSE = ReRMSE, RootDiff = RootDiff, 
                EmpCoverage = EmpCoverage$frac_coverage,
                OneTrans = OneTrans
              ))
  
}

getSlopeError <- function(simulation, estimatedRe, valid_cond){
  date_t3 = simulation$date[valid_cond$t3]
  date_t2 = simulation$date[valid_cond$t2]
  date_first = max(min(estimatedRe$date), date_t2)
  date_int = as.numeric(date_t3 - date_first)
  
  slope_sim = (simulation$Re[simulation$date == date_t3] -
                 simulation$Re[simulation$date == date_first])/date_int
  
  meanEstRe <- estimatedRe %>%
    filter(variable == 'R_mean') %>%
    dplyr::select(date, replicate, value) %>%
    group_by(replicate) %>%
    summarise(slope_est = (value[date == date_t3] -
                             value[date == date_first])/date_int,
              .groups = "drop")
  
  abs_slope_error = slope_sim - meanEstRe$slope_est
  rel_slope_error = (abs_slope_error/slope_sim)
  return(rel_slope_error)
}

getMeanSlopeError <- function(valid_cond_grid){
  SlopeError = data.frame()
  for (row_id in 1:nrow(valid_cond_grid)){
    
    simulation <- read_csv(paste0(valid_cond_grid[row_id, 'simulationDir'], '/', 
                                  valid_cond_grid[row_id, 'filename']))
    
    estimatedInfections <- read_csv(paste0(valid_cond_grid[row_id, 'estimationDir'],
                                           valid_cond_grid[row_id, 'infection_file']))
    
    estimatedRe <- read_csv(paste0(valid_cond_grid[row_id, 'estimationDir'], 
                                   valid_cond_grid[row_id, 're_file']))
    
    new_error <- getSlopeError(simulation, estimatedRe, valid_cond_grid[row_id, ])
    SlopeError <- bind_rows(SlopeError, list(slopeError = mean(new_error)) )
  }
  return(SlopeError)
}
