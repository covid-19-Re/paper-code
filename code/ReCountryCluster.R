
library("lubridate")
library("fitdistrplus")
library("EpiEstim")
library("readxl")
library("tidyverse")


args <- commandArgs(trailingOnly = TRUE)
incidenceDate <- args[1]

popData <- qs::qread(file = "data/popData.qs")
interval_ends <- qs::qread(file = "data/interval_ends.qs")
countryData <- qs::qread(
  file = str_c("data/countryData/", incidenceDate, "_countryData.qs" ))

# get Infection Incidence
# load functions
source("scripts/2_utils_getInfectionIncidence.R")
# load parameter
source("scripts/2_params_InfectionIncidencePars.R")
# load empirical delays
delays_data_path <- str_c("data/delays/", incidenceDate, "_CHE_data_delays.csv")
delays_onset_to_count <- read_csv(delays_data_path,
                                  col_types = cols(
                                    data_type = col_character(),
                                    onset_date = col_date(format = ""),
                                    count_date = col_date(format = ""),
                                    delay = col_number()))
# constant delay distribution
constant_delay_distributions <- list()
for (type_i in unique(names(shape_onset_to_count))) {
  m <- get_vector_constant_waiting_time_distr(
    shape_incubation,
    scale_incubation,
    shape_onset_to_count[[type_i]],
    scale_onset_to_count[[type_i]])

  constant_delay_distributions <- c(constant_delay_distributions, list(m))
}
names(constant_delay_distributions) <- unique(names(shape_onset_to_count))

constant_delay_symptom_to_report_distributions <- list()
for (type_i in unique(names(shape_onset_to_count))) {
  m <- get_vector_constant_waiting_time_distr(
    0,
    0,
    shape_onset_to_count[[type_i]],
    scale_onset_to_count[[type_i]])

  constant_delay_symptom_to_report_distributions <- c(constant_delay_symptom_to_report_distributions, list(m))
}
names(constant_delay_symptom_to_report_distributions) <- paste0('Onset to ',  unique(names(shape_onset_to_count)))

constant_delay_distributions <- c(constant_delay_distributions, constant_delay_symptom_to_report_distributions)

# filter out regions with too few cases for estimation
countryData <- countryData %>%
  filterRegions(thresholdConfirmedCases = 500)
# remove Oxford Stringenxy Index for Re calculation
countryData <- countryData %>%
  filter(data_type != "Stringency Index")

# filter out data_types with 0 total cases
data_type0 <- countryData %>%
  group_by(data_type) %>%
  summarize(total = sum(value), .groups = "drop") %>%
  filter(total == 0) %>%
  .$data_type

countryData <- filter(countryData, !(data_type %in% data_type0))

if (nrow(countryData) > 0) {

  countryData <- countryData %>%
    mutate(
      data_type = fct_drop(data_type)
    )

  right_truncation <- list()

    right_truncation[["Confirmed cases"]] <- 0
    right_truncation[["Confirmed cases / tests"]] <- 0
    right_truncation[["Hospitalized patients"]] <- 0
    right_truncation[["Deaths"]] <- 0
 

  right_truncate <- function(df, data_type, right_truncation) {
      dplyr::filter(df, date <= (max(date) - right_truncation[[unique(data_type)]]))
  }

  countryData <- countryData %>%
    group_by(country, region, source, data_type) %>%
    right_truncate(data_type, right_truncation) %>%
    dplyr::select(-countryIso3, -populationSize) %>%
    ungroup()
  
  # Deconvolution
  deconvolvedData <- list()

  deconvolvedData[[1]] <- get_all_infection_incidence(
    countryData,
    constant_delay_distributions = constant_delay_distributions,
    onset_to_count_empirical_delays = delays_onset_to_count,
    data_types = c("Confirmed cases",
                    "Hospitalized patients",
                    "Deaths"),
    n_bootstrap = 100,
    verbose = FALSE)


    countryDataTests <- countryData %>%
      filter(region == "CHE", data_type == "Confirmed cases / tests")

    deconvolvedData[[2]] <- get_all_infection_incidence(
      countryDataTests,
      constant_delay_distributions = constant_delay_distributions,
      onset_to_count_empirical_delays = delays_onset_to_count,
      data_types = c("Confirmed cases / tests"),
      n_bootstrap = 100,
      verbose = FALSE)

  
  deconvolvedCountryData <- bind_rows(deconvolvedData)
  countryDataPath <- str_c("data/deconvoluted/", incidenceDate,"-CHE", "-DeconvolutedData.rds")
  if (dim(deconvolvedCountryData)[1] == 0) {
    print("no data remaining")
    next
  }
    saveRDS(deconvolvedCountryData, file = countryDataPath)
    # Re Estimation
    source("scripts/3_utils_doReEstimates.R")
    
    swissRegions <- deconvolvedCountryData %>%
      filter(country %in% c("Switzerland", "Liechtenstein")) %>%
      dplyr::select(region) %>%
      distinct() %>%
      .$region
    
    ### Window
    window <- 3
    
    ##TODO this all_delays could be removed because we always deconvolve
    ### Delays applied
    all_delays <- list(
      "infection_Confirmed cases" = c(Cori = 0, WallingaTeunis = -5),
      "infection_Confirmed cases / tests" = c(Cori = 0, WallingaTeunis = -5),
      "infection_Deaths" = c(Cori = 0, WallingaTeunis = -5),
      "infection_Hospitalized patients" = c(Cori = 0, WallingaTeunis = -5),
      "Confirmed cases" = c(Cori = 10, WallingaTeunis = 5),
      "Confirmed cases / tests" = c(Cori = 10, WallingaTeunis = 5),
      "Deaths" = c(Cori = 20, WallingaTeunis = 15),
      "Hospitalized patients" = c(Cori = 8, WallingaTeunis = 3),
      "infection_Excess deaths" = c(Cori = 0, WallingaTeunis = -5),
      "Excess deaths" = c(Cori = 20, WallingaTeunis = 15))
    
    truncations <- list(
      left = c(Cori = 5, WallingaTeunis = 0),
      right = c(Cori = 0, WallingaTeunis = 8))
    
    ### Run EpiEstim
    countryEstimatesRaw <- doAllReEstimations(
      deconvolvedCountryData,
      slidingWindow = window,
      methods = "Cori",
      #variationTypes = c("step", "slidingWindow"),
      variationTypes = c("slidingWindow"),
      all_delays = all_delays,
      truncations = truncations,
      interval_ends = interval_ends,
      swissRegions = swissRegions)

    countryEstimates <- cleanCountryReEstimate(countryEstimatesRaw, method = 'bootstrap') 
    
  countryEstimates <- countryEstimates %>%
      left_join(
        dplyr::select(popData, region, countryIso3),
        by = c("region")
      )
    countryDataPath <- str_c("data/estimates/", incidenceDate, "CHE", "-Estimates.rds")
    saveRDS(countryEstimates, file = countryDataPath)
} else {
  cat(str_c(Sys.time(), " | ", "CHE", ": Not enough cases. Skipping Re calculation.\n"))
}
