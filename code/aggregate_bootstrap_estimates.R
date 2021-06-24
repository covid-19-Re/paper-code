###########################################################
## aggregate_bootstrap_estimates.R
## authors: J. Scire, D. Angst and J.S. Huisman
###########################################################
# Note: this script is included for transparency
# It will not work because we did not upload the 
# individual estimate files
###########################################################

library(tidyverse)
library(lubridate)

WORKDIR <- "some_local_filepath"
ESTDIR <- paste0(WORKDIR, '/estimates')

country_codes <- c("CHE") 

all_estimates <- c()

# Shiny estimates
estimateFiles <- list.files(ESTDIR,
                         pattern = "*.rds",
                         full.names = TRUE,
                         recursive = TRUE)

all_estimates <- c()

for (file in estimateFiles) {
  estimates <- readRDS(file)
  
  date_file <- strptime(
    stringr::str_extract(file, "(\\d{8})"),
    format = "%Y%m%d")
  
  estimates <- estimates %>% 
    filter(region == "CHE", data_type == "Confirmed cases") %>% 
    dplyr::select(-country, - countryIso3, -source) %>% 
    mutate(file_date = as_date(date_file))
  # file_date: date of estimate_file; most recent possible estimate is this date -10
  
  all_estimates <- c(all_estimates, list(estimates))
}

all_estimates <- bind_rows(all_estimates)

write_csv(all_estimates, file.path(WORKDIR, "aggregated_files", "all_bootstrap_estimates.csv"))



