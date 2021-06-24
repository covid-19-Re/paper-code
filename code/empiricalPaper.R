###########################################################
## empiricalPaper.R
## author: J.S. Huisman
###########################################################

library(ggplot2)
library(tidyverse)
library(qs)
library(lubridate)
library(ggrepel)
library(viridis)
library(dplyr)

plot_path = '../figures'
dataDir = "../data"

empiri_countries = c("CHE", "AUT", "BEL", "DNK",
              "FIN", "FRA", "DEU", "ITA", 
              "IRL", "ESP","NOR","POL","PRT",
               "SWE","TUR", "GBR", "NLD",
              "SVN", "ROU", "RUS")

###########################################################
### Calculate R=1 crossings ####

estimates <- qread('../data/allCountryData.qs')
estimates$estimates <- unique(estimates$estimates)

ReSumEstimates <- estimates$estimates %>%
  filter(estimate_type == 'Cori_slidingWindow',
         region == countryIso3,
         date >= as_date("2020-02-01")) %>%
  dplyr::select(-estimate_type, -region, -source) %>%
  group_by(country, countryIso3, data_type) %>%
  mutate(sign_change = sign(median_R_mean - 1) - dplyr::lag(sign(median_R_mean - 1), 1),
    sign_change_low = sign(median_R_lowHPD - 1) - dplyr::lag(sign(median_R_lowHPD - 1), 1),
    sign_change_high = sign(median_R_highHPD - 1) - dplyr::lag(sign(median_R_highHPD - 1), 1)
    ) %>%
  summarise(first_Re_date = min(date),
            first_Re = median_R_mean[which(date == first_Re_date)],
            change_point = ifelse(median_R_mean[which.min(date)] < 1, which.min(date),
                                  which(sign_change < 0)[1] ),
            change_date = date[change_point],
            sign_change_point_low = which(sign_change_low < 0)[which.min(abs(change_point - which(sign_change_low < 0)))],
            sign_change_date_low = date[sign_change_point_low],
            sign_change_point_high = which(sign_change_high < 0)[which.min(abs(which(sign_change_high < 0) - change_point))],
            sign_change_date_high = date[sign_change_point_high],
            .groups = 'drop') %>%
  dplyr::select(-c(change_point, sign_change_point_low, sign_change_point_high))

WWReSumEstimates <- ReSumEstimates %>%
  filter(data_type == 'Confirmed cases') %>%
  dplyr::select(-data_type) %>%
  mutate(first_date = first_Re_date,
         Rcross_date = change_date)

EUReSumEstimates <- ReSumEstimates %>%
  filter(countryIso3 %in% empiri_countries)

#############################################
# Intervention dates #####

pathToInterventionDates <- file.path("../data/interventions.csv")
interventionData <- read_csv(pathToInterventionDates,
                             col_types = cols(
                               date = col_date(format = ""),
                               country = col_character(),
                               measure = col_character(),
                               source = col_character(),
                               name = col_character(),
                               y = col_double(),
                               text = col_character(),
                               tooltip = col_character(),
                               type = col_character(),
                               plotTextPosition = col_character(),
                               region = col_character()))

first_measures <- interventionData %>%
  filter(countryIso3 %in% empiri_countries,
         date > '2020-02-15') %>%
  group_by(country, countryIso3) %>%
  summarise(first_measure = min(date),
            name = paste0(list(name[which(date == first_measure)])),
            .groups = 'drop') %>%
  mutate(first_measure = lubridate::stamp(x = "01-10", orders = "%d-%Om")(first_measure)) %>%
  dplyr::select(-name)

additional_countries = data.frame(country = c("Bosnia_and_Herzegovina",
                                              "Serbia"),
                            first_measure = c("11-03",
                                              "15-03"),
                            countryIso3 = c("BIH",
                                            "SRB"))
first_measures <- bind_rows(first_measures, additional_countries)

lockdown_dates <- interventionData %>%
  filter(countryIso3 %in% empiri_countries,
         name == "lockdown",
         type == 'start') %>%
  mutate(lockdown_date = lubridate::stamp(x = "01-10", orders = "%d-%Om")(date)) %>%
  dplyr::select(lockdown_date, orig_lock_date = date, country, countryIso3)

add_lock_dated = data.frame(country = c("Bosnia and Herzegovina",
                                        "Czechia",
                                        "Serbia",
                                        "Sweden"),
                            lockdown_date = c("18-03",
                                              "12-03",
                                              "15-03",
                                              ""
                                              ),
                            orig_lock_date = as_date(c("2020-03-18",
                                                       "2020-03-12",
                                               "2020-03-15",
                                               ""
                                               )),
                            countryIso3 = c("BIH",
                                        "CZE",
                                        "SRB",
                                        "SWE") )
lockdown_dates <- bind_rows(lockdown_dates, add_lock_dated)

#############################################
# Oxford CGRT data #######

url = 'https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv'

oxfordGRT <- read_csv(url, 
                      col_types = cols('Date' = col_character(),
                                       'RegionName' = col_character(),
                                       'RegionCode'= col_character()))

cleanOxford <- oxfordGRT %>%
  filter(CountryCode %in% empiri_countries,
         is.na(RegionCode)) %>%
  mutate(Date = parse_date(Date, format = "%Y%m%d")) %>%
  dplyr::select(c(country = CountryName, countryIso3 = CountryCode, date = Date, 
                  SI = StringencyIndexForDisplay))

lockdownDF <- cleanOxford %>%
  group_by(country, countryIso3) %>%
  summarise(SI50_date = date[min(which(SI > 50))],
            print_SI_date = lubridate::stamp(x = "01-10", orders = "%d-%Om")(SI50_date),
            .groups = 'drop') 
# this corresponds to known date at least for CH and DE

###############
# Main Manuscript table ########

cc_estimates <- EUReSumEstimates %>%
  group_by(country, countryIso3, data_type) %>%
  filter(data_type == 'Confirmed cases') %>%
  summarise('Rcross' = paste0(lubridate::stamp(x = "01-10", orders = "%d-%Om")(change_date), ' [',
                              lubridate::stamp(x = "01-10", orders = "%d-%Om")(sign_change_date_low), ', ',
                              lubridate::stamp(x = "01-10", orders = "%d-%Om")(sign_change_date_high), ']'),
            first_Re_date = min(first_Re_date),
            change_date,
            .groups = 'drop')

cc_estimates

final_estimates <- cc_estimates %>%
  pivot_wider(names_from = data_type, values_from = Rcross, values_fill = list(Rcross ='')) %>%
  left_join(lockdown_dates, by = 'countryIso3') %>%
  left_join(first_measures, by = 'countryIso3') %>%
  mutate(include = orig_lock_date - first_Re_date,
         time_delay = change_date - orig_lock_date) %>%
  dplyr::select(country, countryIso3, first_measure, lockdown_date, `Confirmed cases`, time_delay, include)#, `Deaths`)

final_estimates


# library('xtable')
# print(xtable(final_estimates, caption = "R = 1 crossings", 
#              type = 'latex'), include.rownames = FALSE, 
#       file = '../figures/R_crossing_cc.tex')

#############################
# supplementary table ######

final_estimates <- EUReSumEstimates %>%
  group_by(country, countryIso3, data_type) %>%
  summarise('Rcross' = paste0(lubridate::stamp(x = "01-10", orders = "%d-%Om")(change_date), ' [',
                              lubridate::stamp(x = "01-10", orders = "%d-%Om")(sign_change_date_low), ', ',
                              lubridate::stamp(x = "01-10", orders = "%d-%Om")(sign_change_date_high), ']'),
            .groups = 'drop') %>%
  pivot_wider(names_from = data_type, values_from = Rcross, values_fill = list(Rcross ='')) %>%
  left_join(lockdown_dates, by = 'countryIso3') %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  dplyr::select(-country.x, - country.y) %>%
  dplyr::select(country, lockdown_date, #first_measure, print_SI_date,
                `Confirmed cases`, `Deaths`, `Hospitalized patients`)

#############################################
## Worldwide analysis ######

##### SI50 analysis #####

lockdownDF <- oxfordGRT %>%
  filter(is.na(RegionCode)) %>%
  mutate(Date = parse_date(Date, format = "%Y%m%d")) %>%
  dplyr::select(c(country = CountryName, countryIso3 = CountryCode, date = Date, 
                  SI = StringencyIndexForDisplay)) %>%
  group_by(country, countryIso3) %>%
  mutate(SIdiff = SI - lag(SI, 7)) %>% 
  filter(!is.na(SIdiff)) %>%
  summarise(SI50_date = date[min(which(SI > 50))],
            .groups = 'drop') 

# if we don't filter first_date < SI50_date; of course we underestimate the fraction 
# where SI50_date < Rcross_date

#Rcross_date is na for India??
test <- WWReSumEstimates %>%
  left_join(lockdownDF, by = 'countryIso3') %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  filter(first_date < SI50_date,
         !is.na(Rcross_date)) %>%
  mutate(lockdown_relevant = Rcross_date >= SI50_date,
         lockdown_signif = sign_change_date_low >= SI50_date,
         qual_check = sign_change_date_high >= sign_change_date_low,
         lock_diff = as.numeric(Rcross_date - SI50_date),
         lock_diff_signif = as.numeric(sign_change_date_low - SI50_date) ) %>%
  dplyr::select(country, countryIso3, lockdown_relevant, 
                lockdown_signif, lock_diff, lock_diff_signif, qual_check) 

nrow(test)
paste(unique(test$country), collapse = ', ')

sum(test$lockdown_relevant, na.rm = T)
sum(!test$lockdown_relevant, na.rm = T) # Andorra

test%>%
  filter(!test$lockdown_relevant)


sum(test$lockdown_signif, na.rm = T)
sum(!test$lockdown_signif, na.rm = T) # Australia, Egypt, Singapore 

test %>%
  filter(!test$lockdown_signif)

test %>%
  filter(lock_diff <= 3)

test %>%
  filter(lock_diff_signif <= 3)


#### SI jump analysis #####


lockdownDF <- oxfordGRT %>%
  filter(is.na(RegionCode)) %>%
  mutate(Date = parse_date(Date, format = "%Y%m%d")) %>%
  dplyr::select(c(country = CountryName, countryIso3 = CountryCode, date = Date, 
                  SI = StringencyIndexForDisplay)) %>%
  group_by(country, countryIso3) %>%
  mutate(SIdiff = SI - lag(SI, 7)) %>% 
  filter(!is.na(SIdiff),
         date < '2020-06-15') %>%
  summarise(SIjump_date = date[which.max(SIdiff)],
            .groups = 'drop') 

test <- WWReSumEstimates %>%
  left_join(lockdownDF, by = 'countryIso3') %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  filter(first_date < SIjump_date,
         !is.na(Rcross_date)) %>%
  mutate(lockdown_relevant = Rcross_date >= SIjump_date,
         lockdown_signif = sign_change_date_low >= SIjump_date,
         lock_diff = as.numeric(Rcross_date - SIjump_date),
         lock_diff_signif = as.numeric(sign_change_date_low - SIjump_date)) %>%
  dplyr::select(country, countryIso3, lockdown_relevant, 
                lockdown_signif, lock_diff, lock_diff_signif)

nrow(test)
paste(unique(test$country), collapse = ', ')

sum(test$lockdown_relevant, na.rm = T)
sum(!test$lockdown_relevant, na.rm = T) # Andorra, Denmark

sum(test$lockdown_signif, na.rm = T)
sum(!test$lockdown_signif, na.rm = T)

test %>%
  filter(!test$lockdown_relevant)

test %>%
  filter(!test$lockdown_signif)

test %>%
  filter(lock_diff <= 3)

test %>%
  filter(lock_diff_signif <= 3)

test %>%
  filter(!test$lockdown_signif)

