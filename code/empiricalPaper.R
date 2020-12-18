###########################################################
## empiricalPaper.R
## author: J.S. Huisman
###########################################################

library(ggplot2)
library(tidyverse)
library(lubridate)
library(ggrepel)
library(viridis)
library(dplyr)

#appDir = "../../covid-19-re-shiny-app"
dataDir = "../data"
plot_path = '../figures'

empiri_countries = c("CHE", "AUT", "BEL", "DNK",
              "FIN", "FRA", "DEU", "ITA", 
              "IRL", "ESP","NOR","POL","PRT",
               "SWE","TUR", "GBR", "NLD",
              "SVN", "ROU", "RUS")

###########################################################
### Calculate R=1 crossings for selected countries ####

estimates <- try(readRDS(paste0(dataDir, '/allCountryData.rds')))
estimates$estimates <- unique(estimates$estimates)

ReSumEstimates <- estimates$estimates %>%
  filter(estimate_type == 'Cori_slidingWindow',
         region == countryIso3,
         countryIso3 %in% empiri_countries) %>%
  dplyr::select(-estimate_type, -region, -source) %>%
  group_by(country, countryIso3, data_type) %>%
  mutate(
    sign_change = sign(median_R_mean - 1) - dplyr::lag(sign(median_R_mean - 1), 1),
    sign_change_low = sign(median_R_lowHPD - 1) - dplyr::lag(sign(median_R_lowHPD - 1), 1),
    sign_change_high = sign(median_R_highHPD - 1) - dplyr::lag(sign(median_R_highHPD - 1), 1)
    ) %>%
  summarise(first_Re_date = min(date),
            first_Re = median_R_mean[which(date == first_Re_date)],
            change_point = ifelse(median_R_mean[which.min(date)] < 1, which.min(date),
                                  which(sign_change < 0)[1] ),
            change_date = date[change_point],
            sign_change_point_low = ifelse(median_R_lowHPD[which.min(date)] < 1, which.min(date),
                                           which(sign_change_low < 0)[1] ),
            sign_change_date_low = date[sign_change_point_low],
            sign_change_point_high = ifelse(median_R_highHPD[which.min(date)] < 1, which.min(date),
                                            which(sign_change_high < 0)[1] ),
            sign_change_date_high = date[sign_change_point_high],
            .groups = 'drop') %>%
  dplyr::select(-c(change_point, sign_change_point_low, sign_change_point_high))

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

#############################################
# Intervention dates #####

# This is the same data as displayed in the 3rd panel on the shiny app
pathToInterventionDates <- file.path(dataDir, "/interventions.csv")
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

lockdown_dates <- interventionData %>%
  filter(countryIso3 %in% empiri_countries,
         name == "lockdown",
         type == 'start') %>%
  mutate(lockdown_date = lubridate::stamp(x = "01-10", orders = "%d-%Om")(date)) %>%
  dplyr::select(lockdown_date, orig_lock_date = date, country, countryIso3)

add_lock_dated = data.frame(country = "Sweden",
                            lockdown_date = "",
                            orig_lock_date = as_date(""),
                            countryIso3 = "SWE" )
lockdown_dates <- bind_rows(lockdown_dates, add_lock_dated)

###############
# Main Manuscript table 1 ########

cc_estimates <- ReSumEstimates %>%
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

#############################
# supplementary table S1 ######

final_estimates <- ReSumEstimates %>%
  group_by(country, countryIso3, data_type) %>%
  summarise('Rcross' = paste0(lubridate::stamp(x = "01-10", orders = "%d-%Om")(change_date), ' [',
                              lubridate::stamp(x = "01-10", orders = "%d-%Om")(sign_change_date_low), ', ',
                              lubridate::stamp(x = "01-10", orders = "%d-%Om")(sign_change_date_high), ']'),
            .groups = 'drop') %>%
  pivot_wider(names_from = data_type, values_from = Rcross, values_fill = list(Rcross ='')) %>%
  left_join(lockdown_dates, by = 'countryIso3') %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  select(-country.x, - country.y) %>%
  dplyr::select(country, lockdown_date, #first_measure, print_SI_date,
                `Confirmed cases`, `Deaths`, `Hospitalized patients`)

#############################################
## Do the worldwide analysis ######

estimates <- try(readRDS(paste0(dataDir, '/allCountryData.rds')))

ReSumEstimates <- estimates$estimates %>%
  ungroup() %>%
    filter(estimate_type == 'Cori_slidingWindow',
           region == countryIso3,
           data_type == 'Confirmed cases') %>%
    dplyr::select(-estimate_type, -data_type, -region, -source) %>%
    group_by(country, countryIso3) %>%  
  mutate(#Rdiff_mean = median_R_mean - dplyr::lag(median_R_mean, 1),
      sign_change = sign(median_R_mean - 1) - dplyr::lag(sign(median_R_mean - 1), 1),
      sign_change_low = sign(median_R_lowHPD - 1) - dplyr::lag(sign(median_R_lowHPD - 1), 1),
      sign_change_high = sign(median_R_highHPD - 1) - dplyr::lag(sign(median_R_highHPD - 1), 1)
    ) %>%
  #filter(!is.na(Rdiff_mean)) %>%
  summarise(first_date = min(date),
            change_point = ifelse(median_R_mean[which.min(date)] < 1, which.min(date),
                                  which(sign_change < 0)[1] ),
            Rcross_date = date[change_point],
            #Rdiff_change = (Rdiff_mean[change_point]/median_R_mean[change_point])*100,
            sign_change_point_low = ifelse(median_R_lowHPD[which.min(date)] < 1, which.min(date),
                                           which(sign_change_low < 0)[1] ),
            sign_change_date_low = date[sign_change_point_low],
            sign_change_point_high = ifelse(median_R_highHPD[which.min(date)] < 1, which.min(date),
                                            which(sign_change_high < 0)[1] ),
            sign_change_date_high = date[sign_change_point_high],
            .groups = 'drop') %>%
  dplyr::select(-c(change_point, sign_change_point_low, sign_change_point_high))

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

# if we would not filter first_date < SI50_date; 
# we would underestimate the fraction 
# where SI50_date < Rcross_date

test <- ReSumEstimates %>%
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
            #print_SI_date = lubridate::stamp(x = "01-10", orders = "%d-%Om")(SIjump_date),
            .groups = 'drop') 

test <- ReSumEstimates %>%
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

#############################################
# Supplementary figures S10 ######

plot_estimates <- ReSumEstimates %>%
  filter(data_type == 'Confirmed cases') %>%
  left_join(lockdown_dates, by = 'countryIso3') %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  dplyr::select(country, countryIso3, orig_lock_date, first_Re, change_date,
                sign_change_date_low, sign_change_date_high) %>%
  mutate(time_delay = as.numeric(change_date - orig_lock_date),
         time_delay_lowHPD = as.numeric(sign_change_date_low - orig_lock_date),
         time_delay_highHPD = as.numeric(sign_change_date_high - orig_lock_date)) %>%
  dplyr::select(-sign_change_date_low, - sign_change_date_high) %>%
  rowwise() %>%
  mutate(SI_lock = ifelse(is.na(orig_lock_date), NA,
                          as.numeric(cleanOxford[cleanOxford$countryIso3 == countryIso3 &
                                                   cleanOxford$date == orig_lock_date, 'SI'])),
         SI_Rcross = ifelse(is.na(change_date), NA,
                            as.numeric(cleanOxford[cleanOxford$countryIso3 == countryIso3 &
                                                     cleanOxford$date == change_date, 'SI'])),
         SI_average = ifelse(is.na(orig_lock_date), NA,
                             (SI_lock+SI_Rcross)/2),
         lock_Re = ifelse(is.na(orig_lock_date), NA, estimates$estimates[
           estimates$estimates$estimate_type == 'Cori_slidingWindow' & 
             estimates$estimates$region == countryIso3 &
             estimates$estimates$date == orig_lock_date &
             estimates$estimates$data_type == 'Confirmed cases', 'median_R_mean'] %>% as.numeric()),
         lock_Re_lowHPD = ifelse(is.na(orig_lock_date), NA, estimates$estimates[
           estimates$estimates$estimate_type == 'Cori_slidingWindow' & 
             estimates$estimates$region == countryIso3 &
             estimates$estimates$date == orig_lock_date &
             estimates$estimates$data_type == 'Confirmed cases', 'median_R_lowHPD'] %>% as.numeric()),
         lock_Re_highHPD = ifelse(is.na(orig_lock_date), NA, estimates$estimates[
           estimates$estimates$estimate_type == 'Cori_slidingWindow' & 
             estimates$estimates$region == countryIso3 &
             estimates$estimates$date == orig_lock_date &
             estimates$estimates$data_type == 'Confirmed cases', 'median_R_highHPD'] %>% as.numeric())
  )

bar.lm <- lm(time_delay ~  lock_Re , plot_estimates  )
summary(bar.lm)

# The very early lower bound for Spain is an artifact 
# (because I pick the 1st date the trace drops below 1)
# as there was a day below 1 at the very start of the timeseries
# but that one is not related to the actual R crossing we're interested in
# this is the hacky solution to circumvent that problem
plot_estimates[plot_estimates$countryIso3 == 'ESP', 'time_delay_lowHPD'] = 4

ggplot(plot_estimates,  
       aes(y = time_delay, x = SI_lock  )) +
  #aes(y = time_delay, x = lock_Re  )) +
  geom_smooth(method='lm', formula= y~x, colour = 'black') +
  geom_point(aes(colour = countryIso3), size = 3, show.legend = FALSE) +
  geom_errorbar(aes(ymin=time_delay_lowHPD, ymax=time_delay_highHPD, 
                    colour = countryIso3), size = 1, show.legend = FALSE) +
  #geom_errorbarh(aes(xmin=lock_Re_lowHPD, xmax=lock_Re_highHPD, 
  #            colour = countryIso3), size = 1, show.legend = FALSE) +
  geom_text_repel(aes(label= countryIso3, colour = countryIso3), 
                  show.legend = FALSE, size = 7) +
  scale_colour_viridis(discrete = T) + 
  theme_minimal() +
  labs(x = 'SI on day of lockdown', y = 'Time between lockdown and R<1') +
  #labs(x = 'Re on day of lockdown', y = 'Time between lockdown and R<1') +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    axis.text.y= element_text(size=14),
    axis.text.x= element_text(size=14),
    axis.title.y =  element_text(size=17),
    axis.title.x =  element_text(size=17),
    legend.title = element_text(size=17),
    legend.text = element_text(size=15)
  )

plotPath <- paste0(plot_path_spec, '/SI_timetoR1.pdf')
#plotPath <- paste0(plot_path_spec, '/Re_timetoR1.pdf')
ggsave(plotPath, width = 12, height = 10)


# Supplementary figure S11 ######
# the minimum R attained during lockdown was read off
# from our shiny app, which is why it's called "byHand"
byHand_df <- read_csv('../data/Rmin_SIlock.csv')

bar.lm <- lm(R_min ~ SI_lock, byHand_df)
summary(bar.lm)

byHand_df %>%
  mutate(time_of_lock = date_end - date_start,
         time_to_min = date_R_min - date_start) %>%
  dplyr::select(country, SI_lock, time_of_lock, time_to_min)


ggplot(byHand_df ,  
       aes(y = R_min, x = SI_lock  )) +
  #aes(y = R_pre_end, x = SI_lock  )) +
  geom_smooth(method='lm', formula= y~x, colour = 'black') +
  geom_point(aes(colour = countryIso3), size = 3, show.legend = FALSE) +
  geom_errorbar(aes(ymin=R_min_lowHPD, ymax=R_min_highHPD, 
                    #geom_errorbar(aes(ymin=R_pre_end_lowHPD, ymax=R_pre_end_highHPD, 
                    colour = countryIso3), size = 1,  show.legend = FALSE) +
  geom_text_repel(aes(label= countryIso3, colour = countryIso3), 
                  show.legend = FALSE, size = 7) +
  scale_colour_viridis(discrete = T) + 
  theme_minimal() +
  labs(x = 'SI during lockdown', y = 'Minimum Re during lockdown') +
  #labs(x = 'SI during lockdown', y = 'R 7-days prior to lockdown end') +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    axis.text.y= element_text(size=14),
    axis.text.x= element_text(size=14),
    axis.title.y =  element_text(size=17),
    axis.title.x =  element_text(size=17),
    legend.title = element_text(size=17),
    legend.text = element_text(size=15)
  )

plotPath <- paste0(plot_path_spec, '/Rmin_SIlock.pdf')
#plotPath <- paste0(plot_path_spec, '/R_pre_end_SIlock.pdf')
ggsave(plotPath, width = 12, height = 10)








