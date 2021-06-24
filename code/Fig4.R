###########################################################
## Fig4.R
## author: J.S. Huisman
###########################################################

library(ggplot2)
library(ggrepel)
library(qs)
library(tidyverse)
library(lubridate)
library(viridis)

plot_path = '../figures'
dataDir = "../data"

###########################################################

ReEstimates <- qread('../data/allCountryData.qs')

continentData <- read_csv(paste0(dataDir, '/continents.csv'))
newContinentData <- continentData %>%
  distinct(countryIso3, .keep_all = T)

countries = unique(ReEstimates$estimates$countryIso3)
length(countries)

###########################################################
url = 'https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv'

oxfordGRT <- read_csv(url, 
                      col_types = cols('Date' = col_character(),
                                       'RegionName' = col_character(),
                                       'RegionCode'= col_character()))
#write_csv(oxfordGRT, '../data/oxfordGRT.csv')

cleanOxford <- oxfordGRT %>%
  filter(is.na(RegionCode)) %>%
  mutate(Date = parse_date(Date, format = "%Y%m%d")) %>%
  dplyr::select(c(country = CountryName, countryIso3 = CountryCode, date = Date, 
                  GRI = GovernmentResponseIndexForDisplay, 
                  SI = StringencyIndexForDisplay, 
                  CHI = ContainmentHealthIndexForDisplay))

###########################################################
cleanRe <- ReEstimates$estimates %>%
  ungroup() %>%
  filter(region == countryIso3,
         data_type == 'Confirmed cases',
         estimate_type == 'Cori_slidingWindow') %>%
  dplyr::select(date, country, countryIso3, R = median_R_mean) %>%
  left_join(newContinentData, by = 'countryIso3') %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  dplyr::select(-country.x, - country.y)

## Individual oxford indices #####
new_colnames <- c("country", "countryIso3", "date", "C1", "C1_Flag", "C2", "C2_Flag", "C3", "C3_Flag", 
                  "C4", "C4_Flag", "C5", "C5_Flag", "C6", "C6_Flag",
                  "C7", "C7_Flag", "C8")

rowAny <- function(x) rowSums(x) > 0

# Calculate individual indices according to the codebook
# keep all weeks
# diff is daily change; then summed over week
index_df <- oxfordGRT %>%
  filter(is.na(RegionCode)) %>%
  mutate(date = parse_date(Date, format = "%Y%m%d")) %>%
  dplyr::select(country = CountryName, 
                countryIso3 = CountryCode,
                date, matches('^C._', perl = TRUE)) %>%
  rename_all(~new_colnames) %>%
  group_by(country, countryIso3) %>%
  mutate(week = isoweek(date),
         year = year(date),
         C1_index = ifelse(is.na(C1_Flag), 0, 100*(C1 - 0.5*(1-C1_Flag))/3),
         C2_index = ifelse(is.na(C2_Flag), 0, 100*(C2 - 0.5*(1-C2_Flag))/3),
         C3_index = ifelse(is.na(C3_Flag), 0, 100*(C3 - 0.5*(1-C3_Flag))/2),
         C4_index = ifelse(is.na(C4_Flag), 0, 100*(C4 - 0.5*(1-C4_Flag))/4),
         C5_index = ifelse(is.na(C5_Flag), 0, 100*(C5 - 0.5*(1-C5_Flag))/2),
         C6_index = ifelse(is.na(C6_Flag), 0, 100*(C6 - 0.5*(1-C6_Flag))/3),
         C7_index = ifelse(is.na(C7_Flag), 0, 100*(C7 - 0.5*(1-C7_Flag))/2),
         C8_index = 100*(C8/4),
  ) %>%
  mutate(across(ends_with('_index'), ~ .x - lag(.x, 1),
                .names = '{col}_diff')) %>%
  group_by(country, countryIso3, week, year) %>%
  summarise(across(ends_with('_index_diff'), ~ sum(.x, na.rm = T)),
            across(ends_with('_index'), ~ .x[which.max(date)]),
            date = max(date),
            .groups = 'drop') %>%
  filter(rowAny(across(ends_with('_index_diff'), ~ !is.na(.x))) )

index_effects <- cleanRe %>%
  mutate(Rlead7 = lead(R, 7) - R,
         Rlead10 = lead(R, 10) - R,
         Rlead14 = lead(R, 14) - R) %>%
  right_join(index_df, by = c('countryIso3', 'date')) %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  dplyr::select(-country.x, -country.y) %>%
  mutate(roundR7 = round(Rlead7, digits = 1), 
         roundR10 = round(Rlead10, digits = 1), 
         roundR14 = round(Rlead14, digits = 1), 
         across(ends_with('_index_diff'),
                ~ factor(sign(.x), levels = c(-1, 0, 1)),
                .names = 'sign_{col}'),
         dummy = 1) %>% 
  filter(!is.na(countryIso3),
         !is.na(continent))

## THE SI index #####
SI_df <- oxfordGRT %>%
  filter(is.na(RegionCode)) %>%
  mutate(date = parse_date(Date, format = "%Y%m%d"),
         week = isoweek(date),
         year = year(date)) %>%
  dplyr::select(c(country = CountryName, 
                  countryIso3 = CountryCode,
                  date, week, year, 
                  SI_index = StringencyIndexForDisplay)) %>%
  group_by(country, countryIso3) %>%
  mutate(SI_diff = SI_index - lag(SI_index, 1)) %>%
  group_by(country, countryIso3, week, year) %>%
  summarise(SI_diff = sum(SI_diff, na.rm = T),
            SI_value = SI_index[which.max(date)],
            date = max(date),
            .groups = 'drop')


SI_effects <- cleanRe %>%
  mutate(Rlead7 = lead(R, 7) - R,
         Rlead10 = lead(R, 10) - R,
         Rlead14 = lead(R, 14) - R) %>%
  right_join(SI_df, by = c('countryIso3', 'date')) %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  dplyr::select(-country.x, - country.y) %>%
  mutate(roundR7 = round(Rlead7, digits = 1),
         roundR10 = round(Rlead10, digits = 1),
         roundR14 = round(Rlead14, digits = 1),
         signCdiff = factor(sign(SI_diff), levels = c(-1, 0, 1))) %>%
  filter(!is.na(continent))

plotSI <- SI_effects %>% 
  filter(!is.na(countryIso3)) %>%
  mutate(dummy = 1,
         roundSI = round(SI_diff, digits = -1))

###########################################################
## NPI abolishment/implementation #####

SI_Rediff_plot <- function(plot_data, SI_var = 'SI_diff',
                           group_var = 'roundSI',
                           sign_var = 'signCdiff', R_var = 'Rlead7'){
  
  #plot_data <- plotSI
  
  p <- ggplot(plot_data ) +
    geom_violin(aes(x = get(R_var), y = get(group_var), group = get(group_var)),
                # outlier.shape = NA, 
                draw_quantiles = c(0.5),
                show.legend = T, orientation = 'y') +
    #geom_point(aes(x = get(R_var), y = get(SI_var)), colour = 'black',
    #         alpha = 0.2, show.legend = T) +
    scale_colour_continuous(limits = c(-1.5, 1.5)) +
    coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-50, 50)) +
    geom_vline(aes(xintercept = 0), colour = 'black', size = 0.5) +
    geom_hline(aes(yintercept = 0), colour = 'black', size = 0.5) +
    facet_grid(cols = vars(continent), scale = 'free') +
    labs(y = 'Change in SI', x = 'Change in Re over 7 days', 
         #fill = 'Change in \nStringency\nIndex') +
         colour = 'Change \nin Re') +
    theme_minimal() +
    theme(axis.text.y= element_text(size=17),
          axis.text.x= element_text(size=17),
          strip.text =  element_text(size=20), #, angle = 0
          panel.spacing.y = unit(0.75, "lines"),
          axis.title.y =  element_text(size=20),
          axis.title.x =  element_text(size=20),
          legend.title = element_text(size=20),
          legend.text = element_text(size=17)
    )
  
  plotPath <- paste0(plot_path, '/SIdiff_Rediff_distribution_country_', SI_var, '.pdf')
  ggsave(plotPath, plot = p, width = 17, height = 10)
  #ggsave(plotPath, plot = p, width = 12, height = 8)
  return(p)
}

SI_Rediff_plot(plotSI, SI_var = 'SI_diff', R_var = 'Rlead7')

###########################################################
## Direct relations ####
plot_effects <- index_effects %>%
  dplyr::select(-paste0('C', 1:8, '_index'), -paste0('sign_C', 1:8, '_index_diff'),
                #-discreteR, -signR, 
                -roundR7,-roundR10,-roundR14, -dummy
                ) %>%
  pivot_longer(cols = paste0('C', 1:8, '_index_diff'), names_to = 'index', 
               values_to = 'index_change') %>%
  mutate(index_change = factor(round(index_change)))

index_names <- setNames(paste0('C', 1:8), sort(unique(plot_effects$index)))


ggplot(plot_effects, aes(y = Rlead7, x = index_change, 
                         group = index)) +
  geom_boxplot(aes(group = index_change, fill = continent), outlier.shape = NA, show.legend = F)  +
  #geom_violin(aes(group = index_change, fill = continent), show.legend = F)  +
  #geom_point() +
  geom_smooth(aes(colour = continent), method='lm', formula = y~x, show.legend = F) +
  geom_hline(aes(yintercept = 0), colour = 'black', size = 0.5) +
  geom_vline(aes(xintercept = 0), colour = 'black', size = 0.5) +
  facet_grid(rows = vars(index), cols = vars(continent),
             labeller = labeller(index = index_names), scale = 'free_y') +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(fill = 'Continent', y = 'Change in Re over 7 days', 
       x = 'Change in Index') +
  scale_x_discrete(breaks = seq(-100, 100, 50)) +
  theme_minimal() +
  theme(axis.text.y= element_text(size=17),
        axis.text.x= element_text(size=17, angle = 45),
        strip.text.x =  element_text(size=20),
        strip.text.y =  element_text(size=20),
        strip.placement = "outside",
        panel.spacing.x = unit(1.5, "lines"),
        panel.spacing.y = unit(1, "lines"),
        axis.title.y =  element_text(size=20),
        axis.title.x =  element_text(size=20)
  )

plotPath <- paste0(plot_path, '/All_index_lms.pdf')
ggsave(plotPath, width = 14, height = 18)


