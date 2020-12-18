###########################################################
## interventionImpact.R
## author: J.S. Huisman
###########################################################

library(ggplot2)
library(ggrepel)

library(tidyverse)
library(lubridate)
library(viridis)

#library(dplyr)

appDir = "../../covid-19-re-shiny-app"
dataDir = "../data"
plot_path = '../figures'

###########################################################

ReEstimates <- readRDS(paste0(dataDir, '/allCountryData.rds'))

#####
continentData <- read_csv(paste0(dataDir, '/continents.csv'))
# some countries are in 2 continents, here we keep only one
newContinentData <- continentData %>%
  distinct(countryIso3, .keep_all = T)

countries = unique(ReEstimates$estimates$countryIso3)
length(countries)

empiri_countries = c("CHE", "AUT", "BEL", "DNK",
              "FIN", "FRA", "DEU", "ITA", 
              "IRL", "ESP","NOR","POL","PRT",
              "SWE","TUR", "GBR", "NLD",
              "SVN", "ROU", "RUS")

##### Used for Fig 4 #####
cleanRe <- ReEstimates$estimates %>%
  ungroup() %>%
  filter(region == countryIso3,
         data_type == 'Confirmed cases',
         estimate_type == 'Cori_slidingWindow') %>%
  dplyr::select(date, country, countryIso3, R = median_R_mean) %>%
  left_join(newContinentData, by = 'countryIso3') %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  dplyr::select(-country.x, - country.y)

###########################################################
url = 'https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv'

oxfordGRT <- read_csv(url, 
                      col_types = cols('Date' = col_character(),
                                       'RegionName' = col_character(),
                                       'RegionCode'= col_character()))

# RegionCode is NA if the numbers refer to the whole country
cleanOxford <- oxfordGRT %>%
  filter(is.na(RegionCode)) %>%
  mutate(Date = parse_date(Date, format = "%Y%m%d")) %>%
  dplyr::select(c(country = CountryName, countryIso3 = CountryCode, date = Date, 
                  GRI = GovernmentResponseIndexForDisplay, 
                  SI = StringencyIndexForDisplay, 
                  CHI = ContainmentHealthIndexForDisplay))

## Individual oxford indices #####
new_colnames <- c("country", "countryIso3", "date", "C1", "C1_Flag", "C2", "C2_Flag", "C3", "C3_Flag", 
                  "C4", "C4_Flag", "C5", "C5_Flag", "C6", "C6_Flag",
                  "C7", "C7_Flag", "C8")
#https://github.com/OxCGRT/covid-policy-tracker/blob/master/documentation/index_methodology.md

rowAny <- function(x) rowSums(x) > 0

ox_index_df <- oxfordGRT %>%
  filter(is.na(RegionCode)) %>%
  mutate(date = parse_date(Date, format = "%Y%m%d")) %>%
  dplyr::select(country = CountryName, 
                countryIso3 = CountryCode,
                date, matches('^C._', perl = TRUE)) %>%
  rename_all(~new_colnames) %>%
  group_by(country, countryIso3) %>%
  mutate(week = isoweek(date),
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
  group_by(country, countryIso3, week) %>%
  summarise(across(ends_with('_index_diff'), ~ sum(.x, na.rm = T)),
            across(ends_with('_index'), ~ .x[which.max(date)]),
            date = max(date),
            .groups = 'drop') %>%
  filter( rowAny(across(ends_with('_index_diff'), ~ .x != 0)) ) 



###########################################################
## SI50 Regression - Fig. 3b ####
lockdownDF <- cleanOxford %>%
  group_by(country, countryIso3) %>%
  mutate(diff = SI - lag(SI, 7)) %>% 
  filter(!is.na(diff)) %>%
  summarise(cross50_date = date[min(which(SI > 50))],
            SI50_diff = diff[which(date == cross50_date)],
            .groups = 'drop') 

slopesDF <- ReEstimates$estimates %>%
  ungroup() %>%
  filter(region == countryIso3,
         data_type == 'Confirmed cases',
         estimate_type == 'Cori_slidingWindow') %>%
  dplyr::select(date, country, countryIso3, R = median_R_mean) %>%
  left_join(lockdownDF, by = 'countryIso3') %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  group_by(country, countryIso3, cross50_date, SI50_diff) %>%
  mutate(Rdiff = R - lag(R, 1)) %>% 
  filter(!is.na(Rdiff),
         !is.na(cross50_date)) %>%
  summarise(first_R_date = min(date),
            R_SI50 = R[date == cross50_date],
            Rdiff_SI50 = Rdiff[date == cross50_date],
           .groups = 'drop')

slopesDF <- slopesDF %>%
  left_join(newContinentData, by = 'countryIso3') %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) )
## the country names do not entirely match between Re estimates and 
# Oxford stringency index

nrow(slopesDF)
paste(unique(slopesDF$country), collapse = ', ')

for (option in c('all', 'eu_select')){
  print('##############################')
  print(option)
  
  if (option == 'all'){
    finalDF <- slopesDF
    colourvar = 'continent'
  } else if (option == 'eu_select'){
    finalDF <- slopesDF %>%
      filter(continent == 'Europe')
    colourvar = 'countryIso3'
  } 
  
  print(nrow(finalDF))
  print(paste(unique(finalDF$country), collapse = ', '))
  
  bar.lm <- lm(Rdiff_SI50 ~ SI50_diff, finalDF)
  print(summary(bar.lm))

  # this if-else loop is just so that the EU-only plots
  # do not have coloured points
  
  if (option == 'all'){
    ggplot(finalDF, aes(y = Rdiff_SI50, x = SI50_diff)) +
      geom_smooth(method='lm', formula= y~x, colour = 'black') + 
      geom_point(aes(colour = get(colourvar)), show.legend = FALSE) +
      geom_text_repel(aes(label= countryIso3, colour = get(colourvar)), 
                      show.legend = FALSE, size = 7) +
      scale_colour_viridis(direction = -1, discrete = T) + 
      theme_minimal() + 
      labs(x = '7-day SI change prior to SI>50', y = '1-day Re change')+
      theme_minimal() +
      theme(
        axis.text.y= element_text(size=20),
        axis.text.x= element_text(size=20),
        axis.title.y =  element_text(size=25),
        axis.title.x =  element_text(size=25)
      )
  
  } else if (option == 'eu_select'){
    ggplot(finalDF, aes(y = Rdiff_SI50, x = SI50_diff)) +
      geom_smooth(method='lm', formula= y~x, colour = 'black') + 
      geom_point(show.legend = FALSE) +
      geom_text_repel(aes(label= countryIso3), 
                      show.legend = FALSE, size = 7) +
      theme_minimal() + 
      labs(x = '7-day SI change prior to SI>50', y = '1-day Re change')+
      theme_minimal() +
      theme(
        axis.text.y= element_text(size=20),
        axis.text.x= element_text(size=20),
        axis.title.y =  element_text(size=25),
        axis.title.x =  element_text(size=25)
      )
    
  }
  
  plotPath <- paste0(plot_path, '/SI50_regression_', option, '.pdf')
  ggsave(plotPath, width = 12, height = 9)
  
}

###########################################################
## NPI abolishment/implementation - Fig 4, S12; Table S3, S4 #####

index_effects <- cleanRe %>%
  mutate(Rlead7 = lead(R, 7) - R) %>%
  left_join(ox_index_df, by = c('countryIso3', 'date')) %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  dplyr::select(-country.x, -country.y) %>%
  filter(rowAny(across(ends_with('_index_diff'), ~ !is.na(.x))) ) %>%
  mutate(discreteR = cut(Rlead7, 
                         breaks=seq(-2, 2, 0.1)),
         signR = factor(sign(Rlead7), levels = c(-1, 0, 1)),
         roundR = round(Rlead7, digits = 1), 
         across(ends_with('_index_diff'),
                ~ factor(sign(.x), levels = c(-1, 0, 1)),
                .names = 'sign_{col}'),
         dummy = 1) %>% 
  filter(!is.na(discreteR), 
         !is.na(countryIso3))

## THE SI index #####
SI_df <- oxfordGRT %>%
  filter(is.na(RegionCode)) %>%
  mutate(date = parse_date(Date, format = "%Y%m%d"),
         week = isoweek(date)) %>%
  dplyr::select(c(country = CountryName, 
                  countryIso3 = CountryCode,
                  date, week,
                  SI_index = StringencyIndexForDisplay)) %>%
  group_by(country, countryIso3) %>%
  mutate(SI_diff = SI_index - lag(SI_index, 1)) %>%
  group_by(country, countryIso3, week) %>%
  summarise(SI_diff = sum(SI_diff, na.rm = T),
            SI_value = SI_index[which.max(date)],
            date = max(date),
            .groups = 'drop') %>%
  filter(SI_diff != 0)


SI_effects <- cleanRe %>%
  mutate(Rlead7 = lead(R, 7) - R,
         normRlead7 = Rlead7/R) %>%
  left_join(SI_df, by = c('countryIso3', 'date')) %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  dplyr::select(-country.x, - country.y) %>%
  filter(!is.na(SI_diff)) %>%
  mutate(discreteR = cut(Rlead7, breaks=seq(-2, 2, 0.1)),
         signR = factor(sign(Rlead7), levels = c(-1, 0, 1)),
         roundR = round(Rlead7, digits = 1),
         signCdiff = factor(sign(SI_diff), levels = c(-1, 0, 1)))

plotSI <- SI_effects %>% 
  filter(!is.na(discreteR), 
         !is.na(countryIso3)) %>%
  mutate(dummy = 1,
         roundSI = round(SI_diff))

## Compare distributions - Table S3 #######

test_SI_RE_dist <- function(data, method = 'mean',
                            sign_var = 'signCdiff'){
  p_vec <- vector(mode = 'numeric',
                  length = length(unique(data$continent)))
  count = 1
  for(cont in unique(data$continent)){
    # after the continent implemented measures
    cont_implement <- data %>% 
      filter(get(sign_var) == 1,
             continent == cont) %>%
      pull(Rlead7)
    
    # after the measures were abolished
    cont_decrease <- data %>% 
      filter(get(sign_var) == -1,
             continent == cont) %>%
      pull(Rlead7)
    
    if (method == 'mean'){
      # for diff in means:
      obs =  mean(cont_decrease) - mean(cont_implement)
    } else {
      # to check whether neg mean further from 0 than pos
      obs = abs(mean(cont_implement)) - abs(mean(cont_decrease))
    }

    n_inc = length(cont_implement)
    cont_all = c(cont_implement, cont_decrease)

    tstat = vector(mode = 'numeric', length = 10000)
    for(i in 1:10000){
      new_inc = sample(cont_all, n_inc)
      new_dec = setdiff(cont_all, new_inc)

      if (method == 'mean'){
        tstat[i] = mean(new_inc) - mean(new_dec)
      } else {
        tstat[i] = abs(mean(new_inc)) - abs(mean(new_dec))
      }

    }

    p_val = sum(tstat > obs)/length(tstat)
    
    p_vec[count] = p_val
    count = count+1
    
    if (p_val < (0.05/6)){
      print('##################')
      print(cont)
      print(length(cont_implement))
      print(length(cont_decrease))
      print(p_val)
    }
  }
  return(p_vec)
}

test_SI_RE_dist(SI_effects, method = 'mean')
test_SI_RE_dist(SI_effects, method = 'abs')

# the same but for individual indices - Table S3 ####
p_df = data.frame(matrix(nrow = 8, ncol = 6))
for (i in 1:8){
  print(paste0('sign_C', i, '_index_diff'))
  p_df[i, ] = test_SI_RE_dist(index_effects, method = 'mean',
                  sign_var = paste0('sign_C', i, '_index_diff'))
}

colnames(p_df) <- unique(SI_effects$continent)

p_df %>% 
  mutate(index = 1:8) %>%
  pivot_longer(cols = -index, names_to = 'continent') %>%
  arrange(value) %>%
  mutate(alpha_factor = 0.05/(48:1),
         signif_b = value < 0.05/48,
         signif_bh = value < alpha_factor,
         signif = value < 0.05) %>%
  filter(signif)
# bh stands for bonferroni-holm



## Distribution of SI_diff and Re diff  - count countries - Fig 4A ####
#(i.e. coloured by implementation/release of measures)

SI_Rediff_plot <- function(plot_data, SI_var = 'SI_diff',
                           sign_var = 'signCdiff'){
  
  p <- ggplot(plot_data) +
    geom_bar(data = plot_data %>% filter(get(sign_var) == 1),
             aes(x = roundR, group = -get(SI_var),  y = dummy, fill = get(SI_var)),
             stat = 'identity',
             position = 'stack', show.legend = T) +
    geom_bar(data = plot_data %>% filter(get(sign_var) == -1),
             aes(x = roundR, group = get(SI_var),  y = -dummy, fill = get(SI_var)),
             stat = 'identity',
             position = 'stack', show.legend = T) +
    scale_fill_gradient2(mid = 'white') +
    xlim(c(-1.5, 1.5)) +
    geom_vline(aes(xintercept = 0), colour = 'black', size = 0.5) +
    geom_hline(aes(yintercept = 0), colour = 'black', size = 0.5) +
    facet_grid(rows = vars(continent), scale = 'free') +
    labs(y = 'Number of countries', x = 'Change in Re over 7 days', 
         fill = 'Change in \nStringency\nIndex') +
    theme_minimal() +
    theme(axis.text.y= element_text(size=17),
          axis.text.x= element_text(size=17),
          strip.text.y =  element_text(size=20), #, angle = 0
          panel.spacing.y = unit(0.75, "lines"),
          axis.title.y =  element_text(size=20),
          axis.title.x =  element_text(size=20),
          legend.title = element_text(size=20),
          legend.text = element_text(size=17)
    )
  
  plotPath <- paste0(plot_path, '/SIdiff_Rediff_distribution_country_', SI_var, '.pdf')
  ggsave(plotPath, plot = p, width = 12, height = 18)
  return(p)
}

SI_Rediff_plot(plotSI)


for (i in 1:8){
  SI_Rediff_plot(index_effects, SI_var = paste0('C', i, '_index_diff'),
                 sign_var = paste0('sign_C', i, '_index_diff'))
}

## SI-diff through time (i.e. implementation/release of measures) - Fig. 4B ####
# coloured by R diff over following week

SI_time_plot <- function(plot_data, SI_var = 'SI_diff'){
  
    q <- ggplot(plot_data) +
    geom_bar(aes(x = date, fill = roundR, group = roundR, y = get(SI_var)),
             stat = 'identity', position = 'stack') +
    scale_fill_gradient2(mid = 'grey', limits = c(-1.5, 1.5)) +
    geom_hline(aes(yintercept = 0), colour = 'black', size = 0.5) +
    facet_grid(rows = vars(continent), scale = 'free_y') +
    labs(x = 'Date', fill = 'Change in \nRe over \n7 days', 
         y = 'Change in Stringency Index') +
    theme_minimal() +
    theme(
      axis.text.y= element_text(size=17),
      axis.text.x= element_text(size=17),
      strip.text.y =  element_text(size=20),
      panel.spacing.y = unit(0.75, "lines"),
      axis.title.y =  element_text(size=20),
      axis.title.x =  element_text(size=20),
      legend.title = element_text(size=20),
      legend.text = element_text(size=17)
    )
  
  plotPath <- paste0(plot_path, '/SIdiff_through_time_', SI_var, '.pdf')
  ggsave(plotPath, plot = q, width = 12, height = 18)
  return(q)
}

SI_time_plot(SI_effects)

for (i in 1:8){
  SI_time_plot(index_effects, SI_var = paste0('C', i, '_index_diff'))
}



## Linear regressions for all oxford indices - Fig. S12, Table S4 ####

plot_effects <- index_effects %>%
  dplyr::select(-paste0('C', 1:8, '_index'), -paste0('sign_C', 1:8, '_index_diff'),
         -discreteR, -signR, -roundR, -dummy) %>%
  pivot_longer(cols = paste0('C', 1:8, '_index_diff'), names_to = 'index', 
               values_to = 'index_value') %>%
  mutate(index_value = factor(round(index_value)))

index_names <- setNames(paste0('C', 1:8), sort(unique(plot_effects$index)))


ggplot(plot_effects, aes(y = Rlead7, x = index_value, 
                         group = index)) +
  geom_boxplot(aes(group = index_value, fill = continent), show.legend = F)  +
  geom_smooth(aes(colour = continent), method='lm', formula = y~x, show.legend = F) +
  geom_hline(aes(yintercept = 0), colour = 'black', size = 0.5) +
  geom_vline(aes(xintercept = 0), colour = 'black', size = 0.5) +
  facet_grid(rows = vars(index), cols = vars(continent),
             labeller = labeller(index = index_names)) +
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



p_df = data.frame(matrix(nrow = 8, ncol = 6))
adj_r_df = data.frame(matrix(nrow = 8, ncol = 6))
for(j in 1:6){
  cont = unique(index_effects$continent)[j]
    for (i in 1:8){
      
      sub_data <- index_effects %>%
        dplyr::select(-paste0('C', 1:8, '_index'), -paste0('sign_C', 1:8, '_index_diff'),
                      -discreteR, -signR, -roundR, -dummy) %>%
        filter(continent == cont,
               get(paste0('C', i, '_index_diff')) != 0 )
      
      fit <- lm(Rlead7 ~ get(paste0('C', i, '_index_diff')) , sub_data)
      
      
      p_val = summary(fit)$coefficients[2,4]
      p_df[i,j] = p_val
      
      if (p_val < 0.05){
        print('#################')
        print(paste0(cont, ' Index: C', i))
        print(summary(fit))
        adj_r_df[i,j] = summary(fit)$adj.r.squared
      } 
    }
}

colnames(p_df) <- unique(index_effects$continent)
colnames(adj_r_df) <- unique(index_effects$continent)

p_df %>% 
  mutate(index = 1:8) %>%
  pivot_longer(cols = -index, names_to = 'continent') %>%
  arrange(value) %>%
  mutate(alpha_factor = 0.05/(48:1),
         signif_b = value < 0.05/48,
         signif_bh = value < alpha_factor,
         signif = value < 0.05) %>%
  filter(signif)


