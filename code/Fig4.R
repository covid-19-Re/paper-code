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

mobData <- qread('../data/allMobilityDataGoogle.qs')

mob_df <- mobData %>%
  mutate(week = isoweek(date),
         year = year(date)) %>%
  group_by(countryIso3, region, data_type, week, year) %>%
  summarise(mob_value = mean(change),
            date = max(date),
            .groups = 'drop') %>%
  arrange(countryIso3, date, data_type)

# ggplot(mobData %>% filter(countryIso3 == 'NLD')) +
#   geom_point(aes(x = date, y = change)) +
#   facet_wrap(vars(data_type))

###########################################################
#url = 'https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv'

oxfordGRT <- read_csv('../data/oxfordGRT.csv', #url, 
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
         Rpast7 = R - lag(R, 7),
         Rlead14 = lead(R, 14) - R) %>%
  right_join(SI_df, by = c('countryIso3', 'date')) %>%
  mutate(country = ifelse(!is.na(country.y), country.y, country.x) ) %>%
  dplyr::select(-country.x, - country.y) %>%
  mutate(#roundR7 = round(Rlead7, digits = 1),
         #roundRpast7 = round(Rpast7, digits = 1),
         #roundR14 = round(Rlead14, digits = 1),
         #signCdiff = factor(sign(SI_diff), levels = c(-1, 0, 1)),
         dummy = 1,
         roundSI = round(SI_diff, digits = -1)) %>%
  filter(!is.na(continent))

plotSI <- SI_effects %>% 
  filter(!is.na(countryIso3)) %>%
  pivot_longer(c('Rlead7', 'Rpast7', 'Rlead14'), names_to = 'R_var', values_to = 'R_val') %>%
  mutate(R_var = factor(R_var, levels = c('Rpast7', 'Rlead7', 'Rlead14'),
                        labels = c('R(t)', 'R(t+7)', 'R(t+14)'))) %>%
  group_by(continent, roundSI, R_var) %>%
  summarise(Rmed = median(R_val, na.rm = T),
            R_low = quantile(R_val, 0.25, na.rm = T),
            R_high = quantile(R_val, 0.75, na.rm = T))

###########################################################
## NPI abolishment/implementation #####

SI_plot <- ggplot(plotSI) +
  # geom_violin(aes(x = R_val, y = roundSI, group = interaction(roundSI, R_var), colour = R_var),
  #             draw_quantiles = c(0.5),
  #             show.legend = T, orientation = 'y') +
  geom_errorbar(aes(xmin = R_low, xmax = R_high, y = roundSI, colour = R_var),
               show.legend = T, orientation = 'y') +
  geom_point(aes(x = Rmed, y = roundSI, colour = R_var), show.legend = T) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-50, 50)) +
  scale_colour_manual(values = viridis(4)[1:3])+
  geom_vline(aes(xintercept = 0), colour = 'black', size = 0.5) +
  geom_hline(aes(yintercept = 0), colour = 'black', size = 0.5) +
  facet_grid(cols = vars(continent), rows = vars(R_var), scale = 'free') +
  labs(y = 'Change in SI', x = 'Change in Re', 
       colour = 'Change in Re') +
  theme_minimal() +
  theme(legend.position = 'bottom',
    axis.text.y= element_text(size=17),
        axis.text.x= element_text(size=17, angle = 45, hjust = 1, vjust = 1),
        strip.text =  element_text(size=20), #, angle = 0
        panel.spacing.y = unit(0.75, "lines"),
        axis.title.y =  element_text(size=20),
        axis.title.x =  element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=17)
  )
SI_plot

plotPath <- paste0(plot_path, '/SIdiff_Rediff_distribution_country_SI_diff.pdf')
ggsave(plotPath, width = 17, height = 10)
#ggsave(plotPath, plot = p, width = 12, height = 8)

###########################################################
plotMob <- cleanRe %>%
  mutate(Rpast7 = R - lag(R, 7),
         Rlead7 = lead(R, 7) - R,
         Rlead14 = lead(R, 14) - R) %>%
  right_join(mob_df, by = c('countryIso3', 'date')) %>%
  mutate(roundMob = round(mob_value, digits = -1)) %>%
  filter(!is.na(continent), 
         !is.na(countryIso3),
         date <= as_date('2021-05-03')) %>%
  pivot_longer(c('Rlead7', 'Rpast7', 'Rlead14'), names_to = 'R_var', values_to = 'R_val') %>%
  mutate(R_var = factor(R_var, levels = c('Rpast7', 'Rlead7', 'Rlead14'),
                        labels = c('R(t)', 'R(t+7)', 'R(t+14)'))) %>%
  group_by(continent, data_type, roundMob, R_var) %>%
  summarise(Rmed = median(R_val, na.rm = T),
            R_low = quantile(R_val, 0.25, na.rm = T),
            R_high = quantile(R_val, 0.75, na.rm = T)) %>%
  mutate(data_type = factor(data_type, levels = c("Grocery And Pharmacy","Parks", 
                          "Residential", "Retail And Recreation", "Transit Stations","Workplaces" ),
                          labels = c("Grocery", "Parks", "Residential", 
                                     "Retail And Rec.", "Transit","Workplaces" ))) %>%
  filter(R_var == 'R(t)')

Mob_plot <- ggplot(plotMob %>% filter(continent == 'Europe')) +
  #ggplot(plotMob %>% filter(data_type %in% c('Grocery', 'Residential', 'Workplaces'))) +
  geom_errorbar(aes(xmin = R_low, xmax = R_high, y = roundMob, colour = R_var),
                show.legend = F, orientation = 'y') +
  geom_point(aes(x = Rmed, y = roundMob, colour = R_var), show.legend = F) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-100, 100)) +
  scale_colour_manual(values = viridis(4)[1:3])+
  geom_vline(aes(xintercept = 0), colour = 'black', size = 0.5) +
  geom_hline(aes(yintercept = 0), colour = 'black', size = 0.5) +
  facet_grid(cols = vars(data_type), rows = vars(R_var), scale = 'free') +
  #facet_grid(cols = vars(continent), rows = vars(data_type), scale = 'free') +
  labs(x = 'Change in Re', y = 'Mobility wrt. baseline\n(Europe)', 
       colour = 'Change in Re') +
  theme_minimal() +
  theme(axis.text.y= element_text(size=17),
        axis.text.x= element_text(size=17, angle = 45, hjust = 1, vjust = 1),
        strip.text =  element_text(size=20), #, angle = 0
        panel.spacing = unit(0.75, "lines"),
        axis.title.y =  element_text(size=20),
        axis.title.x =  element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=17)
  )
Mob_plot
plotPath <- paste0(plot_path, '/Mobdiff_Rediff_country.pdf')
ggsave(plotPath, width = 17, height = 10)

##
library(patchwork)

SI_plot + Mob_plot +
  plot_layout(ncol =1, heights = c(3, 1), guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 25),
        legend.position = 'bottom')

plotPath <- paste0(plot_path, '/Mob_SI_Rediff.pdf')
ggsave(plotPath, width = 17, height = 16)

######################################################################################################################
library(lme4)
library(lmerTest)

###########################################################

SI_mod <- lmer(Rlead14 ~ SI_diff + (1|countryIso3), SI_effects)
summary(SI_mod)
anova(SI_mod)

###########################################################
mob_stats_df <- cleanRe %>%
  mutate(Rpast7 = R - lag(R, 7),
         Rlead7 = lead(R, 7) - R,
         Rlead14 = lead(R, 14) - R) %>%
  right_join(mob_df, by = c('countryIso3', 'date')) %>%
  mutate(roundMob = round(mob_value, digits = -1)) %>%
  filter(!is.na(continent), 
         !is.na(countryIso3),
         date <= as_date('2021-05-03')) 


lme_mod <- lmer(Rpast7 ~ mob_value + (1|countryIso3) , mob_stats_df %>% filter(data_type == 'Workplaces') )
summary(lme_mod)
anova(lme_mod) # for the F value
confint(lme_mod, oldNames = FALSE) # for the CIs
ranef(lme_mod) # for the random effects (per country)
plot(lme_mod)

qqnorm(ranef(lme_mod)$countryIso3[,"(Intercept)"], 
       main = "Random effects")
qqnorm(resid(lme_mod), main = "Residuals")

# new model
lme_mod <- lmer(Rpast7 ~ mob_value:data_type + (1|continent/countryIso3), mob_stats_df)
summary(lme_mod)
anova(lme_mod) # for the F value
confint(lme_mod, oldNames = FALSE) # for the CIs
ranef(lme_mod) # for the random effects (per country)
continent_df <- rownames_to_column(ranef(lme_mod)$countryIso3) %>%
  separate('rowname', c('countryIso3', 'continent'), sep = ':') %>%
   rename(#countryIso3 = rowname,
          effect = "(Intercept)") #%>%
  # left_join(continentData, by = 'countryIso3')
ggplot(continent_df) +
  geom_point(aes(x = continent, y = effect))

plot(lme_mod)


######################################################################################################################
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



## Direct relations ####

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


