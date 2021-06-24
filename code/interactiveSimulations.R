###########################################################
## interactiveSimulations.R
## author: J.S. Huisman
###########################################################

library(tidyverse)
library(lubridate)
library(viridis)
library(cowplot)

plot_path = '../figures'
simulationDir = '../simulations'

source('generateSimulations.R')
###########################################################

###### Plot Validation ######

# 6 times, 3 plateaus; 9 times, 5 values/4 plateaus
shift_times = c(0, 25, 50, 60, 80, 100, 110, 120, 130, 160)
R_levels = c(2.5, 0.5, 1.2, 0.95, 1.1)

IncubationParams <- getGammaParams(meanParam = 5.3, sdParam = 3.2)
OnsetToCountParams = getGammaParams(4.5, 4.9)

simulation <- simulateTS(shift_times, R_levels,
                         IncubationParams, OnsetToCountParams, init_infection = 10,
                         noise = list(), smooth_R = FALSE,
                         timevarying = FALSE)

longSim <- simulation %>%
  pivot_longer(cols = c(Re, infections, observations),
               names_to = 'type') %>%
  mutate(type = factor(R.utils::capitalize(type), levels = c('Re', 'Infections', 'Observations')),
         source_type = 'simulation',
         plot_row = ifelse(type == 'Re', 'Re', 'Cases'),
         plot_row = factor(plot_row, levels = c('Re', 'Cases')) )

estimatedInfections <- estimateInfectionTS(simulation, IncubationParams, OnsetToCountParams,
                                           smooth_param = TRUE, fixed_shift = FALSE,
                                           timevarying = FALSE, n_boot = 100)
estimatedRe <- estimateReTS(estimatedInfections, delay = 0)

longInfections <- estimatedInfections %>%
  dplyr::select(c(date, replicate, value)) %>% 
  group_by(date) %>%
  summarize(
    median_val = median(value),
    high_quant = quantile(value, probs=0.975, na.rm=T),
    low_quant = quantile(value, probs=0.025, na.rm=T),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  mutate(type = factor('Infections', levels = c('Re', 'Infections', 'Observations'), ordered = TRUE),
         source_type = 'estimation',
         plot_row = factor('Cases', levels = c('Re', 'Cases'), ordered = TRUE))

longRe <- cleanReTSestimate(estimatedRe) %>%
  dplyr::select(-c(country, region, source, data_type,estimate_type)) %>% 
  mutate(type = factor('Re', levels = c('Re', 'Infections', 'Observations'), ordered = TRUE),
         source_type = 'estimation',
         plot_row = factor('Re', levels = c('Re', 'Cases')) )

ggplot() +
  geom_line(data = longSim, aes(x = date, y = value, colour = type), size = 2) +
  geom_ribbon(data = longRe %>% filter(median_R_mean < 4), aes(x = date, ymin = median_R_lowHPD,
                                 ymax = median_R_highHPD), alpha = 0.7) +
  geom_line(data = longRe %>% filter(median_R_mean < 4), aes(x = date, y = median_R_mean), size=1.1) +
  geom_ribbon(data = longInfections, aes(x = date, ymin = low_quant,
                                         ymax = high_quant), alpha = 0.7) +
  geom_line(data = longInfections, aes(x = date, y = median_val), size=1.1) +
  labs( x = 'Date', colour = 'Simulated') +
  scale_colour_viridis(direction = -1, begin = 0.3, discrete = T) + 
  facet_grid(rows = vars(plot_row), scale = 'free', switch = 'y') +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = 25),
    strip.placement = 'outside',
    axis.text.y= element_text(size=20),
    axis.text.x= element_text(size=20),
    axis.title.y =  element_blank(), #element_text(size=17),
    axis.title.x =  element_text(size=25),
    legend.title = element_text(size=25),
    legend.text = element_text(size=20),
    legend.position = c(0.8, 0.8)
  )

plotPath <- paste0(plot_path, "/First_simulations.pdf")
ggsave(plotPath, width = 8, height = 8)

