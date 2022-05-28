###########################################################
## bootstrap_feathers.R
## authors: J. Scire, D. Angst and J.S. Huisman
###########################################################

library(tidyverse)
library(lubridate)

library(ggplot2)
library(viridis)

all_estimates <- read_csv("../data/all_bootstrap_estimates.csv")

theme_set(theme_minimal() + 
              theme(
                strip.text = element_text(size=18),
                axis.text.y= element_text(size=18),
                axis.text.x= element_text(size=16),
                axis.title.y =  element_text(size=20),
                legend.text = element_text(size=19),
                panel.spacing = unit(2, "lines")
              )
)


plot_begin_date <- as_date("2020-09-01")
plot_end_date <- as_date("2021-04-30")

# CHE estimates 
CHE_estimates <- all_estimates %>% 
  filter(date > plot_begin_date,
         !(region %in% c("grR Eastern Switzerland", 
                         "grR Central Switzerland", 
                         "grR Espace Mittelland",
                         "grR Lake Geneva Region",
                         "grR Northwestern Switzerland",
                         "grR Ticino",
                         "grR Zurich"))) %>% 
  mutate(estimate_type = factor(estimate_type, 
      levels = c( "Cori_slidingWindow_EpiEstim", "Cori_slidingWindow"),
    labels = c("EpiEstim", "Combined CIs") ))  %>%
  mutate(truncation = "0",
         week_day = factor((wday(file_date) + 5) %% 7),
         month = month(file_date, label = TRUE))

last_Re <- CHE_estimates %>%
  group_by(estimate_type, file_date) %>%
  summarise(last_Re_from_file = max(date),
            .groups = 'drop')

first_Re <- CHE_estimates %>%
  group_by(estimate_type, date) %>%
  summarise(first_file_for_date = min(file_date),
            .groups = 'drop')

CHE_estimates <- CHE_estimates %>% 
  left_join(first_Re, by = c('estimate_type', 'date')) %>% 
  left_join(last_Re, by = c('estimate_type', 'file_date')) %>%
  filter(!(is.na(median_R_lowHPD) | is.na(median_R_highHPD)) )


max(CHE_estimates$file_date, na.rm = T)
min(CHE_estimates$file_date, na.rm = T)

### First feather plots ####
plot_Re_data <- CHE_estimates %>% 
  filter(date >= last_Re_from_file - 21)

Re_plot <- ggplot(plot_Re_data ,
       aes( x = date) ) +
  facet_wrap(vars(estimate_type)) +
  
  geom_ribbon(aes(ymin = median_R_lowHPD, ymax = median_R_highHPD, fill = week_day, group = file_date),
               alpha=0.1) +
  geom_line(aes(y = median_R_mean, colour = week_day, group = file_date), lwd = 0.5) +
  scale_x_date(date_labels = '%b\n%d',
               limits = c(plot_begin_date, plot_end_date) ) +
  coord_cartesian(ylim = c(0.5, 2.25)) +
  scale_colour_manual(values=rev(viridis(7)),
                      aesthetics = c("colour", "fill")) +
  xlab("") +
  ylab(("Estimated Re")) +
  theme(
    legend.position = "none" #remove legend
  )

Re_plot

ggsave("../figures/feathers_CHE_confirmed.png")

##################################
quantifyReOverlap <- function(bot.x, top.x, bot.y, top.y, method = 'x_in_y'){
  
  if (bot.x >= top.y){
    overlap = 0
  } else if (bot.y >= top.x){
    overlap = 0
  } else{
    
    if(method == 'equal'){
      #union
      total = max(top.x, top.y) - min(bot.x, bot.y)
    } else if (method == 'x_in_y'){
      total = top.x - bot.x
    } else if (method == 'y_in_x'){
      total = top.y - bot.y
    }
    
    intersection = min(top.x, top.y) - max(bot.x, bot.y)
    overlap = (intersection/total)
  }
  return(overlap)
}
##################################
most_recent <-  CHE_estimates %>% filter(!is.na(file_date)) %>% pull(file_date) %>% max()

estimates_firsttime <- CHE_estimates %>%
  group_by(region, date, data_type, estimate_type) %>%
  filter(file_date == min(file_date) ) %>%
  rename(first_median_R_highHPD = median_R_highHPD, first_median_R_lowHPD = median_R_lowHPD)

estimates_mostrecent <- CHE_estimates %>%
  filter(file_date == most_recent) %>%
  dplyr::select(region, date,  data_type, estimate_type, median_R_mean, median_R_highHPD, median_R_lowHPD)

estimates_coverage <- estimates_firsttime %>%
  left_join(estimates_mostrecent, by = c('region', 'date', 'data_type', 'estimate_type')) %>%
  mutate(#covered = (median_R_mean.x > median_R_lowHPD) & (median_R_mean.x < median_R_highHPD),
         covered = (median_R_mean.y > first_median_R_lowHPD) & (median_R_mean.y < first_median_R_highHPD),
         ci_covered = quantifyReOverlap(first_median_R_lowHPD, first_median_R_highHPD, 
                                        median_R_lowHPD, median_R_highHPD, method = 'x_in_y'),
         ci_inv_covered = quantifyReOverlap(first_median_R_lowHPD, first_median_R_highHPD, 
                                        median_R_lowHPD, median_R_highHPD, method = 'y_in_x')) 


#### CI 

CI_subset <- estimates_coverage %>% 
  ungroup() %>%
  filter( region == "CHE", data_type == 'Confirmed cases' ) %>%
  dplyr::select(date, file_date, estimate_type, covered, ci_covered) %>%
  group_by(date) %>% 
  slice_min(order_by = file_date)

CI_plot <- ggplot(CI_subset, aes( x = date) ) +
  facet_wrap(vars(estimate_type), nrow =1) +
  geom_bar(aes(y = ci_covered*100, fill = covered), position =  'identity', stat = 'identity') +
  geom_point(aes(y = ci_covered*100, colour = covered), size = 3, show.legend = F) +
  scale_x_date(date_labels = '%b\n%d',
               limits = c(plot_begin_date, plot_end_date)) +
  scale_fill_manual(values = viridis(2)[c(2, 1)], aesthetics = c('colour', 'fill')) +
  labs(y ="% first CI in stable CI", x = "", fill = 'Stable Re point estimate in first CI') +
  theme(
    strip.text.x = element_blank(),
    legend.title = element_text(size=20),
    legend.text = element_text(size=19),
    legend.position = "bottom" #remove legend
  )

CI_plot
ggsave(file.path(WORKDIR, "plots", "bootstrap_ci_coverage_CHE_sliding.png"), units = "cm", width = 40, height = 25)

###########################################################
library(patchwork)
Re_plot + CI_plot +
  plot_layout(heights = c(1, 1), ncol = 1) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 30, hjust = 1, vjust = 0.5) )

ggsave("../figures/Figure2.pdf", width = 17, height = 12)
