library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(readxl)
library(viridis)
library(lubridate)
library(ggpubr)

DATADIR <- "placeholder"
OUTDIR <- "placeholder"

country_code <- "CHE"

estimates_without_onset_variable_delays <- read_rds(file.path(DATADIR, "variable_delays_without_onset", paste0(country_code, "-Estimates.rds")))
estimates_without_onset_variable_delays <- estimates_without_onset_variable_delays %>%
  filter(estimate_type == "Cori_slidingWindow", data_type == "Confirmed cases", region == country_code) %>%
  mutate(analysis_type = "Without symptom onset data - Variable delays",
         include_onset = F,
         variable_delays = T)

estimates_without_onset_fixed_delays <- read_rds(file.path(DATADIR, "fixed_delays_without_onset", paste0(country_code, "-Estimates.rds")))
estimates_without_onset_fixed_delays <- estimates_without_onset_fixed_delays %>%
  filter(estimate_type == "Cori_slidingWindow", data_type == "Confirmed cases", region == country_code) %>%
  mutate(analysis_type = "Without symptom onset data - Fixed delays",
         include_onset = F,
         variable_delays = F)

estimates_with_onset_fixed_delays <- read_rds(file.path(DATADIR, "fixed_delays_with_onset", paste0(country_code, "-Estimates.rds")))
estimates_with_onset_fixed_delays <- estimates_with_onset_fixed_delays %>%
  filter(estimate_type == "Cori_slidingWindow", data_type == "Confirmed cases", region == country_code) %>%
  mutate(analysis_type = "With symptom onset data - Fixed delays",
         include_onset = T,
         variable_delays = F)

estimates_with_onset_variable_delays <- read_rds(file.path(DATADIR, "variable_delays_with_onset", paste0(country_code, "-Estimates.rds")))
estimates_with_onset_variable_delays <- estimates_with_onset_variable_delays %>%
  filter(estimate_type == "Cori_slidingWindow", data_type == "Confirmed cases", region == country_code) %>%
  mutate(analysis_type = "With symptom onset data - Variable delays",
         include_onset = T,
         variable_delays = T)

estimates <- bind_rows(estimates_with_onset_fixed_delays, estimates_without_onset_fixed_delays, estimates_with_onset_variable_delays, estimates_without_onset_variable_delays) %>%
  dplyr::select(date, median_R_mean, median_R_highHPD, median_R_lowHPD, analysis_type, include_onset, variable_delays)

estimates$analysis_type <- factor(estimates$analysis_type,
                                  c("Without symptom onset data - Fixed delays",
                                    "Without symptom onset data - Variable delays",
                                    "With symptom onset data - Fixed delays",
                                    "With symptom onset data - Variable delays"))
min_plot_date <- as.Date("2020-03-01")
max_plot_date <-  as.Date("2021-04-20")

ggplot(estimates, aes(x = date)) +
  geom_line(aes(y = median_R_mean, colour = analysis_type), lwd = 1.0) +
  geom_ribbon(aes(ymin = median_R_lowHPD, ymax = median_R_highHPD, fill = analysis_type),
              alpha=0.1) +
  scale_colour_manual(values=rev(viridis(4)),
                      breaks = c("With symptom onset data - Fixed delays",
                                 "With symptom onset data - Variable delays",
                                 "Without symptom onset data - Fixed delays",
                                 "Without symptom onset data - Variable delays"),
                      name  ="",
                      aesthetics = c("colour", "fill")) +
  scale_x_date(date_breaks = "1 month",
               date_labels = '%b-%d\n%Y',
               limits = c(min_plot_date, max_plot_date)) +
  ylab("Reproductive number") +
  coord_cartesian(ylim = c(0, 3.5)) +
  xlab("") +
  theme_bw() +
  theme(
    legend.position="top",
    strip.background = element_blank(),
    strip.text.y = element_text(size=18),
    axis.text.y= element_text(size=18),
    axis.text.x= element_text(size=14),
    axis.title.y =  element_text(size=20),
    legend.text = element_text(size=15)
  ) +
  guides(fill=guide_legend(ncol=2)) +
  guides(colour=guide_legend(ncol=2))
