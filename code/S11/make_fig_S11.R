library("lubridate")
library("readr")
library("gridExtra")
library("reshape2")
library("gdata")
library("fitdistrplus")
library("here")
library("viridis")
library("tidyverse")

data_dir <- "placeholder"
plot_dir <- "placeholder"

estimates_with_imports <- read_rds(file.path(data_dir, "CHE-Estimates_imports.rds"))
estimates_without_imports <- read_rds(file.path(data_dir, "CHE-estimates_no_imports.rds"))


estimates_with_imports <- estimates_with_imports %>% mutate(imports_included = "TRUE")
estimates_without_imports <- estimates_without_imports %>% mutate(imports_included = "FALSE")

estimates_CH <- bind_rows(estimates_with_imports, estimates_without_imports) %>% filter(region == "CHE", estimate_type == "Cori_slidingWindow", data_type == "Confirmed cases")

color_scale <- viridis(6)

ggplot(estimates_CH) +
  geom_line(aes(x = date, y = median_R_mean, colour = imports_included), size = 1.1) +
  geom_ribbon(aes(x = date, ymax = median_R_highHPD, ymin = median_R_lowHPD, fill = imports_included), alpha = 0.15, colour = NA) +
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_x_date(date_breaks = "1 month",
               date_labels = '%d %b\n%Y',
               limits = c(as.Date("2020-03-01"), as.Date("2021-04-01"))) +
  coord_cartesian(ylim=c(0,3.3)) +
  theme_bw() +
  xlab("") +
  ylab("Reproductive number") +
  theme(
    axis.text.y= element_text(size=14),
    axis.text.x= element_text(size=14),
    axis.title.y =  element_text(size=17),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14),
    legend.position="top"
  ) +
  scale_colour_manual(values=color_scale[c(1,4)],
                      labels=c("Ignoring import cases", "Accounting for imports"),
                      breaks=c("FALSE", "TRUE"),
                      name  ="",
                      aesthetics = c("fill", "color"))
