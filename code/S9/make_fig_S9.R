library("lubridate")
library("readr")
library("gridExtra")
library("reshape2")
library("gdata")
library("fitdistrplus")
library("here")
library("tidyverse")

get_dates_to_average_over <- function(i, dates, weeks_averaged) {
  if (i < (weeks_averaged * 7)) {
    return(dates[seq_len(weeks_averaged * 7)])
  } else {
    return(dates[(i - (weeks_averaged * 7)):i])
  }
}

#TODO set
delays_data_path <- here::here("placeholder.csv")

delays_onset_to_count <- read_csv(delays_data_path,
                                  col_types = cols(
                                    data_type = col_character(),
                                    onset_date = col_date(format = ""),
                                    count_date = col_date(format = ""),
                                    delay = col_number(),
                                    country = col_character(),
                                    region = col_character(),
                                    countryIso3 = col_character(),
                                    source = col_character()))
latestDate <- Sys.Date()

swiss_delays <- delays_onset_to_count %>%
  filter(country == "Switzerland") %>%
  filter(onset_date < latestDate - 15)



getMeansDelays <- function(onset_to_report_empirical_delays,
                              min_number_cases = 300){

  all_dates <- seq(min(onset_to_report_empirical_delays$onset_date), max(onset_to_report_empirical_delays$onset_date), by = "days")
  mean_gamma_fit <- c()
  raw_mean <- c()
  raw_mean_last_n_cases <- c()

  for(i in 1:length(all_dates)) {

    weeks_averaged <- 0
    repeat{
      weeks_averaged <- weeks_averaged + 1
      recent_counts_distribution <- onset_to_report_empirical_delays %>%
        dplyr::filter( onset_date %in% get_dates_to_average_over(i, all_dates, weeks_averaged) )

      if(nrow( recent_counts_distribution ) >= min_number_cases) {
        break
      }
    }

    recent_delays <- recent_counts_distribution %>% pull(delay)

    gamma_fit <- fitdist(recent_delays + 1, distr = "gamma", method = "mme")

    shape_fit <- gamma_fit$estimate["shape"]
    rate_fit <- gamma_fit$estimate["rate"]

    mean_fit <- shape_fit/rate_fit - 1

    mean_gamma_fit <- c(mean_gamma_fit, mean_fit)
    raw_mean <- c(raw_mean, mean(recent_delays))

    ## now just take last n cases
    if(nrow(onset_to_report_empirical_delays %>% filter( onset_date <= all_dates[i])) >= min_number_cases) {
      last_n_delays <- onset_to_report_empirical_delays %>% filter( onset_date <= all_dates[i]) %>% arrange(onset_date) %>% slice_tail(n = min_number_cases) %>%  pull(delay)
    } else {
      last_n_delays <- onset_to_report_empirical_delays %>% arrange(onset_date) %>% slice_head(n = min_number_cases) %>% pull(delay)
    }

    raw_mean_last_n_cases <- c(raw_mean_last_n_cases, mean(last_n_delays))
  }

  names(mean_gamma_fit) <- c()

  result <- bind_rows(tibble(date = all_dates, value = mean_gamma_fit, variable = "Mean - Gamma fit"),
                      tibble(date = all_dates, value = raw_mean, variable = "Mean - Full weeks"),
                      tibble(date = all_dates, value = raw_mean_last_n_cases, variable = "Mean - last N"))

  return(result)
}


min_number_cases <- 300

onset_to_report_empirical_delays <- swiss_delays %>% filter(data_type == "Confirmed cases")
a <- getMeansDelays(onset_to_report_empirical_delays, min_number_cases = min_number_cases)
confirmed_means <- a %>% mutate(data_type = "Confirmed cases")

onset_to_report_empirical_delays <- swiss_delays %>% filter(data_type == "Deaths")
a <- getMeansDelays(onset_to_report_empirical_delays,  min_number_cases = min_number_cases)
deaths_means <-  a %>% mutate(data_type = "Deaths")

onset_to_report_empirical_delays <- swiss_delays %>% filter(data_type == "Hospitalized patients")
a <- getMeansDelays(onset_to_report_empirical_delays,  min_number_cases = min_number_cases)
hosp_means <- a %>% mutate(data_type = "Hospitalized patients")

mean_fits <- bind_rows(confirmed_means, deaths_means, hosp_means)

#SI paper figure
ggplot(data = filter(mean_fits, variable == "Mean - last N")) +
  geom_line(aes(x = date, y = value, colour = data_type), lwd = 2) +
  ylab("Number of days") +
  xlab("") +
  theme_minimal() +
  coord_cartesian(ylim=c(0,17)) +
  scale_y_continuous(breaks=seq(0,18,by=2)) +
  scale_x_date(date_breaks = "1 month",
               date_labels = '%b\n%d\n%Y') +
  theme(
    axis.text.y= element_text(size=17),
    axis.text.x= element_text(size=17),
    axis.title.y =  element_text(size=18),
    legend.text = element_text(size=16)
  ) +
  labs(colour="") +
  theme(legend.position="top")
