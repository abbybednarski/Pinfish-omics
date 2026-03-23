### code from the one and only Logan F. Turner ###

#### getting buoy data

#let's download the data from NOAA's website
#########
#NOTE: for all files, air temp is in column 14 and water temp is in column 15.
#Therefore, we will focus on column 15 for our records
#########

library(dplyr)
library(lubridate)

#For Orange Beach, AL
#Station 42012 - Orange Beach - 44 nm SE of Mobile, AL
#Data available from 2009-2024 so going to use 2021-24 to get good estimate
#2021 data:
DI21_data = read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=42012h2021.txt.gz&dir=data/historical/stdmet/")
head(DI21_data)

#extract column 15 info with sea surface temps
DI21_wtmp = DI21_data[,15]

#make 999 = NA
DI21_wtmp = replace(DI21_wtmp, DI21_wtmp==999, NA)
head(DI21_wtmp, n=100)

max(DI21_wtmp, na.rm = TRUE)
#33
min(DI21_wtmp, na.rm = TRUE)
#14.6
mean(DI21_wtmp, na.rm = TRUE)
#23.60275

#2022 data:
DI22_data = read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=42012h2022.txt.gz&dir=data/historical/stdmet/")
DI22_wtmp = DI22_data[,15]
DI22_wtmp = replace(DI22_wtmp, DI22_wtmp==999, NA)
head(DI22_wtmp, n=100)
max(DI22_wtmp, na.rm = TRUE)
#31.8
min(DI22_wtmp, na.rm = TRUE)
#13.5
mean(DI22_wtmp, na.rm = TRUE)
#22.60599


#2023 data:
DI23_data = read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=42012h2023.txt.gz&dir=data/historical/stdmet/")
DI23_wtmp = DI23_data[,15]
DI23_wtmp = replace(DI23_wtmp, DI23_wtmp==999, NA)
head(DI23_wtmp, n=100)
max(DI23_wtmp, na.rm = TRUE)
#32.8
min(DI23_wtmp, na.rm = TRUE)
#15.5
mean(DI23_wtmp, na.rm = TRUE)
#24.03651

#2024 data:
DI24_data = read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=42012h2024.txt.gz&dir=data/historical/stdmet/")
DI24_wtmp = DI24_data[,15]
DI24_wtmp = replace(DI24_wtmp, DI24_wtmp==999, NA)
head(DI24_wtmp, n=100)
max(DI24_wtmp, na.rm = TRUE)
#32.7
min(DI24_wtmp, na.rm = TRUE)
#14.9
mean(DI24_wtmp, na.rm = TRUE)
#24.42861

#I want to find average of all max temps for the data I input so here is the code for that - I also want to do this with min so that code will follow
max_21 <- max(DI21_wtmp, na.rm = TRUE)
max_22 <- max(DI22_wtmp, na.rm = TRUE)
max_23 <- max(DI23_wtmp, na.rm = TRUE)
max_24 <- max(DI24_wtmp, na.rm = TRUE)

avg_max <- mean(c(max_21, max_22, max_23, max_24))
avg_max
#32.575

#avg min below
min_21 <- min(DI21_wtmp, na.rm = TRUE)
min_22 <- min(DI22_wtmp, na.rm = TRUE)
min_23 <- min(DI23_wtmp, na.rm = TRUE)
min_24 <- min(DI24_wtmp, na.rm = TRUE)

avg_min <- mean(c(min_21, min_22, min_23, min_24))
avg_min
#14.625

#avg mean below
mean_21 <- mean(DI21_wtmp, na.rm = TRUE)
mean_22 <- mean(DI22_wtmp, na.rm = TRUE)
mean_23 <- mean(DI23_wtmp, na.rm = TRUE)
mean_24 <- mean(DI24_wtmp, na.rm = TRUE)

avg_mean <- mean(c(mean_21, mean_22, mean_23, mean_24))
avg_mean
#23.66847

#find hourly temp fluctuations in one of the years (picked most recent so 2024)

#clean temp in full df
DI24_data$V15[DI24_data$V15 >= 999] <- NA

#make sure data is in time order
DI24_data <- DI24_data[order(DI24_data$V1, DI24_data$V2, DI24_data$V3, DI24_data$V4), ]

#calc hourly change
temp_change <- diff(DI24_data$V15)

#hourly max change
max(abs(temp_change), na.rm = TRUE)
#1.2

#dont need anything beyond here - this is Dr. Logan Turners work after this for CTmax so needed for in depth temp info!
  #keeping code in case the next student needs :) good luck!

#clean the data for violin plots
#DI23_clean = DI23_data %>% filter(DI23_data[[15]] !=999)

#let's merge our cleaned datasets to make one Dauphin Island dataset
#merged_DI <- bind_rows(DI21_clean, DI22_clean, DI23_clean)
#comparison_DI = bind_rows(DI22_clean, DI23_clean)


#Great, now those are all loaded in. 
#Let's look at the range of data for each site and year.

#let's get the min and max values from each site across all sampling times

# Dauphin Island 
#max(merged_DI[15])
#max(DI21_clean[15])
#max(DI22_clean[15])
#max(DI23_clean[15])

#Port Aransas
#max(merged_PA[15])
#max(PA23_clean[15])
#max(PA22_clean[15])
#max(PA21_clean[15])
#max(PA17_clean[15])

#Sant Augustine
#max(merged_SA[15])

#What if we try to make a violin plot of the data across years?

library(ggplot2)

#For Dauphin Island datasets
#ggplot(merged_DI, aes(x = as.factor(merged_DI[[1]]), y = merged_DI[[15]], fill = as.factor(merged_DI[[1]]))) +
#  geom_violin(trim = FALSE) +
#  labs(title = "Dauphin Island Sea Surface Temperature Distribution",
#       x = "Year",
#       y = "Sea Surface Temperature (°C)",
#       fill = "Year") +
#  theme_minimal() +
#  scale_fill_brewer(palette = "Set3")

#For Port A
#ggplot(merged_PA, aes(x = as.factor(merged_PA[[1]]), y = merged_PA[[15]], fill = as.factor(merged_PA[[1]]))) +
#  geom_violin(trim = FALSE) +
#  labs(title = "Port Aransas Sea Surface Temperature Distribution",
#       x = "Year",
#       y = "Sea Surface Temperature (°C)",
#       fill = "Year") +
#  theme_minimal() +
#  scale_fill_brewer(palette = "Set3")

#For Saint Augustine 
#ggplot(merged_SA, aes(x = as.factor(merged_SA[[1]]), y = merged_SA[[15]], fill = as.factor(merged_SA[[1]]))) +
#  geom_violin(trim = FALSE) +
#  labs(title = "Saint Augustine Sea Surface Temperature Distribution",
#       x = "Year",
#       y = "Sea Surface Temperature (°C)",
#       fill = "Year") +
#  theme_minimal() +
#  scale_fill_brewer(palette = "Set3")


########################
# now let's look at everything by hour, day, etc and compare them to one another

# Function to process any buoy dataframe
#process_buoy_data <- function(df, temp_col = "WTMP") {
#  df %>%
 #   mutate(
#      datetime = make_datetime(V1, V2, V3, V4, V5),    # Combine to full datetime
#      WTMP = as.numeric(V15)             # Make sure WTMP is numeric
 #   ) %>%
#    filter(WTMP > 30) %>%                              # Filter to SST > 30°C
#    mutate(
#      date = as_date(datetime),
#      hour = hour(datetime)
#    ) %>%
#    group_by(date, hour) %>%
#    summarize(
#      n_obs = n(),                                     # Number of hot obs in that hour
  #    mean_temp = mean(WTMP, na.rm = TRUE),
 #     .groups = "drop"
 #   )
#}

#PA_hot30 = process_buoy_data(merged_PA)
#summary(PA_hot30)

#PA_hot30 %>%
#  group_by(date) %>%
#  summarize(hot_hours = n())

#PA_hot30 %>%
 # arrange(date, hour)

#SA_hot30 = process_buoy_data(merged_SA)

#SA_hot30 %>%
 # group_by(date) %>%
 # summarize(hot_hours = n())

#SA_hot30 %>%
 # arrange(date, hour)

#DI_hot30 = process_buoy_data(merged_DI)

#DI_hot30 %>%
 # group_by(date) %>%
 # summarize(hot_hours = n())

#DI_hot30 %>%
 # arrange(date, hour)

#summarize site data into days and hours above 30C
#summarize_site <- function(df, site_name) {
 # df %>%
 #   mutate(year = year(date)) %>%
 #   group_by(site = site_name, year) %>%
 #   summarize(
 #     days_above_30 = n_distinct(date),
 #     hours_above_30 = n(),  # each row = 1 hour
 #     .groups = "drop"
 #   )
#}

#PA_summary <- summarize_site(PA_hot30, "TX")
#DI_summary <- summarize_site(DI_hot30, "AL")
#SA_summary <- summarize_site(SA_hot30, "FL")

#combine into one summary
#combined_summary <- bind_rows(PA_summary, SA_summary, DI_summary)

#library(tidyr)

#wide_summary <- combined_summary %>%
#  pivot_longer(cols = c(days_above_30, hours_above_30), names_to = "metric") %>%
#  unite("year_metric", year, metric) %>%
#  pivot_wider(names_from = year_metric, values_from = value)

# First, reshape your combined summary back into long format if needed
#long_summary <- combined_summary %>%
#  pivot_longer(cols = c(days_above_30, hours_above_30),
#               names_to = "metric", values_to = "value")

# bar plots of temp data
#assign colors
#state_colors <- c(
#  "AL" = "#2CAFA2",  # seafoam/teal
#  "TX" = "#E07B3C",  # orange
#  "FL" = "#2C72B0"   # blue
#)

#days above 30 C
#ggplot(combined_summary, aes(x = factor(year), y = days_above_30, fill = site)) +
#  geom_bar(stat = "identity", position = position_dodge()) +
#  labs(
#    title = "Days with Sea Temperature > 30°C",
#    x = "Year",
#    y = "Number of Days",
#    fill = "Site"
#  ) +
#  theme_minimal() +
#  scale_fill_manual(values = state_colors)

#hours above 30 C
#ggplot(combined_summary, aes(x = factor(year), y = hours_above_30, fill = site)) +
#  geom_bar(stat = "identity", position = position_dodge()) +
#  labs(
#    title = "Hours with Sea Temperature > 30°C",
#    x = "Year",
#    y = "Number of Hours",
#    fill = "Site"
#  ) +
#  theme_minimal() +
#  scale_fill_manual(values = state_colors)

# Heatmap for days above 30
#ggplot(PA_summary %>% bind_rows(DI_summary, SA_summary), 
#       aes(x = factor(year), y = site, fill = days_above_30)) +
#  geom_tile(color = "white") +
#  scale_fill_gradient(low = "white", high = "red") +
#  labs(title = "Days with SST > 30°C", x = "Year", y = "Site") +
#  theme_minimal()

#let's compare and get values from each site at each year. 
#we can use the data subsets that we already have and run those. 

#2019
#PA19_hot = process_buoy_data(PA19_clean)
#SA19_hot = process_buoy_data(SA19_clean)

#PA19_summary = summarize_site(PA19_hot, "TX")
#SA19_summary = summarize_site(SA19_hot, "FL")

#summary19 = bind_rows(PA19_summary, SA19_summary)

#print(summary19)

#2020
#PA20_hot = process_buoy_data(PA20_clean)
#SA20_hot = process_buoy_data(SA20_clean)

#PA20_summary = summarize_site(PA20_hot, "TX")
#SA20_summary = summarize_site(SA20_hot, "FL")

#summary20 = bind_rows(PA20_summary, SA20_summary)

#print(summary20)

#2021
#PA21_hot = process_buoy_data(PA21_clean)
#SA21_hot = process_buoy_data(SA21_clean)

#PA21_summary = summarize_site(PA21_hot, "TX")
#SA21_summary = summarize_site(SA21_hot, "FL")

#summary21 = bind_rows(PA21_summary, SA21_summary)

#print(summary21)


#2022
#PA22_hot = process_buoy_data(PA22_clean)
#SA22_hot = process_buoy_data(SA22_clean)
#DI22_hot = process_buoy_data(DI22_clean)

#PA22_summary = summarize_site(PA22_hot, "TX")
#SA22_summary = summarize_site(SA22_hot, "FL")
#DI22_summary = summarize_site(DI22_hot, "AL")

#summary22 = bind_rows(PA22_summary, SA22_summary, DI22_summary)

#print(summary22)

#2023
#PA23_hot = process_buoy_data(PA23_clean)
#SA23_hot = process_buoy_data(SA23_clean)
#DI23_hot = process_buoy_data(DI23_clean)

#PA23_summary = summarize_site(PA23_hot, "TX")
#SA23_summary = summarize_site(SA23_hot, "FL")
#DI23_summary = summarize_site(DI23_hot, "AL")

#summary23 = bind_rows(PA23_summary, SA23_summary, DI23_summary)

#print(summary23)

#summary_all <- bind_rows(
#  PA19_summary, SA19_summary,
#  PA20_summary, SA20_summary,
#  PA21_summary, SA21_summary,
#  PA22_summary, SA22_summary, DI22_summary,
#  PA23_summary, SA23_summary, DI23_summary
#)

#days_model <- lm(days_above_30 ~ site, data = summary_all)
#summary(days_model)

#hours_model <- lm(hours_above_30 ~ site, data = summary_all)
#summary(hours_model)

#ggplot(summary_all, aes(x = year, y = days_above_30, color = site)) +
 # geom_line() +
 # geom_point(size = 3) +
 # theme_minimal() +
 # labs(title = "Days Above 30°C by Site", y = "Days Above 30°C", x = "Year")

#summary_recent <- bind_rows(
#  PA21_summary, SA21_summary,
#  PA22_summary, SA22_summary, DI22_summary,
#  PA23_summary, SA23_summary, DI23_summary
#)

#days_model_recent <- lm(days_above_30 ~ site, data = summary_recent)
#summary(days_model_recent)

#hours_model_recent <- lm(hours_above_30 ~ site, data = summary_recent)
#summary(hours_model_recent)

