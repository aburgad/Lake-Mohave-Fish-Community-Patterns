#--------------------------------------
library(tidyverse) # ggplot and dplyr
library(lubridate) # dates 
#--------------------------------------

# read in csv
df <- read.csv(file = "mohave_roundup_catch_data.csv")

# check structure
head(df)
str(df)
summary(df)

# Function to summarize sampling trips 

sample_trips_summary <- function(df){

# change dates from factor to date format and create months
samp_df <- df %>%
  mutate(date_begin = mdy(date_begin),
         date_end = mdy(date_end),
         month_begin = month(date_begin),
         month_end = month(date_end)) %>%
  select(year, season, 
         date_begin, date_end, 
         net_units, month_begin, 
         month_end)

print(paste("Range of years"))
print(paste(range(samp_df$year)))

# total number of trips 
print(paste("Total number of sampling trips:",nrow(samp_df)))

# total net units 
print(paste("Total net units:", sum(samp_df$net_units)))

# summarise trips by month begin
month_trips <- samp_df %>%
  group_by(month_begin) %>%
  summarise(trips = n(),
            avg_net_units = mean(net_units))

print(paste("Summary of trips by month"))
print(month_trips)

# years that didn't have two trips 
one_trip_year <- df %>%
  group_by(year) %>%
  summarise(trips = n()) %>%
  filter(trips < 2)

print(paste("Years with only one trip"))
print(one_trip_year)

# years that had more than two trips included 
mult_trip_year <- df %>%
  group_by(year) %>%
  summarise(trips = n()) %>%
  filter(trips > 2)

print(paste("Years with multiple trips"))
print(mult_trip_year)

}

sample_trips_summary(df)
