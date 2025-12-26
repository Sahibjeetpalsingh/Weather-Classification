# Renewable Energy Classifier
# This script processes weather data and classifies optimal renewable energy sources

library(jsonlite)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(purrr)

# Load and clean data
df <- read_csv("cleaned_data.csv")

# Show column names
colnames(df)

# Select relevant columns and pivot observations
df <- df |>
  select(c("date", "name", "elevation", "latitude", "longitude", "value", "observation")) |>
  pivot_wider(
    names_from  = observation,
    values_from = value
  ) |>
  na.omit()

# Handle list columns
df <- df |>
  dplyr::mutate(
    dplyr::across(
      where(is.list),
      ~ purrr::map_dbl(
        .x,
        ~ if (length(.x) == 0 || is.null(.x)) NA_real_ else as.numeric(.x[1])
      )
    )
  )

# Classification function
classify_energy <- function(temp_avg, wind_speed, precipitation, elevation, snow_depth = 0) {
  # Bound inputs
  temp_avg <- max(min(temp_avg, 50), -20)
  wind_speed <- max(min(wind_speed, 20), 0)
  precipitation <- max(precipitation, 0)
  
  solar_score <- 0
  wind_score <- 0
  hydro_score <- 0
  
  # Solar Logic
  if (temp_avg > 15) solar_score <- solar_score + 3
  if (temp_avg > 25) solar_score <- solar_score + 2
  if (precipitation < 50) solar_score <- solar_score + 2
  if (precipitation > 150) solar_score <- solar_score - 3  # Penalize for heavy precip
  if (snow_depth < 10) solar_score <- solar_score + 1
  if (elevation < 2000) solar_score <- solar_score + 1
  
  # Wind Logic
  if (wind_speed >= 4) wind_score <- wind_score + 2
  if (wind_speed >= 6) wind_score <- wind_score + 3
  if (elevation > 500 && elevation < 2000) wind_score <- wind_score + 1
  if (precipitation < 100) wind_score <- wind_score + 1
  
  # Hydro Logic
  if (precipitation > 100) hydro_score <- hydro_score + 2
  if (precipitation > 150) hydro_score <- hydro_score + 2
  if (snow_depth > 20) hydro_score <- hydro_score + 2
  if (elevation > 300 && elevation < 2000) hydro_score <- hydro_score + 2
  
  # Determine best resource
  scores <- c(Solar = solar_score, Wind = wind_score, Hydro = hydro_score)
  best_resource <- names(which.max(scores))
  
  return(list(
    solar = solar_score,
    wind = wind_score,
    hydro = hydro_score,
    best = best_resource
  ))
}

# Apply classification to each location
# Note: Adjust column names based on your actual data structure
# df$best_resource <- sapply(1:nrow(df), function(i) {
#   result <- classify_energy(
#     temp_avg = df$TAVG[i] / 10,  # Convert from tenths of degrees
#     wind_speed = df$AWND[i] / 10,
#     precipitation = df$PRCP[i] / 10,
#     elevation = df$elevation[i],
#     snow_depth = df$SNWD[i] / 10
#   )
#   result$best
# })

print(paste("Total rows:", nrow(df)))
