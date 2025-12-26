# Renewable Energy Classifier
# Classifies optimal renewable energy sources based on weather/climate data
# Uses GHCN weather station data to recommend Solar, Wind, or Hydro power
# Load packages
load_packages <- function() {
  required_packages <- c(
    "jsonlite",   # JSON data handling and API integration
    "ggplot2",    # Advanced data visualization
    "tidyr",      # Data tidying and reshaping
    "dplyr",      # Data manipulation and transformation
    "readr",      # Fast CSV reading and writing
    "purrr",      # Functional programming tools
    "scales",     # Scale functions for visualization
    "viridis",    # Color palettes for plots
    "sf",         # Spatial data handling
    "leaflet"     # Interactive mapping
  )
  # Check and install missing packages
  missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, repos = "https://cran.r-project.org")
  }
  # Load packages
  suppressPackageStartupMessages({
    library(jsonlite)
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    library(readr)
    library(purrr)
  })
  message("âœ“ All packages loaded successfully")
}
# Initialize packages
load_packages()
# GLOBAL CONFIGURATION
CONFIG <- list(
  # Temperature thresholds (Celsius)
  temp = list(
    solar_optimal_min = 15,
    solar_optimal_max = 25,
    min_viable = -20,
    max_viable = 50
  ),
  # Wind speed thresholds (m/s)
  wind = list(
    min_viable = 4,
    optimal = 6,
    max_viable = 20
  ),
  # Precipitation thresholds (mm)
  precipitation = list(
    solar_low = 50,
    solar_penalty = 150,
    hydro_min = 100,
    hydro_optimal = 150,
    wind_max = 100
  ),
  # Elevation thresholds (meters)
  elevation = list(
    solar_max = 2000,
    wind_min = 500,
    wind_max = 2000,
    hydro_min = 300,
    hydro_max = 2000
  ),
  # Snow depth thresholds (cm)
  snow = list(
    solar_max = 10,
    hydro_min = 20
  ),
  # Maximum possible scores for normalization
  max_scores = list(
    solar = 9,
    wind = 7,
    hydro = 8
  )
)
# DATA LOADING AND PREPROCESSING
load_weather_data <- function(file_path = "cleaned_data.csv", verbose = TRUE) {
  if (verbose) message("Loading weather data from: ", file_path)
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("Error: Data file not found at ", file_path)
  }
  # Read CSV with appropriate column types
  df <- read_csv(
    file_path,
    col_types = cols(
      .default = col_guess()
    ),
    show_col_types = FALSE
  )
  if (verbose) {
    message("âœ“ Loaded ", nrow(df), " rows and ", ncol(df), " columns")
    message("  Columns: ", paste(colnames(df), collapse = ", "))
  }
  return(df)
}
preprocess_data <- function(df) {
  message("Preprocessing weather data...")
  # Select and pivot relevant columns
  df_processed <- df |>
    select(any_of(c("date", "name", "elevation", "latitude", "longitude", "value", "observation"))) |>
    pivot_wider(
      names_from = observation,
      values_from = value
    ) |>
    # Remove rows with missing critical values
    filter(!is.na(name) & !is.na(latitude) & !is.na(longitude))
  # Handle list columns (from pivot_wider)
  df_processed <- df_processed |>
    mutate(
      across(
        where(is.list),
        ~ map_dbl(
          .x,
          ~ if (length(.x) == 0 || is.null(.x)) NA_real_ else as.numeric(.x[1])
        )
      )
    )
  # Remove rows with all NA observations
  df_processed <- df_processed |>
    filter(rowSums(!is.na(select(., -any_of(c("date", "name", "latitude", "longitude", "elevation"))))) > 0)
  message("âœ“ Preprocessing complete: ", nrow(df_processed), " valid locations")
  return(df_processed)
}
validate_numeric_data <- function(df) {
  message("Validating numeric data ranges...")
  # Define valid ranges for each variable
  valid_ranges <- list(
    TAVG = c(-60, 60),      # Temperature in Celsius (tenths)
    TMAX = c(-60, 60),
    TMIN = c(-60, 60),
    PRCP = c(0, 1000),      # Precipitation in mm (tenths)
    AWND = c(0, 50),        # Wind speed in m/s (tenths)
    SNWD = c(0, 500),       # Snow depth in cm
    elevation = c(-500, 9000)  # Elevation in meters
  )
  # Apply range validation
  for (col_name in names(valid_ranges)) {
    if (col_name %in% colnames(df)) {
      range_vals <- valid_ranges[[col_name]]
      df <- df |>
        mutate(
          !!col_name := ifelse(
            .data[[col_name]] < range_vals[1] | .data[[col_name]] > range_vals[2],
            NA_real_,
            .data[[col_name]]
          )
        )
    }
  }
  message("âœ“ Data validation complete")
  return(df)
}
# SOLAR ENERGY SCORING
calculate_solar_score <- function(temp_avg, precipitation, snow_depth = 0, 
                                   elevation = 0, cloud_cover = NULL) {
  score <- 0
  breakdown <- list()
  # Temperature scoring (max 5 points)
  # Higher temperatures correlate with more sunshine hours
  if (!is.na(temp_avg)) {
    if (temp_avg > CONFIG$temp$solar_optimal_min) {
      score <- score + 3
      breakdown$temp_base <- "+3 (temp > 15Â°C: good solar irradiance)"
    }
    if (temp_avg > CONFIG$temp$solar_optimal_max) {
      score <- score + 2
      breakdown$temp_bonus <- "+2 (temp > 25Â°C: optimal conditions)"
    }
  }
  # Precipitation scoring (max 2 points, potential -3 penalty)
  # Low precipitation indicates clear skies
  if (!is.na(precipitation)) {
    if (precipitation < CONFIG$precipitation$solar_low) {
      score <- score + 2
      breakdown$precip_low <- "+2 (precipitation < 50mm: clear skies likely)"
    }
    if (precipitation > CONFIG$precipitation$solar_penalty) {
      score <- score - 3
      breakdown$precip_penalty <- "-3 (precipitation > 150mm: cloudy conditions)"
    }
  }
  # Snow depth scoring (max 1 point)
  # Less snow means better panel efficiency
  if (!is.na(snow_depth)) {
    if (snow_depth < CONFIG$snow$solar_max) {
      score <- score + 1
      breakdown$snow <- "+1 (snow < 10cm: panels remain efficient)"
    }
  }
  # Elevation scoring (max 1 point)
  # Moderate elevation for accessibility and maintenance
  if (!is.na(elevation)) {
    if (elevation < CONFIG$elevation$solar_max) {
      score <- score + 1
      breakdown$elevation <- "+1 (elevation < 2000m: accessible installation)"
    }
  }
  # Optional cloud cover adjustment
  if (!is.null(cloud_cover) && !is.na(cloud_cover)) {
    cloud_adjustment <- round((100 - cloud_cover) / 25)  # 0-4 points
    score <- score + cloud_adjustment
    breakdown$cloud <- paste0("+", cloud_adjustment, " (cloud cover: ", cloud_cover, "%)")
  }
  return(list(
    score = max(score, 0),  # Ensure non-negative
    max_score = CONFIG$max_scores$solar,
    percentage = round(max(score, 0) / CONFIG$max_scores$solar * 100, 1),
    breakdown = breakdown
  ))
}
# WIND ENERGY SCORING
calculate_wind_score <- function(wind_speed, elevation = 0, precipitation = 0,
                                  terrain_roughness = NULL) {
  score <- 0
  breakdown <- list()
  # Wind speed scoring (max 5 points)
  # Crucial factor for wind energy viability
  if (!is.na(wind_speed)) {
    if (wind_speed >= CONFIG$wind$min_viable) {
      score <- score + 2
      breakdown$wind_base <- "+2 (wind â‰¥ 4 m/s: minimum viable)"
    }
    if (wind_speed >= CONFIG$wind$optimal) {
      score <- score + 3
      breakdown$wind_optimal <- "+3 (wind â‰¥ 6 m/s: optimal generation)"
    }
  }
  # Elevation scoring (max 1 point)
  # Higher plateaus often have consistent wind patterns
  if (!is.na(elevation)) {
    if (elevation > CONFIG$elevation$wind_min && elevation < CONFIG$elevation$wind_max) {
      score <- score + 1
      breakdown$elevation <- "+1 (elevation 500-2000m: consistent wind corridor)"
    }
  }
  # Precipitation scoring (max 1 point)
  # Lower precipitation often indicates stable weather patterns
  if (!is.na(precipitation)) {
    if (precipitation < CONFIG$precipitation$wind_max) {
      score <- score + 1
      breakdown$precipitation <- "+1 (precipitation < 100mm: stable conditions)"
    }
  }
  # Optional terrain roughness adjustment
  if (!is.null(terrain_roughness) && !is.na(terrain_roughness)) {
    # Low roughness (open terrain) is better for wind
    roughness_bonus <- round((1 - terrain_roughness) * 2)
    score <- score + roughness_bonus
    breakdown$terrain <- paste0("+", roughness_bonus, " (terrain roughness: ", 
                                 round(terrain_roughness * 100), "%)")
  }
  return(list(
    score = max(score, 0),
    max_score = CONFIG$max_scores$wind,
    percentage = round(max(score, 0) / CONFIG$max_scores$wind * 100, 1),
    breakdown = breakdown
  ))
}
# HYDRO ENERGY SCORING
calculate_hydro_score <- function(precipitation, snow_depth = 0, elevation = 0,
                                   river_proximity = NULL) {
  score <- 0
  breakdown <- list()
  # Precipitation scoring (max 4 points)
  # Water availability is crucial for hydro power
  if (!is.na(precipitation)) {
    if (precipitation > CONFIG$precipitation$hydro_min) {
      score <- score + 2
      breakdown$precip_base <- "+2 (precipitation > 100mm: adequate water supply)"
    }
    if (precipitation > CONFIG$precipitation$hydro_optimal) {
      score <- score + 2
      breakdown$precip_bonus <- "+2 (precipitation > 150mm: excellent water supply)"
    }
  }
  # Snow depth scoring (max 2 points)
  # Snowmelt provides seasonal water flow
  if (!is.na(snow_depth)) {
    if (snow_depth > CONFIG$snow$hydro_min) {
      score <- score + 2
      breakdown$snow <- "+2 (snow > 20cm: snowmelt contribution)"
    }
  }
  # Elevation scoring (max 2 points)
  # Elevation difference is key for hydroelectric generation
  if (!is.na(elevation)) {
    if (elevation > CONFIG$elevation$hydro_min && elevation < CONFIG$elevation$hydro_max) {
      score <- score + 2
      breakdown$elevation <- "+2 (elevation 300-2000m: suitable head height)"
    }
  }
  # Optional river proximity adjustment
  if (!is.null(river_proximity) && !is.na(river_proximity)) {
    if (river_proximity < 10) {
      proximity_bonus <- 2
    } else if (river_proximity < 50) {
      proximity_bonus <- 1
    } else {
      proximity_bonus <- 0
    }
    score <- score + proximity_bonus
    breakdown$river <- paste0("+", proximity_bonus, " (river proximity: ", 
                               river_proximity, " km)")
  }
  return(list(
    score = max(score, 0),
    max_score = CONFIG$max_scores$hydro,
    percentage = round(max(score, 0) / CONFIG$max_scores$hydro * 100, 1),
    breakdown = breakdown
  ))
}
# MAIN CLASSIFICATION FUNCTION
classify_energy <- function(temp_avg, wind_speed, precipitation, elevation, 
                            snow_depth = 0, return_details = FALSE) {
  # Input validation and bounding
  temp_avg <- if (is.na(temp_avg)) 20 else max(min(temp_avg, CONFIG$temp$max_viable), CONFIG$temp$min_viable)
  wind_speed <- if (is.na(wind_speed)) 3 else max(min(wind_speed, CONFIG$wind$max_viable), 0)
  precipitation <- if (is.na(precipitation)) 75 else max(precipitation, 0)
  elevation <- if (is.na(elevation)) 500 else max(elevation, 0)
  snow_depth <- if (is.na(snow_depth)) 0 else max(snow_depth, 0)
  # Calculate individual scores
  solar_result <- calculate_solar_score(temp_avg, precipitation, snow_depth, elevation)
  wind_result <- calculate_wind_score(wind_speed, elevation, precipitation)
  hydro_result <- calculate_hydro_score(precipitation, snow_depth, elevation)
  # Determine best resource
  scores <- c(
    Solar = solar_result$score,
    Wind = wind_result$score,
    Hydro = hydro_result$score
  )
  # Handle ties by prioritizing based on sustainability factors
  max_score <- max(scores)
  best_options <- names(scores[scores == max_score])
  # Priority order for ties: Hydro > Wind > Solar (based on reliability)
  priority <- c("Hydro", "Wind", "Solar")
  best_resource <- priority[priority %in% best_options][1]
  # Calculate confidence based on score difference
  sorted_scores <- sort(scores, decreasing = TRUE)
  if (length(sorted_scores) >= 2) {
    confidence <- round((sorted_scores[1] - sorted_scores[2]) / sorted_scores[1] * 100, 1)
  } else {
    confidence <- 100
  }
  # Build result
  result <- list(
    # Raw scores
    solar_score = solar_result$score,
    wind_score = wind_result$score,
    hydro_score = hydro_result$score,
    # Percentages
    solar_pct = solar_result$percentage,
    wind_pct = wind_result$percentage,
    hydro_pct = hydro_result$percentage,
    # Recommendation
    best_resource = best_resource,
    confidence = confidence,
    # All scores vector
    all_scores = scores
  )
  # Add detailed breakdown if requested
  if (return_details) {
    result$details <- list(
      solar = solar_result,
      wind = wind_result,
      hydro = hydro_result,
      input_params = list(
        temp_avg = temp_avg,
        wind_speed = wind_speed,
        precipitation = precipitation,
        elevation = elevation,
        snow_depth = snow_depth
      )
    )
  }
  return(result)
}
classify_dataframe <- function(df, temp_col = "TAVG", wind_col = "AWND",
                                precip_col = "PRCP", elev_col = "elevation",
                                snow_col = "SNWD") {
  message("Classifying ", nrow(df), " locations...")
  # Apply classification to each row
  results <- pmap(
    list(
      temp = if (temp_col %in% names(df)) df[[temp_col]] / 10 else rep(NA, nrow(df)),
      wind = if (wind_col %in% names(df)) df[[wind_col]] / 10 else rep(NA, nrow(df)),
      precip = if (precip_col %in% names(df)) df[[precip_col]] / 10 else rep(NA, nrow(df)),
      elev = if (elev_col %in% names(df)) df[[elev_col]] else rep(NA, nrow(df)),
      snow = if (snow_col %in% names(df)) df[[snow_col]] / 10 else rep(0, nrow(df))
    ),
    function(temp, wind, precip, elev, snow) {
      classify_energy(temp, wind, precip, elev, snow)
    }
  )
  # Extract results into columns
  df$solar_score <- sapply(results, function(x) x$solar_score)
  df$wind_score <- sapply(results, function(x) x$wind_score)
  df$hydro_score <- sapply(results, function(x) x$hydro_score)
  df$best_resource <- sapply(results, function(x) x$best_resource)
  df$confidence <- sapply(results, function(x) x$confidence)
  message("âœ“ Classification complete")
  return(df)
}
# DATA ANALYSIS AND STATISTICS
generate_summary_stats <- function(df) {
  if (!"best_resource" %in% names(df)) {
    stop("Dataframe must be classified first. Use classify_dataframe()")
  }
  # Overall distribution
  resource_distribution <- df |>
    count(best_resource) |>
    mutate(percentage = round(n / sum(n) * 100, 1))
  # Score statistics
  score_stats <- df |>
    summarise(
      solar_mean = mean(solar_score, na.rm = TRUE),
      solar_sd = sd(solar_score, na.rm = TRUE),
      wind_mean = mean(wind_score, na.rm = TRUE),
      wind_sd = sd(wind_score, na.rm = TRUE),
      hydro_mean = mean(hydro_score, na.rm = TRUE),
      hydro_sd = sd(hydro_score, na.rm = TRUE),
      avg_confidence = mean(confidence, na.rm = TRUE)
    )
  # Geographic distribution
  if ("latitude" %in% names(df)) {
    geo_stats <- df |>
      mutate(
        hemisphere = ifelse(latitude >= 0, "Northern", "Southern"),
        climate_zone = case_when(
          abs(latitude) < 23.5 ~ "Tropical",
          abs(latitude) < 35 ~ "Subtropical",
          abs(latitude) < 55 ~ "Temperate",
          TRUE ~ "Polar"
        )
      ) |>
      group_by(climate_zone) |>
      count(best_resource) |>
      pivot_wider(names_from = best_resource, values_from = n, values_fill = 0)
  } else {
    geo_stats <- NULL
  }
  return(list(
    distribution = resource_distribution,
    scores = score_stats,
    geographic = geo_stats,
    total_locations = nrow(df)
  ))
}
print_summary_report <- function(stats) {
  cat("\n")
  cat("â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n")
  cat("           RENEWABLE ENERGY CLASSIFICATION REPORT              \n")
  cat("â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n")
  cat("ðŸ“Š OVERALL DISTRIBUTION\n")
  cat("â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€\n")
  for (i in 1:nrow(stats$distribution)) {
    resource <- stats$distribution$best_resource[i]
    count <- stats$distribution$n[i]
    pct <- stats$distribution$percentage[i]
    icon <- switch(resource, "Solar" = "â˜€ï¸", "Wind" = "ðŸ’¨", "Hydro" = "ðŸ’§", "â“")
    cat(sprintf("  %s %-6s: %5d locations (%5.1f%%)\n", icon, resource, count, pct))
  }
  cat("\nðŸ“ˆ SCORE STATISTICS\n")
  cat("â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€\n")
  cat(sprintf("  â˜€ï¸ Solar:  Mean = %.2f, SD = %.2f\n", 
              stats$scores$solar_mean, stats$scores$solar_sd))
  cat(sprintf("  ðŸ’¨ Wind:   Mean = %.2f, SD = %.2f\n", 
              stats$scores$wind_mean, stats$scores$wind_sd))
  cat(sprintf("  ðŸ’§ Hydro:  Mean = %.2f, SD = %.2f\n", 
              stats$scores$hydro_mean, stats$scores$hydro_sd))
  cat(sprintf("  ðŸŽ¯ Average Confidence: %.1f%%\n", stats$scores$avg_confidence))
  if (!is.null(stats$geographic)) {
    cat("\nðŸŒ GEOGRAPHIC DISTRIBUTION\n")
    cat("â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€\n")
    print(stats$geographic, n = Inf)
  }
  cat("\nâ•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n")
  cat(sprintf("  Total Locations Analyzed: %d\n", stats$total_locations))
  cat("â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n")
}
# VISUALIZATION FUNCTIONS
plot_distribution <- function(df) {
  if (!"best_resource" %in% names(df)) {
    stop("Dataframe must be classified first")
  }
  # Define colors
  energy_colors <- c(
    "Solar" = "#ea580c",
    "Wind" = "#0891b2",
    "Hydro" = "#2563eb"
  )
  # Create plot
  p <- df |>
    count(best_resource) |>
    mutate(percentage = n / sum(n) * 100) |>
    ggplot(aes(x = reorder(best_resource, -n), y = n, fill = best_resource)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, percentage)), 
              vjust = -0.5, size = 4) +
    scale_fill_manual(values = energy_colors) +
    labs(
      title = "Renewable Energy Classification Distribution",
      subtitle = paste("Based on", nrow(df), "weather stations"),
      x = "Energy Type",
      y = "Number of Locations"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray50"),
      legend.position = "none",
      axis.text = element_text(size = 11),
      panel.grid.major.x = element_blank()
    ) +
    expand_limits(y = max(table(df$best_resource)) * 1.2)
  return(p)
}
plot_score_comparison <- function(df) {
  required_cols <- c("solar_score", "wind_score", "hydro_score")
  if (!all(required_cols %in% names(df))) {
    stop("Dataframe must contain score columns")
  }
  # Reshape data for plotting
  score_data <- df |>
    select(all_of(required_cols)) |>
    pivot_longer(
      cols = everything(),
      names_to = "energy_type",
      values_to = "score"
    ) |>
    mutate(
      energy_type = case_when(
        energy_type == "solar_score" ~ "Solar",
        energy_type == "wind_score" ~ "Wind",
        energy_type == "hydro_score" ~ "Hydro"
      )
    )
  # Define colors
  energy_colors <- c(
    "Solar" = "#ea580c",
    "Wind" = "#0891b2",
    "Hydro" = "#2563eb"
  )
  # Create boxplot
  p <- ggplot(score_data, aes(x = energy_type, y = score, fill = energy_type)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 21) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = energy_colors) +
    labs(
      title = "Score Distribution by Energy Type",
      subtitle = "Boxplot with individual data points",
      x = "Energy Type",
      y = "Score"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray50"),
      legend.position = "none"
    )
  return(p)
}
plot_geographic <- function(df) {
  if (!all(c("latitude", "longitude", "best_resource") %in% names(df))) {
    stop("Dataframe must contain latitude, longitude, and best_resource columns")
  }
  energy_colors <- c(
    "Solar" = "#ea580c",
    "Wind" = "#0891b2",
    "Hydro" = "#2563eb"
  )
  p <- ggplot(df, aes(x = longitude, y = latitude, color = best_resource)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = energy_colors) +
    labs(
      title = "Geographic Distribution of Optimal Energy Sources",
      subtitle = "Each point represents a weather station",
      x = "Longitude",
      y = "Latitude",
      color = "Best Resource"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray50"),
      legend.position = "bottom"
    ) +
    coord_fixed(ratio = 1.3)
  return(p)
}
# EXPORT FUNCTIONS
export_to_csv <- function(df, output_path = "data.csv", include_scores = TRUE) {
  message("Exporting results to: ", output_path)
  # Select columns to export
  if (include_scores) {
    export_cols <- c("name", "longitude", "latitude", "elevation",
                     "solar_score", "wind_score", "hydro_score",
                     "best_resource", "confidence")
  } else {
    export_cols <- c("name", "longitude", "latitude", "best_resource")
  }
  # Filter to available columns
  export_cols <- export_cols[export_cols %in% names(df)]
  # Export
  df |>
    select(all_of(export_cols)) |>
    write_csv(output_path)
  message("âœ“ Exported ", nrow(df), " rows to ", output_path)
}
export_to_json <- function(df, output_path = "data.json") {
  message("Exporting results to JSON: ", output_path)
  # Convert to JSON
  json_data <- toJSON(df, pretty = TRUE, auto_unbox = TRUE)
  # Write to file
  write(json_data, output_path)
  message("âœ“ Exported to ", output_path)
}
# MAIN EXECUTION
main <- function(input_file = "cleaned_data.csv", 
                 output_file = "data.csv",
                 generate_plots = FALSE) {
  cat("\n")
  cat("â•”â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•—\n")
  cat("â•‘       RENEWABLE ENERGY CLASSIFIER - BATCH PROCESSING         â•‘\n")
  cat("â•šâ•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n")
  tryCatch({
    # Step 1: Load data
    df <- load_weather_data(input_file)
    # Step 2: Preprocess
    df <- preprocess_data(df)
    # Step 3: Validate
    df <- validate_numeric_data(df)
    # Step 4: Classify
    df <- classify_dataframe(df)
    # Step 5: Generate summary
    stats <- generate_summary_stats(df)
    print_summary_report(stats)
    # Step 6: Export results
    export_to_csv(df, output_file)
    # Step 7: Generate plots if requested
    if (generate_plots) {
      message("\nGenerating visualization plots...")
      p1 <- plot_distribution(df)
      ggsave("distribution_plot.png", p1, width = 10, height = 6, dpi = 300)
      message("âœ“ Saved distribution_plot.png")
      p2 <- plot_score_comparison(df)
      ggsave("score_comparison.png", p2, width = 10, height = 6, dpi = 300)
      message("âœ“ Saved score_comparison.png")
      if (all(c("latitude", "longitude") %in% names(df))) {
        p3 <- plot_geographic(df)
        ggsave("geographic_map.png", p3, width = 12, height = 8, dpi = 300)
        message("âœ“ Saved geographic_map.png")
      }
    }
    message("\nâœ… Processing complete!")
    return(invisible(df))
  }, error = function(e) {
    message("\nâŒ Error: ", e$message)
    return(NULL)
  })
}
# SCRIPT EXECUTION
# Run main function if script is executed directly
if (interactive()) {
  message("Renewable Energy Classifier loaded. Use main() to run batch processing.")
  message("Or use classify_energy() for individual classifications.")
} else {
  # Attempt to run main if data file exists
  if (file.exists("cleaned_data.csv")) {
    result <- main()
  } else {
    message("No cleaned_data.csv found. Run main('your_data.csv') with your data file.")
  }
}
# Script ready
