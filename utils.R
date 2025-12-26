# =============================================================================
# UTILITY FUNCTIONS FOR RENEWABLE ENERGY CLASSIFIER
# =============================================================================
# Author: Sahibjeet Pal Singh
# Date: December 2025
# Description: Helper functions and utilities for data processing, validation,
#              and common operations used throughout the classifier
# =============================================================================

# =============================================================================
# INPUT VALIDATION UTILITIES
# =============================================================================

#' Validate numeric input within specified bounds
#' @param value Numeric value to validate
#' @param min_val Minimum allowed value
#' @param max_val Maximum allowed value
#' @param default Default value if NA or out of bounds
#' @param name Name of the parameter (for error messages)
#' @return Validated numeric value
validate_numeric <- function(value, min_val = -Inf, max_val = Inf, 
                              default = NA, name = "value") {
  # Handle NA
  if (is.na(value) || is.null(value)) {
    if (!is.na(default)) {
      return(default)
    }
    warning(paste(name, "is NA, using NA"))
    return(NA_real_)
  }
  
  # Check if numeric
  if (!is.numeric(value)) {
    value <- tryCatch(
      as.numeric(value),
      error = function(e) {
        warning(paste(name, "could not be converted to numeric"))
        return(default)
      }
    )
  }
  
  # Bound to range
  if (value < min_val) {
    value <- min_val
  } else if (value > max_val) {
    value <- max_val
  }
  
  return(value)
}

#' Validate temperature value
#' @param temp Temperature in Celsius
#' @return Validated temperature
validate_temperature <- function(temp) {
  validate_numeric(temp, min_val = -60, max_val = 60, default = 20, name = "temperature")
}

#' Validate wind speed value
#' @param wind_speed Wind speed in m/s
#' @return Validated wind speed
validate_wind_speed <- function(wind_speed) {
  validate_numeric(wind_speed, min_val = 0, max_val = 50, default = 3, name = "wind_speed")
}

#' Validate precipitation value
#' @param precip Precipitation in mm
#' @return Validated precipitation
validate_precipitation <- function(precip) {
  validate_numeric(precip, min_val = 0, max_val = 1000, default = 75, name = "precipitation")
}

#' Validate elevation value
#' @param elevation Elevation in meters
#' @return Validated elevation
validate_elevation <- function(elevation) {
  validate_numeric(elevation, min_val = -500, max_val = 9000, default = 500, name = "elevation")
}

#' Validate snow depth value
#' @param snow_depth Snow depth in cm
#' @return Validated snow depth
validate_snow_depth <- function(snow_depth) {
  validate_numeric(snow_depth, min_val = 0, max_val = 500, default = 0, name = "snow_depth")
}

#' Validate geographic coordinates
#' @param latitude Latitude value
#' @param longitude Longitude value
#' @return List with validated lat/lon
validate_coordinates <- function(latitude, longitude) {
  lat <- validate_numeric(latitude, min_val = -90, max_val = 90, 
                          default = 0, name = "latitude")
  lon <- validate_numeric(longitude, min_val = -180, max_val = 180, 
                          default = 0, name = "longitude")
  return(list(latitude = lat, longitude = lon))
}

# =============================================================================
# DATA TRANSFORMATION UTILITIES
# =============================================================================

#' Convert GHCN temperature values (tenths of degrees) to Celsius
#' @param temp_tenths Temperature in tenths of degrees Celsius
#' @return Temperature in Celsius
ghcn_temp_to_celsius <- function(temp_tenths) {
  if (is.na(temp_tenths)) return(NA_real_)
  return(temp_tenths / 10)
}

#' Convert GHCN precipitation values (tenths of mm) to mm
#' @param precip_tenths Precipitation in tenths of mm
#' @return Precipitation in mm
ghcn_precip_to_mm <- function(precip_tenths) {
  if (is.na(precip_tenths)) return(NA_real_)
  return(precip_tenths / 10)
}

#' Convert GHCN wind speed values (tenths of m/s) to m/s
#' @param wind_tenths Wind speed in tenths of m/s
#' @return Wind speed in m/s
ghcn_wind_to_ms <- function(wind_tenths) {
  if (is.na(wind_tenths)) return(NA_real_)
  return(wind_tenths / 10)
}

#' Convert Fahrenheit to Celsius
#' @param fahrenheit Temperature in Fahrenheit
#' @return Temperature in Celsius
fahrenheit_to_celsius <- function(fahrenheit) {
  if (is.na(fahrenheit)) return(NA_real_)
  return((fahrenheit - 32) * 5 / 9)
}

#' Convert Celsius to Fahrenheit
#' @param celsius Temperature in Celsius
#' @return Temperature in Fahrenheit
celsius_to_fahrenheit <- function(celsius) {
  if (is.na(celsius)) return(NA_real_)
  return(celsius * 9 / 5 + 32)
}

#' Convert inches to mm
#' @param inches Measurement in inches
#' @return Measurement in mm
inches_to_mm <- function(inches) {
  if (is.na(inches)) return(NA_real_)
  return(inches * 25.4)
}

#' Convert mm to inches
#' @param mm Measurement in mm
#' @return Measurement in inches
mm_to_inches <- function(mm) {
  if (is.na(mm)) return(NA_real_)
  return(mm / 25.4)
}

#' Convert feet to meters
#' @param feet Measurement in feet
#' @return Measurement in meters
feet_to_meters <- function(feet) {
  if (is.na(feet)) return(NA_real_)
  return(feet * 0.3048)
}

#' Convert meters to feet
#' @param meters Measurement in meters
#' @return Measurement in feet
meters_to_feet <- function(meters) {
  if (is.na(meters)) return(NA_real_)
  return(meters / 0.3048)
}

#' Convert km/h to m/s
#' @param kmh Speed in km/h
#' @return Speed in m/s
kmh_to_ms <- function(kmh) {
  if (is.na(kmh)) return(NA_real_)
  return(kmh / 3.6)
}

#' Convert m/s to km/h
#' @param ms Speed in m/s
#' @return Speed in km/h
ms_to_kmh <- function(ms) {
  if (is.na(ms)) return(NA_real_)
  return(ms * 3.6)
}

#' Convert mph to m/s
#' @param mph Speed in mph
#' @return Speed in m/s
mph_to_ms <- function(mph) {
  if (is.na(mph)) return(NA_real_)
  return(mph * 0.44704)
}

#' Convert m/s to mph
#' @param ms Speed in m/s
#' @return Speed in mph
ms_to_mph <- function(ms) {
  if (is.na(ms)) return(NA_real_)
  return(ms / 0.44704)
}

# =============================================================================
# CLIMATE ZONE UTILITIES
# =============================================================================

#' Determine climate zone based on latitude
#' @param latitude Geographic latitude
#' @return Climate zone classification
get_climate_zone <- function(latitude) {
  abs_lat <- abs(latitude)
  
  if (is.na(abs_lat)) return("Unknown")
  
  if (abs_lat < 10) {
    return("Equatorial")
  } else if (abs_lat < 23.5) {
    return("Tropical")
  } else if (abs_lat < 35) {
    return("Subtropical")
  } else if (abs_lat < 55) {
    return("Temperate")
  } else if (abs_lat < 66.5) {
    return("Subarctic/Subantarctic")
  } else {
    return("Polar")
  }
}

#' Determine hemisphere based on latitude
#' @param latitude Geographic latitude
#' @return Hemisphere (Northern/Southern)
get_hemisphere <- function(latitude) {
  if (is.na(latitude)) return("Unknown")
  return(ifelse(latitude >= 0, "Northern", "Southern"))
}

#' Get season based on month and hemisphere
#' @param month Month number (1-12)
#' @param hemisphere "Northern" or "Southern"
#' @return Season name
get_season <- function(month, hemisphere = "Northern") {
  if (is.na(month) || month < 1 || month > 12) return("Unknown")
  
  # Northern hemisphere seasons
  north_seasons <- c(
    "Winter", "Winter",   # Jan, Feb
    "Spring", "Spring", "Spring",  # Mar, Apr, May
    "Summer", "Summer", "Summer",  # Jun, Jul, Aug
    "Fall", "Fall", "Fall",  # Sep, Oct, Nov
    "Winter"  # Dec
  )
  
  season <- north_seasons[month]
  
  # Flip for southern hemisphere
  if (hemisphere == "Southern") {
    season <- switch(season,
                     "Winter" = "Summer",
                     "Spring" = "Fall",
                     "Summer" = "Winter",
                     "Fall" = "Spring")
  }
  
  return(season)
}

#' Estimate daylight hours based on latitude and day of year
#' @param latitude Geographic latitude
#' @param day_of_year Day of year (1-365)
#' @return Estimated daylight hours
estimate_daylight_hours <- function(latitude, day_of_year) {
  if (is.na(latitude) || is.na(day_of_year)) return(12)
  
  # Convert to radians
  lat_rad <- latitude * pi / 180
  
  # Calculate declination angle
  declination <- 23.45 * sin(2 * pi * (284 + day_of_year) / 365) * pi / 180
  
  # Calculate hour angle
  cos_hour_angle <- -tan(lat_rad) * tan(declination)
  
  # Handle polar day/night
  if (cos_hour_angle > 1) return(0)   # Polar night
  if (cos_hour_angle < -1) return(24)  # Midnight sun
  
  hour_angle <- acos(cos_hour_angle)
  daylight_hours <- 2 * hour_angle * 180 / (15 * pi)
  
  return(round(daylight_hours, 1))
}

# =============================================================================
# STATISTICAL UTILITIES
# =============================================================================

#' Calculate confidence interval
#' @param x Numeric vector
#' @param confidence Confidence level (default 0.95)
#' @return Named vector with lower, mean, upper bounds
calculate_ci <- function(x, confidence = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(c(lower = NA, mean = NA, upper = NA))
  
  mean_x <- mean(x)
  se <- sd(x) / sqrt(length(x))
  alpha <- 1 - confidence
  t_value <- qt(1 - alpha / 2, df = length(x) - 1)
  margin <- t_value * se
  
  return(c(
    lower = mean_x - margin,
    mean = mean_x,
    upper = mean_x + margin
  ))
}

#' Calculate coefficient of variation
#' @param x Numeric vector
#' @return Coefficient of variation (as percentage)
calculate_cv <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0 || mean(x) == 0) return(NA_real_)
  return(sd(x) / mean(x) * 100)
}

#' Normalize values to 0-1 range
#' @param x Numeric vector
#' @return Normalized vector
normalize <- function(x) {
  x <- as.numeric(x)
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)
  
  if (max_x == min_x) return(rep(0.5, length(x)))
  
  return((x - min_x) / (max_x - min_x))
}

#' Standardize values (z-score)
#' @param x Numeric vector
#' @return Standardized vector
standardize <- function(x) {
  x <- as.numeric(x)
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  
  if (sd_x == 0) return(rep(0, length(x)))
  
  return((x - mean_x) / sd_x)
}

#' Calculate mode (most frequent value)
#' @param x Vector of values
#' @return Mode value
calculate_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  
  uniq <- unique(x)
  uniq[which.max(tabulate(match(x, uniq)))]
}

#' Calculate weighted mean
#' @param x Numeric vector of values
#' @param weights Numeric vector of weights
#' @return Weighted mean
weighted_mean <- function(x, weights) {
  if (length(x) != length(weights)) {
    stop("x and weights must have same length")
  }
  
  valid <- !is.na(x) & !is.na(weights)
  if (sum(valid) == 0) return(NA_real_)
  
  return(sum(x[valid] * weights[valid]) / sum(weights[valid]))
}

# =============================================================================
# STRING UTILITIES
# =============================================================================

#' Create formatted percentage string
#' @param value Numeric value (0-1 or 0-100)
#' @param is_proportion Whether value is proportion (0-1) or percentage (0-100)
#' @param digits Number of decimal places
#' @return Formatted percentage string
format_percentage <- function(value, is_proportion = TRUE, digits = 1) {
  if (is.na(value)) return("N/A")
  
  if (is_proportion) {
    value <- value * 100
  }
  
  return(paste0(round(value, digits), "%"))
}

#' Create formatted score string
#' @param score Current score
#' @param max_score Maximum possible score
#' @return Formatted score string
format_score <- function(score, max_score) {
  if (is.na(score)) return("N/A")
  
  percentage <- round(score / max_score * 100, 1)
  return(sprintf("%d/%d (%.1f%%)", score, max_score, percentage))
}

#' Clean and standardize location names
#' @param name Location name string
#' @return Cleaned name
clean_location_name <- function(name) {
  if (is.na(name) || is.null(name)) return("Unknown")
  
  name <- as.character(name)
  name <- trimws(name)
  name <- gsub("\\s+", " ", name)  # Remove extra spaces
  name <- tools::toTitleCase(tolower(name))
  
  return(name)
}

#' Generate unique identifier
#' @param prefix Optional prefix for the ID
#' @return Unique identifier string
generate_id <- function(prefix = "LOC") {
  timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
  random <- sprintf("%04d", sample(1:9999, 1))
  return(paste0(prefix, "_", timestamp, "_", random))
}

# =============================================================================
# DATE/TIME UTILITIES
# =============================================================================

#' Get day of year from date
#' @param date Date object or string
#' @return Day of year (1-365/366)
get_day_of_year <- function(date) {
  if (is.character(date)) {
    date <- as.Date(date)
  }
  return(as.numeric(format(date, "%j")))
}

#' Check if year is leap year
#' @param year Year as numeric
#' @return Logical indicating leap year
is_leap_year <- function(year) {
  return((year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0))
}

#' Get number of days in month
#' @param month Month number (1-12)
#' @param year Year (for leap year calculation)
#' @return Number of days
days_in_month <- function(month, year = 2024) {
  days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  if (month == 2 && is_leap_year(year)) {
    return(29)
  }
  
  return(days[month])
}

# =============================================================================
# FILE UTILITIES
# =============================================================================

#' Check if file exists and is readable
#' @param file_path Path to file
#' @return Logical indicating file is valid
is_valid_file <- function(file_path) {
  return(file.exists(file_path) && file.access(file_path, mode = 4) == 0)
}

#' Create directory if it doesn't exist
#' @param dir_path Path to directory
#' @return Logical indicating success
ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created directory: ", dir_path)
  }
  return(dir.exists(dir_path))
}

#' Get file extension
#' @param file_path Path to file
#' @return File extension (lowercase, without dot)
get_file_extension <- function(file_path) {
  ext <- tools::file_ext(file_path)
  return(tolower(ext))
}

#' Generate timestamped filename
#' @param base_name Base name for the file
#' @param extension File extension
#' @return Timestamped filename
timestamped_filename <- function(base_name, extension = "csv") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  return(paste0(base_name, "_", timestamp, ".", extension))
}

# =============================================================================
# LOGGING UTILITIES
# =============================================================================

#' Log message with timestamp
#' @param message Message to log
#' @param level Log level (INFO, WARNING, ERROR)
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
}

#' Log info message
#' @param message Message to log
log_info <- function(message) {
  log_message(message, "INFO")
}

#' Log warning message
#' @param message Message to log
log_warning <- function(message) {
  log_message(message, "WARNING")
}
  
#' Log error message
#' @param message Message to log
log_error <- function(message) {
  log_message(message, "ERROR")
}

# =============================================================================
# PRINT VERSION INFO
# =============================================================================

cat("📦 Utility Functions v1.0.0 loaded\n")
cat("   ", length(ls()), " functions available\n\n")
