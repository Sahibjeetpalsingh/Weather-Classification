# =============================================================================
# DATA FETCHING AND API INTEGRATION FOR RENEWABLE ENERGY CLASSIFIER
# =============================================================================
# Author: Sahibjeet Pal Singh
# Date: December 2025
# Description: Functions for fetching weather data from NOAA GHCN API,
#              processing responses, and managing data downloads
# =============================================================================

# =============================================================================
# PACKAGE DEPENDENCIES
# =============================================================================

#' Load required packages for API operations
load_api_packages <- function() {
  required <- c("httr", "jsonlite", "curl", "readr", "dplyr", "purrr")
  
  for (pkg in required) {
    if (pkg %in% installed.packages()[, "Package"]) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
}

load_api_packages()

# =============================================================================
# API CONFIGURATION
# =============================================================================

#' NOAA API Configuration
NOAA_CONFIG <- list(
  # Base URL for NOAA CDO API
  base_url = "https://www.ncei.noaa.gov/cdo-web/api/v2/",
  
  # GHCN Daily dataset ID
  dataset_id = "GHCND",
  
  # Data types we're interested in
  datatypes = list(
    TAVG = "Average Temperature",
    TMAX = "Maximum Temperature",
    TMIN = "Minimum Temperature",
    PRCP = "Precipitation",
    AWND = "Average Wind Speed",
    SNWD = "Snow Depth",
    SNOW = "Snowfall"
  ),
  
  # Rate limiting
  requests_per_second = 5,
  max_retries = 3,
  retry_delay = 2,
  
  # Pagination
  max_results_per_request = 1000
)

# =============================================================================
# API TOKEN MANAGEMENT
# =============================================================================

#' Get NOAA API token from environment or config file
#' @param config_file Path to config file
#' @return API token string
get_api_token <- function(config_file = ".noaa_token") {
  # Try environment variable first
  token <- Sys.getenv("NOAA_API_TOKEN")
  
  if (token != "") {
    return(token)
  }
  
  # Try config file
  if (file.exists(config_file)) {
    token <- trimws(readLines(config_file, n = 1, warn = FALSE))
    if (nchar(token) > 0) {
      return(token)
    }
  }
  
  # Prompt for token
  message("NOAA API token not found.")
  message("Get your free token at: https://www.ncdc.noaa.gov/cdo-web/token")
  message("Set it as environment variable NOAA_API_TOKEN or save to .noaa_token file")
  
  return(NULL)
}

#' Validate API token
#' @param token API token to validate
#' @return Logical indicating valid token
validate_token <- function(token) {
  if (is.null(token) || nchar(token) < 10) {
    return(FALSE)
  }
  
  # Test API call
  tryCatch({
    response <- httr::GET(
      url = paste0(NOAA_CONFIG$base_url, "datasets"),
      httr::add_headers(token = token),
      httr::timeout(10)
    )
    
    return(httr::status_code(response) == 200)
  }, error = function(e) {
    return(FALSE)
  })
}

# =============================================================================
# API REQUEST FUNCTIONS
# =============================================================================

#' Make authenticated API request to NOAA
#' @param endpoint API endpoint
#' @param params Query parameters
#' @param token API token
#' @return Parsed JSON response or NULL on error
make_api_request <- function(endpoint, params = list(), token = NULL) {
  if (is.null(token)) {
    token <- get_api_token()
    if (is.null(token)) {
      stop("No valid API token available")
    }
  }
  
  # Build URL
  url <- paste0(NOAA_CONFIG$base_url, endpoint)
  
  # Make request with retries
  for (attempt in 1:NOAA_CONFIG$max_retries) {
    tryCatch({
      response <- httr::GET(
        url = url,
        query = params,
        httr::add_headers(token = token),
        httr::timeout(30)
      )
      
      # Check status
      if (httr::status_code(response) == 200) {
        content <- httr::content(response, as = "text", encoding = "UTF-8")
        return(jsonlite::fromJSON(content))
      } else if (httr::status_code(response) == 429) {
        # Rate limited - wait and retry
        message("Rate limited. Waiting ", NOAA_CONFIG$retry_delay * attempt, " seconds...")
        Sys.sleep(NOAA_CONFIG$retry_delay * attempt)
      } else {
        warning("API request failed with status: ", httr::status_code(response))
        return(NULL)
      }
    }, error = function(e) {
      warning("Request attempt ", attempt, " failed: ", e$message)
      Sys.sleep(NOAA_CONFIG$retry_delay)
    })
  }
  
  warning("All retry attempts failed")
  return(NULL)
}

#' Paginated API request for large datasets
#' @param endpoint API endpoint
#' @param params Base query parameters
#' @param token API token
#' @param max_pages Maximum number of pages to fetch
#' @return Combined results from all pages
make_paginated_request <- function(endpoint, params = list(), token = NULL, max_pages = 10) {
  all_results <- list()
  offset <- 1
  page <- 1
  
  params$limit <- NOAA_CONFIG$max_results_per_request
  
  repeat {
    if (page > max_pages) {
      message("Reached maximum page limit (", max_pages, ")")
      break
    }
    
    params$offset <- offset
    
    message("Fetching page ", page, "...")
    response <- make_api_request(endpoint, params, token)
    
    if (is.null(response) || is.null(response$results) || length(response$results) == 0) {
      break
    }
    
    all_results <- c(all_results, list(response$results))
    
    # Check if more pages available
    if (length(response$results) < NOAA_CONFIG$max_results_per_request) {
      break
    }
    
    offset <- offset + NOAA_CONFIG$max_results_per_request
    page <- page + 1
    
    # Rate limiting pause
    Sys.sleep(1 / NOAA_CONFIG$requests_per_second)
  }
  
  # Combine all results
  if (length(all_results) > 0) {
    return(dplyr::bind_rows(all_results))
  }
  
  return(NULL)
}

# =============================================================================
# STATION DATA FUNCTIONS
# =============================================================================

#' Fetch list of weather stations
#' @param location_id Location ID (country, state, etc.)
#' @param extent Geographic bounding box (minLat, minLon, maxLat, maxLon)
#' @param token API token
#' @return Dataframe of stations
fetch_stations <- function(location_id = NULL, extent = NULL, token = NULL) {
  params <- list(
    datasetid = NOAA_CONFIG$dataset_id,
    limit = 1000
  )
  
  if (!is.null(location_id)) {
    params$locationid <- location_id
  }
  
  if (!is.null(extent)) {
    params$extent <- paste(extent, collapse = ",")
  }
  
  message("Fetching weather stations...")
  response <- make_paginated_request("stations", params, token)
  
  if (!is.null(response)) {
    message("✓ Found ", nrow(response), " stations")
  }
  
  return(response)
}

#' Get station details by ID
#' @param station_id Station identifier
#' @param token API token
#' @return Station details
get_station_details <- function(station_id, token = NULL) {
  endpoint <- paste0("stations/", station_id)
  make_api_request(endpoint, list(), token)
}

#' Search stations by name
#' @param name_pattern Name pattern to search
#' @param stations Dataframe of stations (optional, will fetch if NULL)
#' @return Matching stations
search_stations <- function(name_pattern, stations = NULL) {
  if (is.null(stations)) {
    stations <- fetch_stations()
  }
  
  if (is.null(stations)) {
    return(NULL)
  }
  
  stations %>%
    filter(grepl(name_pattern, name, ignore.case = TRUE))
}

# =============================================================================
# WEATHER DATA FUNCTIONS
# =============================================================================

#' Fetch weather data for a station
#' @param station_id Station identifier
#' @param start_date Start date (YYYY-MM-DD)
#' @param end_date End date (YYYY-MM-DD)
#' @param datatypes Vector of data types to fetch
#' @param token API token
#' @return Dataframe of weather observations
fetch_station_data <- function(station_id, start_date, end_date, 
                                datatypes = NULL, token = NULL) {
  params <- list(
    datasetid = NOAA_CONFIG$dataset_id,
    stationid = station_id,
    startdate = start_date,
    enddate = end_date,
    units = "metric"
  )
  
  if (!is.null(datatypes)) {
    params$datatypeid <- paste(datatypes, collapse = ",")
  }
  
  message("Fetching data for station: ", station_id)
  response <- make_paginated_request("data", params, token)
  
  if (!is.null(response)) {
    message("✓ Retrieved ", nrow(response), " observations")
  }
  
  return(response)
}

#' Fetch weather data for multiple stations
#' @param station_ids Vector of station IDs
#' @param start_date Start date
#' @param end_date End date
#' @param datatypes Data types to fetch
#' @param token API token
#' @param progress Show progress bar
#' @return Combined dataframe of all station data
fetch_multiple_stations <- function(station_ids, start_date, end_date,
                                     datatypes = NULL, token = NULL,
                                     progress = TRUE) {
  all_data <- list()
  n <- length(station_ids)
  
  for (i in seq_along(station_ids)) {
    if (progress) {
      message(sprintf("[%d/%d] Processing station: %s", i, n, station_ids[i]))
    }
    
    data <- fetch_station_data(
      station_id = station_ids[i],
      start_date = start_date,
      end_date = end_date,
      datatypes = datatypes,
      token = token
    )
    
    if (!is.null(data) && nrow(data) > 0) {
      all_data[[station_ids[i]]] <- data
    }
    
    # Rate limiting
    Sys.sleep(0.5)
  }
  
  # Combine all data
  if (length(all_data) > 0) {
    combined <- bind_rows(all_data, .id = "station_id")
    message("✓ Total observations: ", nrow(combined))
    return(combined)
  }
  
  return(NULL)
}

#' Fetch data for a geographic region
#' @param min_lat Minimum latitude
#' @param max_lat Maximum latitude
#' @param min_lon Minimum longitude
#' @param max_lon Maximum longitude
#' @param start_date Start date
#' @param end_date End date
#' @param token API token
#' @return Dataframe of weather data
fetch_region_data <- function(min_lat, max_lat, min_lon, max_lon,
                               start_date, end_date, token = NULL) {
  # First get stations in region
  extent <- c(min_lat, min_lon, max_lat, max_lon)
  stations <- fetch_stations(extent = extent, token = token)
  
  if (is.null(stations) || nrow(stations) == 0) {
    message("No stations found in region")
    return(NULL)
  }
  
  # Limit to first 50 stations to avoid excessive API calls
  if (nrow(stations) > 50) {
    message("Limiting to first 50 stations")
    stations <- stations[1:50, ]
  }
  
  # Fetch data for all stations
  data <- fetch_multiple_stations(
    station_ids = stations$id,
    start_date = start_date,
    end_date = end_date,
    token = token
  )
  
  # Add station metadata
  if (!is.null(data)) {
    data <- data %>%
      left_join(
        stations %>% select(id, name, latitude, longitude, elevation),
        by = c("station_id" = "id")
      )
  }
  
  return(data)
}

# =============================================================================
# DATA PROCESSING FUNCTIONS
# =============================================================================

#' Process raw NOAA data into analysis-ready format
#' @param raw_data Raw dataframe from API
#' @return Processed dataframe
process_noaa_data <- function(raw_data) {
  if (is.null(raw_data) || nrow(raw_data) == 0) {
    return(NULL)
  }
  
  message("Processing raw data...")
  
  # Pivot data types to columns
  processed <- raw_data %>%
    mutate(date = as.Date(date)) %>%
    select(station_id, date, name, latitude, longitude, elevation, datatype, value) %>%
    pivot_wider(
      names_from = datatype,
      values_from = value,
      values_fn = mean  # Average if multiple values per day
    )
  
  # Convert units (GHCN uses tenths for many values)
  if ("TAVG" %in% names(processed)) {
    processed$TAVG <- processed$TAVG / 10  # tenths of Celsius to Celsius
  }
  if ("TMAX" %in% names(processed)) {
    processed$TMAX <- processed$TMAX / 10
  }
  if ("TMIN" %in% names(processed)) {
    processed$TMIN <- processed$TMIN / 10
  }
  if ("PRCP" %in% names(processed)) {
    processed$PRCP <- processed$PRCP / 10  # tenths of mm to mm
  }
  if ("AWND" %in% names(processed)) {
    processed$AWND <- processed$AWND / 10  # tenths of m/s to m/s
  }
  
  message("✓ Processed ", nrow(processed), " records")
  
  return(processed)
}

#' Aggregate daily data to monthly averages
#' @param daily_data Daily observation dataframe
#' @return Monthly aggregated dataframe
aggregate_to_monthly <- function(daily_data) {
  if (is.null(daily_data)) return(NULL)
  
  daily_data %>%
    mutate(
      year = lubridate::year(date),
      month = lubridate::month(date)
    ) %>%
    group_by(station_id, name, latitude, longitude, elevation, year, month) %>%
    summarise(
      across(where(is.numeric) & !c(year, month), 
             ~ mean(.x, na.rm = TRUE)),
      n_observations = n(),
      .groups = "drop"
    )
}

#' Aggregate to station-level climatology
#' @param data Weather dataframe
#' @return Station-level averages
aggregate_to_climatology <- function(data) {
  if (is.null(data)) return(NULL)
  
  data %>%
    group_by(station_id, name, latitude, longitude, elevation) %>%
    summarise(
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      n_observations = n(),
      date_range = paste(min(date), "to", max(date)),
      .groups = "drop"
    )
}

# =============================================================================
# CACHING FUNCTIONS
# =============================================================================

#' Cache directory management
CACHE_DIR <- ".cache"

#' Initialize cache directory
init_cache <- function() {
  if (!dir.exists(CACHE_DIR)) {
    dir.create(CACHE_DIR, recursive = TRUE)
    message("Created cache directory: ", CACHE_DIR)
  }
}

#' Generate cache key from parameters
#' @param ... Parameters to hash
#' @return Cache key string
generate_cache_key <- function(...) {
  params <- list(...)
  hash <- digest::digest(params, algo = "md5")
  return(hash)
}

#' Check if cached data exists and is fresh
#' @param cache_key Cache key
#' @param max_age_days Maximum age in days
#' @return Logical indicating cache validity
is_cache_valid <- function(cache_key, max_age_days = 7) {
  cache_file <- file.path(CACHE_DIR, paste0(cache_key, ".rds"))
  
  if (!file.exists(cache_file)) {
    return(FALSE)
  }
  
  file_age <- difftime(Sys.time(), file.mtime(cache_file), units = "days")
  return(as.numeric(file_age) <= max_age_days)
}

#' Read from cache
#' @param cache_key Cache key
#' @return Cached data or NULL
read_cache <- function(cache_key) {
  cache_file <- file.path(CACHE_DIR, paste0(cache_key, ".rds"))
  
  if (file.exists(cache_file)) {
    message("Reading from cache...")
    return(readRDS(cache_file))
  }
  
  return(NULL)
}

#' Write to cache
#' @param cache_key Cache key
#' @param data Data to cache
write_cache <- function(cache_key, data) {
  init_cache()
  cache_file <- file.path(CACHE_DIR, paste0(cache_key, ".rds"))
  saveRDS(data, cache_file)
  message("Data cached to: ", cache_file)
}

#' Fetch data with caching
#' @param station_id Station ID
#' @param start_date Start date
#' @param end_date End date
#' @param use_cache Whether to use caching
#' @param token API token
#' @return Weather data
fetch_with_cache <- function(station_id, start_date, end_date, 
                              use_cache = TRUE, token = NULL) {
  cache_key <- generate_cache_key(station_id, start_date, end_date)
  
  if (use_cache && is_cache_valid(cache_key)) {
    return(read_cache(cache_key))
  }
  
  data <- fetch_station_data(station_id, start_date, end_date, token = token)
  
  if (!is.null(data) && use_cache) {
    write_cache(cache_key, data)
  }
  
  return(data)
}

#' Clear cache
#' @param older_than_days Only clear files older than this many days (NULL = all)
clear_cache <- function(older_than_days = NULL) {
  if (!dir.exists(CACHE_DIR)) {
    message("No cache directory exists")
    return(invisible(NULL))
  }
  
  files <- list.files(CACHE_DIR, full.names = TRUE)
  
  if (!is.null(older_than_days)) {
    file_ages <- sapply(files, function(f) {
      as.numeric(difftime(Sys.time(), file.mtime(f), units = "days"))
    })
    files <- files[file_ages > older_than_days]
  }
  
  if (length(files) > 0) {
    file.remove(files)
    message("Removed ", length(files), " cached files")
  } else {
    message("No files to remove")
  }
}

# =============================================================================
# SAMPLE DATA GENERATION
# =============================================================================

#' Generate sample weather data for testing
#' @param n_locations Number of locations
#' @param seed Random seed
#' @return Dataframe with sample data
generate_sample_data <- function(n_locations = 100, seed = 42) {
  set.seed(seed)
  
  # Generate random locations worldwide
  data <- data.frame(
    name = paste0("Station_", sprintf("%03d", 1:n_locations)),
    latitude = runif(n_locations, -60, 70),
    longitude = runif(n_locations, -170, 170)
  )
  
  # Generate climate-appropriate values based on latitude
  data <- data %>%
    mutate(
      # Temperature varies with latitude
      TAVG = 30 - abs(latitude) * 0.5 + rnorm(n_locations, 0, 5),
      
      # Elevation (higher near equator for mountains)
      elevation = pmax(0, 
                       ifelse(abs(latitude) < 30, 
                              rexp(n_locations, 1/500) + 200,
                              rexp(n_locations, 1/200))),
      
      # Wind speed (higher at coast and higher latitudes)
      AWND = 3 + abs(latitude) * 0.03 + rexp(n_locations, 0.5),
      
      # Precipitation (higher in tropics and temperate regions)
      PRCP = ifelse(
        abs(latitude) < 15,
        rnorm(n_locations, 150, 50),
        ifelse(abs(latitude) < 45,
               rnorm(n_locations, 80, 30),
               rnorm(n_locations, 40, 20))
      ),
      
      # Snow depth (only in cold regions)
      SNWD = ifelse(
        TAVG < 0,
        rexp(n_locations, 1/30),
        0
      )
    ) %>%
    mutate(
      # Ensure reasonable ranges
      TAVG = pmax(-40, pmin(45, TAVG)),
      AWND = pmax(0, pmin(20, AWND)),
      PRCP = pmax(0, pmin(500, PRCP)),
      SNWD = pmax(0, pmin(200, SNWD)),
      elevation = round(elevation, 0)
    )
  
  message("✓ Generated sample data for ", n_locations, " locations")
  
  return(data)
}

#' Save sample data to CSV
#' @param data Sample data
#' @param filename Output filename
save_sample_data <- function(data, filename = "sample_weather_data.csv") {
  readr::write_csv(data, filename)
  message("✓ Saved sample data to: ", filename)
}

# =============================================================================
# INITIALIZATION
# =============================================================================

cat("🌐 Data Fetching Module v1.0.0 loaded\n")
cat("   API: NOAA Climate Data Online (CDO)\n")
cat("   Endpoints: stations, data\n")
cat("   Features: caching, pagination, rate limiting\n\n")
