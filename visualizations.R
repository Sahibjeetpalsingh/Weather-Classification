# =============================================================================
# ADVANCED VISUALIZATION FUNCTIONS FOR RENEWABLE ENERGY CLASSIFIER
# =============================================================================
# Author: Sahibjeet Pal Singh
# Date: December 2025
# Description: Comprehensive visualization library for creating publication-
#              quality plots, interactive charts, and geographic maps
# =============================================================================

# =============================================================================
# PACKAGE DEPENDENCIES
# =============================================================================

#' Load visualization packages
load_viz_packages <- function() {
  required_packages <- c(
    "ggplot2",      # Core plotting
    "scales",       # Scale functions
    "viridis",      # Color palettes
    "RColorBrewer", # Additional color palettes
    "gridExtra",    # Arrange multiple plots
    "patchwork",    # Compose ggplots
    "ggthemes",     # Additional themes
    "ggrepel"       # Text label repulsion
  )
  
  # Load available packages
  for (pkg in required_packages) {
    if (pkg %in% installed.packages()[, "Package"]) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
}

# Load packages
load_viz_packages()

# =============================================================================
# THEME DEFINITIONS
# =============================================================================

#' Custom theme for renewable energy visualizations
#' @param base_size Base font size
#' @param base_family Base font family
#' @return ggplot2 theme
theme_energy <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Title styling
      plot.title = element_text(
        size = base_size * 1.4,
        face = "bold",
        hjust = 0,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        size = base_size * 1.0,
        color = "gray50",
        hjust = 0,
        margin = margin(b = 15)
      ),
      plot.caption = element_text(
        size = base_size * 0.8,
        color = "gray60",
        hjust = 1,
        margin = margin(t = 10)
      ),
      
      # Axis styling
      axis.title = element_text(size = base_size * 0.9, face = "bold"),
      axis.text = element_text(size = base_size * 0.85),
      axis.line = element_line(color = "gray80", linewidth = 0.5),
      
      # Legend styling
      legend.title = element_text(size = base_size * 0.9, face = "bold"),
      legend.text = element_text(size = base_size * 0.85),
      legend.position = "bottom",
      legend.box = "horizontal",
      
      # Panel styling
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Plot margins
      plot.margin = margin(20, 20, 20, 20)
    )
}

#' Dark theme variant
#' @param base_size Base font size
#' @return ggplot2 theme
theme_energy_dark <- function(base_size = 12) {
  theme_energy(base_size) +
    theme(
      plot.background = element_rect(fill = "#1a1a1a", color = NA),
      panel.background = element_rect(fill = "#1a1a1a", color = NA),
      text = element_text(color = "white"),
      axis.text = element_text(color = "gray80"),
      axis.line = element_line(color = "gray60"),
      panel.grid.major = element_line(color = "gray30"),
      legend.background = element_rect(fill = "#1a1a1a", color = NA)
    )
}

# =============================================================================
# COLOR PALETTES
# =============================================================================

#' Get energy type color palette
#' @param type Palette type: "main", "light", "dark", or "gradient"
#' @return Named vector of colors
energy_colors <- function(type = "main") {
  palettes <- list(
    main = c(
      "Solar" = "#ea580c",
      "Wind" = "#0891b2",
      "Hydro" = "#2563eb"
    ),
    light = c(
      "Solar" = "#fed7aa",
      "Wind" = "#a5f3fc",
      "Hydro" = "#bfdbfe"
    ),
    dark = c(
      "Solar" = "#9a3412",
      "Wind" = "#0e7490",
      "Hydro" = "#1d4ed8"
    ),
    gradient = list(
      solar = c("#fff7ed", "#ea580c", "#7c2d12"),
      wind = c("#ecfeff", "#0891b2", "#164e63"),
      hydro = c("#eff6ff", "#2563eb", "#1e3a8a")
    )
  )
  
  return(palettes[[type]])
}

#' Get sequential color palette for heatmaps
#' @param n Number of colors
#' @param palette Palette name
#' @return Vector of colors
sequential_palette <- function(n = 9, palette = "viridis") {
  switch(palette,
         viridis = viridis::viridis(n),
         magma = viridis::magma(n),
         plasma = viridis::plasma(n),
         inferno = viridis::inferno(n),
         blues = colorRampPalette(c("#eff6ff", "#1e40af"))(n),
         greens = colorRampPalette(c("#f0fdf4", "#166534"))(n),
         oranges = colorRampPalette(c("#fff7ed", "#9a3412"))(n)
  )
}

# =============================================================================
# BAR CHARTS
# =============================================================================

#' Create distribution bar chart
#' @param df Classified dataframe
#' @param show_percentages Whether to show percentage labels
#' @param horizontal Whether to make horizontal bars
#' @return ggplot object
create_distribution_chart <- function(df, show_percentages = TRUE, horizontal = FALSE) {
  # Prepare data
  plot_data <- df %>%
    count(best_resource, name = "count") %>%
    mutate(
      percentage = count / sum(count) * 100,
      label = if (show_percentages) {
        sprintf("%d\n(%.1f%%)", count, percentage)
      } else {
        as.character(count)
      }
    )
  
  # Create base plot
  p <- ggplot(plot_data, aes(
    x = reorder(best_resource, -count),
    y = count,
    fill = best_resource
  )) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(
      aes(label = label),
      vjust = if (horizontal) 0.5 else -0.5,
      hjust = if (horizontal) -0.1 else 0.5,
      size = 4,
      fontface = "bold"
    ) +
    scale_fill_manual(values = energy_colors("main")) +
    labs(
      title = "Renewable Energy Classification Distribution",
      subtitle = sprintf("Analysis of %s locations", format(nrow(df), big.mark = ",")),
      x = "Energy Type",
      y = "Number of Locations",
      caption = "Source: GHCN Weather Data"
    ) +
    theme_energy()
  
  # Adjust for orientation
  if (horizontal) {
    p <- p + coord_flip()
  } else {
    p <- p + expand_limits(y = max(plot_data$count) * 1.2)
  }
  
  return(p)
}

#' Create grouped bar chart comparing scores
#' @param df Classified dataframe
#' @return ggplot object
create_score_bars <- function(df) {
  # Calculate mean scores
  score_data <- df %>%
    summarise(
      Solar = mean(solar_score, na.rm = TRUE),
      Wind = mean(wind_score, na.rm = TRUE),
      Hydro = mean(hydro_score, na.rm = TRUE)
    ) %>%
    pivot_longer(everything(), names_to = "energy_type", values_to = "mean_score")
  
  # Add max scores
  score_data <- score_data %>%
    mutate(
      max_score = case_when(
        energy_type == "Solar" ~ 9,
        energy_type == "Wind" ~ 7,
        energy_type == "Hydro" ~ 8
      ),
      percentage = mean_score / max_score * 100
    )
  
  # Create plot
  ggplot(score_data, aes(x = energy_type, y = mean_score, fill = energy_type)) +
    geom_col(width = 0.7) +
    geom_text(
      aes(label = sprintf("%.1f / %d\n(%.0f%%)", mean_score, max_score, percentage)),
      vjust = -0.5,
      size = 3.5,
      fontface = "bold"
    ) +
    scale_fill_manual(values = energy_colors("main")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Average Scores by Energy Type",
      subtitle = "Mean score achieved out of maximum possible",
      x = NULL,
      y = "Average Score"
    ) +
    theme_energy() +
    theme(legend.position = "none")
}

# =============================================================================
# BOXPLOTS AND DISTRIBUTIONS
# =============================================================================

#' Create score distribution boxplot
#' @param df Classified dataframe
#' @param show_points Whether to show individual data points
#' @return ggplot object
create_score_boxplot <- function(df, show_points = TRUE) {
  # Reshape data
  score_data <- df %>%
    select(solar_score, wind_score, hydro_score) %>%
    pivot_longer(
      everything(),
      names_to = "energy_type",
      values_to = "score"
    ) %>%
    mutate(
      energy_type = case_when(
        energy_type == "solar_score" ~ "Solar",
        energy_type == "wind_score" ~ "Wind",
        energy_type == "hydro_score" ~ "Hydro"
      ),
      energy_type = factor(energy_type, levels = c("Solar", "Wind", "Hydro"))
    )
  
  # Create plot
  p <- ggplot(score_data, aes(x = energy_type, y = score, fill = energy_type)) +
    geom_boxplot(
      alpha = 0.8,
      outlier.shape = if (show_points) NA else 21,
      width = 0.6
    )
  
  if (show_points) {
    p <- p + geom_jitter(
      width = 0.15,
      alpha = 0.4,
      size = 1.5,
      color = "gray30"
    )
  }
  
  p <- p +
    scale_fill_manual(values = energy_colors("main")) +
    labs(
      title = "Score Distribution by Energy Type",
      subtitle = "Boxplot showing median, quartiles, and individual scores",
      x = NULL,
      y = "Score"
    ) +
    theme_energy() +
    theme(legend.position = "none")
  
  return(p)
}

#' Create violin plot for score distributions
#' @param df Classified dataframe
#' @return ggplot object
create_score_violin <- function(df) {
  # Reshape data
  score_data <- df %>%
    select(solar_score, wind_score, hydro_score) %>%
    pivot_longer(
      everything(),
      names_to = "energy_type",
      values_to = "score"
    ) %>%
    mutate(
      energy_type = case_when(
        energy_type == "solar_score" ~ "Solar",
        energy_type == "wind_score" ~ "Wind",
        energy_type == "hydro_score" ~ "Hydro"
      )
    )
  
  ggplot(score_data, aes(x = energy_type, y = score, fill = energy_type)) +
    geom_violin(alpha = 0.8, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 18,
      size = 3,
      color = "darkred"
    ) +
    scale_fill_manual(values = energy_colors("main")) +
    labs(
      title = "Score Distribution - Violin Plot",
      subtitle = "Showing probability density with embedded boxplot",
      x = NULL,
      y = "Score"
    ) +
    theme_energy() +
    theme(legend.position = "none")
}

# =============================================================================
# HISTOGRAMS AND DENSITY PLOTS
# =============================================================================

#' Create histogram of a specific variable
#' @param df Dataframe with data
#' @param column Column name to plot
#' @param bins Number of bins
#' @param fill_color Fill color
#' @return ggplot object
create_histogram <- function(df, column, bins = 30, fill_color = "#3b82f6") {
  # Check column exists
  if (!column %in% names(df)) {
    stop(paste("Column", column, "not found in dataframe"))
  }
  
  ggplot(df, aes(x = .data[[column]])) +
    geom_histogram(
      bins = bins,
      fill = fill_color,
      color = "white",
      alpha = 0.8
    ) +
    geom_density(
      aes(y = after_stat(count)),
      color = "darkred",
      linewidth = 1
    ) +
    labs(
      title = paste("Distribution of", column),
      subtitle = paste("N =", nrow(df), "observations"),
      x = column,
      y = "Frequency"
    ) +
    theme_energy()
}

#' Create faceted histogram by energy type
#' @param df Classified dataframe
#' @param column Column to plot
#' @return ggplot object
create_faceted_histogram <- function(df, column) {
  if (!column %in% names(df) || !"best_resource" %in% names(df)) {
    stop("Required columns not found")
  }
  
  ggplot(df, aes(x = .data[[column]], fill = best_resource)) +
    geom_histogram(bins = 20, alpha = 0.8, color = "white") +
    facet_wrap(~best_resource, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = energy_colors("main")) +
    labs(
      title = paste(column, "Distribution by Recommended Energy Type"),
      x = column,
      y = "Count"
    ) +
    theme_energy() +
    theme(legend.position = "none")
}

#' Create density comparison plot
#' @param df Classified dataframe
#' @param column Column to compare
#' @return ggplot object
create_density_comparison <- function(df, column) {
  ggplot(df, aes(x = .data[[column]], fill = best_resource, color = best_resource)) +
    geom_density(alpha = 0.4, linewidth = 1) +
    scale_fill_manual(values = energy_colors("main")) +
    scale_color_manual(values = energy_colors("dark")) +
    labs(
      title = paste(column, "Density by Energy Type"),
      subtitle = "Overlapping density curves for comparison",
      x = column,
      y = "Density",
      fill = "Energy Type",
      color = "Energy Type"
    ) +
    theme_energy()
}

# =============================================================================
# SCATTER PLOTS
# =============================================================================

#' Create scatter plot of two variables
#' @param df Classified dataframe
#' @param x_var X-axis variable
#' @param y_var Y-axis variable
#' @param color_by Variable to color by (default: best_resource)
#' @return ggplot object
create_scatter <- function(df, x_var, y_var, color_by = "best_resource") {
  p <- ggplot(df, aes(
    x = .data[[x_var]],
    y = .data[[y_var]],
    color = .data[[color_by]]
  )) +
    geom_point(alpha = 0.6, size = 2) +
    labs(
      title = paste(y_var, "vs", x_var),
      subtitle = paste("Colored by", color_by),
      x = x_var,
      y = y_var
    ) +
    theme_energy()
  
  if (color_by == "best_resource") {
    p <- p + scale_color_manual(values = energy_colors("main"), name = "Energy Type")
  }
  
  return(p)
}

#' Create scatter plot with regression lines
#' @param df Classified dataframe
#' @param x_var X-axis variable
#' @param y_var Y-axis variable
#' @return ggplot object
create_scatter_with_regression <- function(df, x_var, y_var) {
  ggplot(df, aes(
    x = .data[[x_var]],
    y = .data[[y_var]],
    color = best_resource
  )) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    scale_color_manual(values = energy_colors("main")) +
    labs(
      title = paste(y_var, "vs", x_var, "with Regression Lines"),
      subtitle = "Linear regression by energy type",
      x = x_var,
      y = y_var,
      color = "Energy Type"
    ) +
    theme_energy()
}

# =============================================================================
# GEOGRAPHIC PLOTS
# =============================================================================

#' Create geographic scatter map
#' @param df Classified dataframe with lat/lon
#' @param point_size Size of points
#' @return ggplot object
create_geo_scatter <- function(df, point_size = 2) {
  if (!all(c("latitude", "longitude", "best_resource") %in% names(df))) {
    stop("Dataframe must contain latitude, longitude, and best_resource columns")
  }
  
  ggplot(df, aes(x = longitude, y = latitude, color = best_resource)) +
    geom_point(alpha = 0.6, size = point_size) +
    scale_color_manual(values = energy_colors("main")) +
    labs(
      title = "Geographic Distribution of Optimal Energy Sources",
      subtitle = "Each point represents a classified weather station",
      x = "Longitude",
      y = "Latitude",
      color = "Recommended\nEnergy Source"
    ) +
    theme_energy() +
    coord_fixed(ratio = 1.3) +
    theme(legend.position = "right")
}

#' Create world map with classifications
#' @param df Classified dataframe with lat/lon
#' @return ggplot object
create_world_map <- function(df) {
  # Get world map data
  if (!requireNamespace("maps", quietly = TRUE)) {
    message("Package 'maps' not installed. Using simple scatter plot.")
    return(create_geo_scatter(df))
  }
  
  world <- map_data("world")
  
  ggplot() +
    geom_polygon(
      data = world,
      aes(x = long, y = lat, group = group),
      fill = "gray95",
      color = "gray80",
      linewidth = 0.2
    ) +
    geom_point(
      data = df,
      aes(x = longitude, y = latitude, color = best_resource),
      alpha = 0.7,
      size = 2
    ) +
    scale_color_manual(values = energy_colors("main")) +
    labs(
      title = "Global Renewable Energy Classification Map",
      subtitle = "Optimal energy source by geographic location",
      color = "Energy Type"
    ) +
    theme_energy() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-60, 85))
}

#' Create latitude band analysis chart
#' @param df Classified dataframe with latitude
#' @return ggplot object
create_latitude_analysis <- function(df) {
  if (!"latitude" %in% names(df)) {
    stop("Dataframe must contain latitude column")
  }
  
  # Create latitude bands
  lat_data <- df %>%
    mutate(
      lat_band = cut(
        latitude,
        breaks = seq(-90, 90, by = 15),
        labels = paste0(seq(-90, 75, by = 15), "° to ", seq(-75, 90, by = 15), "°"),
        include.lowest = TRUE
      )
    ) %>%
    filter(!is.na(lat_band)) %>%
    count(lat_band, best_resource) %>%
    group_by(lat_band) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ungroup()
  
  ggplot(lat_data, aes(x = lat_band, y = percentage, fill = best_resource)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = energy_colors("main")) +
    labs(
      title = "Energy Type Distribution by Latitude Band",
      subtitle = "Percentage of each energy type across latitude zones",
      x = "Latitude Band",
      y = "Percentage (%)",
      fill = "Energy Type"
    ) +
    theme_energy() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# =============================================================================
# HEATMAPS
# =============================================================================

#' Create correlation heatmap
#' @param df Dataframe with numeric columns
#' @param variables Vector of variable names to include
#' @return ggplot object
create_correlation_heatmap <- function(df, variables = NULL) {
  # Select numeric columns
  if (is.null(variables)) {
    numeric_df <- df %>% select(where(is.numeric))
  } else {
    numeric_df <- df %>% select(all_of(variables))
  }
  
  # Calculate correlation matrix
  cor_matrix <- cor(numeric_df, use = "pairwise.complete.obs")
  
  # Reshape for plotting
  cor_data <- as.data.frame(as.table(cor_matrix))
  names(cor_data) <- c("Var1", "Var2", "Correlation")
  
  ggplot(cor_data, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = sprintf("%.2f", Correlation)),
      color = ifelse(abs(cor_data$Correlation) > 0.5, "white", "black"),
      size = 3
    ) +
    scale_fill_gradient2(
      low = "#2563eb",
      mid = "white",
      high = "#ea580c",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    labs(
      title = "Correlation Heatmap",
      subtitle = "Pairwise correlations between variables",
      x = NULL,
      y = NULL
    ) +
    theme_energy() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    coord_fixed()
}

#' Create score heatmap by climate zone
#' @param df Classified dataframe with latitude
#' @return ggplot object
create_climate_heatmap <- function(df) {
  # Add climate zone
  heatmap_data <- df %>%
    mutate(
      climate_zone = case_when(
        abs(latitude) < 23.5 ~ "Tropical",
        abs(latitude) < 35 ~ "Subtropical",
        abs(latitude) < 55 ~ "Temperate",
        TRUE ~ "Polar/Subarctic"
      )
    ) %>%
    group_by(climate_zone) %>%
    summarise(
      Solar = mean(solar_score, na.rm = TRUE),
      Wind = mean(wind_score, na.rm = TRUE),
      Hydro = mean(hydro_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(Solar, Wind, Hydro),
      names_to = "energy_type",
      values_to = "avg_score"
    )
  
  ggplot(heatmap_data, aes(x = energy_type, y = climate_zone, fill = avg_score)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(
      aes(label = sprintf("%.1f", avg_score)),
      color = "white",
      fontface = "bold",
      size = 5
    ) +
    scale_fill_viridis_c(option = "plasma", name = "Avg Score") +
    labs(
      title = "Average Scores by Climate Zone",
      subtitle = "Heatmap of mean scores across different climate zones",
      x = "Energy Type",
      y = "Climate Zone"
    ) +
    theme_energy() +
    theme(panel.grid = element_blank()) +
    coord_fixed()
}

# =============================================================================
# RADAR/SPIDER CHARTS
# =============================================================================

#' Create radar chart for single location
#' @param solar_score Solar score
#' @param wind_score Wind score
#' @param hydro_score Hydro score
#' @param title Chart title
#' @return ggplot object
create_radar_chart <- function(solar_score, wind_score, hydro_score, 
                                title = "Energy Potential Profile") {
  # Normalize scores to percentage
  data <- data.frame(
    energy = c("Solar", "Wind", "Hydro"),
    score = c(solar_score / 9, wind_score / 7, hydro_score / 8) * 100,
    max_score = c(9, 7, 8)
  )
  
  # Create circular coordinates
  n <- nrow(data)
  angles <- seq(0, 2 * pi, length.out = n + 1)[1:n]
  
  data$x <- data$score * cos(angles - pi/2)
  data$y <- data$score * sin(angles - pi/2)
  
  # Create polygon data
  polygon_data <- rbind(data, data[1, ])
  
  # Grid circles
  grid_levels <- c(25, 50, 75, 100)
  
  ggplot() +
    # Grid circles
    lapply(grid_levels, function(r) {
      geom_path(
        data = data.frame(
          x = r * cos(seq(0, 2*pi, length.out = 100)),
          y = r * sin(seq(0, 2*pi, length.out = 100))
        ),
        aes(x = x, y = y),
        color = "gray80",
        linewidth = 0.3
      )
    }) +
    # Axis lines
    geom_segment(
      data = data,
      aes(x = 0, y = 0, xend = 100 * cos(angles - pi/2), yend = 100 * sin(angles - pi/2)),
      color = "gray70",
      linewidth = 0.5
    ) +
    # Score polygon
    geom_polygon(
      data = polygon_data,
      aes(x = x, y = y),
      fill = "#3b82f6",
      alpha = 0.4
    ) +
    geom_path(
      data = polygon_data,
      aes(x = x, y = y),
      color = "#1d4ed8",
      linewidth = 1.5
    ) +
    # Points
    geom_point(
      data = data,
      aes(x = x, y = y),
      color = "#1d4ed8",
      size = 4
    ) +
    # Labels
    geom_text(
      data = data,
      aes(
        x = 115 * cos(angles - pi/2),
        y = 115 * sin(angles - pi/2),
        label = paste0(energy, "\n", round(score), "%")
      ),
      size = 4,
      fontface = "bold"
    ) +
    coord_fixed() +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
}

# =============================================================================
# COMPOSITE VISUALIZATIONS
# =============================================================================

#' Create comprehensive dashboard of all visualizations
#' @param df Classified dataframe
#' @return Combined ggplot object (requires patchwork package)
create_dashboard <- function(df) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' required for dashboard. Install with install.packages('patchwork')")
  }
  
  library(patchwork)
  
  # Create individual plots
  p1 <- create_distribution_chart(df)
  p2 <- create_score_boxplot(df, show_points = FALSE)
  
  if (all(c("latitude", "longitude") %in% names(df))) {
    p3 <- create_geo_scatter(df, point_size = 1)
  } else {
    p3 <- create_score_bars(df)
  }
  
  p4 <- create_score_violin(df)
  
  # Combine plots
  dashboard <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = "Renewable Energy Classification Dashboard",
      subtitle = paste("Comprehensive analysis of", nrow(df), "locations"),
      caption = paste("Generated:", Sys.time()),
      theme = theme(
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 14, color = "gray50")
      )
    )
  
  return(dashboard)
}

# =============================================================================
# EXPORT FUNCTIONS
# =============================================================================

#' Save plot with consistent settings
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution
save_plot <- function(plot, filename, width = 10, height = 6, dpi = 300) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  message("✓ Saved plot to: ", filename)
}

#' Generate and save all standard visualizations
#' @param df Classified dataframe
#' @param output_dir Output directory for plots
generate_all_plots <- function(df, output_dir = "plots") {
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("Generating visualizations...")
  
  # Distribution chart
  p <- create_distribution_chart(df)
  save_plot(p, file.path(output_dir, "distribution.png"))
  
  # Score boxplot
  p <- create_score_boxplot(df)
  save_plot(p, file.path(output_dir, "score_boxplot.png"))
  
  # Score violin
  p <- create_score_violin(df)
  save_plot(p, file.path(output_dir, "score_violin.png"))
  
  # Score bars
  p <- create_score_bars(df)
  save_plot(p, file.path(output_dir, "score_bars.png"))
  
  # Geographic plots (if coordinates available)
  if (all(c("latitude", "longitude") %in% names(df))) {
    p <- create_geo_scatter(df)
    save_plot(p, file.path(output_dir, "geographic_scatter.png"), width = 12, height = 8)
    
    p <- create_latitude_analysis(df)
    save_plot(p, file.path(output_dir, "latitude_analysis.png"))
  }
  
  # Climate heatmap (if latitude available)
  if ("latitude" %in% names(df)) {
    p <- create_climate_heatmap(df)
    save_plot(p, file.path(output_dir, "climate_heatmap.png"), width = 8, height = 6)
  }
  
  message("✓ All visualizations saved to: ", output_dir)
}

# =============================================================================
# PRINT VERSION INFO
# =============================================================================

cat("📊 Visualization Functions v1.0.0 loaded\n")
cat("   Themes: theme_energy(), theme_energy_dark()\n")
cat("   Palettes: energy_colors(), sequential_palette()\n")
cat("   Charts: 15+ visualization functions available\n\n")
