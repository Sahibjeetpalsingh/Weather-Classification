# Renewable Energy Classifier

A web-based tool that determines the optimal renewable energy source (Solar, Wind, or Hydro) based on climate and geographic conditions.

## Project Structure

```
renewable-energy-classifier/
├── index.html          # Interactive web interface
├── classify_energy.R   # R script for data processing and classification
├── data.csv            # Output data with location classifications
└── README.md           # This file
```

## Files Description

### index.html
Interactive web interface that allows users to:
- Input climate parameters (temperature, wind speed, precipitation, elevation, snow depth)
- See real-time classification results
- View score breakdowns for each energy type
- Explore preset climate profiles (Desert, Coastal, Mountain, etc.)
- Visualize how parameters affect energy scores

### classify_energy.R
R script that:
- Loads and processes GHCN weather data
- Applies scoring algorithms for each energy type
- Classifies optimal renewable energy sources per location

### data.csv
Contains classification results with columns:
- `name`: Location name
- `longitude`, `latitude`: Coordinates
- `best_resource`: Recommended energy source (Solar, Wind, or Hydro)

## Scoring Algorithm

### Solar Energy (max 9 points)
- Temperature > 15°C: +3 pts
- Temperature > 25°C: +2 pts
- Precipitation < 50mm: +2 pts
- Precipitation > 150mm: -3 pts (penalty)
- Snow depth < 10cm: +1 pt
- Elevation < 2000m: +1 pt

### Wind Energy (max 7 points)
- Wind speed ≥ 4 m/s: +2 pts
- Wind speed ≥ 6 m/s: +3 pts
- Elevation 500-2000m: +1 pt
- Precipitation < 100mm: +1 pt

### Hydro Energy (max 8 points)
- Precipitation > 100mm: +2 pts
- Precipitation > 150mm: +2 pts
- Snow depth > 20cm: +2 pts
- Elevation 300-2000m: +2 pts

## Usage

### Web Interface
Simply open `index.html` in a web browser. No server required.

### R Script
```r
source("classify_energy.R")
```

## Data Source
GHCN (Global Historical Climatology Network) weather data containing temperature, precipitation, wind, snow, and elevation measurements from weather stations worldwide.
