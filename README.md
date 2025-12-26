# 🌍 Renewable Energy Classifier

[![Live Demo](https://img.shields.io/badge/demo-live-brightgreen?style=for-the-badge&logo=netlify)](https://weather-energy.netlify.app/)
[![HTML](https://img.shields.io/badge/HTML-95.2%25-orange?style=for-the-badge&logo=html5)](https://github.com/Sahibjeetpalsingh/Weather-Classification)
[![R](https://img.shields.io/badge/R-4.8%25-blue?style=for-the-badge&logo=r)](https://github.com/Sahibjeetpalsingh/Weather-Classification)
[![License](https://img.shields.io/badge/license-MIT-green?style=for-the-badge)](LICENSE)

> A smart web-based tool that determines the **optimal renewable energy source** (☀️ Solar, 💨 Wind, or 💧 Hydro) based on climate and geographic conditions using an intelligent scoring algorithm.

![Renewable Energy Classifier Preview](https://img.shields.io/badge/Preview-Interactive%20Web%20App-blueviolet?style=for-the-badge)

---

## 🚀 Live Demo

**[👉 Try the Live App](https://weather-energy.netlify.app/)**

Experience the classifier in action! Input your climate parameters and instantly see which renewable energy source is best suited for your location.

---

## ✨ Features

| Feature | Description |
|---------|-------------|
| 🎛️ **Interactive Sliders** | Adjust temperature, wind speed, precipitation, elevation, and snow depth in real-time |
| 📊 **Live Score Visualization** | See score breakdowns for Solar, Wind, and Hydro energy |
| 🌄 **Climate Presets** | Quick-select profiles for Desert, Coastal, Mountain, Tropical, and more |
| 📱 **Responsive Design** | Works seamlessly on desktop, tablet, and mobile devices |
| ⚡ **No Server Required** | Pure client-side application - just open and use |
| 🎨 **Modern UI/UX** | Clean, intuitive interface with smooth animations |

---

## 📁 Project Structure

```
renewable-energy-classifier/
│
├── 📄 index.html          # Interactive web interface (main application)
├── 📜 classify_energy.R   # R script for data processing & classification
├── 📊 data.csv            # Output data with location classifications
└── 📖 README.md           # Project documentation
```

---

## 📋 Files Description

### 🌐 `index.html`
The main interactive web interface that provides:
- **Real-time Classification**: Input climate parameters and instantly see results
- **Parameter Controls**: Sliders for temperature, wind speed, precipitation, elevation, snow depth
- **Score Breakdown**: Visual representation of scores for each energy type
- **Preset Profiles**: Pre-configured climate scenarios for quick testing
  - 🏜️ Desert Climate
  - 🏖️ Coastal Region
  - ⛰️ Mountain Area
  - 🌴 Tropical Zone
  - ❄️ Arctic/Cold Region
- **Responsive Layout**: Optimized for all screen sizes

### 📊 `classify_energy.R`
R script for batch processing weather data:
```r
# Key Libraries Used
library(jsonlite)    # JSON handling
library(ggplot2)     # Data visualization
library(tidyr)       # Data tidying
library(dplyr)       # Data manipulation
library(readr)       # CSV reading
library(purrr)       # Functional programming
```

**Functionality:**
- Loads and cleans GHCN weather station data
- Pivots observation data for analysis
- Applies scoring algorithms to determine optimal energy sources
- Outputs classification results to CSV

### 📈 `data.csv`
Classification results containing:

| Column | Description |
|--------|-------------|
| `name` | Weather station/location name |
| `longitude` | Geographic longitude coordinate |
| `latitude` | Geographic latitude coordinate |
| `best_resource` | Recommended energy source (Solar, Wind, or Hydro) |

---

## 🧮 Scoring Algorithm

The classifier uses a weighted scoring system based on scientific research into optimal conditions for each energy type.

### ☀️ Solar Energy (Maximum: 9 points)

| Condition | Points | Rationale |
|-----------|--------|-----------|
| Temperature > 15°C | +3 | Warmer regions typically have more sunshine |
| Temperature > 25°C | +2 | Hot climates correlate with high solar irradiance |
| Precipitation < 50mm | +2 | Low rainfall indicates clearer skies |
| Precipitation > 150mm | -3 | Heavy rain means cloud cover (penalty) |
| Snow depth < 10cm | +1 | Less snow = better panel efficiency |
| Elevation < 2000m | +1 | Moderate elevation for accessibility |

### 💨 Wind Energy (Maximum: 7 points)

| Condition | Points | Rationale |
|-----------|--------|-----------|
| Wind speed ≥ 4 m/s | +2 | Minimum viable wind for turbines |
| Wind speed ≥ 6 m/s | +3 | Optimal wind speed range |
| Elevation 500-2000m | +1 | High plateaus often have consistent winds |
| Precipitation < 100mm | +1 | Drier conditions for maintenance |

### 💧 Hydro Energy (Maximum: 8 points)

| Condition | Points | Rationale |
|-----------|--------|-----------|
| Precipitation > 100mm | +2 | Sufficient water supply |
| Precipitation > 150mm | +2 | High rainfall for reservoir filling |
| Snow depth > 20cm | +2 | Snowmelt provides seasonal water |
| Elevation 300-2000m | +2 | Elevation drop for hydroelectric generation |

### 🏆 Winner Determination
The energy source with the **highest score** is recommended. In case of a tie, the system prioritizes based on sustainability factors.

---

## 🛠️ Installation & Usage

### Option 1: Web Interface (Recommended)

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Sahibjeetpalsingh/Weather-Classification.git
   cd Weather-Classification
   ```

2. **Open the application:**
   ```bash
   # Simply open index.html in your browser
   # On Windows:
   start index.html
   
   # On macOS:
   open index.html
   
   # On Linux:
   xdg-open index.html
   ```

3. **Or visit the live demo:** [weather-energy.netlify.app](https://weather-energy.netlify.app/)

### Option 2: R Script (Batch Processing)

1. **Install required R packages:**
   ```r
   install.packages(c("jsonlite", "ggplot2", "tidyr", "dplyr", "readr", "purrr"))
   ```

2. **Run the classification script:**
   ```r
   source("classify_energy.R")
   ```

3. **View results:**
   ```r
   # Results will be processed and can be exported to data.csv
   ```

---

## 📊 Data Source

**GHCN (Global Historical Climatology Network)**

The project uses data from NOAA's GHCN database, which includes:
- 🌡️ **Temperature** - Average, min, max temperatures
- 🌧️ **Precipitation** - Rainfall measurements
- 💨 **Wind Speed** - Average wind velocities
- ❄️ **Snow Depth** - Snow accumulation data
- ⛰️ **Elevation** - Station altitude above sea level

> Data collected from weather stations worldwide, providing comprehensive climate metrics for accurate energy recommendations.

---

## 🎯 Use Cases

| Scenario | Application |
|----------|-------------|
| 🏗️ **Urban Planning** | Determine best renewable energy for new developments |
| 🏠 **Homeowners** | Decide between solar panels or small wind turbines |
| 🏢 **Businesses** | Plan sustainable energy investments |
| 📚 **Education** | Learn about renewable energy and climate factors |
| 🔬 **Research** | Analyze global renewable energy potential |

---

## 🤝 Contributing

Contributions are welcome! Here's how you can help:

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/AmazingFeature`)
3. **Commit** your changes (`git commit -m 'Add AmazingFeature'`)
4. **Push** to the branch (`git push origin feature/AmazingFeature`)
5. **Open** a Pull Request

### Ideas for Contribution:
- 🌍 Add more climate presets for different regions
- 📊 Implement data visualization charts
- 🌐 Add multilingual support
- 📱 Create a mobile app version
- 🔌 Add geothermal energy classification

---

## 📜 License

This project is open source and available under the [MIT License](LICENSE).

---

## 👨‍💻 Author

**Sahibjeet Pal Singh**

[![GitHub](https://img.shields.io/badge/GitHub-@Sahibjeetpalsingh-181717?style=flat-square&logo=github)](https://github.com/Sahibjeetpalsingh)

---

## 🙏 Acknowledgments

- [NOAA GHCN](https://www.ncdc.noaa.gov/ghcn-daily-description) for climate data
- [Netlify](https://netlify.com) for hosting
- The open-source community for inspiration

---

<p align="center">
  <strong>⭐ Star this repo if you found it helpful!</strong>
</p>

<p align="center">
  Made with ❤️ for a sustainable future 🌱
</p>
