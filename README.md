<div align="center">

# 🌍 Renewable Energy Classifier

[![Live App](https://img.shields.io/badge/Try%20Live-weather--energy.netlify.app-brightgreen?style=for-the-badge&logo=netlify)](https://weather-energy.netlify.app/)
[![R](https://img.shields.io/badge/R-3000%2B%20lines-276DC3?style=for-the-badge&logo=r)](https://github.com/Sahibjeetpalsingh/Weather-Classification)
[![HTML](https://img.shields.io/badge/HTML%2FJS-Interactive%20Web%20App-E34F26?style=for-the-badge&logo=html5)](https://weather-energy.netlify.app/)
[![NOAA](https://img.shields.io/badge/Data-NOAA%20GHCN-0072C6?style=for-the-badge)](https://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](LICENSE)

*Solar, Wind, or Hydro — what does the climate of your location actually support?*

</div>

---

## A Question Nobody Had a Good Answer To

A few years ago, the conversation around renewable energy shifted. It stopped being about whether to transition away from fossil fuels and started being about how and where to do it. Governments are setting timelines. Developers are planning projects. Homeowners are installing rooftop panels or wondering whether they should. And underneath all of that activity is a question that sounds simple but almost never gets answered clearly: given the actual climate of a specific place, which renewable energy source is the most viable?

Not in general. Not "solar is good in sunny places." Specifically. With numbers. With a justification you can point to and argue with.

The frustrating reality is that most people who need this answer — a homeowner in Vancouver, a city planner in rural Alberta, a student writing a sustainability report — have almost no good tools to get it. Energy consulting firms exist but charge thousands of dollars. Academic energy siting literature exists but requires a PhD to navigate. Government planning databases exist but are designed for engineers, not for curious non-specialists who just want a clear answer.

This project was built to fill that gap. The goal was a tool that anyone could open in a browser, input the climate conditions of their location, and immediately get a scored, justified, transparent recommendation: Solar, Wind, or Hydro — and here is exactly why.

---

## Going Deep Into the Problem: What Does Each Energy Type Actually Need?

Before writing a single line of code, the first task was understanding the domain. What does solar energy actually need to be viable? What are the hard requirements for wind? What makes hydro work in one watershed and fail in another? The answers to these questions come from the physics and economics of energy generation, and they are well-established in the engineering literature — but they are not the kind of thing most people learn unless they specifically go looking.

**Solar photovoltaics** convert sunlight into electricity. The efficiency of that conversion depends on solar irradiance — how much sunlight reaches the panel surface per unit time. Temperature correlates with irradiance in most climates, not because heat itself powers the panels (in fact, very high temperatures reduce panel efficiency slightly), but because warm, dry climates tend to have clearer skies and longer effective daylight hours. Cloud cover is the enemy: every cloud reduces the energy reaching the panel. Precipitation is therefore a strong negative signal for solar viability — regions that receive more than 150mm of monthly precipitation are persistently cloudy. Snow is doubly problematic: it blocks panels and reduces output. Elevation below 2000m matters for installation feasibility and infrastructure access.

**Wind turbines** convert kinetic energy in moving air into electricity. The power output of a wind turbine scales with the cube of wind speed — which means a small increase in average wind speed produces a disproportionately large increase in power generation. This cubic relationship also creates a hard economic threshold: at average wind speeds below roughly 4 meters per second, most commercial turbines do not generate electricity at economically viable rates. Below that threshold, the turbine spins but the energy is barely worth the capital and maintenance cost. Above 6 m/s, turbines enter their optimal operating range. Elevation matters enormously for wind because terrain channels and accelerates airflow — high plateaus and ridgelines reliably see faster, steadier winds than low valleys. Low precipitation is a secondary positive factor because corrosion from moisture is a significant maintenance cost driver for turbine components.

**Hydroelectric generation** converts the gravitational potential energy of water into electricity. It fundamentally needs two things: water, and height. Water comes from precipitation and snowmelt. Height comes from terrain — water must fall to spin a turbine, and the greater the vertical drop, the more power is generated per unit volume of flow. This is why countries like Norway, Switzerland, Canada, and Brazil — all of which have significant mountain ranges and high precipitation or snowmelt — are the world's dominant hydroelectric producers. Snowpack is particularly important because it stores precipitation from the winter months and releases it gradually through spring and early summer, providing a sustained water supply for reservoirs. Regions with heavy snowfall above 20cm consistently exhibit this pattern. Elevation between 300m and 2000m is the productive range for hydro — too low and there is not enough elevation drop; too high and the terrain is often too steep and remote for practical installation.

Understanding these requirements clearly and precisely — not as vague intuitions but as specific quantitative thresholds — was the foundation that everything else in the project was built on.

---

## The Design Decision That Shaped Everything

Once the physical requirements were mapped out, the question became: what kind of tool should encode them?

The first instinct, for anyone with a machine learning background, is to train a classifier. Find a dataset of locations labeled with their optimal energy type, train a model on it, deploy it. But this approach runs into a fundamental problem: that dataset does not exist in any reliable form. There is no public, well-labeled dataset saying "this weather station should use solar, that one should use wind." You would have to build the labels yourself — which means you are already making the expert judgments you were hoping the model would make for you.

More critically: even if you could build such a dataset, a trained classifier would be a black box. It would give you a prediction with no explanation. For an energy planning tool aimed at non-specialists — people who are going to make real, expensive infrastructure decisions based on the output — a recommendation with no reasoning is not just unhelpful. It is actively dangerous. If the tool says "Wind" and you cannot see why, you have no way to evaluate whether that recommendation makes sense for your specific context, or whether the model has learned a spurious correlation from its training data.

The decision was therefore to build a **rule-based weighted scoring system** that directly encodes the physical science into explicit, readable logic. Each energy type is scored independently across the five input dimensions. The scores are grounded in the physical thresholds from the research phase — not arbitrary numbers, but thresholds that correspond to real engineering viability limits. The full scoring breakdown is visible to the user. You can see exactly which conditions fired, exactly how many points each one contributed, and exactly why the winner won. If you disagree with the recommendation, you know precisely where to push back.

This design choice also made a secondary observation possible: a rule-based system lets you reason about uncertainty. When the winner's score is far ahead of the runners-up, the recommendation is high confidence. When the scores are close, the tool says so explicitly — "Low confidence, consider a hybrid approach" — rather than pretending to more certainty than the data supports. That kind of calibrated honesty is something trained classifiers rarely provide without significant extra engineering.

---

## Building the Scoring Algorithm From First Principles

The scoring algorithm assigns points to each energy type based on whether specific climate conditions are met. The thresholds are grounded in the domain research described above.

Solar energy earns up to 9 points. A location with average temperature above 15°C earns 3 points — this threshold captures the bulk of the warm-climate, high-irradiance zone. If the temperature exceeds 25°C, it earns 2 additional points because the correlation with clear, intense sunlight strengthens further. Monthly precipitation below 50mm earns 2 points, reflecting low cloud cover and high sunshine hours. Heavy precipitation above 150mm subtracts 3 points — this is a penalty, not just an absence of points, because persistent cloud cover actively suppresses solar output below the baseline of mild conditions. Snow depth below 10cm earns 1 point for keeping panels clear and accessible. Elevation below 2000m earns 1 point for infrastructure viability.

Wind energy earns up to 7 points. Wind speed at or above 4 m/s earns 2 points — this is the minimum viability threshold for commercial turbines. Wind speed at or above 6 m/s earns 3 additional points, putting the location in the optimal operating range. Elevation between 500m and 2000m earns 1 point for the terrain advantage in wind channeling. Precipitation below 100mm earns 1 point for reduced corrosion maintenance.

Hydro energy earns up to 8 points. Monthly precipitation above 100mm earns 2 points for baseline water supply. Precipitation above 150mm earns 2 more points for reservoir-grade rainfall. Snow depth above 20cm earns 2 points for the snowmelt contribution — this single feature captures much of the "mountain hydro" signal. Elevation between 300m and 2000m earns 2 points for the gravitational head required to generate power.

The winner is the energy type with the highest total score. The confidence signal is derived from the margin between winner and runner-up: a margin of 6 or more points is high confidence, 3 to 5 points is moderate confidence, and 2 or fewer points is low confidence with a hybrid approach suggested. This confidence signal was not an afterthought — it directly addresses the most common failure mode of simple decision tools, which is presenting every recommendation with equal certainty regardless of how close the call actually was.

---

## The Web Application: Built for Accessibility, Not Complexity

The front-facing part of this project is a pure HTML, CSS, and JavaScript application deployed on Netlify at weather-energy.netlify.app. The choice of vanilla web technologies — no React, no Vue, no Node backend, no build pipeline — was deliberate. The application needs to be instantly accessible to anyone with a browser, including people in environments where installing software is not straightforward. A single HTML file that opens in any browser and works immediately is the most accessible possible delivery mechanism.

The interface is structured around five climate parameter sliders. Average Temperature ranges from minus 20 to plus 50 degrees Celsius, covering virtually all inhabited climates on Earth. Wind Speed runs from 0 to 15 m/s, capturing the range from dead calm to storm-force. Monthly Precipitation runs from 0 to 300mm, spanning desert to monsoon. Elevation runs from sea level to 4000m. Snow Depth runs from 0 to 200cm. As the user moves any slider, the scoring algorithm runs synchronously in JavaScript and the results panel updates immediately — no network call, no server, no latency. The recommendation changes in real time as the parameters change.

Eight climate preset buttons — Desert, Coastal, Mountain, Tropical, Arctic, Plains, Mediterranean, and Monsoon — provide one-click population of all five sliders with representative values for each climate archetype. These presets serve two purposes: they let first-time users immediately see how the tool behaves across diverse climate types without having to know what values to enter, and they provide a quick reference for common scenarios that planners and researchers encounter repeatedly.

The results panel shows the recommended energy source with a large, unambiguous visual label, the confidence level, the numerical scores for all three energy types displayed against their respective maximums, and a full condition-by-condition breakdown showing which scoring rules fired and how many points each contributed. A user can read the breakdown and understand not just what the recommendation is, but why every single point was awarded — which conditions were met and which were not.

The design is responsive and tested on desktop, tablet, and mobile viewports. This matters because planners in the field who are assessing a location on a phone should be able to get the same quality of recommendation as someone sitting at a desktop workstation.

---

## The R Analysis Suite: 3,000+ Lines of Research Infrastructure

While the web interface is designed for accessibility, the R suite underneath it is designed for rigor. It is the kind of analytical infrastructure that would sit behind a professional energy assessment study or an academic publication.

The main classification engine in `classify_energy.R` — roughly 800 lines — implements the same scoring algorithm in R, but optimized for batch processing. Rather than scoring a single location at a time, it processes entire dataframes of locations in a vectorized fashion, which is orders of magnitude faster for large-scale geographic analysis. The scoring functions — `calculate_solar_score()`, `calculate_wind_score()`, and `calculate_hydro_score()` — are written as independent, testable units that each return not just a total score but a complete breakdown of which conditions fired and how many points each contributed. The overall `classify_energy()` function orchestrates these three functions and adds confidence calculation and tie-breaking logic. A global `CONFIG` object holds all scoring thresholds, which means every threshold in the system can be adjusted in one place without touching any of the scoring functions themselves. The module also includes `generate_summary_stats()`, which produces a statistical summary of the classification results across a full dataset, and `main()`, which executes the full pipeline from loading raw data through to exporting classified results to both CSV and JSON.

The statistical analysis module in `analysis.R` — roughly 600 lines — is where the empirical validation of the scoring algorithm happens. It runs one-way ANOVA tests to check whether the climate variables (temperature, wind speed, precipitation, elevation, snow depth) differ significantly across the three energy recommendation categories. If the algorithm is working correctly, you should see highly significant ANOVA results — the climate conditions in Solar-recommended locations should be measurably different from Wind-recommended and Hydro-recommended locations. It computes full correlation matrices using Pearson, Spearman, and Kendall methods with associated p-values, revealing which input variables are most strongly related to each other. It builds multivariate linear regression models with variance inflation factor analysis to check for multicollinearity — if two input variables are so closely correlated that they carry redundant information, this shows up in elevated VIF scores and suggests the scoring algorithm might be double-counting. It generates confusion matrix metrics, F1 scores, precision, recall, and Cohen's Kappa to measure how consistently the classification algorithm performs across different data subsets. And it runs permutation-based feature importance analysis to empirically determine which input variables have the most influence on the final classification — confirming or challenging the assumptions built into the scoring weights.

The visualization module in `visualizations.R` — roughly 700 lines — produces publication-quality figures using ggplot2. The library includes distribution bar charts showing how energy recommendations are distributed across a dataset, score boxplots and violin plots comparing the score distributions for each energy type across different climate zones, histograms with density overlays for individual climate variables, scatter plots with regression lines for pairwise variable relationships, geographic scatter maps that plot energy recommendations spatially using latitude and longitude coordinates, correlation heatmaps with color-coded significance levels, radar charts that profile a single location across all five climate dimensions, and multi-panel dashboard layouts that combine multiple chart types into a single comprehensive view. Two custom ggplot2 themes are defined — `theme_energy()` for light-background publications and `theme_energy_dark()` for dark-background presentations — along with `energy_colors()`, a consistent color palette function that maps the three energy types to the same orange, blue, and green across every chart type.

The data fetching module in `data_fetch.R` — roughly 500 lines — connects to NOAA's Climate Data Online (CDO) API, which is the authoritative source for global historical climate records. The module handles the full complexity of working with a rate-limited, paginated REST API: it manages API key authentication via HTTP headers, implements exponential backoff and retry logic for rate limit responses, handles pagination to retrieve datasets that span thousands of records across multiple API pages, converts GHCN raw unit codes to standard meteorological values (GHCN stores precipitation in tenths of millimeters, temperatures in tenths of degrees Celsius, wind speeds in tenths of meters per second), aggregates raw daily observations to monthly climatological averages, and implements a file-based caching system that saves API responses to disk so that repeated analysis runs do not re-fetch data that has already been downloaded. The module can fetch weather station metadata for any geographic bounding box, retrieve observations from specific stations or entire regions, and process the raw GHCN format into clean, analysis-ready dataframes.

The utility module in `utils.R` — roughly 500 lines — is the shared infrastructure layer that the other four modules draw on. It includes input validation functions that check numeric bounds and geographic coordinate ranges, unit conversion functions for temperature, precipitation, wind speed, and elevation between imperial and metric systems, climate zone classification from latitude and altitude (tropical, subtropical, temperate, subarctic, arctic), daylight hour calculations by hemisphere, latitude, and season, statistical utilities including confidence interval computation, coefficient of variation, normalization methods, and weighted mean, string formatting helpers for producing human-readable percentage and score strings, timestamped logging functions for pipeline execution tracking, and file system utilities for creating output directories and generating timestamped filenames.

The data that all of this infrastructure operates on comes from NOAA's Global Historical Climatology Network — GHCN — which integrates climate records from weather stations worldwide into a consistent, quality-controlled database. The project includes `data.csv`, which contains classified location data generated by running the R pipeline on a set of NOAA station records.

---

## What the Analysis Actually Found

Running the full R analysis suite on real NOAA climate data produced results that were both validating and genuinely informative.

The ANOVA tests confirmed strong statistical separation between the climate conditions in Solar-recommended, Wind-recommended, and Hydro-recommended locations. Temperature showed the largest F-statistic among all variables — the difference in average temperature between Solar and Hydro locations was highly significant, consistent with the physical expectation that hot dry climates favor solar while wet climates favor hydro. Precipitation showed the second-largest separation. Wind speed showed the most distinctive separation for Wind-recommended locations, with a clear bimodal pattern: Wind recommendations clustered at wind speeds well above the global median, while Solar and Hydro recommendations were much more evenly distributed at lower wind speeds.

The correlation analysis revealed that temperature and precipitation are negatively correlated across real-world climate data — hotter places tend to be drier, colder places tend to be wetter. This structural relationship in real climate data has a direct consequence for the scoring algorithm: Solar and Hydro recommendations tend to be naturally separated by this correlation, making the solar-versus-hydro distinction easier to make than the wind-versus-either distinction. Wind speed has a weak and inconsistent correlation with both temperature and precipitation, which means Wind recommendations can appear across a wide range of temperature and precipitation values — they are primarily defined by the elevation and wind speed signals rather than temperature or precipitation.

The feature importance analysis, using permutation-based assessment where each input variable is randomly shuffled and the resulting decrease in classification consistency is measured, ranked temperature as the most important variable overall, followed closely by precipitation, then wind speed, then elevation, then snow depth. This ordering is intuitive: temperature and precipitation between them define most of the solar-versus-hydro boundary, wind speed defines most of the wind-versus-neither signal, and elevation and snow depth contribute important refinements — elevation particularly for the wind signal and snow depth particularly for the mountain-hydro signal.

The demo presets encode the most instructive results in an accessible form. A Desert preset — 35°C temperature, 4 m/s wind, 10mm precipitation, 300m elevation, no snow — produces Solar at 9 out of 9 possible points with high confidence. Every solar scoring condition fires: high temperature earns 5 points, low precipitation earns 2 points, low snow earns 1 point, moderate elevation earns 1 point. Wind gets 3 points (wind speed just barely clears the 4 m/s threshold) and Hydro gets 0 (no water supply). This is exactly the expected result for Sahara-type conditions. A Mountain preset — 5°C, 7 m/s wind, 120mm precipitation, 2000m elevation, 30cm snow — produces Wind at 5 out of 7 with low confidence, with Hydro close behind at 4 out of 8. The tool correctly flags this as a close call. This is exactly the expected result for alpine conditions where wind turbines and run-of-river hydro both have genuine merit. A Monsoon preset — 26°C, 5 m/s wind, 200mm precipitation, 500m elevation, no snow — produces Hydro at 6 out of 8 with moderate confidence, while Solar is penalized for the heavy precipitation that implies persistent cloud cover.

The finding that matters most from a design perspective: the scoring algorithm produced results that are not just statistically separable but also consistently intuitive. Users who try the presets and think through why each recommendation makes sense invariably find the reasoning transparent and convincing. That combination — statistically grounded AND intuitively comprehensible — is what makes the tool genuinely useful rather than just formally correct.

---

## What This Project Ultimately Demonstrates

This project is simultaneously a practical tool and an argument about methodology.

The practical dimension: it is a real, working, deployed web application that gives anyone a transparent, scored, evidence-based renewable energy recommendation for any climate on Earth, for free, in under ten seconds. That is something that did not exist before this project was built, and it fills a genuine gap between academic energy research and public accessibility.

The methodological argument: in domains where the decision logic is already well-understood from scientific research, a rule-based system with explicit, interpretable scoring can be more appropriate and more trustworthy than a machine learning model trained on labels. The ML approach sounds more sophisticated. But "more sophisticated" and "more reliable" are not the same thing, especially when the training data is scarce, the domain knowledge is rich, and the people using the tool need to understand and trust the reasoning behind the recommendation. The 3,000 lines of R behind this project were written in part to validate that argument empirically — the ANOVA and feature importance results confirm that the scoring thresholds the algorithm uses genuinely correspond to statistically meaningful separations in real-world climate data.

The two layers of the project — the accessible web app and the rigorous R analysis suite — reflect a belief that good data science should serve both audiences: the non-specialist who needs a clear answer and the expert who needs to see the methodology and verify the results.

---

## Running It

To use the web app, visit [weather-energy.netlify.app](https://weather-energy.netlify.app/) or clone the repository and open `index.html` in any browser. No installation required.

To run the R analysis suite:

```r
# Install dependencies
install.packages(c("jsonlite", "ggplot2", "tidyr", "dplyr", "readr", "purrr"))

# Run the classification engine on data.csv
source("classify_energy.R")

# Run statistical analysis
source("analysis.R")

# Generate all visualizations
source("visualizations.R")
```

To fetch real NOAA climate data for a specific region, set your NOAA CDO API token in `data_fetch.R` and call `fetch_region_data()` with a geographic bounding box.

The project structure:

```
Weather-Classification/
├── index.html           # Complete web application — works in any browser
├── classify_energy.R    # Main R engine: scoring, batch classification, export (~800 lines)
├── analysis.R           # Statistical analysis: ANOVA, correlation, regression (~600 lines)
├── visualizations.R     # ggplot2 visualization library (~700 lines)
├── data_fetch.R         # NOAA CDO API integration with caching (~500 lines)
├── utils.R              # Shared utilities: validation, conversion, logging (~500 lines)
├── data.csv             # Real NOAA station data with classifications
└── README.md            # This document
```

---

*Built by Sahibjeet Pal Singh — [GitHub](https://github.com/Sahibjeetpalsingh) · [Live App](https://weather-energy.netlify.app/) · [LinkedIn](https://linkedin.com/in/sahibjeet-pal-singh-418824333)*
