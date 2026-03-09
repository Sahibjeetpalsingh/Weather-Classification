<div align="center">

# 🌍 Renewable Energy Classifier

[![Live Demo](https://img.shields.io/badge/Try%20It%20Live-weather--energy.netlify.app-brightgreen?style=for-the-badge&logo=netlify)](https://weather-energy.netlify.app/)
[![R](https://img.shields.io/badge/R-Analysis%20Suite-276DC3?style=for-the-badge&logo=r)](https://github.com/Sahibjeetpalsingh/Weather-Classification)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](LICENSE)

</div>

---

## The Question That Started Everything

It began with a surprisingly simple question that turned out to have no simple answer: if you were dropped anywhere on Earth and told you could only build one type of renewable energy source, how would you figure out which one to pick?

Not theoretically. Not "it depends." Actually, concretely, using the real climate of that place — the sun it gets, the wind that moves through it, the rain that falls on it, the elevation it sits at. The answer should be different for the Sahara versus the Scottish Highlands versus the Amazon basin versus rural British Columbia. And yet most people, including most non-specialist planners and even students studying sustainability, do not have a fast, honest, transparent way to get that answer.

That gap is what this project was built to close.

---

## Understanding the Problem

The deeper you look into renewable energy planning, the more you realize the decision is not "solar vs wind vs hydro" in the abstract — it is "solar vs wind vs hydro *here*, given *these* conditions, for *this* use case." Each energy type has a completely different set of environmental requirements, and getting those requirements wrong is expensive. Building solar infrastructure in a region that is cloudy 200 days a year is a waste. Putting wind turbines where average wind speeds never cross the viability threshold is money left on the ground. Designing hydro capacity where rainfall is seasonal and unpredictable creates a system that underperforms half the year.

The research phase of this project involved going deep into what those requirements actually are — not from opinion, but from the physical and economic science of energy generation. Solar panels need heat and clear skies. Photovoltaic efficiency drops with cloud cover, and panels covered in snow or in persistently rainy climates significantly underperform their rated output. Wind turbines have a hard minimum: below roughly 4 meters per second of average wind speed, commercial turbines do not generate economical power. Above 6 m/s, they move into their optimal range. Elevation matters because high plateaus channel faster, more consistent winds. Hydro is different again — it needs water, which means both precipitation and elevation change. Water has to fall to generate power, and the greater the fall the more power you get. Snowmelt is a huge part of this equation for mountain regions, because snow accumulation in winter creates a powerful seasonal water pulse in spring.

Once those physical rules were mapped out clearly, the question became: how do you turn them into a decision tool?

---

## The Thinking Behind the Solution

The first instinct for anyone with a machine learning background is to train a classifier. Feed it labeled examples — this location is solar, that one is wind — and let the model learn the boundaries. But that instinct runs into a wall quickly: there is no reliable public dataset of "correct" energy recommendations for thousands of global locations. And even if you built one, you would be training a black box that makes decisions you cannot explain or audit. For an energy planning tool used by non-specialists, explainability is not optional — if someone is about to make a real infrastructure investment, they need to understand *why* the recommendation says solar and not wind.

So instead of machine learning, a rule-based weighted scoring system was built, directly encoding the physical science into logic. Each of the three energy types gets scored independently across five dimensions: temperature, wind speed, precipitation, elevation, and snow depth. The scores are grounded in real thresholds from energy engineering research — not made up, not arbitrary. A region with temperatures consistently above 15°C starts earning solar points. A region above 25°C earns more. Heavy rainfall above 150mm penalizes solar because persistent cloud cover directly reduces panel output. Wind above 4 m/s earns wind points; above 6 m/s earns more. High precipitation and snowpack earn hydro points. Meaningful elevation — between 300m and 2000m — earns hydro points because water needs height to fall.

The energy type with the highest total score wins. The margin between the winner and the runner-up determines how confident that recommendation is. A narrow margin means the tool tells you honestly: this is a close call, consider a hybrid approach. A large margin means high confidence.

This design choice — explicit rules over learned models — turned out to be the right one, because it produces recommendations you can actually argue with. You can look at the score breakdown and say "I think precipitation should matter more here" and understand what to change. You cannot do that with a neural network.

---

## Building It

The project was built in two layers that reflect two different audiences. The front layer is a live interactive web application. The back layer is a rigorous R analysis suite that could sit behind a professional energy feasibility study.

The web app — available at weather-energy.netlify.app — is a pure HTML, CSS, and JavaScript application. No server. No frameworks. No dependency chain to manage. Open it in a browser and it works. The interface gives you five sliders corresponding to the five climate inputs, and eight preset buttons for common climate archetypes: Desert, Coastal, Mountain, Tropical, Arctic, Plains, Mediterranean, and Monsoon. As you move any slider, the scoring algorithm runs instantly in JavaScript and the results panel updates in real time — the recommended energy source, the confidence level, the scores for all three energy types displayed against their maximum possible scores, and a full breakdown of which conditions fired and how many points each one contributed. You can see exactly why the tool said Solar and not Wind.

The R analysis suite underneath is a different kind of tool entirely. It is over 3,000 lines of code organized into four modules. The classification engine in classify_energy.R implements the same scoring logic but at scale — it can process dataframes of thousands of locations, generate statistical summary reports, and export results to CSV and JSON. The analysis module in analysis.R runs the kind of statistical machinery you would find in a peer-reviewed study: ANOVA tests to check whether climate variables differ significantly across energy recommendation categories, correlation matrices across all input dimensions, regression models with multicollinearity checks, confusion matrix metrics, and permutation-based feature importance analysis to understand which inputs are actually driving the classifications. The visualizations module produces publication-quality charts in ggplot2 — boxplots, scatter maps, radar charts, correlation heatmaps, full dashboards. And the data fetching module connects directly to NOAA's Climate Data Online API to pull real historical climate data for any geographic region on Earth, with rate limiting, caching, and unit conversion handled automatically.

The data used in the project comes from NOAA's GHCN — the Global Historical Climatology Network — which aggregates climate records from weather stations worldwide.

---

## What the Data Actually Showed

Running the R analysis suite on real NOAA station data confirmed something that the scoring design assumed but had not empirically verified: the physical thresholds embedded in the scoring rules do in fact correspond to meaningful separations in real-world climate data.

When ANOVA tests were run across the three energy recommendation categories, temperature and precipitation came out as the dominant discriminators between Solar and Hydro regions. That makes physical sense — hot dry places score for solar, wet rainy places score for hydro, and the climate variables that distinguish those two are precisely temperature and precipitation. Wind speed and elevation came out as the primary drivers of Wind recommendations, again exactly consistent with the physical reality that wind energy is fundamentally about how fast air is moving at hub height, and elevation is one of the most reliable predictors of that.

The correlation analysis showed that temperature and precipitation are negatively correlated in the dataset — hotter places tend to be drier — which is why Solar and Hydro recommendations rarely compete closely with each other. They tend to be geographically separated. Wind competes more often with both, because wind speed is less correlated with temperature or precipitation and can appear in almost any climate type.

The preset demo results tell the story clearly across the eight climate archetypes. Desert conditions — high temperature, minimal rain, low elevation, negligible snow — push Solar to a perfect 9 out of 9 possible points, while Hydro earns zero because there is no water to work with. Mountain conditions — cold, windy, moderate precipitation, high elevation, deep snow — produce a close contest between Wind at 5/7 and Hydro at 4/8, with Wind narrowly winning but the tool correctly flagging it as low confidence and suggesting a hybrid approach makes sense. Monsoon conditions — warm, heavy rainfall, moderate elevation, no snow — push Hydro to 6/8 while Solar is penalized for the cloud cover that monsoon rainfall implies.

The finding that matters most is not any individual result — it is that a transparent rule-based system built from physical science can produce recommendations that are both intuitive to a non-specialist and defensible to an engineer. You do not need machine learning to make a good decision tool. You need a clear model of the problem.

---

## What Is in the Repository

The project lives at github.com/Sahibjeetpalsingh/Weather-Classification. The index.html file is the complete web application — download it and open it in any browser, or use the live version. The classify_energy.R file is the main R classification engine, roughly 800 lines. The analysis.R file handles statistical testing and model evaluation, roughly 600 lines. The visualizations.R file is the plotting library, roughly 700 lines. The data_fetch.R file is the NOAA API integration layer, roughly 500 lines. The utils.R file contains all helper functions — input validation, unit conversions, logging — roughly 500 lines. The data.csv file contains real classified location data generated from NOAA station records.

To run the web app you just need a browser. To run the R suite you need R installed with the packages jsonlite, ggplot2, tidyr, dplyr, readr, and purrr, then call source on the relevant file.

---

## The Bigger Picture

This project ended up being about more than renewable energy. It was about the question of when rules beat models, when explainability matters more than raw accuracy, and how to take domain knowledge that already exists in scientific literature and turn it into something a regular person can actually use. The scoring algorithm is not sophisticated by machine learning standards — it is deliberately simple. But that simplicity is the point. A homeowner can look at the score breakdown and immediately understand what it is saying. A city planner can explain it to their council. A student can see the logic and learn from it.

The live app has been visited by people using it to think about solar installations, compare climates across different cities, and explore how changing a single variable — adding more precipitation, reducing temperature — shifts the recommendation. That interactive quality, where cause and effect are instantly visible, turned out to be the most valuable thing about the project.

---

*Built by Sahibjeet Pal Singh — [GitHub](https://github.com/Sahibjeetpalsingh) · [Live App](https://weather-energy.netlify.app/)*
