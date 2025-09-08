# Scenario-Based Projection Tool for Pandemic Spread

## Overview

This is the repository for the scenario-based projection tool for the spread of the pandemic. This tool simulates COVID-19 infections in Germany (2020–2021) using weekly age- and sex-specific incidence data from RKI SurvStat and a simple microsimulation model (`Susceptible → Infected → Recovered`). The tool is developed within the Leibniz Lab for Pandemic Preparedness by the team at DIW Berlin (Alex Lepe and Sabine Zinn).

The app is also available as an online demo via [GitHub Pages](https://bieneschwarze.github.io/ProjectionPandemicIncidence/). Runtime may be slow, especially with larger sample sizes, so for optimal performance we recommend running the app locally by following the instructions below.

## Repository contents

-   app.R - Shiny tool for the simulation

-   00_prepping_rki_data.qmd - Code to load and clean the RKI data

-   00_prepping_rki_data.html - Rendered HTML of the quarto file

-   docs/ - Folder for shinylive deployment

-   data/

    -   SurvStat_complete_data - All the processed data in both CSV and RDS formats
    -   SurvStat_female_data - Processed data for females in both CSV and RDS formats
    -   SurvStat_male_data - Processed data for males in both CSV and RDS formats
    -   weekly - Folder with all the raw weekly data files from SurvStat\@RKI 2.0

-   shinycovid-simulation.png – Example output (screenshot)

-   README.md - Description of the project

-   ProjectionPandemicIncidence.Rproj - R project

    **Note:** The app reads the `.rds` files, so ensure these are available

## Running the tool

### Prerequisites

Ensure the following packages are installed: `shiny`, `shinyBS`, `MicSim`, `data.table`, and `ggplot2`

### Running the tool

From the project root, run:

``` r
shiny::runApp("app.R")
```

Adjust inputs in the left-hand panel and click **Run simulation**. **Note:** larger values of **N** (cohort size) increase runtime.

If you run the tool with **all default settings**, the output should match the screenshot ![Simulation output](shinycovid-simulation.png)

## How the projection tool works

1.  Import infection rates: Weekly SARS‑CoV‑2 infection rates for 2020–2021 are read from `data/weekly/` (source: Robert Koch Institute – SurvStat\@RKI 2.0).
2.  Generate a synthetic cohort: `N` individuals are sampled with uniform ages between the user‑chosen age range and the sample is evenly split by sex.
3.  Microsimulation (using MicSim): Each individual can move through the states `Susceptible → Infected → Recovered`
    -   Transition rates can be scaled by `1 – (intervention effectiveness / 100)` starting on the selected `intervention date`.
    -   Illness duration determines how long infections remain ‘active’.
    -   For more details on how MicSim works, see the [CRAN documentation](https://cran.r-project.org/web/packages/MicSim/index.html).
4.  Post‑processing: Calculate daily counts for
    -   `Infected (7-day total)`: 7-day rolling sum
    -   `Infected (active, user-set duration)`: rolling sum over your chosen illness days
    -   `Recovered (cumulative)`: total number of people recovered from illness
    -   `Susceptible (remaining)`: remaining number of people susceptible to infection
5.  Visualization: The Shiny UI plots the percentage of people actively infected and recovered per day.

## Inputs (left panel)

-   **Sample size (N):** 1,000–30,000

    *Note: Values above 30,000 can also be used, but the simulation will run more slowly.*

-   **Minimum/Maximum age at start:** age bounds for the initial cohort

-   **Duration of illness (days):** used to define recovery timing, which controls “active infection” window

-   **Intervention effectiveness (%):** percentage reduction in incidence

-   **Intervention starts:** date when the intervention begins

-   **Compare intervention against baseline:** run both a baseline scenario (without intervention) and an intervention scenario, and compare the results

-   **Break down by sex:** Display plots and indicators by sex (Female/Male)

-   **Show series** — select which time series to plot:

    -   *Infected (7-day total)*
    -   *Infected (active, user-set duration)*
    -   *Recovered (cumulative)*
    -   *Susceptible*

-   **Run simulation:** executes the model(s) with current settings

### Outputs

-   **Plot:** daily percentages according to the selected options
-   **Indicators:** peak infected, peak date, recovered/susceptible by end
-   **Downloads** - **CSV**: the daily values used to generate the plots - **PNG**: current plot
