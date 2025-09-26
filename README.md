# DSAA811 – Clustering Visualisation Dashboard with R Shiny

**Subject:** DSAA811 – Data Analytics and Visualisation (Enhanced)  
**Author:** Yeongjin Yu  
**Date:** June 2025  

This repository presents an interactive R Shiny dashboard built for the DSAA811 final exam. It demonstrates hierarchical and K-means clustering techniques using two real-world biomedical datasets: the NCI60 cancer gene expression dataset** and the WCGS heart disease dataset.

The app allows users to experiment with clustering parameters, visualise clustering steps, and interpret the results interactively — combining both educational content and analytical functionality.

---

## Project Overview

### Task 1: Hierarchical Clustering (NCI60 Dataset)
- **Dataset:** Gene expression data from 60 human cancer cell lines (6,830 genes per line)
- **Goal:** Explore similarity among cancer types based on molecular profiles
- **Visualisations:** Dendrogram, heatmap, and cancer-type summary table
- **Controls:** Distance metric, linkage method, gene variance filtering, cut height, palette

### Task 2: K-means Clustering (WCGS Dataset)
- **Dataset:** Clinical and behavioural data from 3,154 men (Western Collaborative Group Study)
- **Goal:** Identify patient groups with distinct heart disease risk profiles
- **Features:** Age, blood pressure, cholesterol, smoking, etc.
- **Visualisations:** Step-by-step K-means animation, WCSS table, elbow plot
- **Algorithm:** Implemented K-means from scratch with interactive iteration steps

---

## Files Included

| File               | Description                                      |
|--------------------|--------------------------------------------------|
| `global.R`         | Loads libraries, datasets, and constants         |
| `ui.R`             | R Shiny user interface definition                |
| `server.R`         | Server-side logic for clustering and visualisation |
| `NCI60_data.csv`   | Gene expression dataset (NCI60 cancer cell lines) |
| `NCI60_labs.csv`   | Cancer type labels for NCI60 data                |
| `wcgs.csv`         | WCGS health and behavioural dataset              |
| `README.md`        | This summary file                                |

---

## Technologies & Packages

- **Language:** R (v4.3+)
- **Framework:** R Shiny (with `shinydashboard`)
- **Key Packages:**
  - `shiny`, `shinydashboard`, `ggplot2`, `dendextend`, `heatmaply`, `factoextra`, `cluster`, `DT`, `RColorBrewer`, `tidyverse`, `animation`

---

## Features

- Dynamic control of clustering parameters (distance metric, linkage, number of genes/clusters)
- Interactive dendrogram and heatmap generation
- Step-by-step K-means visualisation with iteration control
- Educational explanations for clustering theory and algorithms
- Colour-coded cancer types and summary tables
- Elbow plot and WCSS table for selecting optimal K

---

## Running the App

To run this dashboard locally:

```r
# Clone the repository
git clone https://github.com/yjyuwisely/dsaa811-clustering-shiny-app.git
setwd("dsaa811-clustering-shiny-app")

# Launch the Shiny app
shiny::runApp()
```

Install required packages if not already installed:

```r
install.packages(c(
  "shiny", "shinydashboard", "DT", "ggplot2", "dendextend",
  "heatmaply", "dplyr", "cluster", "RColorBrewer", 
  "tidyr", "factoextra", "animation"
))
```

---

## Application Use Cases

| Task    | Use Case                                                                 |
|---------|--------------------------------------------------------------------------|
| Task 1  | Discover cancer subtypes with similar gene expression for targeted therapy |
| Task 2  | Cluster patient profiles based on cardiovascular risk indicators         |


---
