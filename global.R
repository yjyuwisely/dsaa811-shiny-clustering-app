# DSAA811 Final Exam Task 1: Hierarchical Clustering with NCI60 Data
# Author: Yeongjin Yu
# global.R - Loads libraries and data

# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(dendextend)
library(heatmaply)
library(dplyr)
library(cluster)
library(RColorBrewer)
library(tidyr)  # For pivot_wider function

# Load the data (at startup)
# Read the gene expression data - rows are cell lines, columns are genes
nci60_data <- read.csv("NCI60_data.csv", row.names = 1)

# Read the cancer type labels - no header in the file
nci60_labels <- read.csv("NCI60_labs.csv", header = TRUE, stringsAsFactors = TRUE)

# Verify dimensions match
if(nrow(nci60_data) != nrow(nci60_labels)) {
  warning("Warning: Number of rows in data (", nrow(nci60_data), 
          ") doesn't match number of labels (", nrow(nci60_labels), ")")
}

# Define color palettes for use in both UI and server
color_palette_options <- c(
  "Spectral" = "Spectral",
  "RdYlBu" = "RdYlBu", 
  "YlOrRd" = "YlOrRd", 
  "YlGnBu" = "YlGnBu"
)

# Define distance metrics and their descriptions
distance_metrics <- c(
  "Euclidean" = "euclidean", 
  "Maximum" = "maximum", 
  "Manhattan" = "manhattan", 
  "Canberra" = "canberra"
)

# Define linkage methods and their descriptions
linkage_methods <- c(
  "Complete" = "complete", 
  "Average" = "average", 
  "Single" = "single", 
  "Ward's" = "ward.D2"
)