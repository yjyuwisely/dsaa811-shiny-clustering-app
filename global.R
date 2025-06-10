
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

# Read the cancer type labels file - skip the header row
nci60_labels <- read.csv("NCI60_labs.csv", skip = 1)

# Extract just the cancer types (column 3)
cancer_types <- nci60_labels[, 3]

# Create a clean labels data frame
nci60_labels <- data.frame(
  cancer_type = cancer_types,
  stringsAsFactors = TRUE
)

# Make sure we have exactly the right number of labels
if(nrow(nci60_labels) != nrow(nci60_data)) {
  warning(paste("Number of labels (", nrow(nci60_labels), 
                ") doesn't match number of data rows (", nrow(nci60_data), ")"))
  
  # Trim to match data dimensions
  nci60_labels <- nci60_labels[1:nrow(nci60_data), , drop = FALSE]
}

# Print confirmation of data dimensions
message("NCI60 data: ", nrow(nci60_data), " rows × ", ncol(nci60_data), " columns")
message("Cancer labels: ", nrow(nci60_labels), " rows")

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
