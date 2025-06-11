# DSAA811 Final Exam (Autumn 2025)
# Task 1: Hierarchical Clustering (NCI60 Gene Expression Data)
# Task 2: K-means Clustering (WCGS Heart Disease Data)
# Author: Yeongjin Yu
# global.R - Loads libraries and data

# Task 1 ----
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

library(factoextra)  # For cluster validation plots
library(animation)   # For step-by-step visualization

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

# Task 2 ----
# Load WCGS (Heart) dataset
tryCatch({
  heart_data <- read.csv("wcgs.csv", stringsAsFactors = TRUE)
  cat("Heart dataset loaded successfully:", nrow(heart_data), "rows\n")
}, error = function(e) {
  warning("Could not load wcgs.csv: ", e$message)
  # Create dummy data if file doesn't exist
  heart_data <- data.frame(
    age = rnorm(100, 45, 10),
    sbp = rnorm(100, 140, 20),
    dbp = rnorm(100, 85, 15),
    chol = rnorm(100, 200, 40),
    weight = rnorm(100, 170, 30)
  )
})

# Select only numeric columns for K-means clustering
heart_numeric <- heart_data[, sapply(heart_data, is.numeric)]
heart_numeric <- na.omit(heart_numeric)  # Remove rows with NA

# Scale the data for fair clustering (K-means is sensitive to scale)
heart_numeric_scaled <- scale(heart_numeric)

# K-means clustering method option for UI
clustering_methods <- c(
  "K-means (from scratch)" = "kmeans_scratch"
)
k_values <- 2:10  # Allowed values for K

# Data quality checks
if (nrow(heart_data) != nrow(heart_numeric)) {
  cat("Note:", nrow(heart_data) - nrow(heart_numeric), "rows with missing values were removed\n")
}

cat("Data loading complete:\n")
cat("- NCI60 data:", nrow(nci60_data), "cell lines,", ncol(nci60_data), "genes\n")
cat("- Heart data:", nrow(heart_numeric), "observations,", ncol(heart_numeric), "variables\n")