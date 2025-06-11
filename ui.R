
# DSAA811 Final Exam Task 1: Hierarchical Clustering with NCI60 Data
# Author: Yeongjin Yu
# ui.R - User Interface

# UI Component
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "DSAA811 Final Exam: Clustering Analysis"),
  
  dashboardSidebar(
    width = 300,
    
    sidebarMenu(
      id = "tabs",
      menuItem("Task 1: Hierarchical Clustering", tabName = "task1", icon = icon("sitemap")),
      menuItem("Task 2: K-means Clustering", tabName = "task2", icon = icon("chart-scatter")),
      menuItem("Educational Content", tabName = "education", icon = icon("book")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    
    # Task 1 conditionalPanel ---- 
    conditionalPanel(
      condition = "input.tabs == 'task1'",
      
      # Clustering method controls
      h4("Clustering Parameters", style = "padding-left: 15px;"),
      
      selectInput("distance", "Distance Method:", 
                  choices = distance_metrics,
                  selected = "euclidean"),
      
      selectInput("linkage", "Linkage Method:", 
                  choices = linkage_methods,
                  selected = "complete"),
      
      # Data filtering
      h4("Data Selection", style = "padding-left: 15px;"),
      
      sliderInput("num_genes", "Number of Most Variable Genes:", 
                  min = 50, max = 1000, value = 100, step = 50),
      
      # Dendrogram controls
      h4("Dendrogram Options", style = "padding-left: 15px;"),
      
      sliderInput("cut_height", "Cut Height for Clusters:", 
                  min = 0, max = 100, value = 80),
      
      checkboxInput("color_labels", "Color Labels by Cancer Type", value = TRUE),
      
      # Display options
      h4("Visualization Options", style = "padding-left: 15px;"),
      
      selectInput("color_palette", "Color Palette:", 
                  choices = color_palette_options,
                  selected = "Spectral"),
      
      # Reset button
      hr(),
      actionButton("reset", "Reset to Defaults", 
                   icon = icon("sync"), 
                   style = "color: #fff; background-color: #337ab7; border-color: #2e6da4; width: 90%;")
    ),
    
    # Task 2 conditionalPanel ---- 
    conditionalPanel(
      condition = "input.tabs == 'task2'",
      
      h4("K-means Parameters", style = "padding-left: 15px;"),
      
      sliderInput("k_clusters", "Number of Clusters (K):", 
                  min = 2, max = 8, value = 3, step = 1),
      
      sliderInput("max_iterations", "Maximum Iterations:", 
                  min = 10, max = 100, value = 50, step = 10),
      
      h4("Visualization Controls", style = "padding-left: 15px;"),
      
      sliderInput("iteration_step", "Show Iteration Step:", 
                  min = 0, max = 50, value = 0, step = 1),
      
      checkboxInput("show_centers", "Show Cluster Centers", value = TRUE),
      
      actionButton("run_kmeans", "Run K-means Algorithm", 
                   icon = icon("play"), 
                   style = "color: #fff; background-color: #28a745; border-color: #28a745; width: 90%; margin-bottom: 10px;"),
      
      hr(),
      actionButton("reset_task2", "Reset Task 2", 
                   icon = icon("sync"), 
                   style = "color: #fff; background-color: #337ab7; border-color: #2e6da4; width: 90%;")
    )
  ),
  
  dashboardBody(
    tabItems(
      # Task 1 tab content ---- 
      tabItem(tabName = "task1",
              fluidRow(
                box(
                  title = "Dendrogram Visualization", status = "primary", solidHeader = TRUE,
                  width = 12, height = 580,
                  plotOutput("dendrogram", height = "450px"),
                  footer = "Cut height shown as red dashed line. Colors indicate cancer types."
                )
              ),
              
              fluidRow(
                box(
                  title = "Heatmap of Gene Expression", status = "primary", solidHeader = TRUE,
                  width = 12, height = 500,
                  plotOutput("heatmap", height = "450px")
                )
              ),
              
              fluidRow(
                box(
                  title = "Number of Clusters", status = "warning", solidHeader = TRUE,
                  width = 12,
                  textOutput("num_clusters")
                )
              ),
              
              fluidRow(
                box(
                  title = "Cluster Information", status = "info", solidHeader = TRUE,
                  width = 12, height = "auto",
                  DTOutput("cluster_info", width = "100%")
                )
              )
      ),
      
      # Task 2 tab content ---- 
      tabItem(tabName = "task2",
              fluidRow(
                box(
                  title = "K-means Clustering Visualization", status = "primary", solidHeader = TRUE,
                  width = 8, height = 620,
                  plotOutput("kmeans_plot", height = "500px"),
                  footer = "Interactive K-means clustering showing step-by-step algorithm execution"
                ),
                
                box(
                  title = "Algorithm Progress", status = "info", solidHeader = TRUE,
                  width = 4, height = 620,
                  verbatimTextOutput("kmeans_progress"),
                  hr(),
                  textOutput("current_iteration"),
                  br(),
                  textOutput("convergence_status")
                )
              ),
              
              fluidRow(
                box(
                  title = "Elbow Method for Optimal K", status = "warning", solidHeader = TRUE,
                  width = 6, height = 400,
                  plotOutput("elbow_plot", height = "350px")
                ),
                
                box(
                  title = "Within-Cluster Sum of Squares", status = "success", solidHeader = TRUE,
                  width = 6, height = 400,
                  DTOutput("wcss_table")
                )
              )
      ),
      
      # Educational Content tab ---- 
      tabItem(tabName = "education",
              fluidRow(
                box(
                  title = "Understanding Hierarchical Clustering", status = "primary", solidHeader = TRUE,
                  width = 12,
                  
                  h3("What is Hierarchical Clustering?"),
                  p("Hierarchical clustering is an algorithm that groups similar objects into groups called clusters. The endpoint is a set of clusters where each cluster is distinct from the other clusters, and the objects within each cluster are broadly similar to each other."),
                  
                  h3("Types of Hierarchical Clustering"),
                  p(strong("Agglomerative (bottom-up):"), " Starts with each observation as its own cluster and merges them until all are in a single cluster."),
                  p(strong("Divisive (top-down):"), " Starts with all observations in one cluster and recursively splits them."),
                  
                  h3("Key Parameters"),
                  
                  h4("Distance Metrics:"),
                  tags$ul(
                    tags$li(strong("Euclidean:"), " Straight-line distance between two points. Good for continuous data."),
                    tags$li(strong("Manhattan:"), " Sum of absolute differences. Less sensitive to outliers."),
                    tags$li(strong("Maximum:"), " Maximum difference across dimensions. Sensitive to outliers."),
                    tags$li(strong("Canberra:"), " Sum of fractions of differences. Good for sparse data with non-negative values.")
                  ),
                  
                  h4("Linkage Methods:"),
                  tags$ul(
                    tags$li(strong("Complete:"), " Maximum distance between any two points in the clusters."),
                    tags$li(strong("Average:"), " Average distance between all pairs of points."),
                    tags$li(strong("Single:"), " Minimum distance between any two points."),
                    tags$li(strong("Ward's:"), " Minimizes variance within clusters. Often produces compact clusters.")
                  ),
                  
                  h3("Interpreting a Dendrogram"),
                  p("A dendrogram is a tree diagram that shows the hierarchical relationship between objects:"),
                  tags$ul(
                    tags$li("The y-axis represents the distance or dissimilarity between clusters."),
                    tags$li("Each vertical line represents a cluster or observation."),
                    tags$li("Horizontal lines connect clusters that are merged."),
                    tags$li("The height of horizontal lines indicates the distance at which clusters are merged.")
                  ),
                  p("By cutting the dendrogram at a specific height, we can determine the number of clusters.")
                )
              ),
              
              fluidRow(
                box(
                  title = "Hierarchical Clustering Step-by-Step", status = "info", solidHeader = TRUE,
                  width = 12,
                  
                  h3("Algorithm Steps (Agglomerative Clustering)"),
                  tags$ol(
                    tags$li("Start with n observations as n individual clusters."),
                    tags$li("Calculate the distance matrix between all pairs of clusters."),
                    tags$li("Find the closest pair of clusters and merge them."),
                    tags$li("Update the distance matrix to reflect the merged cluster."),
                    tags$li("Repeat steps 3-4 until all observations are in one cluster.")
                  ),
                  
                  h3("Advantages and Limitations"),
                  h4("Advantages:"),
                  tags$ul(
                    tags$li("No need to specify the number of clusters in advance."),
                    tags$li("The dendrogram provides a visually interpretable result."),
                    tags$li("Can capture nested cluster structure.")
                  ),
                  
                  h4("Limitations:"),
                  tags$ul(
                    tags$li("Computationally intensive for large datasets (O(n²) or O(n³) complexity)."),
                    tags$li("Sensitive to noise and outliers."),
                    tags$li("Cannot revisit once a merge/split decision is made (greedy algorithm).")
                  )
                )
              )
      ),
      
      # About tab content ---- 
      tabItem(tabName = "about",
              fluidRow(
                box(
                  title = "About the NCI60 Dataset", status = "info", solidHeader = TRUE,
                  width = 12,
                  
                  h3("NCI60 Cancer Cell Line Panel"),
                  p("The NCI60 dataset contains gene expression data from 60 diverse human cancer cell lines used by the National Cancer Institute to screen for new anti-cancer compounds."),
                  
                  h4("Dataset Composition:"),
                  tags$ul(
                    tags$li(strong("Samples:"), " 64 cancer cell lines from 9 different tissue origins"),
                    tags$li(strong("Features:"), " Expression levels of 6,830 genes per cell line"),
                    tags$li(strong("Cancer Types:"), " Includes melanoma, leukemia, lung, colon, central nervous system (CNS), ovarian, renal, prostate, and breast cancer cell lines")
                  ),
                  
                  h3("Application in This Shiny App"),
                  p("In this application, we use hierarchical clustering to:"),
                  tags$ul(
                    tags$li("Group similar cancer cell lines based on their gene expression profiles"),
                    tags$li("Identify patterns that might correlate with cancer types"),
                    tags$li("Explore how different clustering parameters affect these groupings"),
                    tags$li("Visualize gene expression patterns across different cancer types")
                  ),
                  
                  h3("About This Application"),
                  p("This Shiny app was developed as part of the DSAA811 final exam (Autumn 2025)."),
                  p("Author: Yeongjin Yu")
                )
              )
      )
    )
  )
)
