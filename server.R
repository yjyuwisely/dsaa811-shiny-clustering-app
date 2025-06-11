
# DSAA811 Final Exam Task 1: Hierarchical Clustering with NCI60 Data
# Author: Yeongjin Yu
# server.R - Server Logic

# Server logic
server <- function(input, output, session) {
  
  # Task 1 server logic ----
  
  # Process data to select most variable genes
  processed_data <- reactive({
    # Transpose data for clustering by samples (cell lines)
    transposed_data <- t(nci60_data)
    
    # Calculate variance for each gene
    gene_variance <- apply(transposed_data, 2, var)
    
    # Select top n most variable genes
    top_genes <- order(gene_variance, decreasing = TRUE)[1:input$num_genes]
    
    # Return the selected genes
    return(transposed_data[, top_genes])
  })
  
  # Perform hierarchical clustering
  hclust_result <- reactive({
    # Calculate distance matrix
    dist_matrix <- dist(processed_data(), method = input$distance)
    
    # Perform hierarchical clustering
    hclust(dist_matrix, method = input$linkage)
  })
  
  # Create dendrogram
  dendrogram <- reactive({
    # Convert hclust to dendrogram
    dend <- as.dendrogram(hclust_result())
    
    if(input$color_labels) {
      # Get labels
      labs <- nci60_labels$x
      
      # Create a color palette based on cancer types
      cancer_types <- unique(labs)
      
      # Use a predefined color palette for better color distinction
      cancer_colors <- brewer.pal(min(9, length(cancer_types)), input$color_palette)
      
      # If we have more than 9 cancer types, recycle colors
      if(length(cancer_types) > 9) {
        cancer_colors <- colorRampPalette(cancer_colors)(length(cancer_types))
      }
      
      names(cancer_colors) <- cancer_types
      
      # Color the labels based on cancer type
      label_colors <- cancer_colors[labs]
      
      # Apply colors to dendrogram
      dend <- set(dend, "labels_col", label_colors)
      dend <- set(dend, "labels_cex", 0.8)  # Adjust label size
    }
    
    return(dend)
  })
  
  # Cut the dendrogram to get clusters
  clusters <- reactive({
    # Calculate the actual height to cut based on slider percentage
    max_height <- max(hclust_result()$height)
    cut_height <- (input$cut_height / 100) * max_height
    
    # Cut the dendrogram at this height
    cutree(hclust_result(), h = cut_height)
  })
  
  # Render number of clusters (moved here to avoid reactive context issues)
  output$num_clusters <- renderText({
    paste("Number of clusters:", length(unique(clusters())))
  })
  
  # 1. Render dendrogram - MODIFIED to put legend outside
  output$dendrogram <- renderPlot({
    # Create layout: 75% for plot, 25% for legend
    layout(matrix(c(1,2), nrow=1), widths=c(0.75, 0.25))
    
    # Plot dendrogram in first panel
    par(mar = c(5, 4, 2, 1))  # Reduce right margin since legend is separate
    
    # Plot the dendrogram
    plot(dendrogram(), 
         main = "Hierarchical Clustering Dendrogram",
         xlab = "Cancer Cell Lines", 
         ylab = "Height (Dissimilarity)",
         horiz = FALSE)
    
    # Add a horizontal line at the cut height
    max_height <- max(hclust_result()$height)
    cut_height <- (input$cut_height / 100) * max_height
    abline(h = cut_height, col = "red", lty = 2, lwd = 2)
    
    # Create legend in second panel (completely separate)
    if(input$color_labels) {
      par(mar = c(5, 0, 2, 2))  # No left margin for legend panel
      plot.new()
      
      labs <- nci60_labels$x
      cancer_types <- unique(labs)
      cancer_colors <- brewer.pal(min(9, length(cancer_types)), input$color_palette)
      if(length(cancer_types) > 9) {
        cancer_colors <- colorRampPalette(cancer_colors)(length(cancer_types))
      }
      names(cancer_colors) <- cancer_types
      
      legend("left", 
             legend = cancer_types,
             fill = cancer_colors,
             title = "Cancer Types",
             cex = 0.8,
             bty = "n")  # No box around legend
    } else {
      # If colors are disabled, create empty plot for consistent layout
      par(mar = c(5, 0, 2, 2))
      plot.new()
    }
  })
  
  # Render heatmap
  output$heatmap <- renderPlot({
    par(mar = c(5, 4, 4, 8))
    # Get processed data
    data_matrix <- as.matrix(processed_data())
    
    # Get cluster order from hierarchical clustering
    cluster_order <- hclust_result()$order
    
    # Reorder data according to clustering
    ordered_data <- data_matrix[cluster_order, ]
    
    # Scale the data for better visualization (by row)
    scaled_data <- t(scale(t(ordered_data)))
    
    # Create a color palette for the heatmap
    heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
    
    # Get cancer labels for row annotation
    cancer_labels <- nci60_labels$x
    
    # Create a color palette for cancer types
    cancer_types <- unique(cancer_labels)
    cancer_colors <- brewer.pal(min(9, length(cancer_types)), input$color_palette)
    if(length(cancer_types) > 9) {
      cancer_colors <- colorRampPalette(cancer_colors)(length(cancer_types))
    }
    names(cancer_colors) <- cancer_types
    
    # Create row annotation colors
    row_colors <- cancer_colors[cancer_labels[cluster_order]]
    
    # 2. Plot the heatmap
    heatmap(scaled_data, 
            Rowv = NA,
            Colv = NA,
            col = heatmap_colors,
            scale = "none",
            labRow = NA,
            labCol = NA,
            margins = c(5, 2),  # Reduce left margin from 5 to 2
            main = "Gene Expression Heatmap",
            xlab = paste("Top", input$num_genes, "Variable Genes"),
            ylab = "Cancer Cell Lines (ordered by clustering)",
            RowSideColors = row_colors)
    
    legend("right", 
           legend = cancer_types,
           fill = cancer_colors,
           title = "Cancer Types",
           cex = 0.7,
           xpd = TRUE
    )
  })
  
  # Render cluster information
  output$cluster_info <- renderDT({
    # Get cluster assignments
    cluster_assignments <- clusters()
    
    # Get cell line names from the processed data
    cell_lines <- rownames(processed_data())
    
    # Get cancer types from the labels file - ensure it's the right length
    cancer_types <- nci60_labels$x[1:length(cell_lines)]
    
    # Create data frame with cell line, cancer type, and cluster
    cluster_data <- data.frame(
      CellLine = cell_lines,
      CancerType = cancer_types,
      Cluster = cluster_assignments,
      stringsAsFactors = FALSE
    )
    
    # 3. Count cancer types per cluster
    cluster_summary <- cluster_data %>%
      group_by(Cluster, CancerType) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      pivot_wider(names_from = CancerType, values_from = Count, values_fill = list(Count = 0)) %>%
      mutate(TotalSize = rowSums(across(where(is.numeric)))) %>%
      arrange(Cluster)
    
    # 4. Format for display
    datatable(cluster_summary, 
              options = list(
                pageLength = 10,
                autoWidth = TRUE,
                dom = 'tip',  # table, information, pagination
                scrollX = TRUE
              ),
              rownames = FALSE,
              caption = "Distribution of Cancer Types Across Clusters") %>%
      formatStyle('TotalSize', 
                  background = styleColorBar(range(cluster_summary$TotalSize), 'lightblue'),
                  fontWeight = 'bold')
  })
  
  # Reset button logic
  observeEvent(input$reset, {
    # Reset all inputs to their default values
    updateSelectInput(session, "distance", selected = "euclidean")
    updateSelectInput(session, "linkage", selected = "complete")
    updateSliderInput(session, "num_genes", value = 100)
    updateSliderInput(session, "cut_height", value = 80)  # Fixed: was 50, should match UI default
    updateCheckboxInput(session, "color_labels", value = TRUE)
    updateSelectInput(session, "color_palette", selected = "Spectral")
  })
  
  
  # Task 2 server logic ----
  # K-means from scratch implementation
  kmeans_scratch <- reactive({
    input$run_kmeans  # Dependency on button click
    
    # Your K-means algorithm implementation here
    # This will be the main algorithm you need to write
  })
  
  # Task 2 outputs
  output$kmeans_plot <- renderPlot({
    # K-means visualization plot
    plot(1:10, 1:10, main = "K-means Plot - To be implemented")
  })
  
  output$kmeans_progress <- renderText({
    "Algorithm progress will be shown here"
  })
  
  output$current_iteration <- renderText({
    "Current iteration: 0"
  })
  
  output$convergence_status <- renderText({
    "Status: Ready to run"
  })
  
  output$elbow_plot <- renderPlot({
    # Elbow method plot
    plot(1:10, 10:1, main = "Elbow Method - To be implemented")
  })
  
  output$wcss_table <- renderDT({
    # WCSS table
    data.frame(K = 2:8, WCSS = c(100, 80, 60, 50, 45, 42, 40))
  })
  
  # Reset button for Task 2
  observeEvent(input$reset_task2, {
    updateSliderInput(session, "k_clusters", value = 3)
    updateSliderInput(session, "max_iterations", value = 50)
    updateSliderInput(session, "iteration_step", value = 0)
    updateCheckboxInput(session, "show_centers", value = TRUE)
  })
}
