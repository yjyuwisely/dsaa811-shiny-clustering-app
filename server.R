# DSAA811 Final Exam Task 1: Hierarchical Clustering with NCI60 Data
# Author: Yeongjin Yu
# server.R - Server Logic

# Server logic
server <- function(input, output, session) {
  # Process data to select most variable genes
  processed_data <- reactive({
    # Transpose data for clustering by samples (cell lines)
    transposed_data <- t(nci60_data)
    
    # Calculate variance for each gene
    gene_variance <- apply(transposed_data, 2, var)
    
    # Select top n most variable genes
    top_genes <- order(gene_variance, decreasing = TRUE)[1:input$num_genes]
    
    output$num_clusters <- renderText({
      length(unique(clusters()))
    })
    
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
  
  # 1. Render dendrogram
  output$dendrogram <- renderPlot({
    # Set plotting parameters: increase right margin to 14
    par(mar = c(5, 4, 2, 14))
    
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
    
    # Add a legend for cancer types if colors are enabled
    if(input$color_labels) {
      labs <- nci60_labels$x
      cancer_types <- unique(labs)
      cancer_colors <- brewer.pal(min(9, length(cancer_types)), input$color_palette)
      if(length(cancer_types) > 9) {
        cancer_colors <- colorRampPalette(cancer_colors)(length(cancer_types))
      }
      names(cancer_colors) <- cancer_types
      
      legend("bottomright",
             legend = cancer_types,
             fill = cancer_colors,
             title = "Cancer Types",
             cex = 0.7,
             xpd = TRUE,
             inset = c(-0.11, 0)) # Adjust the first value as needed for your screen
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
           xpd = TRUE,
           #inset = c(-0.11, 0)
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
    updateSliderInput(session, "cut_height", value = 50)
    updateCheckboxInput(session, "color_labels", value = TRUE)
    updateSelectInput(session, "color_palette", selected = "Spectral")
  })
}