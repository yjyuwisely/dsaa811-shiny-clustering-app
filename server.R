# DSAA811 Final Exam Task 1: Hierarchical Clustering with NCI60 Data
# and Task 2: K-means Clustering with wcgs Data
# Author: Yeongjin Yu
# server.R - Server Logic

server <- function(input, output, session) {
  
  ## ========== Task 1: Hierarchical Clustering ==========
  processed_data <- reactive({
    transposed_data <- t(nci60_data)
    gene_variance <- apply(transposed_data, 2, var)
    top_genes <- order(gene_variance, decreasing = TRUE)[1:input$num_genes]
    transposed_data[, top_genes]
  })
  
  hclust_result <- reactive({
    dist_matrix <- dist(processed_data(), method = input$distance)
    hclust(dist_matrix, method = input$linkage)
  })
  
  dendrogram <- reactive({
    dend <- as.dendrogram(hclust_result())
    if (input$color_labels) {
      labs <- nci60_labels$x
      cancer_types <- unique(labs)
      cancer_colors <- brewer.pal(min(9, length(cancer_types)), input$color_palette)
      if (length(cancer_types) > 9) {
        cancer_colors <- colorRampPalette(cancer_colors)(length(cancer_types))
      }
      names(cancer_colors) <- cancer_types
      label_colors <- cancer_colors[labs]
      dend <- set(dend, "labels_col", label_colors)
      dend <- set(dend, "labels_cex", 0.8)
    }
    dend
  })
  
  clusters <- reactive({
    max_height <- max(hclust_result()$height)
    cut_height <- (input$cut_height / 100) * max_height
    cutree(hclust_result(), h = cut_height)
  })
  
  output$num_clusters <- renderText({
    paste("Number of clusters:", length(unique(clusters())))
  })
  
  # Dendrogram with separate legend panel
  output$dendrogram <- renderPlot({
    layout(matrix(c(1,2), nrow=1), widths=c(0.75, 0.25))
    par(mar = c(5, 4, 2, 1))
    plot(dendrogram(), 
         main = "Hierarchical Clustering Dendrogram",
         xlab = "Cancer Cell Lines", 
         ylab = "Height (Dissimilarity)",
         horiz = FALSE)
    max_height <- max(hclust_result()$height)
    cut_height <- (input$cut_height / 100) * max_height
    abline(h = cut_height, col = "red", lty = 2, lwd = 2)
    # Separate legend panel
    if (input$color_labels) {
      par(mar = c(5, 0, 2, 2))
      plot.new()
      labs <- nci60_labels$x
      cancer_types <- unique(labs)
      cancer_colors <- brewer.pal(min(9, length(cancer_types)), input$color_palette)
      if (length(cancer_types) > 9) {
        cancer_colors <- colorRampPalette(cancer_colors)(length(cancer_types))
      }
      names(cancer_colors) <- cancer_types
      legend("left", 
             legend = cancer_types,
             fill = cancer_colors,
             title = "Cancer Types",
             cex = 0.8,
             bty = "n")
    } else {
      par(mar = c(5, 0, 2, 2))
      plot.new()
    }
  })
  
  # Heatmap (legend to the right)
  output$heatmap <- renderPlot({
    par(mar = c(5, 4, 4, 8))
    data_matrix <- as.matrix(processed_data())
    cluster_order <- hclust_result()$order
    ordered_data <- data_matrix[cluster_order, ]
    scaled_data <- t(scale(t(ordered_data)))
    heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
    cancer_labels <- nci60_labels$x
    cancer_types <- unique(cancer_labels)
    cancer_colors <- brewer.pal(min(9, length(cancer_types)), input$color_palette)
    if (length(cancer_types) > 9) {
      cancer_colors <- colorRampPalette(cancer_colors)(length(cancer_types))
    }
    names(cancer_colors) <- cancer_types
    row_colors <- cancer_colors[cancer_labels[cluster_order]]
    heatmap(scaled_data, 
            Rowv = NA,
            Colv = NA,
            col = heatmap_colors,
            scale = "none",
            labRow = NA,
            labCol = NA,
            margins = c(5, 2),
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
  
  output$cluster_info <- renderDT({
    cluster_assignments <- clusters()
    cell_lines <- rownames(processed_data())
    cancer_types <- nci60_labels$x[1:length(cell_lines)]
    cluster_data <- data.frame(
      CellLine = cell_lines,
      CancerType = cancer_types,
      Cluster = cluster_assignments,
      stringsAsFactors = FALSE
    )
    cluster_summary <- cluster_data %>%
      group_by(Cluster, CancerType) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      pivot_wider(names_from = CancerType, values_from = Count, values_fill = list(Count = 0)) %>%
      mutate(TotalSize = rowSums(across(where(is.numeric)))) %>%
      arrange(Cluster)
    datatable(cluster_summary, 
              options = list(
                pageLength = 10,
                autoWidth = TRUE,
                dom = 'tip',
                scrollX = TRUE
              ),
              rownames = FALSE,
              caption = "Distribution of Cancer Types Across Clusters") %>%
      formatStyle('TotalSize', 
                  background = styleColorBar(range(cluster_summary$TotalSize), 'lightblue'),
                  fontWeight = 'bold')
  })
  
  observeEvent(input$reset, {
    updateSelectInput(session, "distance", selected = "euclidean")
    updateSelectInput(session, "linkage", selected = "complete")
    updateSliderInput(session, "num_genes", value = 100)
    updateSliderInput(session, "cut_height", value = 80)
    updateCheckboxInput(session, "color_labels", value = TRUE)
    updateSelectInput(session, "color_palette", selected = "Spectral")
  })
  
  ## ========== Task 2: K-means (from scratch) ==========
  
  kmeans_result <- eventReactive(input$run_kmeans, {
    set.seed(123)
    K <- input$k_clusters
    max_iter <- input$max_iterations
    dat <- as.matrix(heart_numeric_scaled)
    n <- nrow(dat)
    p <- ncol(dat)
    
    # Randomly initialise centroids
    centroids <- dat[sample(1:n, K), , drop = FALSE]
    
    cluster_history <- list()
    center_history <- list()
    for (iter in 1:max_iter) {
      # Assign clusters
      dists <- as.matrix(dist(rbind(centroids, dat)))[1:K, (K+1):(K+n)]
      clusters <- apply(dists, 2, which.min)
      cluster_history[[iter]] <- clusters
      center_history[[iter]] <- centroids
      # Update centroids
      new_centroids <- t(sapply(1:K, function(k) {
        if (sum(clusters == k) == 0) {
          centroids[k, ]  # Prevent empty clusters
        } else {
          colMeans(dat[clusters == k, , drop = FALSE])
        }
      }))
      # Check for convergence
      if (all(abs(new_centroids - centroids) < 1e-6)) break
      centroids <- new_centroids
    }
    list(
      clusters = cluster_history,
      centers = center_history,
      n_iter = length(cluster_history),
      K = K,
      n = n
    )
  })
  
  output$kmeans_plot <- renderPlot({
    res <- kmeans_result()
    if (is.null(res)) return()
    iter <- input$iteration_step
    if (iter < 1) iter <- 1
    if (iter > res$n_iter) iter <- res$n_iter
    df <- as.data.frame(heart_numeric_scaled)
    xcol <- 1
    ycol <- 2
    plot(df[, xcol], df[, ycol], col = res$clusters[[iter]], pch = 19,
         xlab = colnames(df)[xcol], ylab = colnames(df)[ycol],
         main = paste("K-means Clustering (Iteration", iter, ")"))
    if (input$show_centers) {
      points(res$centers[[iter]][, xcol], res$centers[[iter]][, ycol],
             col = 1:res$K, pch = 8, cex = 2, lwd = 2)
    }
  })
  
  output$kmeans_progress <- renderText({
    res <- kmeans_result()
    if (is.null(res)) return("Not run yet.")
    paste("Total iterations completed:", res$n_iter)
  })
  
  output$current_iteration <- renderText({
    paste("Current iteration:", input$iteration_step)
  })
  
  output$convergence_status <- renderText({
    res <- kmeans_result()
    if (is.null(res)) return("Status: Not run yet.")
    if (input$iteration_step == res$n_iter)
      "Status: Converged"
    else
      "Status: Still running"
  })
  
  # Elbow plot
  output$elbow_plot <- renderPlot({
    dat <- as.matrix(heart_numeric_scaled)
    wss <- numeric()
    for (k in 2:8) {
      set.seed(123)
      clust <- kmeans(dat, centers = k, nstart = 10, iter.max = 100)
      wss[k] <- clust$tot.withinss
    }
    plot(2:8, wss[2:8], type = "b", xlab = "Number of Clusters (K)",
         ylab = "Within-Cluster Sum of Squares (WCSS)",
         main = "Elbow Method for Choosing K")
  })
  
  # WCSS table
  output$wcss_table <- renderDT({
    dat <- as.matrix(heart_numeric_scaled)
    WCSS <- sapply(2:8, function(k) {
      set.seed(123)
      clust <- kmeans(dat, centers = k, nstart = 10, iter.max = 100)
      clust$tot.withinss
    })
    data.frame(K = 2:8, WCSS = WCSS)
  })
  
  # Reset button for Task 2
  observeEvent(input$reset_task2, {
    updateSliderInput(session, "k_clusters", value = 3)
    updateSliderInput(session, "max_iterations", value = 50)
    updateSliderInput(session, "iteration_step", value = 0)
    updateCheckboxInput(session, "show_centers", value = TRUE)
  })
}