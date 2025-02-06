create_pca_plot <- function(se_perseus, n_top_loadings = 5) {
  # Load necessary libraries
  library(ggplot2)
  library(ggfortify)
  library(dplyr)
  library(ggrepel)
  
  # Perform PCA
  pca_result <- prcomp(t(assay(se_perseus)), scale. = TRUE)
  pca_scores <- as.data.frame(pca_result$x) %>% cbind(colData(se_perseus))
  
  # Get top contributing loadings
  loadings <- as.data.frame(pca_result$rotation)
  top_loadings <- loadings %>% mutate(contrib = abs(PC1) + abs(PC2)) %>% 
    top_n(n_top_loadings, contrib)
  
  # Scale loadings for better visualization
  arrow_scale <- max(abs(pca_scores$PC1)) / max(abs(loadings$PC1)) * 0.2  
  
  # Create the PCA plot
  p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = frailty)) +
    geom_point(size = 2, alpha = 0.5) +
    stat_ellipse(aes(fill = frailty), geom = "polygon", alpha = 0.2, color = NA) +
    geom_segment(data = top_loadings, aes(x = 0, y = 0, xend = PC1 * arrow_scale, 
                                          yend = PC2 * arrow_scale),
                 arrow = arrow(length = unit(0.2, "cm")), color = "gray50", 
                 size = 0.6, alpha = 0.5) +
    geom_text_repel(data = top_loadings, aes(x = PC1 * arrow_scale, 
                                             y = PC2 * arrow_scale, 
                                             label = rownames(top_loadings)),
                    size = 3, color = "black", max.overlaps = Inf) +
    labs(x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
         y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)"),
         title = paste0("PCA with Top ", n_top_loadings, " Contributing Loadings")) +
    theme_minimal()
  
  # Return the plot
  return(p)
}

