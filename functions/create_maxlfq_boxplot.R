library(ggplot2)
library(reshape2)

# Define function
create_maxlfq_boxplot <- function(se, col_data, group_col = "frailty",
                                  sample_col = "sample", 
                                  fill_colors = c("FT" = "red", 
                                                  "NFT" = "purple")) {
  # Transform se matrix to long dataframe
  se_assay <- assay(se)
  se_df <- as.data.frame(se_assay)
  se_df$protein_id <- rownames(se_df)
  rownames(se_df) <- NULL
  se_long <- melt(se_df)
  colnames(se_long) <- c("protein_id", "Samples", "MaxLFQ")
  
  # Extract first 7 characters of Samples for matching
  se_long$Samples <- substr(se_long$Samples, 1, 7)
  
  # Merge with col_data to add sample metadata
  se_long <- as.data.frame(merge(se_long, col_data[, c(sample_col,
                                                       group_col)],
                                 by.x = "Samples", by.y = sample_col))
  # Summary of MaxLFQ
  maxlfq_summary <- summary(se_long$MaxLFQ)
  
  # Plot MaxLFQ vs samples
  plot_1 <- ggplot(se_long, aes(x = Samples, y = MaxLFQ, fill = .data[[group_col]])) + 
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
    theme(
      axis.text.x = element_text(angle = 90, size = 6),
      legend.position = c(0.95, 0.95),
      legend.background = element_rect(fill = alpha("white", 0.5))) + 
    scale_fill_manual(values = fill_colors) +
    facet_wrap(~.data[[group_col]], scales = "free_x", ncol = 1)
  
  # Plot MaxLFQ vs groups
  plot_2 <- ggplot(se_long, aes(x = .data[[group_col]], y = MaxLFQ, fill = .data[[group_col]])) + 
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
    theme(
      axis.text.x = element_text(angle = 90, size = 6),
      legend.position = c(0.95, 0.95),
      legend.background = element_rect(fill = alpha("white", 0.5))) + 
    scale_fill_manual(values = fill_colors)
  
  # Return plots
  return(list(se_long = se_long, maxlfq_summary = maxlfq_summary,
              plot_samples = plot_1, plot_groups = plot_2))
}

