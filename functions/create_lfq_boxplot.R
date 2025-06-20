library(ggplot2)
library(reshape2)

# Define LFQ boxplot function
create_lfq_boxplot <- function(se_assay, col_data, group_col = "frailty",
                               sample_col = "sample", 
                               fill_colors = c("FT" = "red", 
                                               "NFT" = "purple"),
                               fill_labels = c("FT" = "Frail",
                                               "NFT" = "Non-Frail"),
                               lfq_type = "maxlfq",
                               y_lab = "Protein MaxLFQ intensities") {
  
  # Transform se matrix to long dataframe
  se_df <- as.data.frame(se_assay)
  se_df$protein_id <- rownames(se_df)
  rownames(se_df) <- NULL
  se_long <- melt(se_df)
  colnames(se_long) <- c("protein_id", "sample", lfq_type)
  
  # Merge with col_data to add sample metadata
  se_long <- as.data.frame(merge(se_long, col_data[, c(sample_col,
                                                       group_col)],
                                 by = sample_col))
  # Summary
  lfq_summary <- summary(se_long[[lfq_type]])
  
  # SD
  lfq_sd <- sd(se_long[[lfq_type]], na.rm = TRUE)
  
  # Plot intensity vs samples
  plot_1 <- ggplot(se_long, aes(x = .data[[sample_col]],
                                y = .data[[lfq_type]],
                                fill = .data[[group_col]])) + 
    geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
    labs(x = "Samples", y = y_lab) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(), #element_text(angle = 90, size = 6),
      legend.position = "top",
      legend.title = element_blank()) + 
    scale_fill_manual(values = fill_colors,
                      labels = fill_labels) +
    facet_wrap(~.data[[group_col]], scales = "free_x", ncol = 1)
  
  # Plot intensity vs groups
  plot_2 <- ggplot(se_long, aes(x = .data[[group_col]],
                                y = .data[[lfq_type]],
                                fill = .data[[group_col]])) + 
    geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
    labs(x = "Frailty group", y = y_lab) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 8),
      legend.position = "top",
      legend.title = element_blank()) + 
    scale_fill_manual(values = fill_colors,
                      labels = fill_labels)
  
  # Return
  return(list(se_long = se_long, lfq_summary = lfq_summary, lfq_sd = lfq_sd,
              plot_samples = plot_1, plot_groups = plot_2))
}
