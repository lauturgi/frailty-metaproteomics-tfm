library(ggplot2)
library(ggrepel)

create_volcano_plot <- function(diff, sig = 0.05, fold = 0.7, 
                                title = "Volcano plot", 
                                save_path = "volcano_plot.png") {
  # Add a column for differential expression labels
  diff$diff.exp <- rep("Not significant", nrow(diff))
  diff$diff.exp[(diff$logFC > fold & diff$adj.P.Val < sig)] <- "Up"
  diff$diff.exp[(diff$logFC < -fold & diff$adj.P.Val < sig)] <- "Down"
  diff$diff.exp[(abs(diff$logFC) <= fold &
                   diff$adj.P.Val < sig)] <- "Significant"
  
  diff$diff.exp <- as.factor(diff$diff.exp)
  
  # Create a column for labeling significant proteins
  diff$label <- ifelse(diff$adj.P.Val < sig, rownames(diff), NA)
  
  # Generate the plot
  g <- ggplot(data = diff, aes(x = logFC, y = -log10(adj.P.Val),
                               colour = diff.exp)) +
    geom_point(alpha = 0.8, size = 1) +
    scale_colour_manual(values = c("Not significant" = "gray", 
                                   "Significant" = "green", 
                                   "Up" = "red", 
                                   "Down" = "blue")) +
    coord_cartesian(
      xlim = c(-4, 4),  #c(-1, 1) * (max(abs(diff$logFC), na.rm = TRUE) + 1),
      ylim = c(0, 2)  #c(0, max(-log10(diff$adj.P.Val), na.rm = TRUE) + 0.5)
    ) +
    xlab("logFC") +
    ylab("-log10 adjusted p-value") +
    labs(title = title) +
    theme_minimal(base_size = 12) +
    theme(legend.title = element_blank()) + 
    geom_vline(xintercept = c(-fold, fold), colour = "black", size = 0.2) +
    geom_hline(yintercept = -log10(sig), colour = "black", size = 0.2) +
    geom_text_repel(aes(label = label), size = 2.5, na.rm = TRUE, color = "black")
  
  # Print the plot
  print(g)
  
  # Save the plot
  ggsave(filename = save_path, plot = g, width = 8, height = 6, dpi = 300)
}
