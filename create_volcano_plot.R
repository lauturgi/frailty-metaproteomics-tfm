library(ggplot2)
library(ggrepel)

create_volcano_plot <- function(diff, sig = 0.05, fold = 0.7, 
                                title = "Volcano plot", 
                                save_path = "volcano_plot.png") {
  # Add a column for differential expression labels
  diff$diff.exp <- rep("not diff.", nrow(diff))
  diff$diff.exp[diff$adj.P.Val < sig] <- "diff.FDR"   
  diff$diff.exp[(diff$logFC > fold & diff$adj.P.Val < sig) |
                  (diff$logFC < -fold & diff$adj.P.Val < sig)] <- "diff.logFC.FDR"
  
  diff$diff.exp <- as.factor(diff$diff.exp)
  
  # Create a column for labeling significant proteins
  diff$label <- ifelse(diff$adj.P.Val < sig, rownames(diff), NA)
  
  # Generate the plot
  g <- ggplot(data = diff, aes(x = logFC, y = -log10(adj.P.Val), colour = diff.exp)) +
    geom_point(alpha = 0.8, size = 0.8) +
    xlim(c(-8, 8)) + ylim(c(0, max(-log10(diff$adj.P.Val), na.rm = TRUE) + 1)) +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    labs(title = title) +
    theme(
      axis.title.x = element_text(colour = "black", size = 18),
      axis.title.y = element_text(colour = "black", size = 18),
      title = element_text(colour = "black", size = 18),
      axis.text.x = element_text(colour = "black", size = 18),
      axis.text.y = element_text(colour = "black", size = 18),
      legend.position = "none"
    ) +
    geom_vline(xintercept = c(-fold, fold), colour = "black", size = 0.2) +
    geom_hline(yintercept = -log10(sig), colour = "black", size = 0.2) +
    geom_text_repel(aes(label = label), size = 3, na.rm = TRUE, color = "black")
  
  # Print the plot
  print(g)
  
  # Save the plot
  ggsave(filename = save_path, plot = g, width = 8, height = 6, dpi = 300)
}
