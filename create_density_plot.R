# Load ggplot2
library(ggplot2)

# Density plot function 
create_density_plot <- function(data, variable, title = "Density plot", save_path = paste0("density_plot_", variable, ".png")) {
  # Create the density plot
  p <- ggplot(data, aes(x = .data[[variable]])) +
    geom_density(fill = "blue", alpha = 0.4) +
    labs(title = title, x = variable, y = "Density") +
    theme_minimal()
  
  # Print the plot
  print(p)
  
  # Save the plot
  ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)
}