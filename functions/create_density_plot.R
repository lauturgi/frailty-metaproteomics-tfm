# Load ggplot2
library(ggplot2)

# Density plot function 
create_density_plot <- function(data, variable, x_lab = variable, 
                                title = "Density plot") {
  # Create the density plot
  p <- ggplot(data, aes(x = .data[[variable]])) +
    geom_density(fill = "blue", alpha = 0.3) +
    labs(title = title, x = x_lab, y = "Density") +
    theme_minimal(base_size = 12)
  
  # Print the plot
  print(p)
}
