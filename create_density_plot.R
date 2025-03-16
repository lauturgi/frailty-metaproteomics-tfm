# Load ggplot2
library(ggplot2)

# Density plot function 
create_density_plot <- function(data, variable, title = "Density plot") {
  # Create the density plot
  p <- ggplot(data, aes(x = .data[[variable]])) +
    geom_density(fill = "blue", alpha = 0.4) +
    labs(title = title, x = variable, y = "Density") +
    theme_minimal()
  
  # Print the plot
  print(p)
}