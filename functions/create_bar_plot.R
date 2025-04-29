# Load ggplot2
library(ggplot2)

# Bar plot function
create_bar_plot <- function(data, variable, title = "Bar Plot") {
  # Create the bar plot
  p <- ggplot(data, aes(x = .data[[variable]], fill = .data[[variable]])) +
    geom_bar() +
    labs(title = title, x = variable, y = "Count") +
    theme_minimal() +
    geom_text(stat = "count", aes(label = after_stat(count)),
              vjust = -0.25)
  
  # Print the plot
  print(p)
  
  return(p)

}
