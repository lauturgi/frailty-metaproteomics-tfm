# Load libraries
library(tibble)
library(dplyr)

# Filter function
filter_valids = function(df, conditions, min_count, at_least_one = TRUE) {
  # df = dataframe containing LOG2 for filtering and organized by data type
  # conditions =character vector dictating the grouping
  # min_count = numeric vector (same length conditions) indicating the minimum 
  #             number of valid values for each condition for retention
  # at_least_one = TRUE to keep the row if min_count is met for at least 
  #                one condition
  #                FALSE to keep the row if min_count is met across all
  #                conditions
  df[df==0] <- NA
  col_names <- colnames(df) 
  cond_names = lapply(conditions,
                      function(x) grep(paste0("^", x, "_"),
                                       col_names, value = TRUE))
  cond_filter = sapply(1:length(cond_names), function(i) {
    df2 = df[cond_names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for following command
    sums = rowSums(is.finite(df2)) # Count number of valid values per condition
    sums >= min_count[i] # Calculates whether min_count is met
  })
  if (at_least_one) {
    df$KEEP = apply(cond_filter, 1, any)
  } else {
    df$KEEP = apply(cond_filter, 1, all)
  }
  
  df[is.na(df)] <- 0
  return(df %>%
    rownames_to_column(., var='peptides') %>%
    filter(KEEP) %>%
    select(-KEEP) %>%
    column_to_rownames(., var='peptides') # only keep rows that meet criteria
  )
}