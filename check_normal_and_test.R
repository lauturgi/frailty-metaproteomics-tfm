check_normal_and_test <- function(df, var_col, group_col) {
  t_ft <- shapiro.test(df[[var_col]][df[[group_col]] == "FT"])
  t_nft <- shapiro.test(df[[var_col]][df[[group_col]] == "NFT"])
  
  p_val_ft <- t_ft$p.value
  p_val_nft <- t_nft$p.value
  
  if (p_val_nft > 0.05 & p_val_ft > 0.05) {
    t_res <- t.test(df[[var_col]] ~ df[[group_col]], data = df,
                    var.equal = FALSE)
    return(t_res)
  } else {
    w_res <- wilcox.test(df[[var_col]] ~ df[[group_col]], data = df)
    return(w_res)
  }
}
