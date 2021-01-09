#' @export
plot_dec_space <- function(mogsa_trace, lower, upper) {
  dfs <- lapply(seq_along(mogsa_trace$sets), function(i) {
    data_dec <- as.data.frame(mogsa_trace$sets[[i]]$dec_space)
    colnames(data_dec) <- c("x1", "x2")
    
    data_dec[,"set_id"] <- i
    
    data_dec
  })
  
  df_combined <- do.call(rbind, dfs)
  
  ggplot(data = df_combined, aes(x = x1, y = x2, color = factor(set_id))) +
    geom_point() +
    xlim(lower[1], upper[1]) +
    ylim(lower[2], upper[2])
}

#' @export
plot_obj_space <- function(mogsa_trace) {
  dfs <- lapply(seq_along(mogsa_trace$sets), function(i) {
    data_obj <- as.data.frame(mogsa_trace$sets[[i]]$obj_space)
    colnames(data_obj) <- c("y1", "y2")
    
    data_obj[,"set_id"] <- i
    
    data_obj
  })
  
  df_combined <- do.call(rbind, dfs)
  
  ggplot(data = df_combined, aes(x = y1, y = y2, color = factor(set_id))) +
    geom_point()
}
