#' @export
plot_dec_space <- function(mogsa_trace, lower, upper, color = "set_id") {
  dfs <- lapply(seq_along(mogsa_trace$sets), function(i) {
    data_dec <- as.data.frame(mogsa_trace$sets[[i]]$dec_space)
    colnames(data_dec) <- c("x1", "x2")
    
    data_dec[,"set_id"] <- i
    
    data_dec
  })
  
  dfs_obj <- lapply(seq_along(mogsa_trace$sets), function(i) {
    data_obj <- as.data.frame(mogsa_trace$sets[[i]]$obj_space)
    colnames(data_obj) <- c("y1", "y2")

    data_obj
  })
  
  df_combined <- do.call(rbind, dfs)
  df_combined_obj <- do.call(rbind, dfs_obj)
  
  if (color == "set_id") {
    data_color <- factor(df_combined$set_id)
  } else if (color == "domcount") {
    nds <- ecr::doNondominatedSorting(t(df_combined_obj[,1:2]))
    data_color <- log(nds$dom.counter + 1)
  }

  g <- ggplot(data = df_combined, aes(x = x1, y = x2, color = data_color)) +
    geom_point() +
    xlim(lower[1], upper[1]) +
    ylim(lower[2], upper[2])
  
  if (color == "domcount") {
    g + scale_color_gradientn(colors = fields::tim.colors(500L))
  } else {
    g
  }
}

#' @export
plot_obj_space <- function(mogsa_trace, color = "set_id") {
  dfs <- lapply(seq_along(mogsa_trace$sets), function(i) {
    data_obj <- as.data.frame(mogsa_trace$sets[[i]]$obj_space)
    colnames(data_obj) <- c("y1", "y2")
    
    data_obj[,"set_id"] <- i
    
    data_obj
  })
  
  df_combined <- do.call(rbind, dfs)
  
  if (color == "set_id") {
    data_color <- factor(df_combined$set_id)
  } else if (color == "domcount") {
    nds <- ecr::doNondominatedSorting(t(df_combined[,1:2]))
    data_color <- log(nds$dom.counter + 1)
  }
  
  g <- ggplot(data = df_combined, aes(x = y1, y = y2, color = data_color)) +
    geom_point()
  
  if (color == "domcount") {
    g + scale_color_gradientn(colors = fields::tim.colors(500L))
  } else {
    g
  }
}
