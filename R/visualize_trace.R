#' @export
plot_dec_space <- function(mole_trace, lower, upper, color = "set_id") {
  dfs <- lapply(seq_along(mole_trace$sets), function(i) {
    data_dec <- as.data.frame(mole_trace$sets[[i]]$dec_space)
    colnames(data_dec) <- c("x1", "x2")
    
    data_dec[,"set_id"] <- i
    
    data_dec
  })
  
  dfs_obj <- lapply(seq_along(mole_trace$sets), function(i) {
    data_obj <- as.data.frame(mole_trace$sets[[i]]$obj_space)
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
  
  df_combined <- df_combined[order(data_color, decreasing = TRUE),]

  g <- ggplot(data = df_combined, aes(x = x1, y = x2, color = sort(data_color, decreasing = TRUE))) +
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
plot_obj_space <- function(mole_trace, color = "set_id") {
  dfs <- lapply(seq_along(mole_trace$sets), function(i) {
    data_obj <- as.data.frame(mole_trace$sets[[i]]$obj_space)
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
  
  df_combined <- df_combined[order(data_color, decreasing = TRUE),]
  
  g <- ggplot(data = df_combined, aes(x = y1, y = y2, color = sort(data_color, decreasing = TRUE))) +
    geom_point()
  
  if (color == "domcount") {
    g + scale_color_gradientn(colors = fields::tim.colors(500L))
  } else {
    g
  }
}

#' @export
plot_set_interactions <- function(mole_trace, layout_2d = FALSE, lower = NULL, upper = NULL) {
  set_transitions <- mole_trace$transitions[apply(mole_trace$transitions, 1, function(x) all(x >= 0)),,drop=FALSE] + 1
  nodes <- mole_trace$transitions[,2] %>% unique + 1
  
  colnames(set_transitions) <- c("from", "to")
  tbl_transitions <- tbl_graph(edges = as.data.frame(set_transitions), nodes = data.frame(name = nodes))
  
  weights <- (mole_trace$transitions[mole_trace$transitions[,1]==-1, 2] + 1) %>% tabulate(nbins = max(mole_trace$transitions + 1))
  prop <- compute_reach_proportions(tbl_transitions, weights)
  
  set_nd_counts <- compute_nondominated_sets(mole_trace$sets)
  node_color <- ifelse(set_nd_counts > 0, "darkgreen", "magenta")
  
  if (layout_2d) {
    node_pos <- set_medians(mole_trace$sets)
    
    ggraph(tbl_transitions, layout = "manual", x = node_pos[,1], y = node_pos[,2]) +
      geom_node_point(aes(size = 1), color = "black", shape = 21) +
      geom_node_point(aes(size = prop), color = node_color) +
      geom_edge_fan(arrow = arrow(length = unit(4, "mm")),
                    end_cap = circle(4, "mm")) +
      scale_size_area(limits = c(0,1)) +
      theme_minimal() +
      coord_fixed(xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2])) +
      theme(legend.position = "none",
            panel.background = element_rect(fill = NA, size = 0),
            plot.background = element_rect(fill = NA, size = 0)) +
      labs(x = expression(x[1]),
           y = expression(x[2]))
  } else {
    ggraph(tbl_transitions, layout = "stress") +
      geom_node_point(aes(size = 1), color = "black", shape = 21) +
      geom_node_point(aes(size = prop), color = node_color) +
      geom_edge_fan(arrow = arrow(length = unit(4, "mm")),
                    end_cap = circle(4, "mm")) +
      scale_size_area(limits = c(0,1)) +
      theme(#legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))
  }
}


