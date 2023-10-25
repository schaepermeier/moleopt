compute_reach_proportions <- function(graph, weights) {
  sapply(igraph::V(graph), function(v) {
    reachable_v <- igraph::subcomponent(graph, v, "in")
    sum(weights[reachable_v]) / sum(weights)
  })
}

compute_nondominated_sets <- function(sets) {
  obj_space <- Reduce(rbind, lapply(sets, function(set) set$obj_space))
  
  # When do we change between different sets of the trace?
  set_change <- cumsum(sapply(sets, function(set) nrow(set$obj_space)))
  
  nd <- obj_space %>% nondominated()
  
  sapply(seq_along(set_change), function(i) {
    if (i == 1) {
      lower = 1
    } else {
      lower = set_change[i - 1] + 1
    }
    
    sum(nd[lower:set_change[i]])
  })
}

set_medians <- function(sets) {
  sapply(sets, function(set) {
    dec_space <- set$dec_space
    dec_space[ceiling(nrow(dec_space) / 2),]
  }) %>% t
}

compute_statistics <- function(mole_trace) {
  stat_list <- list()
  
  set_transitions <- mole_trace$transitions[apply(mole_trace$transitions, 1, function(x) all(x >= 0)),,drop=FALSE] + 1
  nodes <- mole_trace$transitions[,2] %>% unique + 1
  
  set_sizes <- sapply(mole_trace$sets, function(s) nrow(s$obj_space))

  colnames(set_transitions) <- c("from", "to")
  
  tbl_transitions <- tbl_graph(edges = as.data.frame(set_transitions), nodes = data.frame(name = nodes))
  
  weights <- (mole_trace$transitions[mole_trace$transitions[,1]==-1, 2] + 1) %>% tabulate(nbins = max(mole_trace$transitions + 1))
  prop <- compute_reach_proportions(tbl_transitions, weights)
  
  accepted_sets <- which(prop > 1 / sum(weights))
  
  set_nd_counts <- compute_nondominated_sets(mole_trace$sets)
  nondom_set_ids <- which(set_nd_counts[accepted_sets] > 0)
  
  stat_list$nsets <- length(accepted_sets)
  stat_list$globsets <- length(nondom_set_ids)
  stat_list$reachability <- min(prop[accepted_sets][nondom_set_ids])
  
  stat_list
}

#' Get reachability for set containing xopt
#'
#' @param mole_trace 
#' @param xopt 
#'
#' @return
#' @export
#'
#' @examples
compute_set_reachability <- function(mole_trace, xopt) {
  set_transitions <- mole_trace$transitions[apply(mole_trace$transitions, 1, function(x) all(x >= 0)),,drop=FALSE] + 1
  nodes <- mole_trace$transitions[,2] %>% unique + 1
  
  distance_to_xopt <- sapply(mole_trace$sets, function(s) {
    sqrt(min(apply(s$dec_space, 1, function(p) sum((p - xopt)**2))))
  })
  
  xopt_sets <- distance_to_xopt < 1e-2
  
  colnames(set_transitions) <- c("from", "to")
  
  tbl_transitions <- tbl_graph(edges = as.data.frame(set_transitions), nodes = data.frame(name = nodes))
  
  weights <- (mole_trace$transitions[mole_trace$transitions[,1]==-1, 2] + 1) %>% tabulate(nbins = max(mole_trace$transitions + 1))
  prop <- compute_reach_proportions(tbl_transitions, weights)
  
  accepted_sets <- which(prop > 1 / sum(weights))
  
  set_nd_counts <- compute_nondominated_sets(mole_trace$sets)
  nondom_set_ids <- which(set_nd_counts[accepted_sets] > 0)
  
  list(
    set_reachability = prop,
    set_nd_counts = set_nd_counts,
    nondom_set_ids = nondom_set_ids,
    xopt_sets = xopt_sets
  )
}
