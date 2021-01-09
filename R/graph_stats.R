compute_reach_proportions <- function(graph, weights) {
  sapply(igraph::V(graph), function(v) {
    reachable_v <- igraph::subcomponent(graph, v, "in")
    sum(weights[reachable_v]) / sum(weights)
  })
}