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
  lapply(sets, function(set) {
    dec_space <- set$dec_space
    dec_space[ceiling(nrow(dec_space) / 2),]
  })
}
