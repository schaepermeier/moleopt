library(tidyverse)

# ========= moPLOT =========

design <- moPLOT::generateDesign(fn, points.per.dimension = 1000L)
design$obj.space <- moPLOT::calculateObjectiveValues(design$dec.space, fn, parallelize = TRUE)

gradients <- moPLOT::computeGradientFieldGrid(design, normalized.scale = FALSE)
divergence <- moPLOT::computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)

less <- moPLOT::localEfficientSetSkeleton(design, gradients, divergence, integration = "fast")

g <- moPLOT::ggplotPLOT(design$dec.space, design$obj.space, less$sinks, less$height) +
  coord_fixed()
g.obj <- moPLOT::ggplotPLOTObjSpace(design$obj.space, less$sinks, less$height)

# ========= MOGSA =========

d <- 2
fid <- 1
iid <- 1
biobj_bbob_data <- generateBiObjBBOBData(d, fid, iid)
fn <- biobj_bbob_data$fn

# fn <- smoof::makeWFG5Function(2, 1, 1)
# fn <- smoof::makeDTLZ1Function(dimensions = 2, n.objectives = 2)
# fn <- makeAsparFunction(2, 2)
# fn <- smoof::makeMultiObjectiveFunction("test", fn = function(x) c(sum(x ** 2), sum(x ** 2)),
#                                         par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-5, -5), upper = c(5, 5)))
# fn <- smoof::makeMMF4Function()
# fn <- makeBiObjMPM2Function()

lower <- smoof::getLowerBoxConstraints(fn)
upper <- smoof::getUpperBoxConstraints(fn)

nstarts <- 100
starting_points <- lapply(1:nstarts, function(x) runif_box(lower, upper))
starting_points <- do.call(rbind, starting_points)

obj.opt.path <- NULL
opt.path <- NULL

log_y <- FALSE
log_x <- FALSE
count_evals <- TRUE

nruns <- 1
run_counter <- 0

while (run_counter < nruns) {
  cat(paste0("New run, ID: ", run_counter, ", currently at: ", ncol(obj.opt.path), "\n"))
  
  f <- fn
  
  if (log_x || log_y) {
    f <- smoof::addLoggingWrapper(f, logg.x = log_x, logg.y = log_y)
  }
  
  if (count_evals) {
    f <- smoof::addCountingWrapper(f)
  }
  
  if (run_counter > 0) {
    starting_points <- lapply(1:nstarts, function(x) runif_box(lower, upper))
    starting_points <- do.call(rbind, starting_points)
  }
  
  mogsa_trace <- run_mogsa(f, starting_points,
                              eps_gradient = 1e-8,
                              eps_initial_step_size = 1e-8,
                              eps_explore_set = 1e-4,
                              # max_explore_set = 1e-3,
                              max_explore_set = sqrt(sum((upper - lower) ** 2)) / 100,
                              # custom_descent_fn = create_lbfgsb_descent(f, lower, upper),
                              # lower = rep(-Inf, length(lower)),
                              # upper = rep(Inf, length(upper)),
                              max_budget = Inf
                           )

  if (log_y) {
    if (is.null(obj.opt.path)) {
      obj.opt.path <- smoof::getLoggedValues(f)$obj.vals
    } else {
      obj.opt.path <- cbind(obj.opt.path, smoof::getLoggedValues(f)$obj.vals)
    }
  }
  
  if (log_x) {
    if (is.null(opt.path)) {
      opt.path <- smoof::getLoggedValues(f)$pars
    } else {
      opt.path <- rbind(opt.path, smoof::getLoggedValues(f)$pars)
    }
  }

  if (count_evals) {
    cat("No. Evals: ", paste0(smoof::getNumberOfEvaluations(f), "\n"))
  }
  
  run_counter <- run_counter + 1
}

# any(is.na(obj.opt.path))

mogsa_trace$sets %>% length

set_sizes <- sapply(mogsa_trace$sets, function(s) nrow(s$obj_space))
sum(set_sizes > 2)

sum(mogsa_trace$transitions[,1] != -1)

plot_dec_space(mogsa_trace, lower, upper, color = "set_id") +
  coord_fixed() +
  theme(legend.position = "none")

plot_obj_space(mogsa_trace) +
  theme(legend.position = "none")

plot_dec_space(mogsa_trace, lower, upper, color = "domcount") +
  coord_fixed() +
  theme(legend.position = "none")

# ggplot() +
#   geom_point(aes(x1, x2), data = as.data.frame(starting_points), shape = "+", color = "black", size = 10) +
#   coord_fixed()

# ggplot() +
#   geom_point(aes(x = x1, y = x2, color = 1:nrow(opt.path)), data = opt.path) +
#   xlim(lower[1], upper[1]) +
#   ylim(lower[2], upper[2]) +
#   coord_fixed()

# ggplot() +
#   geom_point(aes(x = V1, y = V2, color = 1:ncol(obj.opt.path)), data = as.data.frame(t(obj.opt.path)))

# ggplot() +
#   geom_path(aes(x = x1, y = x2), data = opt.path)

if (is.null(obj.opt.path)) {
  obj.opt.path <- t(do.call(rbind, lapply(mogsa_trace$sets, function(l) l$obj_space)))
}

ncol(obj.opt.path)

nondom <- nondominated(t(obj.opt.path))

obj.opt.path <- as.data.frame(t(obj.opt.path))

colnames(obj.opt.path) <- c("y1", "y2")

nrow(obj.opt.path) / d

sum(nondom) / nrow(obj.opt.path)
sum(duplicated(obj.opt.path)) / nrow(obj.opt.path)
# sum(duplicated(opt.path)) / nrow(opt.path)

# ref_point <- apply(obj.opt.path[nondom,1:2], 2, max)
# ideal_point <- apply(obj.opt.path[nondom,1:2], 2, min)

ref_point <- biobj_bbob_data$ref_point
ideal_point <- biobj_bbob_data$ideal_point

hv.norm <- prod(ref_point - ideal_point)
hv <- ecr::computeHV(t(as.matrix(obj.opt.path[nondom,1:2])), ref.point = ref_point)
hv / hv.norm

# log10(5/6 - hv / hv.norm)
log10(0.922987888165046 - hv / hv.norm)

ggplot() +
  geom_point(aes(x = y1, y = y2), data = obj.opt.path[!nondom,], color = "red") +
  geom_point(aes(x = y1, y = y2), data = obj.opt.path[nondom,], color = "green") +
  geom_point(aes(x = ref_point[1], y = ref_point[2]), shape = "+", size = 5) +
  geom_point(aes(x = ideal_point[1], y = ideal_point[2]), shape = "+", size = 5) +
  xlim(ideal_point[1], ref_point[1]) +
  ylim(ideal_point[2], ref_point[2])

# ggplot() +
#   geom_point(aes(x = y1, y = y2, color = 1:nrow(obj.opt.path)), data = obj.opt.path)

# Set transition graph

library(tidygraph)
library(ggraph)

set_transitions <- mogsa_trace$transitions[apply(mogsa_trace$transitions, 1, function(x) all(x >= 0)),,drop=FALSE] + 1
nodes <- mogsa_trace$transitions[,2] %>% unique + 1

colnames(set_transitions) <- c("from", "to")
tbl_transitions <- tbl_graph(edges = as.data.frame(set_transitions), nodes = data.frame(name = nodes))

weights <- (mogsa_trace$transitions[mogsa_trace$transitions[,1]==-1, 2] + 1) %>% tabulate(nbins = max(mogsa_trace$transitions + 1))
prop <- compute_reach_proportions(tbl_transitions, weights)

set_nd_counts <- compute_nondominated_sets(mogsa_trace$sets)
node_color <- ifelse(set_nd_counts > 0, "green", "red")

(1 / prop[set_nd_counts > 0]) %>% sort(decreasing = TRUE)

ggraph(tbl_transitions, layout = "stress") +
  geom_edge_fan(arrow = arrow(length = unit(4, "mm")),
                end_cap = circle(4, "mm")) +
  geom_node_point(aes(size = prop), color = node_color) +
  scale_size_area(limits = c(0,1))

## Find Local Search Traps

condensation <- igraph::contract(tbl_transitions, igraph::components(tbl_transitions, mode = "strong")$membership, vertex.attr.comb = "ignore")
condensation <- igraph::simplify(condensation) # in particular: remove self-loops
sum(igraph::degree(condensation, mode = "out") == 0)

### 2D Decision Spaces ###

node_pos <- set_medians(mogsa_trace$sets)

ggraph(tbl_transitions, layout = "manual", x = node_pos[,1], y = node_pos[,2]) +
  geom_edge_fan(arrow = arrow(length = unit(4, "mm")),
                end_cap = circle(4, "mm")) +
  geom_node_point(aes(size = prop), color = node_color) +
  scale_size_area(limits = c(0,1)) +
  coord_fixed() +
  theme(legend.position = "none")



### Plotly 3D Sets ###

dfs <- lapply(seq_along(mogsa_trace$sets), function(i) {
  data_dec <- as.data.frame(mogsa_trace$sets[[i]]$dec_space)
  colnames(data_dec) <- c("x1", "x2", "x3")
  
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

plotly::plot_ly(data = df_combined,
                type = "scatter3d", mode = "markers",
                x = ~x1, y = ~x2, z = ~x3,
                color = ~factor(set_id))

nds <- ecr::doNondominatedSorting(t(df_combined_obj[,1:2]))
plotly::plot_ly(data = df_combined,
                type = "scatter3d", mode = "markers",
                x = ~x1, y = ~x2, z = ~x3,
                color = log(nds$dom.counter + 1))
