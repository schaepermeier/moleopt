library(tidyverse)

# ========= moPLOT =========

design <- moPLOT::generateDesign(fn, points.per.dimension = 1001L)
design$obj.space <- moPLOT::calculateObjectiveValues(design$dec.space, fn, parallelize = TRUE)

gradients <- moPLOT::computeGradientFieldGrid(design)
divergence <- moPLOT::computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)

less <- moPLOT::localEfficientSetSkeleton(design, gradients, divergence, integration = "fast")

g <- moPLOT::ggplotPLOT(design$dec.space, design$obj.space, less$sinks, less$height)
g.obj <- moPLOT::ggplotPLOTObjSpace(design$obj.space, less$sinks, less$height)

rescaled <- sapply(1:nrow(gradients$multi.objective), function(i) {
  gradients$multi.objective[i,] *
    sqrt(moPLOT:::computeVectorLengthCPP(gradients$single.objective[[1]][i,])) *
    sqrt(moPLOT:::computeVectorLengthCPP(gradients$single.objective[[2]][i,]))
}) %>% t

# grad <- gradients
# grad$multi.objective <- rescaled
# less <- moPLOT::localEfficientSetSkeleton(design, grad, divergence, integration = "fast")
# g <- moPLOT::ggplotPLOT(design$dec.space, design$obj.space, less$sinks, less$height)
# g.obj <- moPLOT::ggplotPLOTObjSpace(design$obj.space, less$sinks, less$height)

# moPLOT::addGGArrows(g, design$dec.space, rescaled, fac = 0.0002, nRows = 25, nColumns = 25)

vnorm <- function(x) sqrt(sum(x ** 2))

original_length <- apply(gradients$multi.objective, 1, function(r) moPLOT:::computeVectorLengthCPP(r))
rescaled_length <- apply(rescaled, 1, function(r) moPLOT:::computeVectorLengthCPP(r))

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = original_length))
moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = rescaled_length))

moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = original_length))
moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = rescaled_length))

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = apply(gradients$multi.objective, 1, vnorm)), log.scale = T)

# norm.1 = apply(gradients$single.objective[[1]], 1, vnorm)
# norm.2 = apply(gradients$single.objective[[2]], 1, vnorm)
# 
# norm.height = 1 / 2 *
#   (sqrt(norm.2 / norm.1) * (design$obj.space[,1] - min(design$obj.space[,1])) +
#    sqrt(norm.1 / norm.2) * (design$obj.space[,2] - min(design$obj.space[,2])))
# norm.height[norm.height == Inf] = 0
# 
# moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = norm.height - min(norm.height)), log.scale = T)
# moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = norm.height - min(norm.height)), log.scale = T)
# 
# moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = sqrt(norm.1 * norm.2) * apply(gradients$multi.objective, 1, vnorm)), log.scale = T)
# moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = apply(gradients$multi.objective, 1, vnorm)), log.scale = T)

# ========= MOGSA =========

d <- 2
fid <- 10
iid <- 1
fn <- smoof::makeBiObjBBOBFunction(d, fid, iid)

# fn <- smoof::makeWFG5Function(2, 1, 1)
# fn <- smoof::makeDTLZ1Function(dimensions = 3, n.objectives = 2)
# fn <- makeAsparFunction(2, 2)
# fn <- smoof::makeMultiObjectiveFunction("test", fn = function(x) c(sum(x ** 2), sum(x ** 2)),
#                                         par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-5, -5), upper = c(5, 5)))
# fn <- smoof::makeMMF4Function()

lower <- smoof::getLowerBoxConstraints(fn)
upper <- smoof::getUpperBoxConstraints(fn)

nstarts <- 101
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
                              eps_initial_step_size = 1e-7,
                              eps_explore_set = 1e-6,
                              max_explore_set = sqrt(sum((upper - lower) ** 2)) / 100,
                              # custom_descent_fn = create_lbfgsb_descent(f, lower, upper)
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

plot_dec_space(mogsa_trace, lower, upper) +
  coord_fixed() +
  theme(legend.position = "none")

plot_obj_space(mogsa_trace) +
  theme(legend.position = "none")

ggplot() +
  geom_point(aes(x = x1, y = x2, color = 1:nrow(opt.path)), data = opt.path) +
  xlim(lower[1], upper[1]) +
  ylim(lower[2], upper[2]) +
  coord_fixed()

ggplot() +
  geom_point(aes(x = V1, y = V2, color = 1:ncol(obj.opt.path)), data = as.data.frame(t(obj.opt.path)))

# ggplot() +
#   geom_path(aes(x = x1, y = x2), data = opt.path)

# obj.opt.path <- t(do.call(rbind, lapply(mogsa_trace$sets, function(l) l$obj_space)))

ncol(obj.opt.path)

nondom <- nondominated(t(obj.opt.path))

obj.opt.path <- as.data.frame(t(obj.opt.path))

colnames(obj.opt.path) <- c("y1", "y2")

nrow(obj.opt.path) / d

sum(nondom) / nrow(obj.opt.path)
sum(duplicated(obj.opt.path)) / nrow(obj.opt.path)
# sum(duplicated(opt.path)) / nrow(opt.path)

# ref.point <- apply(obj.opt.path[nondom,1:2], 2, max)
# ideal.point <- apply(obj.opt.path[nondom,1:2], 2, min)

fn1 <- smoof::makeBBOBFunction(d, 1, 7)
fn2 <- smoof::makeBBOBFunction(d, 1, 8)

ref.point <- c(fn1(smoof::getGlobalOptimum(fn2)$param),
               fn2(smoof::getGlobalOptimum(fn1)$param))
ideal.point <- c(smoof::getGlobalOptimum(fn1)$value,
                 smoof::getGlobalOptimum(fn2)$value)

hv.norm <- prod(ref.point - ideal.point)
hv = ecr::computeHV(t(as.matrix(obj.opt.path[nondom,1:2])), ref.point = ref.point)
hv / hv.norm

log10(5/6 - hv / hv.norm)
# log10(0.890601509335293 - hv / hv.norm)

ggplot() +
  geom_point(aes(x = y1, y = y2), data = obj.opt.path[!nondom,], color = "red") +
  geom_point(aes(x = y1, y = y2), data = obj.opt.path[nondom,], color = "green") +
  geom_point(aes(x = ref.point[1], y = ref.point[2]), shape = "+", size = 5) +
  geom_point(aes(x = ideal.point[1], y = ideal.point[2]), shape = "+", size = 5) +
  xlim(ideal.point[1], ref.point[1]) +
  ylim(ideal.point[2], ref.point[2])

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

# placement <- Reduce(rbind, set_medians(mogsa_trace$sets))
# colnames(placement) = c("x", "y")
# placement <- as.data.frame(placement)

(1 / prop[set_nd_counts > 0]) %>% sort(decreasing = TRUE)

ggraph(tbl_transitions, layout = "stress") +
  geom_edge_fan(arrow = arrow(length = unit(4, "mm")),
                end_cap = circle(4, "mm")) +
  geom_node_point(aes(size = prop), color = node_color) +
  scale_size_area(limits = c(0,1))

# HV Stuff

# ========= Visualize =========

ref_point <- apply(design$obj.space, 2, max)
ref_point <- design$obj.space[12345,]
# ref_point <- fn(c(0.5, 1))

obj_space_hv <- hv_from_obj_space(design$obj.space, ref_point = ref_point, geom_mean = TRUE)
moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = -obj_space_hv + max(obj_space_hv)), log.scale = FALSE)
moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = -obj_space_hv + max(obj_space_hv)), log.scale = FALSE)

# ========= Optimize =========

f <- fn
f <- smoof::addLoggingWrapper(f, logg.x = FALSE, logg.y = TRUE)

gf <- function(ref_point) {
  function(x) {
    y <- f(x)
    G <- moPLOT::estimateGradientBothDirections(f, x, prec.grad = 1e-08, lower = lower, upper = upper)
    
    g1 <- G[1,]
    g2 <- G[2,]
    
    norm_g1 <- moPLOT:::computeVectorLengthCPP(g1)
    norm_g2 <- moPLOT:::computeVectorLengthCPP(g2)
    
    if (ecr::dominates(y, ref_point)) {
      hv_grad <- -((ref_point[2] - y[2]) * g1 +
                   (ref_point[1] - y[1]) * g2)
      

      hv <- (ref_point[1] - y[1]) * (ref_point[2] - y[2])
      
      hv_grad / (2 * sqrt(hv))
    } else {
      - 1 / 2 * (g1 / norm_g1 + g2 / norm_g2) * sqrt(norm_g1) * sqrt(norm_g2)
    }
    
  }
}

results <- lapply(1:11, function(i) {
  x_start <- runif_box(lower, upper)
  # starting_point <- 1 / 3 * smoof::getGlobalOptimum(fn1)$param + 2 / 3 * smoof::getGlobalOptimum(fn2)$param
  # starting_point <- starting_point + c(0.01, -0.01)
  y_start <- f(x_start)
  
  hv_function <- create_hv_function(f, ref_point = y_start, geom_mean = TRUE)
  
  result <- optim(x_start, function(x) hv_function(x), gr = gf(y_start),
                  method = "L-BFGS-B", control = list(fnscale = -1))
  # result <- optim(x_start, function(x) hv_function(x), gr = gf(y_start),
  #                 method = "L-BFGS-B", lower = lower, upper = upper,
  #                 control = list(fnscale = -1))
  # result <- optim(starting_point, function(x) hv_function(x), method = "L-BFGS-B", lower = lower, upper = upper, control = list(fnscale = -1))
  print(result$message)
  result$par
})

obj.opt.path <- smoof::getLoggedValues(f)$obj.vals
ncol(obj.opt.path)

sqrt_hv_grad <- gf(f(starting_point))

obj.opt.path <- Reduce(cbind, lapply(results, fn))
opt.path <- Reduce(rbind, results) %>% as.data.frame
