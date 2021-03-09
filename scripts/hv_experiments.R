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
