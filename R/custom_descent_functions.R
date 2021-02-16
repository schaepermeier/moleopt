#' @export
create_lbfgsb_descent <- function(fn, lower, upper) {
  function(x_start, ref_point) {
    
    # TODO: Overwrite ref_point currently, because
    # the input is not really a ref_point right now
    ref_point <- fn(x_start)
    
    hv_grad <- create_hv_gradient_function(fn, ref_point, lower, upper)
    hv_function <- create_hv_function(fn, ref_point = ref_point, geom_mean = TRUE)
    
    tryCatch({
      result <- optim(x_start, hv_function, gr = hv_grad,
            method = "L-BFGS-B", lower = lower, upper = upper,
            control = list(fnscale = -1))
      
      # print(result)
      
      list(dec_space = result$par,
           obj_space = fn(result$par))
    }, error = function(e) {
      list(dec_space = x_start,
           obj_space = fn(x_start))
    })
    
  }
}

#' @export
create_hv_gradient_function <- function(f, ref_point, lower, upper) {
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
