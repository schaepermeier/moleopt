all_stats <- data.frame()

for (iid in 1:10) {
  d <- 2
  fid <- 10
  # iid <- 10
  biobj_bbob_data <- generateBiObjBBOBData(d, fid, iid)
  fn <- biobj_bbob_data$fn
  
  lower <- smoof::getLowerBoxConstraints(fn)
  upper <- smoof::getUpperBoxConstraints(fn)
  
  nstarts <- 100
  starting_points <- lapply(1:nstarts, function(x) runif_box(lower, upper))
  starting_points <- do.call(rbind, starting_points)
  
  # starting_points <- design$dec.space[less$sinks,]
  
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
    
    mole_trace <- run_mole(f, starting_points,
                           max_local_sets = 1000,
                           epsilon_gradient = 1e-8,
                           descent_direction_min = 1e-6,
                           descent_step_min = 1e-6,
                           descent_step_max = sqrt(sum((upper - lower) ** 2)) / 100,
                           descent_scale_factor = 2,
                           descent_armijo_factor = 1e-4,
                           descent_history_size = 100,
                           descent_max_iter = 1000,
                           explore_step_min = 1e-4,
                           explore_step_max = 1e-2,
                           # explore_step_max = sqrt(sum((upper - lower) ** 2)) / 100,
                           explore_angle_max = 20,
                           explore_scale_factor = 2,
                           refine_after_nstarts = 100,
                           refine_hv_target = 2e-5,
                           # custom_descent_fn = create_lbfgsb_descent(f, lower, upper),
                           # lower = rep(-100, length(lower)),
                           # upper = rep(100, length(upper)),
                           max_budget = Inf,
                           logging = "info"
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
  
  all_stats <- rbind(all_stats, unlist(compute_statistics(mole_trace)))
  
  p_interact <- plot_set_interactions(mole_trace, layout_2d = TRUE)
  p_sets <- plot_dec_space(mole_trace, lower, upper) +
    coord_fixed() +
    theme_minimal() +
    labs(x = expression(x[1]),
         y = expression(x[2])) +
    theme(legend.position = "none")
  
  ggsave(paste0("~/Desktop/thesis-pics/mole/interact-f10-i", iid, ".png"), width = unit(3, "in"), height = unit(3, "in"), p_interact)
  ggsave(paste0("~/Desktop/thesis-pics/mole/sets-f10-i", iid, ".png"), width = unit(3, "in"), height = unit(3, "in"), p_sets)
}

colnames(all_stats) <- c("nsets", "globsets", "reachability")

all_stats
