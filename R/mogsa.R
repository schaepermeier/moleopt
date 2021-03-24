#' @export
run_mogsa <- function(fn,
                      starting_points,
                      epsilon_gradient = 1e-8,
                      descent_direction_min = 1e-8,
                      descent_step_min = 1e-6,
                      descent_step_max = 1e-1,
                      descent_scale_factor = 2,
                      descent_armijo_factor = 1e-4,
                      descent_history_size = 100,
                      descent_max_iter = 1000,
                      explore_step_min = 1e-4,
                      explore_step_max = 1e-1,
                      explore_angle_max = 45,
                      explore_scale_factor = 2,
                      refine_after_nstarts = 10,
                      refine_hv_target = 2e-5,
                      custom_descent_fn = NULL,
                      lower = NULL,
                      upper = NULL,
                      max_budget = Inf,
                      logging = "info") {
  if (is.null(lower)) {
    lower <- smoof::getLowerBoxConstraints(fn)
  }
  
  if (is.null(upper)) {
    upper <- smoof::getUpperBoxConstraints(fn)
  }

  run_mogsa_cpp(fn = fn,
                starting_points = starting_points,
                lower = lower,
                upper = upper, 
                epsilon_gradient = epsilon_gradient,
                descent_direction_min = descent_direction_min,
                descent_step_min = descent_step_min,
                descent_step_max = descent_step_max,
                descent_scale_factor = descent_scale_factor,
                descent_armijo_factor = descent_armijo_factor,
                descent_history_size = descent_history_size,
                descent_max_iter = descent_max_iter,
                explore_step_min = explore_step_min,
                explore_step_max = explore_step_max,
                explore_angle_max = explore_angle_max,
                explore_scale_factor = explore_scale_factor,
                refine_after_nstarts = refine_after_nstarts,
                refine_hv_target = refine_hv_target,
                custom_descent_fn = custom_descent_fn,
                max_budget = max_budget,
                logging = logging)
}
