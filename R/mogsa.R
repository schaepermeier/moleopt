#' @export
run_mogsa <- function(fn, starting_points,
                      eps_gradient = 1e-8, eps_initial_step_size = 1e-6, eps_explore_set = 1e-4, max_explore_set = 1e-2,
                      custom_descent_fn = NULL, lower = NULL, upper = NULL, max_budget = Inf, logging = TRUE) {
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
                epsilon_gradient = eps_gradient,
                epsilon_initial_step_size = eps_initial_step_size,
                epsilon_explore_set = eps_explore_set,
                max_explore_set = max_explore_set,
                custom_descent_fn = custom_descent_fn,
                max_budget = max_budget,
                logging = logging)
}
