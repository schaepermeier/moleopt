#' Execute MOLE algorithm
#'
#' @param fn Bi-objective \code{smoof} function to optimize
#' @param starting_points Matrix of starting points to pass to MOLE
#' @param max_local_sets Maximum number of local sets MOLE may track
#' @param epsilon_gradient Epsilon used for gradient estimation
#' @param descent_direction_min Minimal admissible length of MO gradient
#' @param descent_step_min Minimal admissible step size for MO descent
#' @param descent_step_max Maximal admissible step size for MO descent
#' @param descent_scale_factor Multiplicative factor by which to scale consecutive descent steps
#' @param descent_armijo_factor Armijo factor utilized in line search
#' @param descent_history_size Size of descent history tracked for nonmonotone descent
#' @param descent_max_iter Maximum number of iterations per descent call
#' @param explore_step_min Minimal admissible step size for exploration along locally efficient set
#' @param explore_step_max Maximal admissible step size for exploration along locally efficient set
#' @param explore_angle_max Maximal admissible angle between consecutive steps for exploration along locally efficient set
#' @param explore_scale_factor Multiplicative factor by which to scale consecutive descent steps
#' @param refine_after_nstarts Number of starting points to fully evaluate before starting refinement
#' @param refine_hv_target Target normalized Hypervolume gap
#' @param custom_descent_fn Custom descent function implemented in R (optional)
#' @param lower Lower box constraints, if different than specified in \code{fn} (optional)
#' @param upper Lower box constraints, if different than specified in \code{fn} (optional)
#' @param max_budget Maximum budget, if budget is constrained
#' @param logging Logging level: \code{"none"}, \code{"debug"} or \code{"info"}
#'
#' @return List containing:
#' \itemize{
#'   \item \code{sets}: Discovered locally efficient sets
#'   \item \code{transitions}: Transitions between locally efficient sets
#' }
#' @export
#'
#' @examples
#' 
run_mole <- function(fn,
                      starting_points,
                      max_local_sets = 1000,
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

  run_mole_cpp(fn = fn,
                starting_points = starting_points,
                lower = lower,
                upper = upper, 
                max_local_sets = max_local_sets,
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
