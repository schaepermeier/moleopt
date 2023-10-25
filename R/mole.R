#' Execute MOLE algorithm
#' 
#' Calling `run_mole` executes a run of the MOLE algorithm with the given parameter setting.
#' 
#' The most important parameters to consider are:
#' 
#' - `fn`: The function to optimize
#' - `starting_points`: Pre-defined starting points for MOLE
#' - `descent_step_max`, `explore_step_max`: Here it may be beneficial to deviate from the
#'   default setting of `1e-1`. A reasonable choice in many situations is 1/100 of the diagonal 
#'   of the search space.
#' - `refine_hv_target`: If post-processing for HV optimization should be used,
#'   set this to the target value for the normalized hypervolume.
#'
#' @param fn Bi-objective `smoof` function to optimize
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
#' @param refine_hv_target Normalized Hypervolume gap target
#' @param custom_descent_fn Custom descent function implemented in R (optional)
#' @param lower Lower box constraints, if different than specified in `fn` (optional)
#' @param upper Lower box constraints, if different than specified in `fn` (optional)
#' @param max_budget Maximum budget, if budget is constrained
#' @param logging Logging level: `"none"`, `"debug"` or `"info"`
#'
#' @return
#' 
#' List containing:
#' 
#' - `sets`: Discovered locally efficient sets
#' - `transitions`: Transitions between locally efficient sets
#' - `budget_at_starting_points`: used budget before starting point was used
#' 
#' @export
#'
#' @examples
#' fn <- makeAsparFunction()
#' starting_points <- rbind(c(1,1))
#' mole_trace <- run_mole(fn, starting_points)
#' @useDynLib moleopt
run_mole <- function(fn,
                      starting_points,
                      max_local_sets = 1000,
                      epsilon_gradient = 1e-8,
                      descent_direction_min = 1e-8,
                      descent_step_min = 1e-6,
                      descent_step_max = 1e-1,
                      descent_scale_factor = 2,
                      descent_armijo_factor = 1e-4,
                      descent_history_size = 100L,
                      descent_max_iter = 1000L,
                      explore_step_min = 1e-4,
                      explore_step_max = 1e-1,
                      explore_angle_max = 45,
                      explore_scale_factor = 2,
                      refine_after_nstarts = 10L,
                      refine_hv_target = 2e-5,
                      custom_descent_fn = NULL,
                      lower = NULL,
                      upper = NULL,
                      max_budget = Inf,
                      logging = "info") {
  
  # ==== Configure defaults if values are missing ====
  
  if (is.null(lower)) {
    lower <- smoof::getLowerBoxConstraints(fn)
  }
  
  if (is.null(upper)) {
    upper <- smoof::getUpperBoxConstraints(fn)
  }
  
  # ==== Checkmate ====
  
  checkmate::assert_class(fn, "smoof_function")
  n_dimensions <- smoof::getNumberOfParameters(fn)
  checkmate::assert_true(smoof::getNumberOfObjectives(fn) == 2L)
  
  checkmate::assert_matrix(starting_points, any.missing = FALSE, ncols = n_dimensions)
  
  checkmate::assert_numeric(epsilon_gradient, lower = .Machine$double.eps, finite = TRUE)
  
  checkmate::assert_numeric(descent_direction_min, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_numeric(descent_step_min, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_numeric(descent_step_max, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_true(descent_step_min <= descent_step_max)
  checkmate::assert_numeric(descent_scale_factor, lower = 1 + .Machine$double.eps, finite = TRUE)
  checkmate::assert_numeric(descent_armijo_factor, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_integerish(descent_history_size, lower = 0)
  checkmate::assert_integerish(descent_max_iter, lower = 0)
  
  checkmate::assert_numeric(explore_step_min, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_numeric(explore_step_max, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_true(explore_step_min <= explore_step_max)
  checkmate::assert_numeric(explore_angle_max, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_numeric(explore_scale_factor, lower = 1 + .Machine$double.eps, finite = TRUE)
  checkmate::assert_integerish(refine_after_nstarts, lower = 1)
  checkmate::assert_numeric(refine_hv_target, upper = 1)
  checkmate::assert_function(custom_descent_fn, null.ok = TRUE)

  checkmate::assert_vector(lower, len = n_dimensions, any.missing = FALSE)
  checkmate::assert_vector(upper, len = n_dimensions, any.missing = FALSE)
  checkmate::assert_true(all(lower <= upper))

  checkmate::assert_numeric(max_budget, lower = 1)
  checkmate::assert_choice(logging, c("none", "info", "debug"))
  
  # ==== Call MOLE implementation ====
  
  run_mole_cpp(fn = fn,
                starting_points = starting_points,
                lower_bounds = lower,
                upper_bounds = upper, 
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
