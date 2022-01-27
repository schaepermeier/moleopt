f1_1 <- function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1)
f2_1 <- function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2)
f3_1 <- function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2)
f_2d2d <- function(x) c(f1_1(x), f2_1(x))
f_2d3d <- function(x) c(f1_1(x), f2_1(x), f3_1(x))

f1_2 <- function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1 + x[3] ** 2)
f2_2 <- function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2 + (x[3] - 1) ** 4)
f3_2 <- function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2 + (x[3] - 1) ** 2)
f_3d2d <- function(x) c(f1_2(x), f2_2(x))
f_3d3d <- function(x) c(f1_2(x), f2_2(x), f3_2(x))

#' `smoof` generator for the Aspar problem with 2-3 dimensions and objectives.
#'
#' @return \code{smoof} function of the configured Aspar problem
#' @export
#'
#' @examples
#' fn <- makeAsparFunction(dimensions = 2L, n.objectives = 2L)
makeAsparFunction <- function(dimensions = 2L, n.objectives = 2L) {
  if (dimensions == 2L && n.objectives == 2L) {
    smoof::makeMultiObjectiveFunction(name = "2D->2D Aspar Function", id = "test_2d2d", description = "", fn = f_2d2d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 2L, lower = c(-2,-1), upper = c(2,3)))
  } else if (dimensions == 2L && n.objectives == 3L) {
    smoof::makeMultiObjectiveFunction(name = "2D->3D Aspar Function", id = "test_2d3d", description = "", fn = f_2d3d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 2L, lower = c(-2,-1), upper = c(2,3)))
  } else if (dimensions == 3L && n.objectives == 2L) {
    smoof::makeMultiObjectiveFunction(name = "3D->2D Aspar Function", id = "test_3d2d", description = "", fn = f_3d2d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 3L, lower = c(-2,-1,-2), upper = c(2,3,2)))
  } else if (dimensions == 3L && n.objectives == 3L) {
    smoof::makeMultiObjectiveFunction(name = "3D->3D Aspar Function", id = "test_3d3d", description = "", fn = f_3d3d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 3L, lower = c(-2,-1,-2), upper = c(2,3,2)))
  }
}

#' `smoof` generator for the SGK problem
#'
#' @return \code{smoof} function containing the bi-objective SGK problem
#' @export
#'
#' @examples
#' fn <- makeSGKFunction()
makeSGKFunction <- function() {
  g = function(x, h, c1, c2) {
    h / (1 + 4 * ((x[1] - c1) ** 2 + (x[2] - c2) ** 2))
  }
  
  f = function(x) c(
    1 - 1 / (1 + 4 * ((x[1] - 2 / 3) ** 2 + (x[2] - 1) ** 2)),
    1 - max(
      g(x, 1.5, 0.5, 0),
      g(x, 2, 0.25, 2 / 3),
      g(x, 3, 1, 1)
    )
  )
  
  lower = c(-0.25, -0.25)
  upper = c(1.25, 1.25)
  
  smoof::makeMultiObjectiveFunction(
    name = "SGK Function", id = "sgk_function", description = "", fn = f,
    par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = lower, upper = upper))
}

#' `smoof` generator for bi-objective MPM2 problems
#'
#' @param dimensions Number of search space dimensions
#' @param n.peaks.1 Number of peaks for first MPM2 problem
#' @param topology.1 Type of topology for first MPM2 problem
#' @param seed.1 Random seed for first MPM2 problem
#' @param n.peaks.2 Number of peaks for second MPM2 problem
#' @param topology.2 Type of topology for second MPM2 problem
#' @param seed.2 Random seed for second MPM2 problem
#'
#' @return \code{smoof} function of the configured bi-objective MPM2 function
#' @export
#'
#' @examples
#' fn <- makeBiObjMPM2Function(dimensions = 2,
#'   n.peaks.1 = 3, topology.1 = "random", seed.1 = 4,
#'   n.peaks.2 = 3, topology.2 = "random", seed.2 = 8
#' )
makeBiObjMPM2Function <- function(dimensions = 2, n.peaks.1 = 3, topology.1 = "random", seed.1 = 4,
                                 n.peaks.2 = 3, topology.2 = "random", seed.2 = 8) {
  f1 <- smoof::makeMPM2Function(n.peaks.1, dimensions, topology.1, seed.1)
  f2 <- smoof::makeMPM2Function(n.peaks.2, dimensions, topology.2, seed.2)
  
  smoof::makeGOMOPFunction(dimensions = dimensions, funs = list(f1, f2))
}

#' `smoof` generator for simple Bi-Sphere problem
#'
#' @return \code{smoof} function containing a simple Bi-Sphere problem with
#' centers \code{c(-1, -1)} and \code{c(1, 1)}
#' @export
#'
#' @examples
#' fn <- makeSimpleBiSphereFunction()
makeSimpleBiSphereFunction <- function() {
  f1 <- function(x) (sum((x - c(1,1))**2))
  f2 <- function(x) (sum((x + c(1,1))**2))
  
  f <- function(x) c(f1(x), f2(x))
  
  smoof::makeMultiObjectiveFunction(name = "Simple Bi-Sphere Function", id = "simple_bisphere", description = "", fn = f,
                                    par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-2), upper = c(2,2)))
}

#' `smoof` generator for Bi-Rosenbrock function
#'
#' @return \code{smoof} function containing the Bi-Rosenbrock function
#' @export
#'
#' @examples
#' fn <- makeBiRosenbrockFunction()
makeBiRosenbrockFunction <- function() {
  f1 <- function(x) {
    (1 - x[1]) ** 2 + 1 * (x[2] - x[1] ** 2) ** 2
  }
  
  f2 <- function(x) {
    (1 + x[1]) ** 2 + 1 * (-(x[2] - 3) - x[1] ** 2) ** 2
  }
  
  f <- function(x) {
    c(f1(x), f2(x))
  }
  
  smoof::makeMultiObjectiveFunction(
    name = "Bi-Rosenbrock", id = "bi_rosenbrock_function", description = "", fn = f,
    par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2, 0), upper = c(2, 3)))
}


#' Generates a vector distributed uniformly at random between \code{lower} and
#' \code{upper}
#'
#' @param lower lower bounds of search space
#' @param upper upper bounds of search space
#'
#' @return Vector of length \code{length(lower)} with entries uniformly at
#' random between \code{lower} and \code{upper} per dimension
#' 
#' @export
#' @examples
#' # Create a vector of length 2 with elements uniformly at random between -5 and 5
#' runif_box(c(-5, -5), c(5, 5))
#' 
#' # Generating a matrix of values
#' points <- lapply(1:nstarts, function(x) runif_box(lower, upper))
#' points_matrix <- do.call(rbind, points)
runif_box <- function(lower, upper) {
  checkmate::assert_true(length(lower) == length(upper))
  checkmate::assert_true(all(lower <= upper))
  
  u <- runif(length(lower))
  
  u * (upper - lower) + lower
}

# f_noglobal = function(i) {
#   x = i[1]
#   y = i[2]
#   
#   x1 = -exp(0.2*(-(x+3)**2-(y-3)**2)) -exp(0.2*(-(x-3)**2-(y+3)**2)) + exp(-(2*x-y)**2)
#   x2 = -exp(0.2*(-(x-3)**2-(y-3)**2)) -exp(0.2*(-(x+3)**2-(y+3)**2)) + exp(-(x+2*y)**2)
#   
#   c(x1, x2)
# }
# 
# test.multitrap = smoof::makeMultiObjectiveFunction(name = "Multi-set Trap", id = "", description = "", fn = f_noglobal,
#                                                    par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-5,-5), upper = c(5,5)))
