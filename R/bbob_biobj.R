# single-objective fids used in bbob_biobj
fids <- c(1L, 2L, 6L, 8L, 13L, 14L, 15L, 17L, 20L, 21L)

fid_mapping <- expand.grid(fids, fids)
fid_mapping <- fid_mapping[fid_mapping[, 1L] <= fid_mapping[, 2L], ]
fid_mapping <- fid_mapping[order(fid_mapping[, 1L]),]
fid_mapping <- as.matrix(fid_mapping)
names(fid_mapping) <- NULL

iid_mapping <- matrix(
  c(
    2, 4,
    3, 5,
    7, 8,
    9, 10,
    11, 12,
    13, 14,
    15, 16,
    17, 18,
    19, 21,
    21, 22,
    23, 24,
    25, 26,
    27, 28,
    29, 30,
    31, 34
  ),
  ncol = 2,
  byrow = TRUE
)

#' Generate data for bi-objective BBOB functions
#'
#' @param dimensions number of dimensions for the test problem
#' @param fid function ID (1-55) of the bi-objective BBOB problem
#' @param iid instance ID (1-15) of the bi-objective BBOB problem
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{fn}: \code{smoof} function of the requested bi-objective BBOB function
#'   \item \code{opt_f1}: Optimal input value for the first constituent BBOB function
#'   \item \code{opt_f2}: Optimal input value for the second constituent BBOB function
#'   \item \code{ref_point}: Reference point for this function as used by thebi-objective BBOB for evaluation
#'   \item \code{ideal_point}: Ideal point for this function as used by the bi-objective BBOB for evaluation
#' }
#' @export
#'
#' @examples
#' # generate data for the 1st instance of the 10th function in 2d
#' bbob_biobj_data <- generateBiObjBBOBData(dimensions = 2L, fid = 10L, iid = 1L)
generateBiObjBBOBData <- function(dimensions, fid, iid) {
  assert_choice(dimensions, c(2L, 3L, 5L, 10L, 20L, 40L))
  assert_choice(fid, 1L:55L)
  assert_choice(iid, 1L:15L)
  
  output <- list()
  
  fid_1 <- fid_mapping[fid,1]
  fid_2 <- fid_mapping[fid,2]
  iid_1 <- iid_mapping[iid,1]
  iid_2 <- iid_mapping[iid,2]
  
  fn_1 <- smoof::makeBBOBFunction(dimensions, fid_1, iid_1)
  fn_2 <- smoof::makeBBOBFunction(dimensions, fid_2, iid_2)
  
  # build parameter set (bounds are [-5, 5] for all BBOB funs)
  par.set <- ParamHelpers::makeNumericParamSet("x", len = dimensions, lower = -5, upper = 5)
  
  output$fn <- smoof::makeMultiObjectiveFunction(
    name = sprintf("Bi-Objective BBOB_%i_%i_%i", dimensions, fid, iid),
    id = paste0("biobj_bbob_", dimensions, "d"),
    description = sprintf("%i-th noiseless Bi-Objective BBOB function\n(FID: %i, IID: %i, DIMENSION: %i)",
                          fid, fid, iid, dimensions),
    fn = function(x) {
      c(fn_1(x), fn_2(x))
    },
    par.set = par.set,
    n.objectives = 2L,
    # the single-objective BBOB functions are vectorized,
    # but not the combined one
    vectorized = FALSE
  )
  
  output$opt_f1 <- smoof::getGlobalOptimum(fn_1)$param
  output$opt_f2 <- smoof::getGlobalOptimum(fn_2)$param
  
  output$ref_point <- c(fn_1(smoof::getGlobalOptimum(fn_2)$param),
                        fn_2(smoof::getGlobalOptimum(fn_1)$param))
  
  output$ideal_point <- c(smoof::getGlobalOptimum(fn_1)$value,
                          smoof::getGlobalOptimum(fn_2)$value)
  
  return(output)
}
