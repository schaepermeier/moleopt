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

#' @export
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
