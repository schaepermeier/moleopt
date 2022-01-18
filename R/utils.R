f1_1 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1)
f2_1 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2)
f3_1 = function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2)
f_2d2d = function(x) c(f1_1(x), f2_1(x))
f_2d3d = function(x) c(f1_1(x), f2_1(x), f3_1(x))

f1_2 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1 + x[3] ** 2)
f2_2 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2 + (x[3] - 1) ** 4)
f3_2 = function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2 + (x[3] - 1) ** 2)
f_3d2d = function(x) c(f1_2(x), f2_2(x))
f_3d3d = function(x) c(f1_2(x), f2_2(x), f3_2(x))

makeAsparFunction <- function(dimensions = 2, n.objectives = 2) {
  
  if (dimensions == 2 && n.objectives == 2) {
    smoof::makeMultiObjectiveFunction(name = "2D->2D Aspar Function", id = "test_2d2d", description = "", fn = f_2d2d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))
  } else if (dimensions == 2 && n.objectives == 3) {
    smoof::makeMultiObjectiveFunction(name = "2D->3D Aspar Function", id = "test_2d3d", description = "", fn = f_2d3d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))
  } else if (dimensions == 3 && n.objectives == 2) {
    smoof::makeMultiObjectiveFunction(name = "3D->2D Aspar Function", id = "test_3d2d", description = "", fn = f_3d2d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))
  } else if (dimensions == 3 && n.objectives == 3) {
    smoof::makeMultiObjectiveFunction(name = "3D->3D Aspar Function", id = "test_3d3d", description = "", fn = f_3d3d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))
  }
}

makeSGKFunction = function() {
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

makeBiObjMPM2Function = function(dimensions = 2, n.peaks.1 = 3, topology.1 = "random", seed.1 = 4,
                                 n.peaks.2 = 3, topology.2 = "random", seed.2 = 8) {
  f1 <- smoof::makeMPM2Function(n.peaks.1, dimensions, topology.1, seed.1)
  f2 <- smoof::makeMPM2Function(n.peaks.2, dimensions, topology.2, seed.2)
  
  smoof::makeGOMOPFunction(dimensions = dimensions, funs = list(f1, f2))
}

makeSimpleBiSphereFunction = function() {
  f1 <- function(x) (sum((x - c(1,1))**2))
  f2 <- function(x) (sum((x + c(1,1))**2))
  
  f <- function(x) c(f1(x), f2(x))
  
  smoof::makeMultiObjectiveFunction(name = "Simple Bi-Sphere Function", id = "simple_bisphere", description = "", fn = f,
                                    par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-2), upper = c(2,2)))
}

makeBiRosenbrockFunction = function() {
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

runif_box <- function(lower, upper) {
  u <- runif(length(lower))
  
  u * (upper - lower) + lower
}

f_noglobal = function(i) {
  x = i[1]
  y = i[2]
  
  x1 = -exp(0.2*(-(x+3)**2-(y-3)**2)) -exp(0.2*(-(x-3)**2-(y+3)**2)) + exp(-(2*x-y)**2)
  x2 = -exp(0.2*(-(x-3)**2-(y-3)**2)) -exp(0.2*(-(x+3)**2-(y+3)**2)) + exp(-(x+2*y)**2)
  
  c(x1, x2)
}

test.multitrap = smoof::makeMultiObjectiveFunction(name = "Multi-set Trap", id = "", description = "", fn = f_noglobal,
                                                   par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-5,-5), upper = c(5,5)))
