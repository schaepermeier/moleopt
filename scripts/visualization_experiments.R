library(tidyverse)

# ========= moPLOT =========

design <- moPLOT::generateDesign(fn, points.per.dimension = 1001L)
design$obj.space <- moPLOT::calculateObjectiveValues(design$dec.space, fn, parallelize = TRUE)

gradients <- moPLOT::computeGradientFieldGrid(design)
divergence <- moPLOT::computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)

less <- moPLOT::localEfficientSetSkeleton(design, gradients, divergence, integration = "fast")

g <- moPLOT::ggplotPLOT(design$dec.space, design$obj.space, less$sinks, less$height)
g.obj <- moPLOT::ggplotPLOTObjSpace(design$obj.space, less$sinks, less$height)

rescaled <- sapply(1:nrow(gradients$multi.objective), function(i) {
  gradients$multi.objective[i,] *
    sqrt(moPLOT:::computeVectorLengthCPP(gradients$single.objective[[1]][i,])) *
    sqrt(moPLOT:::computeVectorLengthCPP(gradients$single.objective[[2]][i,]))
}) %>% t

# grad <- gradients
# grad$multi.objective <- rescaled
# less <- moPLOT::localEfficientSetSkeleton(design, grad, divergence, integration = "fast")
# g <- moPLOT::ggplotPLOT(design$dec.space, design$obj.space, less$sinks, less$height)
# g.obj <- moPLOT::ggplotPLOTObjSpace(design$obj.space, less$sinks, less$height)

# moPLOT::addGGArrows(g, design$dec.space, rescaled, fac = 0.0002, nRows = 25, nColumns = 25)

vnorm <- function(x) sqrt(sum(x ** 2))

original_length <- apply(gradients$multi.objective, 1, function(r) moPLOT:::computeVectorLengthCPP(r))
rescaled_length <- apply(rescaled, 1, function(r) moPLOT:::computeVectorLengthCPP(r))

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = original_length))
moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = rescaled_length))

moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = original_length))
moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = rescaled_length))

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = apply(gradients$multi.objective, 1, vnorm)), log.scale = T)

# norm.1 = apply(gradients$single.objective[[1]], 1, vnorm)
# norm.2 = apply(gradients$single.objective[[2]], 1, vnorm)
# 
# norm.height = 1 / 2 *
#   (sqrt(norm.2 / norm.1) * (design$obj.space[,1] - min(design$obj.space[,1])) +
#    sqrt(norm.1 / norm.2) * (design$obj.space[,2] - min(design$obj.space[,2])))
# norm.height[norm.height == Inf] = 0
# 
# moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = norm.height - min(norm.height)), log.scale = T)
# moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = norm.height - min(norm.height)), log.scale = T)
# 
# moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = sqrt(norm.1 * norm.2) * apply(gradients$multi.objective, 1, vnorm)), log.scale = T)
# moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = apply(gradients$multi.objective, 1, vnorm)), log.scale = T)
