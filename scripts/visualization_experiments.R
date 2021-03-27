library(tidyverse)

# ========= moPLOT =========

design <- moPLOT::generateDesign(fn, points.per.dimension = 500L)
design$obj.space <- moPLOT::calculateObjectiveValues(design$dec.space, fn, parallelize = TRUE)

gradients <- moPLOT::computeGradientFieldGrid(design, normalized.scale = TRUE)
divergence <- moPLOT::computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)

less <- moPLOT::localEfficientSetSkeleton(design, gradients, divergence, integration = "fast")

g <- moPLOT::ggplotPLOT(design$dec.space, design$obj.space, less$sinks, less$height) +
  coord_fixed()
g.obj <- moPLOT::ggplotPLOTObjSpace(design$obj.space, less$sinks, less$height)

# === CONTOUR PLOTS ===

contour_height <- apply(design$obj.space, 2, function(c) log(c - min(c)))
contour_height <- apply(contour_height, 2, function(c) {
  c[c == -Inf] <- min(c[c != -Inf]) / 2
  c
})

g +
  geom_contour(aes(x1, x2, z = contour_height[,1]), data = as.data.frame(design$dec.space), bins = 20, color = "red") +
  geom_contour(aes(x1, x2, z = contour_height[,2]), data = as.data.frame(design$dec.space), bins = 20, color = "blue")

opt_f1 <- biobj_bbob_data$opt_f1
opt_f2 <- biobj_bbob_data$opt_f2

ggplot(data = as.data.frame(design$dec.space)) +
  geom_contour(aes(x1, x2, z = contour_height[,1]), color = "red", bins = 25) +
  geom_contour(aes(x1, x2, z = contour_height[,2]), color = "blue", bins = 25) +
  annotate(x = unlist(opt_f1[1]), y = unlist(opt_f1[2]), shape = "+", size = 10, color = "red", geom = "point") +
  annotate(x = unlist(opt_f2[1]), y = unlist(opt_f2[2]), shape = "+", size = 10, color = "blue", geom = "point") +
  coord_fixed() +
  theme_minimal()

ggplot(data = as.data.frame(design$dec.space)) +
  geom_contour(aes(x1, x2, z = contour_height[,2], color = ..level..), bins = 50) +
  coord_fixed() +
  scale_color_gradientn(colors = fields::tim.colors(500)) +
  theme_minimal() +
  labs(x = expression(x[1]), y = expression(x[2])) +
  theme(legend.position = "none")

ggsave("~/Desktop/thesis-pics/peaks-log-contour.pdf", width = unit(3, "in"), height = unit(3, "in"))

ggplot(data = as.data.frame(design$dec.space)) +
  geom_contour(aes(x1, x2, z = exp(contour_height[,2]), color = ..level..), bins = 20) +
  coord_fixed() +
  scale_color_gradientn(colors = fields::tim.colors(500)) +
  theme_minimal()

# === HEATMAPS ===

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height=design$obj.space[,2]-min(design$obj.space[,2])), 
                      log.scale = F) +
  coord_fixed() +
  labs(fill = expression(y - y[opt])) +
  theme(legend.position = "none")
  
ggsave("~/Desktop/thesis-pics/peaks-heatmap.pdf", width = unit(3, "in"), height = unit(3, "in"))

# === 3D ===

plotly::plot_ly(data = cbind.data.frame(design$dec.space),
                x=~unique(x1), y=~unique(x2), z=matrix(design$obj.space[,2], nrow = sqrt(nrow(design$obj.space))),
                type = "surface") %>% plotly::layout(
                  scene = list(
                    xaxis = list(title = "x₁"),
                    yaxis = list(title = "x₂"),
                    zaxis = list(title = "y"),
                    color
                  ))


# === Gradients ===


# norm.1 = apply(gradients$single.objective[[1]], 1, vnorm)
# norm.2 = apply(gradients$single.objective[[2]], 1, vnorm)
# 
# norm.height = 1 / 2 *
#   (sqrt(norm.2 / norm.1) * (design$obj.space[,1] - min(design$obj.space[,1])) +
#    sqrt(norm.1 / norm.2) * (design$obj.space[,2] - min(design$obj.space[,2])))
# norm.height[norm.height == Inf] = 0
# 
moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = norm.height - min(norm.height)), log.scale = T)
# moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = norm.height - min(norm.height)), log.scale = T)
# 
# moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = sqrt(norm.1 * norm.2) * apply(gradients$multi.objective, 1, vnorm)), log.scale = T)
# moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = apply(gradients$multi.objective, 1, vnorm)), log.scale = T)



# === MO VIZ ===

nds <- ecr::doNondominatedSorting(t(design$obj.space))

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = (nds$dom.counter + 1)), log.scale = T) +
  coord_fixed() +
  theme(legend.position = "none")
  
ggsave("~/Desktop/thesis-pics/aspar-cl-dec.png", width = unit(3, "in"), height = unit(3, "in"))

moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = (nds$dom.counter + 1)), log.scale = T) +
  theme(legend.position = "none")

ggsave("~/Desktop/thesis-pics/aspar-cl-obj.png", width = unit(3, "in"), height = unit(3, "in"))

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, less$height)) +
  coord_fixed() +
  theme(legend.position = "none")

ggsave("~/Desktop/thesis-pics/aspar-gfh-dec.png", width = unit(3, "in"), height = unit(3, "in"))

moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, less$height)) +
  theme(legend.position = "none")

ggsave("~/Desktop/thesis-pics/aspar-gfh-obj.png", width = unit(3, "in"), height = unit(3, "in"))

g

ggsave("~/Desktop/thesis-pics/aspar-plot-dec.png", width = unit(3, "in"), height = unit(3, "in"))

g.obj

ggsave("~/Desktop/thesis-pics/aspar-plot-obj.png", width = unit(3, "in"), height = unit(3, "in"))



mog_ch_length <- lapply(1:nrow(gradients$single.objective[[1]]), function(i) {
  g1 <- 100 * gradients$single.objective[[1]][i,]
  g2 <- gradients$single.objective[[2]][i,]
  
  # const float l2 = length_squared(v, w);  // i.e. |w-v|^2 -  avoid a sqrt
  # if (l2 == 0.0) return distance(p, v);   // v == w case
  # const float t = max(0, min(1, dot(p - v, w - v) / l2));
  # const vec2 projection = v + t * (w - v);  // Projection falls on the segment
  # return distance(p, projection);
  
  l2 <- vnorm(g1 - g2) ** 2
  if (l2 == 0) {
    vnorm(g1)
  } else {
    t <- max(0, min(1, sum((-g1) * (g2 - g1)) / l2))
    proj <- g1 + t * (g2 - g1)
    min(vnorm(proj))
  }
  
}) %>% unlist

obj_space_transformed <- t(apply(design$obj.space, 1, function(c) c * c(100, 1)))

vnorm <- function(x) sqrt(sum(x ** 2))

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = mog_ch_length), log.scale = F) +
  coord_fixed()
moPLOT::ggplotObjectiveSpace(cbind.data.frame(obj_space_transformed, height = mog_ch_length), log.scale = F)

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = apply(gradients$multi.objective, 1, vnorm)), log.scale = F) +
  coord_fixed()
moPLOT::ggplotObjectiveSpace(cbind.data.frame(obj_space_transformed, height = apply(gradients$multi.objective, 1, vnorm)), log.scale = F)

rescaled <- sapply(1:nrow(gradients$multi.objective), function(i) {
  gradients$multi.objective[i,] *
    sqrt(vnorm(gradients$single.objective[[1]][i,])) *
    sqrt(vnorm(gradients$single.objective[[2]][i,]))
}) %>% t

original_length <- apply(gradients$multi.objective, 1, vnorm)
rescaled_length <- apply(rescaled, 1, vnorm)

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = original_length), log.scale = F) +
  coord_fixed()
moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height = original_length), log.scale = F)

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, height = rescaled_length), log.scale = F) +
  coord_fixed()
moPLOT::ggplotObjectiveSpace(cbind.data.frame(obj_space_transformed, height = rescaled_length), log.scale = F)

### Original MOGSA ###

library(mogsa)
library(tidyverse)

mogsa.result = mogsa::runMOGSA(c(1,1), fn, scale.step = 0.5, exploration.step = 0.1,
                        lower = smoof::getLowerBoxConstraints(fn), upper = smoof::getUpperBoxConstraints(fn))

g +
  geom_point(data = mogsa.result, mapping = aes(x1, x2, shape = type))

ggsave("~/Desktop/thesis-pics/mogsa-aspar.png", width = unit(3, "in"), height = unit(3, "in"))

yvals <- apply(mogsa.result[,1:2], 1, fn) %>% t
colnames(yvals) <- c("y1", "y2")

g.obj +
  geom_point(data = as.data.frame(yvals), mapping = aes(y1, y2))

### HEATMAPS

moPLOT::ggplotHeatmap(cbind.data.frame(design$dec.space, less$height)) +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("~/Desktop/thesis-pics/bbobbiobj-gfh-dec.png", width = unit(3, "in"), height = unit(3, "in"))

moPLOT::ggplotObjectiveSpace(cbind.data.frame(design$obj.space, less$height)) +
  theme(legend.position = "none")
ggsave("~/Desktop/thesis-pics/bbobbiobj-gfh-obj.png", width = unit(3, "in"), height = unit(3, "in"))
