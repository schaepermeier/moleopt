test_one_archive <- function(archive, ref_point) {
  hv_a <- computeHVArchive(archive, ref_point)
  hv_ecr <- ecr::computeHV(archive, ref_point)
  
  expect_vector(hv_a, size = ncol(archive))
  expect_true(all(diff(hv_a) >= 0))
  expect_equal(tail(hv_a, 1), hv_ecr)
}

test_that("computeArchiveHV is consistent on simple improving problem", {
  ref_point <- c(1, 1)
  
  archive <- rbind(
    c(0.9, 0.8, 0.7, 0.6, 0.5),
    c(0.9, 0.8, 0.7, 0.6, 0.5)
  )
  
  test_one_archive(archive, ref_point)
})

test_that("computeArchiveHV is consistent on simple nondominated set", {
  ref_point <- c(1, 1)
  
  archive <- rbind(
    c(0.2, 0.4, 0.6, 0.8),
    c(0.8, 0.6, 0.4, 0.2)
  )
  
  test_one_archive(archive, ref_point)
})

test_that("computeArchiveHV works with range of different functions from bi-objective BBOB", {
  n_sample <- 1000L
  dims <- c(2L, 10L)
  fids <- 1L:55L
  iids <- 1L
  
  for (dim in dims) {
    for (fid in fids) {
      for (iid in iids) {
        fn_data <- generateBiObjBBOBData(dim, fid, iid)
        
        ref_point <- fn_data$ref_point
        fn <- fn_data$fn
        
        X <- matrix(runif(n_sample * dim, -5, 5), nrow = n_sample, ncol = dim)
        y <- apply(X, 1, fn)
        
        test_one_archive(y, ref_point)
      }
    }
  }
})

test_that("computeArchiveHV works with large sample numbers", {
  n_sample <- 1e7
  
  archive <- rbind(runif(n_sample), runif(n_sample))
  ref_point <- c(1, 1)
  
  test_one_archive(archive, ref_point)
})
