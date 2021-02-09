#' @export
create_hv_function <- function(fn, ref_point, geom_mean = FALSE) {
  function(x) {
    y <- fn(x)
    
    if (geom_mean) {
      prod(sqrt(pmax(ref_point - y, 0)))
    } else {
      prod(pmax(ref_point - y, 0))
    }
  }
}

#' @export
hv_from_obj_space <- function(Y, ref_point, geom_mean = FALSE) {
  if (geom_mean) {
    apply(Y, 1, function(y) {
      prod(sqrt(pmax(ref_point - y, 0)))
    })
  } else {
    apply(Y, 1, function(y) {
      prod(pmax(ref_point - y, 0))
    })
  }
}
