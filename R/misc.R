
remove_small <- function(range_adj, intr, d) {

  a <- terra::expanse(range_adj, unit="km")
  intr <- terra::disagg(intr)
  intr$area <- terra::expanse(intr, unit="km")
  x = intr[intr$area> d,]
  y = intr[intr$area<= d,]
  y = terra::disagg(terra::aggregate(y, dissolve=TRUE))
  terra::values(x) = NULL
  terra::values(y) = NULL
  if (nrow(y) > 0 & nrow(x) > 0) {
    intr <- terra::combineGeoms(x, y, boundary=TRUE, distance=TRUE)
    intr$area <- terra::expanse(intr, unit="km")
  }
  intr
}


adjust_range <- function(x, sp, th, CAmin=100000, CAmax=250000) { 
	r <- x > th
	ca_add <- terra::buffer(sp, CAmin) |> terra::aggregate()
	ca_remove <- terra::buffer(sp, CAmax) |> terra::aggregate()
	m <- terra::mask(r, ca_remove, updatevalue=FALSE)
	r <- (m + r) > 1
	terra::rasterize(ca_add, r, update=TRUE)
}

