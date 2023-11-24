
remove_small <- function(range_adj, intr, d) {

  a <- terra::expanse(range_adj, unit="km")
  #d = p*a$area
  
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

