
get_samplesize <- function(range, fun, nmin=10, minsize=10000) {
	if (inherits(range, "SpatRaster")) {
		range <- terra::as.polygons(range)
		range <- range[range[,1,drop=TRUE]==1]
	}
	a <- terra::expanse(range, unit="km")
	n <- round(fun(a))
	n <- max(nmin, n)
	z <- max(1, min(n, round(a/minsize)))
	return(list(range=range, area=a, n=n, nzones=z))
}


branch_length <- function(tree, samp, adjust, adjfun=log10 ) {

## TODO take into account that you may need multiple obs per leaf

	if (adjust) {
		tab <- table(samp)
		tab <- adjfun(tab) + 1
	}
	sample <- as.character(unique(samp))
    if (is.null(tree$edge.length)) {stop("Tree has no branch lengths, cannot compute pd") }
    
	absent <- tree$tip.label[!(tree$tip.label %in% sample)]
    if (length(sample) == 0) {
        GD <- 0
    } else if (length(sample) == 1) {
        GD <- max(tree$edge.length)
    } else if (length(absent) == 0) {
		if (adjust) {
			i <- match(names(tab), tree$tip.label)
			tree$edge.length[i] <- tree$edge.length[i] * tab
		}
        GD <- sum(tree$edge.length)
    } else {
        tree <- ape::drop.tip(tree, absent)
		if (adjust) {
			i <- match(names(tab), tree$tip.label)
			tree$edge.length[i] <- tree$edge.length[i] * tab
        }
		GD <- sum(tree$edge.length)
    }
    return(GD)
}


get_cover <- function(regions, sample, env=NULL, adjust=TRUE) {

## TODO  RH
# fix the adjust effect such that when you have many observations in one zones
# they can only contribute to their neighbors. Do not increase branch length to avoid that
# one region does not compensate for another
		
	xy <- terra::centroids(regions)
	if (!is.null(env)) {
		# use the ClustGeo approach?
		e <- terra::extract(env, regions, fun=mean, na.rm=TRUE)
		d <- terra::distance(cbind(xy, e))
	} else {
		d <- terra::distance(xy, unit="km")
	}

	x <- stats::hclust(d)
	x <- ape::as.phylo(x)
	terra::values(regions) <- NULL
	
	# if nsamples for a region > threshold (10) one neighbor that is empty can get an observation
	##
	
	s <- terra::extract(regions, sample)

	actual_pd <- branch_length(x, s[,2], adjust=adjust)	
	potential_pd <- branch_length(x, 1:nrow(regions), adjust=adjust)
	min(1, actual_pd/potential_pd)
}


