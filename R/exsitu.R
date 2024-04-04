
n_zones <- function(x, min_area, m=3) {
	round(pmax(1, pmin(x, m*sqrt(x))))
}

get_samplesize <- function(range, fun=n_zones, ...) {
	if (inherits(range, "SpatRaster")) {
		range <- terra::as.polygons(range)
		range <- range[range[,1,drop=TRUE]==1]
	} else if (!inherits(range, "SpatVector")) {
		stop("range should be a SpatVector")
	}
	a <- terra::expanse(range, unit="km")
	n <- fun(a, ...)
	#n <- max(nmin, n)
	#z <- max(1, min(n, round(a/min_area)))
	return(list(range=range, area=a, n=n))
}


branch_length <- function(tree, samp, adjust, adjfun=log10 ) {

## TODO take into account that you may need multiple obs per leaf

	if (adjust) {
		tab <- table(samp)
		tab <- adjfun(tab) + 1
		# or? 
		# tab <- pmax(1, adjfun(tab))

	}
	sample <- as.character(unique(samp))
    if (is.null(tree$edge.length)) {stop("Tree has no branch lengths, cannot compute pd") }
    
	absent <- tree$tip.label[!(tree$tip.label %in% sample)]
    if (length(sample) == 0) {
        GD <- 0
    } else if (length(sample) == 1) {
		# also adjust for sample size
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


small_ssize_penalty <- function(ssize, score, minssize=10) {
# linear, make curvilinear
	ifelse(ssize > minssize, score, score * ssize / minssize)
}


ex_cvs_old <- function(regions, sample, env=NULL, adjust=FALSE, minssize=10) {

## TODO  RH
# fix the adjust effect such that when you have many observations in one zones
# they can only contribute to their neighbors. Do not increase branch length to avoid that
# one region does not compensate for another

	if (minssize <= 0) {
		return(0)
	}
		
	if (nrow(regions) == 1) {
		# cannot make a tree 
		s <- terra::extract(regions, sample)
		if (length(na.omit(terra::extract(regions, s)[,2])) > 1) { 
			return(small_ssize_penalty(length(sample), 1, minssize))
		} else {
			return(0)		
		}
	}
		
	xy <- terra::centroids(regions)
	d <- terra::distance(xy, unit="km")
	if (!is.null(env)) {
		e <- terra::extract(env, regions, fun=mean, na.rm=TRUE)
		ed1 <- dist(e[,2])
		ed2 <- dist(e[,3])		
	} 
	
	a <- cbind(g=as.vector(d), tmp=as.vector(ed1), prc=as.vector(ed2))
	
	
	x <- stats::hclust(d)
	x <- ape::as.phylo(x)
	terra::values(regions) <- NULL
	
	# if nsamples for a region > threshold (10) one neighbor that is empty can get an observation
	##
	
	s <- terra::extract(regions, sample)

	actual_pd <- branch_length(x, s[,2], adjust=adjust)	
	potential_pd <- branch_length(x, 1:nrow(regions), adjust=adjust)
	score <- min(1, actual_pd/potential_pd)
	
	small_ssize_penalty(length(sample), score, minssize)
	
}




ex_cvs <- function(regions, sample, env=NULL, envfun=NULL, adjust=FALSE, minssize=10) {

## TODO  RH
# fix the adjust effect such that when you have many observations in one zones
# they can only contribute to their neighbors. Do not increase branch length to avoid that
# one region does not compensate for another

	stopifnot(minssize > 0)

	if (nrow(sample) <= 0) {
		return(0)
	}
		
	if (nrow(regions) == 1) {
		# cannot make a tree 
		return(small_ssize_penalty(nrow(sample), 1, minssize))
	}
		
	xy <- terra::centroids(regions)
	d <- terra::distance(xy, unit="km")
	if (!is.null(env)) {
		e <- terra::extract(env, regions, fun=mean, na.rm=TRUE, ID=FALSE)
		ed1 <- dist(e[,1])
		ed2 <- dist(e[,2])	
		ed <- data.frame(ed1, ed2)
		names(ed) <- names(env)
		# express environmental distance expressed as geographic distance
		envd <- envfun(ed)
		envd <- structure(envd, class = 'dist', Size=attr(ed1, "Size"))
		# env dist cannot be smaller than geo dist
		envd[envd < d] <- d[envd < d]
		# mean of geo and env dist
		d <- (d + envd) / 2
	} 

	x <- stats::hclust(d)
	x <- ape::as.phylo(x)
	terra::values(regions) <- NULL
	
	# if nsamples for a region > threshold (10) one neighbor that is empty can get an observation
	##
	
	s <- terra::extract(regions, sample)

	actual_pd <- branch_length(x, s[,2], adjust=adjust)	
	potential_pd <- branch_length(x, 1:nrow(regions), adjust=adjust)
	score <- min(1, actual_pd/potential_pd)
	
	small_ssize_penalty(length(sample), score, minssize)
	
}


