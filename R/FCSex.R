
GRSex <- function(srange, ca) {
	bufr <- terra::mask(srange, ca)
	e <- terra::expanse(c(srange, bufr))
	e[2,2] / e[1,2]
}

ERSex <- function(srange, ca, ecoregions) {
	eco <- terra::mask(ecoregions, srange) 
	ueco <- nrow(terra::unique(eco))
	seedeco <- terra::mask(eco, ca)
	seco <- nrow(terra::unique(seedeco))
	if(is.null(seco)) return(0)
	seco / ueco
}

SRSex <- function(s, h) {
	if (s == 0) return(0)
	s / (h + s)
}

FCSex <- function(seeds, herbarium, srange, ecoregions, bsize=50000) {
	srange <- terra::subst(srange, 0, NA)
	ca <- terra::buffer(seeds, bsize) |> terra::aggregate()  
	ca <- terra::rasterize(ca, srange)
	g <- GRSex(srange, ca)
	e <- ERSex(srange, ca, ecoregions)
	nseed <- nrow(seeds)
	nherb <- NROW(herbarium)
	s <- SRSex(nseed, nherb)
	c(nseed=nseed, nherb=nherb, GRSex=g, ERSex=e, SRSex=s, FCS=mean(c(g, e, s)))
}

