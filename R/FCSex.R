#In GapAnalysis::GRSex the CA50 areas outside of the SDM range are considered for the genebank area and in GapAnalysis::ERSex the number of ecoregions considered for the genebank samples includes areas that are outside the SDM. 
# In both cases you could end up with a score > 1. That can be corrected for, but this is inconsistent. Areas outside the range (SDM) should not be included in the analysis.
# They are not included in the analysis below (but can be with "inrange=FALSE").


GRSex <- function(srange, ca, inrange) {
	if (inrange) {
		bufr <- terra::mask(srange, ca)
		e <- terra::expanse(c(srange, bufr), unit="km")
		e[2,2] / e[1,2]
	} else {
		bufr <- terra::expanse(ca, unit="km")
		e <- terra::expanse(srange, unit="km")
		bufr / e[1,2]
		#min(1, bufr / e)
	}
}

ERSex <- function(srange, ca, ecoregions, inrange) {
	eco <- terra::mask(ecoregions, srange) 
	ueco <- nrow(terra::unique(eco))
	if (inrange) {
		seedeco <- terra::mask(eco, ca)
	} else {
		seedeco <- terra::mask(ecoregions, ca)	
	}
	seco <- nrow(terra::unique(seedeco))
	if(is.null(seco)) return(0)
	seco / ueco
	#min(1, seco / ueco)
}

SRSex <- function(s, h) {
	if (s == 0) return(0)
	s / (h + s)
}

FCSex <- function(seeds, herbarium, srange, ecoregions, bsize=50000, inrange=TRUE) {
	srange <- terra::subst(srange, 0, NA)
	ca <- terra::buffer(seeds, bsize) |> terra::aggregate()  
	ca <- terra::rasterize(ca, srange)
	g <- GRSex(srange, ca, inrange)
	e <- ERSex(srange, ca, ecoregions, inrange)
	nseed <- nrow(seeds)
	nherb <- NROW(herbarium)
	s <- SRSex(nseed, nherb)
	c(nseed=nseed, nherb=nherb, GRSex=g, ERSex=e, SRSex=s, FCS=mean(c(g, e, s)))
}

