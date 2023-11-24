\name{FCSex}

\alias{FCSex}

  
\title{FCSex}

\description{
Compute Khoury et al.'s FCSex and its components
}

\usage{
FCSex(seeds, herbarium, srange, ecoregions, bsize=50000)
}

\arguments{
\item{seeds}{SpatVector with seed bank collection locations}
\item{herbarium}{SpatVector with herbarium collection locations}
\item{srange}{SpatRaster with the range (all values > 0 are in the range)}
\item{ecoregions}{SpatRaster with ecoregions (each ecoregion should have a unique ID)}
\item{bsize}{positive number. Buffer size in m}
}

\value{
numeric 
}

\examples{
2 + 2
}

\keyword{spatial}

