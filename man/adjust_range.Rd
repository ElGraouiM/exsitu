\name{adjust_range}

\alias{adjust_range}

  
\title{adjust_range}

\description{
Adjust the range
}

\usage{
adjust_range(x, sp, CAmin=50000, CAmax=250000)
}

\arguments{
\item{x}{SpatRaster with species range predicted with a SDM (NA where absent)}
\item{sp}{SpatVector with species occurrences}
\item{CAmin}{positive number. Buffer size around each occurrence, to be included into the range}
\item{CAmax}{positive number. Buffer size around each occurrence, to be excluded from the range}
}

\value{
numeric 
}

\examples{
2 + 2
#adjust_range( ) 
}

\keyword{spatial}

