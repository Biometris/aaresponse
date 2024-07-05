\name{extractParameters}
\alias{extractParameters}
\alias{getAUC}
\alias{getHeight}
\alias{getTime2Max}
\title{
Extract parameters of interest from Wood curves.
}
\description{
One of the advantages of the Wood curve to describe time courses of
amino-acid levels in blood is that analytical expressions exist for
extracting the parameters of interest: AUC, peak height, and the time to
the peak maximum. These functions take the four curve-fit parameters and
return the parameter of interest.
}
\usage{
getAUC(prs, maxT = 300)
getHeight(prs)
getTime2Max(prs, minutes = TRUE)
extractParameters(prs.df, maxT = 300, minutes = TRUE)
}
\arguments{
  \item{prs}{Either a vector or a one-row data.frame containing the
    names "a", "m", "c" and "d".}
  \item{prs.df}{A multi-row data.frame containing the
    names "a", "m", "c" and "d" in list element "x".}
  \item{maxT}{The maximum for the integration of the curve, expressed in
    minutes.} 
  \item{minutes}{Whether the time to the peak maximum should be
    expressed in minutes. Default is TRUE.}
}
\value{The simple extraction functions return scalars, function
  \code{extractParameters} returns an object of class 'aar'.}
\author{Ron Wehrens}
\seealso{
  \code{\link{checkAAdata}}, \code{\link{extractParameters}}, \code{\link{showCIs}}
}
\examples{
\dontrun{
  ## see example in showCIs
  data(supro)
  supro <- checkAAdata(supro)
  
  result <- fitWoodAll(supro, what = "totals")
  pois <- extractParameters(result)
}
}
\keyword{manip}
