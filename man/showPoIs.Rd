\name{showPoIs}
\alias{showPoIs}
\title{Plot parameters of interest, extracted from Wood curves}
\description{Visualization of estimates for the parameters of interest,
  extracted from Wood curves. Any subset of AUC, Height, and Time2Max
  may be given.}
\usage{
showPoIs(AUC = NULL, Height = NULL, Time2Max = NULL, ...)
}
\arguments{
  \item{AUC, Height, Time2Max}{A data.frame with not only a column
    corresponding to the name of the PoI, but also for the lower and
    upper values of the confidence interval (\code{low} and \code{up}),
    a column \code{Intervention}, and a column \code{AA}, such as give
    by the output of \code{summarizePoIs}.}
  \item{\dots}{Other arguments to the underlying \code{xyplot}
    function. Notable the \code{subset} argument can be useful.}
}
\value{A lattice plot is returned.}
\author{Ron Wehrens}
\seealso{
  \code{\link{summarizePoIs}}
}
\examples{
\dontrun{
  ## see example in showCIs
  data(supro)
  supro <- checkAAdata(supro)
  
  result <- fitWoodAll(supro, what = "totals")
  pois <- extractParameters(result)

  pois.summ <- summarizePoIs(pois, "AUC")
  showPoIs(AUC = pois.summ)
}
}
\keyword{graphics}
