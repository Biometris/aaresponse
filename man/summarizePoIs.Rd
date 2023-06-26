\name{summarizePoIs}
\alias{summarizePoIs}
\title{Summarize parameters of interest from Wood curves.}
\description{Summarize parameters of interest extracted from Wood
  curves, i.e., aggregate them across participants, and present the
  results as average values plus confidence intervals.}
\usage{
summarizePoIs(datf, target = c("AUC", "Height", "Time2Max"), coverage = 0.95)
}
\arguments{
  \item{datf}{A data.frame containing for individual participants
    estimates of the parameters of interest, given Intervention, amino
    acid, and Period. The output of a function like
    \code{extractParameters} or \code{curateParameters}.}
  \item{target}{The name of the PoI.}
  \item{coverage}{The coverage of the confidence interval. Default is
    0.95.}
}
\value{
A data.frame with the following columns: the PoI, \code{low}, \code{up},
\code{Intervention} and \code{AA}.
}
\author{Ron Wehrens}
\seealso{
\code{\link{extractParameters}}, \code{\link{curateParameters}},
\code{\link{showPoIs}}
}
\examples{
\dontrun{
  ## see example in showCIs
  data(supro)
  supro <- checkAAdata(supro)
  
  result <- fitWoodAll(supro, what = "totals")
  pois <- extractParameters(result)

  summarizePoIs(pois, "AUC")
}
}
\keyword{manip}
