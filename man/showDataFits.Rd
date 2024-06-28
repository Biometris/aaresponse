\name{showDataFits}
\alias{showDataFits}
\title{Shows data points and the fitted curve in one plot.}
\description{Plots time-course data as well as fitted curves. Only
  profiles and time courses that are present in both the data and fit
  objects are shown.}
\usage{
showDataFits(datf, pardf, relevantAAs,
             what = c("all", "aas", "essentials", "totals"),
             xlab = "Time (mins)",
             ylab = expression(paste("AA level (", mu, "M)", sep = "")),
             scale = list(y = "free"), baseLineCorr = FALSE,
             points = TRUE, noLegend = FALSE, ...)
}
\arguments{
  \item{datf}{Data frame containing the raw data points. 
    Columns "AA", "Intervention", "Time" and "Participant" must be present.}
  \item{pardf}{Data frame object containing the fit parameters. }
  \item{relevantAAs}{Names of amino acids to be shown
    (optional). Standard sets can be chosen using the \code{what}
    argument.}
  \item{what}{selects individual amino acids, essential amino acids,
    amino acid totals, or all (default).}
  \item{baseLineCorr}{Logical. Setting it to TRUE will have all curves
    start from zero.}
  \item{points}{Logical. Setting it to FALSE shows only the fitted
    curves.}
  \item{noLegend}{Logical, to suppress showing a legend. This can be
    useless in some cases where subsetting is applied.}
  \item{xlab, ylab, scale, \dots}{addditional arguments to the
    underlying plotting functions.}
}
\author{Ron Wehrens}
\seealso{\code{\link{showRawData}}}
\examples{
data(supro)
supro <- checkAAdata(supro)

result <- fitWoodAll(supro[supro$Participant \%in\% 1:3,], what = "totals")
showDataFits(supro, result)
showDataFits(supro, result, baseLineCorr = TRUE, points = FALSE)
}
\keyword{hplot}
