\name{doPCA}
\alias{doPCA}
\alias{plotScores}
\alias{plotLoadings}
\title{PCA for aaresponse-type data.}
\description{
The \code{doPCA} function arranges, in several ways, the original data
in a data matrix and performs PCA. Functions to show loadings and scores
of the resulting object are available. Only columns corresponding to
single amino acids are taken into account - at least three different
ones are required.}
\usage{
doPCA(datf, type = c("par", "pro", "aa", "both"),
      rcenter = TRUE, rscale = FALSE)
plotScores(PCAobj, highlight = c("none", "Intervention", "AA", "Participant"),
           nrowLegend)
plotLoadings(PCAobj, ...)
}
\arguments{
  \item{datf}{Original data frame containing aaresponse data.}
  \item{type}{Determines what the rows of the data matrix are. See
    Details section.}
  \item{rcenter, rscale}{Logicals, related to scaling
    of the time series. Each time series will be scaled individually,
    and only after any scaling the data matrix on which PCA is performed
    will be constructed. The default is to center each time series
    around zero.}
  \item{PCAobj}{PCA object.}
  \item{highlight}{Show grouping in score plots.}
  \item{nrowLegend}{Number of rows in the figure legend (which is
    printed above the figure itself) - can be useful to plot figures
    side-by-side.}
  \item{\dots}{Additional arguments to the underlying \code{xyplot}
    function.}
}
\value{An object of class PCA, with the following elements:
\item{scores }{Scores of the PCA object. See Details section.}
\item{loadings}{Loadings of the PCA object. See Details section.}
\item{var}{Variances of PCs.}
\item{totalvar}{Total variance.}
\item{type}{Type of the data arrangement. See Details section.}
}
\details{The PCA objects that can be generated with \code{doPCA} differ
in how the original data matrices for PCA are generated.
}
\author{Ron Wehrens}
\seealso{\code{\link{checkAAdata}}}
\examples{
data(supro)
supro <- checkAAdata(supro) # mandatory

huhn  <- doPCA(supro, type = "par")
plotScores(huhn) # shows participants
plotLoadings(huhn) # time series for indivual AAs in proteins

huhn  <- doPCA(supro, type = "pro")
plotScores(huhn, highlight = "Intervention") # participant/protein combinations
plotScores(huhn, highlight = "Participant")
plotScores(huhn, highlight = "AA")
plotLoadings(huhn) # time series for individual AAs

huhn  <- doPCA(supro, type = "aa")
plotScores(huhn, highlight = "Intervention") # warning
plotScores(huhn, highlight = "Participant") # participant/AA combinations
plotScores(huhn, highlight = "AA", nrowLegend = 4)
plotLoadings(huhn) # time series for the individual proteins

huhn  <- doPCA(supro, type = "both")
plotScores(huhn, highlight = "Intervention")
plotScores(huhn, highlight = "Participant")
plotScores(huhn, highlight = "AA")
plotLoadings(huhn) # time series
}
\keyword{hplots}
\keyword{multivariate}
