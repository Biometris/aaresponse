\name{imputePoIs}
\alias{imputePoIs}
\title{Impute values for AUC and Height}
\description{If no meaningful Wood curve can be fitted, or if PoI values
  are outside permitted ranges, NA values will be present in the
  data.frame resulting from the analysis. For \code{AUC} and
  \code{Height}, a robust linear regression is employed relating the
  spread of the time course data points to \code{Height} and \code{AUC}
  values. For \code{Time2Max} no meaningful way of doing this exists, so
  there the NA values will remain in the data.
}
\usage{
imputePoIs(df, params, quant = 0.2)
}
\arguments{
  \item{df}{Data frame containing all data.}
  \item{params}{Data frame containing the parameters of interest,
    probably the output of \code{extractParameters} or
    \code{curateParameters}.
}
  \item{quant}{The spread of the time course data is calculated as a
    difference between quantiles - this parameter indicates the
    quantiles used. The default value therefore indicates the difference
    between the quintiles.}
}
\details{
For each amino acid a separate regression line is fitted - in some data
sets we have seen, though, that the lines do not differ substantially.
}
\value{A data fram like the \code{params} input argument. Additional
  columns are present: \code{AUC.i} and \code{Height.i} are the values
  calculated by linear regression (for all records in the data frame), and
  \code{AUC.orig} and \code{Height.orig} are the original results from
  curve fitting, possibly including NA values. In \code{AUC} and
  \code{Height} the \code{NA} values have been replaced by the
  corresponding values in \code{AUC.i} and \code{Height.i}, respectively.}
\author{Ron Wehrens}
\examples{
data(supro)
data(pois)
pois.cur <- curateParameters(pois)
pois.complete <- imputePoIs(supro, pois.cur)
}
\keyword{manip}
