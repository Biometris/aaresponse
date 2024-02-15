\name{trapRule}
\alias{trapRule}
\title{AUC calculation using the trapezoidal rule}
\description{
Calculates the area under the curve (AUC) using the trapezoidal rule. If
argument \code{maxT} is presented only points earlier than this time
point are taken into account.
}
\usage{
trapRule(y, x = 1:length(y), maxT = max(x))
}
\arguments{
  \item{y}{Height of the curve at different points.}
  \item{x}{Location of the curve points.}
  \item{maxT}{The maximal time taken into account when calculating the AUC.}
}
\value{The AUC, a number.}
\author{Ron Wehrens}
\seealso{\code{\link{extractParameters}}}
\examples{
data(supro)
trapRule(supro$Thr[1:10], supro$Time[1:10])
}
\keyword{manip}
