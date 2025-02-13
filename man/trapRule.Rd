\name{trapRule}
\alias{tradHeight}
\alias{tradT2M}
\alias{trapRule}
\title{Traditional calculation of PoIs}
\description{
  Traditional methods to calculate area under the curve (AUC), Height
  and Time2Max. For AUC the trapezoidal rule is used. If
  argument \code{maxT} is presented only points earlier than this time
  point are taken into account. Peak height is the difference between
  the highest point and the first point. Time2Max is the index of the
  time point at which the highest value is recorded. For vectors
  containing \code{NA} values, see the Details section.
}
\usage{
trapRule(y, x = 1:length(y), maxT = max(x))
tradHeight(y)
tradT2M(y)
}
\arguments{
  \item{y}{Height of the curve at different points.}
  \item{x}{Location of the curve points.}
  \item{maxT}{The maximal time taken into account when calculating the AUC.}
}
\value{The PoI, a number. Function \code{tradT2M} returns the index of
  the time point corresponding to the peak maximum, not the time itself.}
\details{If \code{NA} values are present then no AUC can be determined
  and \code{NA} is returned. For the other two PoIs \code{NA} values
  sometimes can be allowed: this is the case when the \code{NA} values
  are not located next to the estimated peak maximum, and for
  \code{Height} not at the first time points.
}
\author{Ron Wehrens}
\seealso{\code{\link{extractParameters}}}
\examples{
data(supro)
trapRule(supro$Thr[1:10], supro$Time[1:10])
tradHeight(supro$Thr[1:10])
tradT2M(supro$Thr[1:10])
}
\keyword{manip}
