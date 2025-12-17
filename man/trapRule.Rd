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
trapRule(y, x = 1:length(y), maxT = max(x), imputePreZero = TRUE)
tradHeight(y, x = 1:length(y), imputePreZero = TRUE)
tradT2M(y, x = 1:length(y), imputePreZero = TRUE)
}
\arguments{
  \item{y}{Height of the curve at different points.}
  \item{x}{Location of the curve points.}
  \item{maxT}{The maximal time taken into account when calculating the
    AUC.}
  \item{imputePreZero}{If multiple time points are taken before or at
    t=0, missing values at these time points are replaced by the average
    of the non-missing values. If all values are missing, NA is
    returned. Note that the default x values will not allow imputation,
    so x must be specified explicitly. }
}
\value{The PoI, a number. Function \code{tradT2M} returns the index of
  the time point corresponding to the peak maximum, not the time itself.}
\details{If \code{NA} values are present (also after imputing preZero
  time points) then no AUC can be determined
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
