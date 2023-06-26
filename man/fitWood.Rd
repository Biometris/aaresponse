\name{fitWood}
\alias{fitWood}
\alias{fitWoodAll}
\title{
Fit a Wood curve through a time series.
}
\description{
A Wood curve is given by y(t) = d + a t^{mc} e^{-ct} where d is the
baseline, m is the time to the maximum, c is related to the rise and
fall of the curve, and a to the height. The fit is done by numerical
optimization using 500 random starts. In normal usage the \code{fitWood}
function will not be called directly by a user.
}
\usage{
fitWood(df)
fitWoodAll(aadata, what = c("all", "aas", "essentials", "totals"))
}
\arguments{
  \item{df}{A two-column data.frame, the first column containing time,
    and the second the level.}
  \item{aadata}{A data frame containing at least the columns
    \code{Participant}, \code{Period}, \code{Intervention}, \code{AA},
    \code{Time} and one or more amino-acid columns.}
  \item{what}{Easy selection of individual amino acids and AA totals.}
}
\value{Function \code{fitWood} returns a vector containing the parameter
  estimates and the final RMS value. \code{fitWoodAll} is mainly a
  wrapper calculating the parameters for all combinations of
  Participant, Period, Intervention and AA (using \code{aggregate}). Its
  result is a data frame where the list element \code{x} contains the
  fitted parameters.
}
\references{
Wood Nature (1967)
}
\author{
Ron Wehrens
}
\seealso{
  \code{\link{showCIs}}, \code{\link{extractParameters}},
  \code{\link{showDataFits}}
}
\examples{
data(supro)
idx <- supro[,"Participant"] == 1 & supro[,"Period"] == 1
df <- supro[idx, c("Time", "His")]
colnames(df) <- c("Time", "value")
result <- fitWood(df)

\dontrun{
  ## see example in showCIs
  supro <- checkAAdata(supro) # to include info on aanames and totalnames
  result.totals <- fitWoodAll(supro, what = "totals")
}
}
\keyword{manip}
