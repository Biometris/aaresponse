\name{PoItradAll}
\alias{PoItradAll}
\title{
Determine PoIs based on traditional methods
}
\description{
The three PoIs, AUC, Height and Time2Max, are determined in this
function using the traditional methods: the trapezoid rule for AUC, the
difference between the highest and lowest values in the time series, and
the time point corresponding to the highest value.
}
\usage{
PoItradAll(aadata, what = c("all", "aas", "essentials", "totals"),
           maxT)
}
\arguments{
  \item{aadata}{A data frame containing at least the columns
    \code{Participant}, \code{Period}, \code{Intervention}, \code{AA},
    \code{Time} and one or more amino-acid columns.}
  \item{what}{Easy selection of individual amino acids and AA totals.}
  \item{maxT}{The maximum for the integration of the curve, expressed in
    minutes. If unspecified the latest time point is used.}
}
\value{An object of class "aar", containing the PoIs for all
  participants in all conditions for the selected amino acids.
}
\details{The value of AUC is 1/15th of the value obtained by function
  \code{trapRule} to be consistent with the current AUC calculations by
  the curve fitting method. This may be changed in the
  future. \code{Time2Max} is given in the units that are the levels of
  \code{aadata$Time}.
}
\author{
Ron Wehrens
}
\seealso{
  \code{\link{trapRule}}
}
\examples{
data(supro)
df <- checkAAdata(supro, quiet = TRUE)
result <- PoItradAll(df, what = "totals")
}
\keyword{manip}
