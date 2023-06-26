\name{curateParameters}
\alias{curateParameters}
\alias{curateFits}
\title{Remove unrealistic peak parameters.}
\description{Some time profiles may not correspond to meaningful peaks
  at all - this function attempts to remove them from the data by
  inspection of the values of the parameters of interest. These may be
  provided by the user. Invoking the function without arguments will
  lead to the defaults being printed. Note that the defaults are
  probably more appropriate for individual amino acids than for
  amino-acid totals (see last example).
}
\usage{
curateFits(fitparams, aa = c(0, 1000), mm = c(0, 100), cc = c(0, 100),
           dd = c(0, 1000), verbose = TRUE)
curateParameters(params, Time2Max = c(15, 200), Height = c(0, 1000),
                 AUC = c(0, Inf), verbose = TRUE)
}
\arguments{
  \item{fitparams, params}{Data frame containing the fitted parameter
    values, or }
  \item{aa, mm, cc, dd, Time2Max, Height, AUC}{Vectors of minimal and
    maximal values.}
  \item{verbose}{Logical, indicating whether feedback is given to the
    terminal about the number of curated records.}
}
\value{A data.frame is returned where the parameter values in all rows
  violating the boundary conditions are set to NA.}
\author{Ron Wehrens}
\examples{
data(pois)
pois.curated <- curateParameters(pois)
}
\keyword{manip}
