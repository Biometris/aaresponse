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
  amino-acid totals (see last example) - if the function can determine
  that both individual AAs and AA totals are processed simultaneously,
  it will give a warning.
}
\usage{
curateFits(fitparams, aa = c(0, 1000), mm = c(0, 100), cc = c(0, 100),
           dd = c(0, 1000), verbose = TRUE)
curateParameters(params, Time2Max = c(15, 200), Height = c(0, 1000),
                 AUC = c(0, Inf), allowNAs = FALSE, verbose = TRUE)
}
\arguments{
  \item{fitparams, params}{Data frame containing the fitted parameter
    values, or }
  \item{aa, mm, cc, dd, Time2Max, Height, AUC}{Vectors of minimal and
    maximal values.}
  \item{allowNAs}{Logical, indicating whether some of the parameters may
    be \code{NA} whereas others have a value. For parameters obtained
    from curve fits, this always should be \code{FALSE}, which is also
    the default. For traditional parameter estimates one can imagine
    that it is possible to estimate peak position and peak height but
    not AUC, so in such a case this should be set to \code{TRUE}.}
  \item{verbose}{Logical, indicating whether feedback is given to the
    terminal about the number of curated entries.}
}
\value{A data.frame is returned where the parameter values in all rows
  violating the boundary conditions are set to NA.}
\details{If curation is attempted to both individual AAs and AA totals
  simultaneously, a warning is issued - almost always these need
  different curation limits.}
\author{Ron Wehrens}
\examples{
data(pois)
pois.curated <- curateParameters(pois)
}
\keyword{manip}
