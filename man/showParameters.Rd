\name{showParameters}
\alias{showParameters}
\title{Provides a dotplot of estimated parameter values.}
\description{The function can be used to visualize parameters of
  interest, like AUC or peak height values, but also curve fit
  parameters used in the Wood curves.}
\usage{
showParameters(params, relevantAAs,
               what = c("all", "aas", "essentials", "totals"), ...)
}
\arguments{
  \item{params}{Data frame of parameter values, for each participant and
    amino acid.}
  \item{relevantAAs}{Names of amino acids to be shown
    (optional). Standard sets can be chosen using the \code{what}
    argument.}
  \item{what}{selects individual amino acids, essential amino acids,
    amino acid totals, or all (default).}
  \item{\dots}{Additional arguments to the underlying (lattice) plotting
    function. Particularly useful when defining subsets, e.g., of
    participants, or of amino acids.}
}
\author{Ron Wehrens}
\seealso{
  \code{\link{showCIs}}, \code{\link{showDataFits}}
}
\examples{
data(supro)
supro <- checkAAdata(supro)

result <- fitWoodAll(supro, what = "totals")
showParameters(result)

data(pois)
showParameters(pois)
}
\keyword{manip}
