\name{showCIs}
\alias{showCIs}
\alias{showCombinedCIs}
\title{Visualize comparisons between protein meals, in the form of
  confidence intervals around estimated differences or ratios.
}
\description{The final outcome of a crossover experiment in protein
  uptake is often a comparison with a reference protein for the
  parameters of interest, AUC, peak height, and time to peak maximum. 
}
\usage{
showCIs(resultsTable, ...)
showCombinedCIs(lAUC = NULL, AUC = NULL, Height = NULL, Time2Max = NULL,
                relevantAAs, what, between, subset, ...)
}
\arguments{
  \item{resultsTable, lAUC, AUC, Height, Time2Max}{An object generated
    by function \code{compareInterventions} and containing, for each
    protein compared to the reference, means, standard deviations, upper
    and lower CI limits, and p values. In addition a p value may be
    present for the Period variable, but this will not be used in the CI
    plots.}
  \item{relevantAAs}{Names of amino acids to be shown
    (optional). Standard sets can be chosen using the \code{what}
    argument.}
  \item{what}{selects individual amino acids, essential amino acids,
    amino acid totals, or all (default).}
  \item{between}{\code{xyplot} argument - default is to separate AUC
    results on a ratio scale from the other PoIs that are always on a
    difference scale.}
  \item{subset}{Optional definition of a subset of amino acids.}
  \item{\dots}{Additional arguments to the underlying plot function,
    most often used to define a subset of the results to be displayed
    (see examples).}
}
\value{A lattice plot object is printed as a side effect, and returned
  invisibly.}
\author{Ron Wehrens}
\examples{
data(pois)

heightComparison <- compareInterventions(pois, target = "Height")
lAUCcomparison <- compareInterventions(pois, target = "AUC",
                                       logTransform = TRUE)

showCIs(lAUCcomparison)
showCombinedCIs(lAUC = lAUCcomparison, Height = heightComparison)
}
\keyword{manip}
