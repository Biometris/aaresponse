\name{testPeriodEffect}
\alias{testPeriodEffect}
\title{Function to test the significance of the fixed factor
  \code{Period} in the default formulation of the mixed models.}
\description{In the experimental data we do not expect \code{Period} to
  be significant, but we should check...}
\usage{
testPeriodEffect(fitModels)
}
\arguments{
  \item{fitModels}{A list of linear (mixed) models, one for each AA of
    interest.}
}
\value{Returns a vector of p values, one for each model. A warning is
  given if any of the p values (after multiple-testing correction using
  the Benjamini-Hochberg procedure) is smaller than 0.05.}
\author{Ron Wehrens}
\seealso{
  \code{\link{compareInterventions}}
}
\examples{
\dontrun{
data(pois)
modList <-
  fitMixedModels(
    pois,
    refIntervention = "WP",
    model.formula = "AUC ~ Period + Intervention + (1 | Participant)",
    lm.alternative = "AUC ~ Period + Intervention",
    mainFun = lmer,
    respondersOnly = TRUE)
pvals <- testPeriodEffect(modList)
}
}
\keyword{manip}
