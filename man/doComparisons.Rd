\name{doComparisons}
\alias{doComparisons}
\title{Calculate contrasts for interventions.}
\description{
Given a list of linear (mixed) models, calculate contrasts and
confidence intervals.
}
\usage{
doComparisons(fitModels, logTransform)
}
\arguments{
  \item{fitModels}{A list of models generated by functions like lme,
    lmer and lm. In principle anything that is supported by emmeans may
    be used.}
  \item{logTransform}{Boolean - whether or not to test ratios (logT =
    TRUE) or differences.}
}
\value{The function returns a data.frame containing estimates and
  confidence intervals for the target PoI. The actual comparison is done
using the 'trt.vs.ctrl' contrasts from the \pkg{emmeans} package, which
employs the Dunnett multiple-testing adjustment.}
\author{Ron Wehrens}
\seealso{
  \code{\link{compareInterventions}}
}
\examples{
data(pois)
modList <- fitMixedModels(pois, target = "AUC")
result <- doComparisons(modList)
}
\keyword{manip}
