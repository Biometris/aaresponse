\name{compareInterventions}
\alias{compareInterventions}
\title{Compare different protein meals (interventions) in terms of amino-acid levels in blood.}
\description{The function fits linear mixed models, tests for the
  significance of the \code{Period} factor, and calculates contrasts for
  different levels of the \code{Intervention} factor, comparing with the
  reference level of the intervention.}
\usage{
compareInterventions(
  dt, refIntervention, target = c("AUC", "Height", "Time2Max"),
  logTransform = FALSE,
  model.formula = "~ Period + Intervention + (1 | Participant)",
  lm.alternative = "~ Period + Intervention",
  mainFun = lmer, singularFun = lm, respondersOnly = TRUE)
}
\arguments{
  \item{dt}{A data frame with at least the columns in the \code{target},
    \code{fixedEffects}, and \code{randomEffects} parameters.}
  \item{refIntervention}{The name of the reference intervention.}
  \item{target}{The variable of interest.}
  \item{logTransform}{Boolean; whether to concentrate on ratios (logT =
    TRUE) or differences. Ratios are common when comparing AUCs.}
  \item{model.formula}{The default definition of the mixed model,
    right-hand side only (text) - the tilde is mandatory.}
  \item{lm.alternative}{The model that will be fit by \code{lm} when the
    mixed model is singular, right-hand side only (text) - the tilde
    is mandatory.}
  \item{mainFun}{The R function used to fit the mixed model.}
  \item{singularFun}{The R function used to fit the alternative model.}
  \item{respondersOnly}{Logical - if TRUE, any imputed values for AUC and
    Height are ignored. If no imputed values for these parameters are
    present, respondersOnly will be set to TRUE, whatever the actual
    value of the argument (with a warning).}
}
\value{A data frame containing estimated differences (or ratios, in the
  case of \code{logTransform == TRUE}) with the reference protein,
  confidence intervals, and p values.
}
\author{Jasper Engel and Ron Wehrens}
\seealso{
\code{\link{showCIs}}, \code{\link{fitMixedModels}},
\code{\link{testPeriodEffect}}, \code{\link{doComparisons}}
}
\examples{
data(pois)  
finalComparison <- compareInterventions(pois, target = "AUC",
                                        logTransform = TRUE)
}
\keyword{manip}
