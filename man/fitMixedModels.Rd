\name{fitMixedModels}
\alias{fitMixedModels}
\title{Fit mixed models to food trial data.}
\description{The function fits a model for the parameter of interest
  (indicated by \code{what}). The default is to use fixed factors for
  \code{Intervention} and \code{Period}, and take \code{Participant} as
  random. If the resulting model is singular, the model is refit using
  only fixed terms.
}
\usage{
fitMixedModels(dt, refIntervention,
               model.formula, lm.alternative,
               mainFun = lmer, respondersOnly = TRUE)
}
\arguments{
  \item{dt}{Data frame containing the relevant information.}
  \item{refIntervention}{Level of the \code{Intervention} factor to be
    used as a reference. If not specified, the first level is used.}
  \item{model.formula}{An object of type \code{formula}, describing the
    mixed model.}
  \item{lm.alternative}{An object of type \code{formula}, describing the
    model using only fixed effects (fit by \code{lm}.} 
  \item{mainFun}{Function used to fit the mixed model.}
  \item{respondersOnly}{Boolean: if TRUE, no imputed values are used.}
}
\value{A list of models, one for each amino acid in the input data frame.}
\author{Ron Wehrens}
\seealso{
  \code{\link{compareInterventions}}
}
\examples{
\dontrun{
  data(pois)
  modList <- fitMixedModels(pois, target = "AUC")
}
}
\keyword{manip}
