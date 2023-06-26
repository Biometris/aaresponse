\name{pois}
\alias{pois}
\docType{data}
\title{Parameters of interest for the supro data.}
\description{The result of applying \code{fitWoodAll} to the supro
  data, for the individual amino acids. Included because it is used in
  several examples.
}
\usage{data("pois")}
\format{
  A data frame with 684 observations on the following 7 variables.
  \describe{
    \item{\code{Participant}}{A factor corresponding to the participants
      in the study.}
    \item{\code{Intervention}}{The three proteins in the study:
      \code{Whey}, \code{BP} and \code{Corn}}
    \item{\code{AA}}{The amino acids present in the data.}
    \item{\code{Period}}{The periods when the meals were distributed - a
    factor with three levels.}
    \item{\code{AUC}, \code{Height}, \code{Time2Max}}{Numeric vectors; the
      parameters of interest in the study.}
  }
}
\examples{
data(pois)
head(pois)
}
\keyword{datasets}
