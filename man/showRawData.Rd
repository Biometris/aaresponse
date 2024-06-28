\name{showRawData}
\alias{showRawData}
\title{
Plot (subsets of) the raw aaresponse data.
}
\description{
Shows raw data in a variety of ways, highlighting particular characteristics.
}
\usage{
showRawData(datf, type = c("aa", "participant", "sequential"),
            relevantAAs,
            what = c("all", "aas", "essentials", "totals"),
            scale = list(y = "free"),
            xlab = "Time (mins)",
            ylab = expression(paste("AA level (", mu, "M)", sep = "")), ...)
}
\arguments{
  \item{datf}{Data frame containing the required columns. These are
    \code{Participant}, \code{Intervention}, \code{Period}, \code{Time}
    and one or more columns corresponding to amino acids or amino-acid
    totals.}
  \item{type}{Determines the character of the plot.}
  \item{relevantAAs}{Names of amino acids to be shown
    (optional). Standard sets can be chosen using the \code{what}
    argument.}
  \item{what}{selects individual amino acids, essential amino acids,
    amino acid totals, or all (default).}
  \item{xlab, ylab}{Default values for axis labels.}
  \item{scale,\dots}{Additional arguments to the internal xyplot
    function.}
}
\author{Ron Wehrens}
\seealso{
  \code{\link{checkAAdata}}
}
\examples{
data(supro)
supro <- checkAAdata(supro)

showRawData(supro, what = "essentials")

showRawData(supro, type = "aa")
showRawData(supro, type = "aa", subset = Participant == 1) 
showRawData(supro, type = "aa", subset = Participant \%in\% 1:3) 
showRawData(supro, type = "participant")
showRawData(supro, type = "participant",
            subset = Participant \%in\% 1:3 & AA \%in\% c("Ile", "Val", "Thr"))
showRawData(supro, type = "sequential")
showRawData(supro, type = "sequential",
            subset = Participant \%in\% 1:3 & AA \%in\% c("Ile", "Val", "Thr"))
}
\keyword{hplot}
