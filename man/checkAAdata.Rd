\name{checkAAdata}
\alias{checkAAdata}
\title{
Function to check the syntactic correctness of data presented to the
functions of this package.
}
\description{
The function checks whether the necessary columns are present in the
data frame. These consist at least of the names \code{Participant},
\code{Period}, \code{Intervention}, \code{Time} and one or more columns
corresponding to amino-acid levels or totals. Running this function
before further analysis with functions in the package is mandatory -
information on the names of the amino acids and amino-acid totals in the
data set is attached to the data as attributes.
}
\usage{
checkAAdata(df, aanames, totalnames, quiet = FALSE)
}
\arguments{
  \item{df}{data frame containing all data.}
  \item{aanames, totalnames}{the names of individual amino acids or
    amino-acid totals, respectively. Vectors of strings.}
  \item{quiet}{Logical - set to TRUE if no text output is required.}
}
\value{The function returns an object of class \code{aar}, which
  inherits from \code{data.frame} but also has attributes \code{aanames}
  and \code{totalnames}. These attributes are preserved upon subsetting.
}
\author{Ron Wehrens}
\examples{
data(supro)
supro <- checkAAdata(supro)
attributes(supro)
}
\keyword{manip}
