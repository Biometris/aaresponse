summary.XOdesign <- function(object, verbose = TRUE) {
  if (!is.null(object$Time))
    object <- object[object$Time == min(object$Time),]

  if (verbose) {
    ## add code here to check how often a particular treatment
    ## preceeds another one, should be a fairly simple tabulation
  } else {
    table(object$Intervention, object$Period)
  }
}

## strategy here: generate a balanced Latin Square design of nInterventions
## times nInterventions, and repeat that as many times as necessary
## (permuting columns) to cover nSubjects. The last design may miss
## one or more lines. In that way no consecutive combination should
## occur more than the other (if the number of subjects is a multiple
## of the number of proteins) or the max diff is 1, if this is not the
## case.
## Lit: Bradley, J. V. Complete counterbalancing of immediate
## sequential effects in a Latin square design.
## J. Amer. Statist. Ass., 1958, 53, 525-528.
## Also E.J. Williams, Australian J. Sci. Res., Series A, 2, 149-168 (1949)
## For odd-numbered Latin Squares this is not possible, and
## V.K. Sharma presents a very simple extension of Bradley's method
## using two Latin Squares

## Code below seems to work, needs testing
generateXOdesign <- function(nInterventions, nSubjects) {
  startfun <- function(start, np) {
    myseq <- rep(1:np, 2)
    myseq[start:(start + np - 1)]
  }
  
  odds <- (1:nInterventions)[c(TRUE, FALSE)]
  evens <- (1:nInterventions)[c(FALSE, TRUE)]
  startingPoint1 <- startingPoint2 <- 1:nInterventions
  startingPoint1[c(odds, rev(evens))] <- 1:nInterventions
  startingPoint2[c(evens, rev(odds))] <- 1:nInterventions
  
  LS1 <- sapply(startingPoint1, startfun, np=nInterventions)
  LS2 <- sapply(startingPoint2, startfun, np=nInterventions)
  LSA <- rbind(LS1, LS2)

  ## determine number of times this matrix needs to be repeated and resampled
  nrep <- ceiling(nSubjects / (2*nInterventions)) - 1
  if (nrep > 0) {
    newLS <- lapply(1:nrep, function(ii) LSA[,sample(1:nInterventions)])
    LSA <- rbind(LSA, do.call(rbind, newLS))
  }
  
  result <- melt(LSA, varnames = c("Participant", "Period"),
                 value.name = "Intervention")

  class(result) <- c("XOdesign", "data.frame")
  result
}
