trapRule <- function(y, x = 1:length(y), maxT = max(x)) {
  y <- y[x <= maxT]
  x <- x[x <= maxT]
  
  dx <- diff(x)
  dy <- abs(diff(y))
  miny <- pmin(head(y, -1), tail(y, -1))

  sum(dx*(miny + dy/2)) - y[1]*diff(range(x))
}

## More elaborate function for obtaining Height in the classical
## way. Basically the difference between the highest point and the
## first point. This function introduces ways to allow for some NAs.
## These are not allowed at the first point or next to the maximum,
## but are allowed elsewhere - in these cases NA is returned. If the
## maximum is at the final time point NA is returned as well.
tradHeight <- function(y) {
    if (is.na(y[1])) return(NA)

    ntime <- length(y)
    
    max.idx <- which.max(y, na.rm = TRUE)
    if (is.na(y[max.idx - 1]) |
        max.idx == ntime |
        (max.idx < ntime & is.na(y[max.idx + 1])))
        return(NA)
    
    y[max.idx] - y[1]
}

tradT2M <- function(y) {
    ntime <- length(y)
    max.idx <- which.max(y, na.rm = TRUE)

    if (is.na(y[max.idx - 1]) |
        max.idx == ntime |
        (max.idx < ntime & is.na(y[max.idx + 1]))) {

        NA
    } else {
        max.idx
    }
}

## Function to extract parameters of interest without using curve
## fitting: for AUC the trapezoid rule is used, for Time2Max the time
## with the highest response, and for Height the highest response
## minus the lowest response
## Based on function fitWoodAll

PoItradAll <- function(aadata,
                       what = c("all", "aas", "essentials", "totals"),
                       maxT) {
    what <- match.arg(what)
    relVars <- switch(what,
                      aas = attr(aadata, "aanames"),
                      totals = attr(aadata, "totalnames"),
                      essentials = intersect(attr(aadata, "aanames"),
                                             aaessentials()),
                      c(attr(aadata, "aanames"), attr(aadata, "totalnames")))
    idvars <- c("Participant", "Period", "Intervention", "Time")
    aadata.df <- melt(aadata[, c(idvars, relVars)], variable.name = "AA",
                      id.vars = c("Participant", "Period",
                                  "Time", "Intervention"))
    ## next lines to remove spurious levels
    aadata.df$AA <- factor(aadata.df$AA) 
    aadata.df$Participant <- factor(aadata.df$Participant)
    aadata.df$Intervention <- factor(aadata.df$Intervention)
    
    if (missing(maxT)) maxT <- max(aadata.df$Time)
    
    result <-
        aggregate(
            1:nrow(aadata.df),
            aadata.df[c("Participant", "AA", "Intervention", "Period")],
            function(ii) # divide by 15 to be consistent with cf approach
                c(AUC = trapRule(y = aadata.df[ii, "value"],
                                 x = aadata.df[ii, "Time"],
                                 maxT = maxT) / 15,
                  Height = tradHeight(aadata.df[ii, "value"]),
                  Time2Max =
                      levels(aadata.df$Time[tradT2M(aadata.df[ii, "value"])))))
    
    finalresult <- cbind(result[,1:4], as.data.frame(result$x))
    attr(finalresult, "aanames") <- intersect(attr(aadata, "aanames"),
                                              levels(finalresult$AA))
    attr(finalresult, "totalnames") <- intersect(attr(aadata, "totalnames"),
                                                 levels(finalresult$AA))
    attr(finalresult, "class") <- attr(aadata, "class")
    
    finalresult
}
