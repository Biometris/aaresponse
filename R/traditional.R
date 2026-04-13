imputePreZero <- function(y, x) {
    ## NA values are allowed at or before t=0, if other points in that
    ## time range are non-NA
    prezero <- which(x <= 0)

    if (!all(is.na(y[prezero]))) {
        if (any(is.na(y[prezero]))) {
            na.idx <- which(is.na(y[prezero]))
            nonna.idx <- prezero[-na.idx]
            y[na.idx] <- mean(y[nonna.idx])
        }
    }
    
    y
}

## As the baseline, the first value in y is taken.
## Update Dec 18 2025: not the first point is taken as the baseline,
## but the average of all points at or before t=0. If x[1] > 0, then
## the first point is taken, indeed.
trapRule <- function(y, x = 1:length(y), maxT = max(x),
                     imputePreZero = TRUE) {
    if (imputePreZero)
        y <- imputePreZero(y, x)
    
    ## if (any(is.na(y))) return(NA)

    y <- y[x <= maxT]
    x <- x[x <= maxT]
    
    dx <- diff(x)
    dy <- abs(diff(y))
    miny <- pmin(head(y, -1), tail(y, -1))

    if (all(x > 0)) {
        zero.idx <- 1
    } else {
        zero.idx <- which(x <= 0)
    }

    sum(dx*(miny + dy/2)) - mean(y[zero.idx])*diff(range(x))
}

## More elaborate function for obtaining Height in the classical
## way. Basically the difference between the highest point and the
## first point. This function introduces ways to allow for some NAs.
## These are not allowed at points before or at zero (unless
## imputePreZero == TRUE), or next to the maximum,
## but are allowed elsewhere. In the cases where NA is not allowed,
## the function returns NA. If the
## maximum is at the final time point NA is returned as well.

## Update Dec 18 2025: not the first point is taken as the baseline,
## but the average of all points at or before t=0. If x[1] > 0, then
## the first point is taken, indeed.
tradHeight <- function(y, x = 1:length(y), imputePreZero = TRUE) {
    if (imputePreZero)
        y <- imputePreZero(y, x)
    
    if (is.na(y[1])) return(NA)

    ntime <- length(y)
    
    max.idx <- which.max(y)
    if ((max.idx > 1 && is.na(y[max.idx - 1])) ||
        (max.idx < ntime && is.na(y[max.idx + 1])))
        return(NA)

    zero.idx <- ifelse(all(x > 0), 1, which(x <= 0))
    
    y[max.idx] - y[zero.idx]
}

tradT2M <- function(y, x = 1:length(y), imputePreZero = TRUE) {
    if (imputePreZero)
        y <- imputePreZero(y, x)

    ntime <- length(y)
    max.idx <- which.max(y)

    if (all(is.na(y)) ||
        (max.idx > 1 && is.na(y[max.idx - 1])) ||
        (max.idx < ntime && is.na(y[max.idx + 1]))) {

        NA
    } else {
        x[max.idx]
    }
}

## Function to extract parameters of interest without using curve
## fitting: for AUC the trapezoid rule is used, for Time2Max the time
## with the highest response, and for Height the highest response
## minus the lowest response
## Based on function fitWoodAll

PoItradAll <- function(aadata,
                       what = c("all", "aas", "essentials", "totals"),
                       maxT, imputePreZero = TRUE) {
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
            function(ii) {
                ## AUC: divide by 15 to be consistent with cf approach
                c(AUC = trapRule(y = aadata.df[ii, "value"],
                                 x = aadata.df[ii, "Time"],
                                 maxT = maxT,
                                 imputePreZero = imputePreZero) / 15,
                  Height = tradHeight(y = aadata.df[ii, "value"],
                                      x = aadata.df[ii, "Time"],
                                      imputePreZero = imputePreZero),
                  Time2Max = tradT2M(y = aadata.df[ii, "value"],
                                     x = aadata.df[ii, "Time"],
                                     imputePreZero = imputePreZero))
            })
    
    finalresult <- cbind(result[,1:4], as.data.frame(result$x))
    attr(finalresult, "aanames") <- intersect(attr(aadata, "aanames"),
                                              levels(finalresult$AA))
    attr(finalresult, "totalnames") <- intersect(attr(aadata, "totalnames"),
                                                 levels(finalresult$AA))
    attr(finalresult, "class") <- attr(aadata, "class")
    
    finalresult
}
