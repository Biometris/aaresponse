woodFun <- function(params, t,
                    woodFor = formula(y ~ a*(t^(m*c))*exp(-c*t) + d)) {
  ##   a*(t^(m*c))*exp(-c*t) + d
  cmd <- tail(as.character(woodFor),1)
  exp <- parse(text=cmd)
  eval(exp, c(as.list(params), list(t = t)))
}

fitWood <- function(df) {
  woodFor <- formula(y ~ a*(t^(m*c))*exp(-c*t) + d)

  df <- df[!is.na(df$value),]
  
  t <- df$Time / 15 ## for better numerical stability
  y <- df$value
  n <- length(y)

  if (length(y) < 8) {
    result <- rep(NA, 5)
  } else {
    huhn <-
      nls_multstart(woodFor,
                    iter = 500,
                    start_lower = c(a = 0, m = 0, c = 0, d = 0),
                    start_upper = c(a = diff(range(y)),
                                    m = max(t), c = 3,
                                    d = max(y)), supp_errors = "Y",
                    lower = c(a = 0, b = 0, c = 0, d = 0))
    
    pars <- coef(huhn)
    
    nullssq <- var(y) * (n-1)
    ssq <- sum((predict(huhn) - y)^2)
    if (ssq > nullssq) {## Fit baseline only
      pars <- c(rep(0, 3), mean(y))
      rms <- sqrt(nullssq / length(y))
    } else {
      rms <- sqrt(ssq / length(y))
    }
    
    result <- c(pars, rms)
  }
  
  names(result) <-
    c(letters[c(1,13,3,4)], "RMS")
  result
}

fitWoodAll <- function(aadata, what = c("all", "aas", "essentials", "totals")) {
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
  
  result <-
    aggregate(1:nrow(aadata.df),
              aadata.df[c("Participant", "AA", "Intervention", "Period")],
              function(ii) fitWood(aadata.df[ii,c("Time", "value")]))

  
  finalresult <- cbind(result[,1:4], as.data.frame(result$x))
  attr(finalresult, "aanames") <- intersect(attr(aadata, "aanames"),
                                        levels(finalresult$AA))
  attr(finalresult, "totalnames") <- intersect(attr(aadata, "totalnames"),
                                           levels(finalresult$AA))
  attr(finalresult, "class") <- attr(aadata, "class")

  finalresult
}

curateFits <- function(fitparams,
                       aa = c(0, 1000),
                       mm = c(0, 100),
                       cc = c(0, 100),
                       dd = c(0, 1000),
                       verbose = TRUE) {
  if (missing(fitparams)) {
    boundaryTable <- rbind(aa, mm, cc, dd)
    dimnames(boundaryTable) <-
      list(c("a", "m", "c", "d"),
           c("Min", "Max"))
    return(boundaryTable)
  }

  if (length(attr(fitparams, "aanames")) > 0 &
      length(attr(fitparams, "totalnames")) > 0)
    warning("Applying the same curation boundaries to individual AAs and AA totals")

  bad.idx <-
    is.na(fitparams[,"a"]) | is.na(fitparams[,"m"]) |
    is.na(fitparams[,"c"]) | is.na(fitparams[,"d"]) |
    fitparams[, "a"] < aa[1] | fitparams[, "a"] > aa[2] |
    fitparams[, "m"] < mm[1] | fitparams[, "m"] > mm[2] |
    fitparams[, "c"] < cc[1] | fitparams[, "c"] > cc[2] |
    fitparams[, "d"] < dd[1] | fitparams[, "d"] > dd[2]
  
  fitparams[bad.idx, "m"] <- NA
  fitparams[bad.idx, "a"] <- NA
  fitparams[bad.idx, "c"] <- NA
  fitparams[bad.idx, "d"] <- NA

  if (verbose) {
    nbad <- sum(bad.idx, na.rm = TRUE)
    if (nbad > 0) {
      message("Removed ", nbad, " fits")
    }
  }

  fitparams
}
