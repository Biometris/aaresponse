getAUC <- function(prs, maxT = 300) {
  required <- c("a", "m", "c", "d")
  if (is.null(names(prs)) & length(prs) == 4)
    names(prs) <- required
  if (!all(required %in% names(prs)))
    stop("Required information absent")
  prs <- prs[required]
  if (is.list(prs)) prs <- unlist(prs)
  if (any(is.na(prs))) return(NA)
      
  unname((prs[1] / prs[3]^(prs[2] * prs[3] + 1)) *
         gammainc(prs[3] * maxT/15, prs[2]*prs[3] + 1)[1])
}

getHeight <- function(prs) {
  required <- c("a", "m", "c", "d")
  if (is.null(names(prs)) & length(prs) == 4)
    names(prs) <- required
  if (!all(required %in% names(prs)))
    stop("Required information absent")
  prs <- prs[required]
  if (is.list(prs)) prs <- unlist(prs)
  if (any(is.na(prs))) return(NA)

  unname((prs[1]*prs[2]^(prs[2]*prs[3]) * exp(-(prs[2] * prs[3]))))
}

getTime2Max <- function(prs, minutes = TRUE) {
  required <- c("a", "m", "c", "d")
  if (is.null(names(prs)) & length(prs) == 4)
    names(prs) <- required
  if (!all(required %in% names(prs)))
    stop("Required information absent")
  prs <- prs[required]
  if (is.list(prs)) prs <- unlist(prs)
  if (any(is.na(prs))) return(NA)
  
  ifelse(minutes, prs[2]*15, prs[2])
}

extractParameters <- function(prs.df, maxT = 300, minutes = TRUE) {
  required <- c("a", "m", "c", "d")
  if (!all(required %in% names(prs.df)))
    stop("Required information absent")
  ddff <- prs.df[required]
    
  result <-
    as.data.frame(
      cbind(apply(ddff, 1, getAUC, maxT = maxT),
            apply(ddff, 1, getHeight),
            apply(ddff, 1, getTime2Max, minutes = minutes)))
  names(result) <- c("AUC", "Height", "Time2Max")

  additional <- c("Participant", "Intervention", "AA", "Period")
  if (all(additional %in% names(prs.df))) {
    cbind(prs.df[additional], as.data.frame(result))
  } else {
    result
  }
}

curateParameters <- function(params,
                             Time2Max = c(15, 200),
                             Height = c(0, 1000),
                             AUC = c(0, Inf),
                             verbose = TRUE) {
  if (missing(params)) {
    boundaryTable <- rbind(Time2Max, Height, AUC)
    dimnames(boundaryTable) <-
      list(c("Time2Max", "Height", "AUC"),
           c("Min", "Max"))
    return(boundaryTable)
  }

  if (any(levels(params$AA) %in% aminoacids()) &
      any(levels(params$AA) %in% aatotals()))
    warning("Curation will probably fail when applied to individual AAs and totals simultaneously")

  if (length(setdiff(levels(params$AA), c(aminoacids(), aatotals()))) > 0)
    warning("Maybe relevant: curation will probably fail when applied to individual AAs and totals simulateously")
  
  bad.idx <-
    is.na(params$Time2Max) | is.na(params$Height) | is.na(params$AUC) |
    params$Time2Max < Time2Max[1] | params$Time2Max > Time2Max[2] |
    params$Height < Height[1] | params$Height > Height[2] |
    params$AUC < AUC[1] | params$AUC > AUC[2]
  
  params$Height[bad.idx] <- NA
  params$AUC[bad.idx] <- NA
  params$Time2Max[bad.idx] <- NA

  if (verbose) {
    nbad <- sum(bad.idx, na.rm = TRUE)
    if (nbad > 0) {
      message("Removed ", nbad, " parameter combinations")
    }
  }

  params
}

summarizePoIs <- function(datf,
                          target = c("AUC", "Height", "Time2Max"),
                          coverage = .95) {
  target <- match.arg(target)
  if (coverage < .5 | coverage > 1)
    stop("Bad value for coverage - choose a value between .5 and 1")
  lowcov <- .5*(1-coverage)
  highcov <- 1 - lowcov

  huhn <-
    by(datf,
       datf["AA"],
       function(xxx) {
         mns <- aggregate(xxx[target], xxx["Intervention"],
                          mean, na.rm = TRUE)[,2]
         sds <- aggregate(xxx[target], xxx["Intervention"],
                          sd, na.rm = TRUE)[,2]
         lowci <-
           aggregate(xxx[target],
                     xxx["Intervention"],
                     function(yy) {
                       if ((nna <- sum(!is.na(yy))) > 1) {
                         qt(lowcov, df = nna - 1) *
                           sd(yy, na.rm = TRUE) / sqrt(nna)
                       } else {
                         NA
                       }
                       })[,2]
         upci <-
           aggregate(xxx[target],
                     xxx["Intervention"],
                     function(yy) {
                       if ((nna <- sum(!is.na(yy))) > 1) {
                         qt(highcov, df = nna - 1) *
                           sd(yy, na.rm = TRUE) / sqrt(nna)
                       } else {
                         NA
                       }
                     })[,2]
         
         res <- as.data.frame(rbind(mns, mns + lowci, mns + upci, sds))
         colnames(res) <- levels(xxx$Intervention)
         rownames(res) <- c(target, "low", "up", "sd")

         t(res)
       })
  
  result <- as.data.frame(do.call(rbind, huhn))
  result$Intervention <- factor(levels(datf$Intervention),
                                levels = levels(datf$Intervention))
  result$AA <- factor(rep(levels(datf$AA),
                          each = nlevels(datf$Intervention)),
                      levels = levels(datf$AA))
  rownames(result) <- NULL

  result[, c(6, 5, 1:4)]
}

## impute NA values for AUC and Height - for Time2Max this does not
## make much sense. The idea is to relate the spread in raw data to
## Height and AUC and then to use that straight line to estimate
## values for NA cases. We add additional columns AUC.i and Height.i.

imputePoIs <- function(df, params, quant = .2) {
  getMD  <- function(subdf, qnt = quant) {
    qs <- c(qnt, 1-qnt)
    diff(quantile(subdf, qs, na.rm = TRUE))
  }

  AA <- NULL # to avoid R CMD check NOTEs

  aas <- levels(params$AA)
  maxdiffs <- aggregate(df[aas],
                        df[c("Participant", "Period")],
                        getMD)
  maxdiffs.df <- melt(maxdiffs,
                      id.vars = c("Participant", "Period"),
                      variable.name = "AA",
                      value.name = "quantDiff")

  paramsNew <- merge(params, maxdiffs.df)
  paramsNew$AUC.orig <- paramsNew$AUC
  paramsNew$Height.orig <- paramsNew$Height
  paramsNew$AUC.i <- paramsNew$Height.i <- NA

  AUCmods <-
    lapply(aas,
           function(aa) {
             rlm(AUC ~ quantDiff, data = paramsNew, subset = AA == aa)
           })

  for (ii in seq(along = aas)) {
    paramsNew$AUC.i[paramsNew$AA == aas[ii] ] <-
      predict(AUCmods[[ii]],
              data.frame(quantDiff =
                           paramsNew$quantDiff[paramsNew$AA == aas[ii] ]))
  }
  
  Heightmods <-
    lapply(aas,
           function(aa) {
             rlm(Height ~ quantDiff, data = paramsNew, subset = AA == aa)
           })

  for (ii in seq(along = aas)) {
    paramsNew$Height.i[paramsNew$AA == aas[ii] ] <-
      predict(Heightmods[[ii]],
              data.frame(quantDiff =
                           paramsNew$quantDiff[paramsNew$AA == aas[ii] ]))
  }

  paramsNew$AUC[is.na(paramsNew$AUC.orig)] <-
    paramsNew$AUC.i[is.na(paramsNew$AUC.orig)]
    
  paramsNew$Height[is.na(paramsNew$Height.orig)] <-
    paramsNew$Height.i[is.na(paramsNew$Height.orig)]
  
  paramsNew                     
}
