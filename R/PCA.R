#############################################################################
## Different forms of PCA for visualizing the behaviour of AA levels
## over time.

## Main function
## We should think about scaling for the different types...
doPCA <- function(datf, type = c("par", "pro", "aa", "both"),
                  rcenter = TRUE, rscale = FALSE) {
  type  <- match.arg(type)

  Xmat <- makeMat(datf, type, rcenter = rcenter, rscale = rscale)
  
  result <- PCA(scale(Xmat, scale = FALSE))
  result$type <- attr(Xmat, "type")
  ## these attributes must be kept for the plotting functions,
  ## otherwise plots will not be consistent with the raw data plots
  result$aaLevels <- attr(datf, "aanames")
  result$intervLevels <- levels(datf$Intervention)
  result$partLevels <- levels(datf$Participant)
    
  result
}

PCA <- function(X, warn = TRUE)
{
  ndf <- nrow(X) - 1
  X.svd <- svd(X)
  varnames <- colnames(X)
  if (is.null(varnames)) varnames <- paste("Var", 1:ncol(X))

  object <- list(scores = X.svd$u %*% diag(X.svd$d),
                 loadings = X.svd$v,
                 var = X.svd$d^2 / ndf,
                 totalvar = sum(X.svd$d^2)/ndf)
  dimnames(object$scores) <- list(rownames(X),
                                  paste("PC", 1:ncol(object$scores)))
  dimnames(object$loadings) <- list(varnames,
                                    paste("PC", 1:ncol(object$loadings)))
  names(object$var) <- paste("PC", 1:length(X.svd$d))

  class(object) <- "PCA"

  object
}

makeMat <- function(datf, type, rcenter, rscale) {
  relVars <- attr(datf, "aanames")
  if (length(relVars) < 3)
    stop("Too few amino acids found for PCA")
  
  idvars <- c("Participant", "Intervention", "Time")
  
  datf <- datf[, c(relVars, idvars)]
  datf.df <- melt(datf, id.vars = idvars)
  datf.array <- acast(datf.df, Time ~ variable ~ Intervention ~ Participant)

  if (rcenter | rscale) {
    timepoints <- dimnames(datf.array)[[1]]
    datf.array <- apply(datf.array, 2:4, scale, center = rcenter,
                        scale = rscale)
    dimnames(datf.array)[[1]] <- timepoints
  }
  
  ## We use the ":" sign to separate dimensions
  result <-
    switch(type,
           both = wrap.array(datf.array, list(2:4, 1)),
           aa = wrap.array(datf.array, list(c(2, 4), c(1, 3))),
           pro = wrap.array(datf.array, list(3:4, 1:2)),
           par = wrap.array(datf.array, list(4, 1:3)))
  attr(result, "type") <- type

  result
}

plotLoadings <- function(PCAobj, ...) {
  Label <- NULL # to avoid R CMD check NOTEs

  if (PCAobj$type == "both") {
    loadings.df <- data.frame(PC1 = PCAobj$loadings[,1],
                              PC2 = PCAobj$loadings[,2])
    x <- as.numeric(rownames(loadings.df))
    xyplot(PC1 + PC2 ~ x, data = loadings.df, type = "l",
           xlab = "Time (min.)", ylab = "Loadings",
           auto.key = list(space = "top", columns = 2,
                           points = FALSE, lines = TRUE),
           panel = function(...) {
             panel.abline(h = 0, v = 0, lty = 2, col = "gray")
             panel.xyplot(...)
           }, ...)
  } else {
    lnames <- strsplit(rownames(PCAobj$loadings), "\\.")
    

    if (PCAobj$type != "par") {
      loadings.df <-
        data.frame(PC1 = PCAobj$loadings[,1],
                   PC2 = PCAobj$loadings[,2],
                   Label = factor(sapply(lnames, "[[", 2)),
                   x = as.numeric(sapply(lnames, "[[", 1)))
      
      xyplot(PC1 + PC2 ~ x, data = loadings.df, type = "b",
             xlab = "Time (min.)", ylab = "Loadings",
             groups = Label, layout = c(1,2), as.table = TRUE,
             scales = list(y = "free"),
             auto.key = list(space = "top",
                             columns = min(6, nlevels(loadings.df$Label))),
             par.settings=list(superpose.symbol=list(pch=1:5)),
             panel = function(...) {
               panel.abline(h = 0, v = 0, lty = 2, col = "gray")
               panel.xyplot(...)
             }, ...)
    } else {
      loadings.df <-
        data.frame(PC = c(PCAobj$loadings[,1],
                          PCAobj$loadings[,2]),
                   what = rep(c("PC1", "PC2"),
                              each = nrow(PCAobj$loadings)),
                   Intervention = factor(sapply(lnames, "[[", 3),
                                         levels = PCAobj$intervLevels),
                   Label = factor(sapply(lnames, "[[", 2)),
                   x = as.numeric(sapply(lnames, "[[", 1)))
      useOuterStrips(
        combineLimits(
          xyplot(PC ~ x | Intervention + what, data = loadings.df, type = "b",
                 xlab = "Time (min.)", ylab = "Loadings",
                 groups = Label, as.table = TRUE,
                 scales = list(y = "free"),
                 auto.key = list(space = "top",
                                 columns = min(6, nlevels(loadings.df$Label))),
                 par.settings=list(superpose.symbol=list(pch=1:5)),
                 panel = function(...) {
                   panel.abline(h = 0, v = 0, lty = 2, col = "gray")
                   panel.xyplot(...)
                 }, ...)
        ))
    }
  }
}

plotScores <-
  function(PCAobj,
           highlight = c("none", "Intervention", "AA", "Participant"),
           nrowLegend) {
  highlight <- match.arg(highlight)
  score.df <- data.frame(PC1 = PCAobj$scores[,1],
                         PC2 = PCAobj$scores[,2])
  names.df <-
    as.data.frame(
      do.call(rbind, sapply(rownames(PCAobj$scores), strsplit, split = "\\.")))
  names(names.df) <- switch(PCAobj$type,
                               par = "Participant",
                               pro = c("Intervention", "Participant"),
                               aa = c("AA", "Participant"),
                               both = c("AA", "Intervention", "Participant"))
  
  if ("Intervention" %in% names(names.df))
    names.df$Intervention <- factor(names.df$Intervention,
                                    levels = PCAobj$intervLevels)
  if ("AA" %in% names(names.df))
    names.df$AA <- factor(names.df$AA,
                          levels = PCAobj$aaLevels)
  if ("Participant" %in% names(names.df))
    names.df$Participant <- factor(names.df$Participant,
                                   levels = PCAobj$partLevels)
  
  ## catch some incompatible argument combinations
  if ((highlight == "AA" & PCAobj$type == "pro") |
      (highlight == "Intervention" & PCAobj$type == "aa")) {
    warning("Incompatible highlight argument")
    highlight <- "none"
  }
      
  Highlight <- NULL # to avoid R CMD check NOTEs

  noLegend <- highlight == "none" | PCAobj$type == "par"
  if (noLegend) {
    score.df$Highlight <- 1
  } else {
    if (highlight == "Participant") {
      score.df$Highlight <-
        factor(names.df[,highlight],
               levels = sort(unique(as.numeric(names.df[,highlight]))))
    } else {
      score.df$Highlight <- names.df[,highlight]
    }
  }
  
  toptext <-
    switch(PCAobj$type,
           par = "Points: participants",
           aa = "Points: AAs in participants",
           pro = "Points: interventions for participants",
           both = "Points: AAs from interventions for participants")
  xlab <- paste("PC1 (",
                formatC(100 * PCAobj$var[1] / PCAobj$totalvar,
                        format = "f", digits = 1), "%)", sep = "")
  ylab <- paste("PC2 (",
                formatC(100 * PCAobj$var[2] / PCAobj$totalvar,
                        format = "f", digits = 1), "%)", sep = "")

  if (noLegend) {
    xyplot(PC2 ~ PC1, data = score.df, main = toptext, groups = Highlight,
           xlab = xlab, ylab = ylab,
           panel = function(x, y, ...) {
             panel.abline(h = 0, v = 0, col = "gray", lty = 2)
             if (PCAobj$type == "par") {
               panel.text(x, y, names.df$Participant)
             } else {
               panel.xyplot(x, y, ...)
             }
           })
  } else {
    if (missing(nrowLegend)) {
      ncolLegend <- min(5, nlevels(score.df$Highlight))
    } else {
      ncolLegend <- nlevels(score.df$Highlight) %/% nrowLegend
      if (nlevels(score.df$Highlight) %% nrowLegend > 0)
        ncolLegend <- ncolLegend + 1
    }
    xyplot(PC2 ~ PC1, data = score.df, main = toptext, groups = Highlight,
           xlab = xlab, ylab = ylab,
           auto.key = list(space = "top",
                           columns = ncolLegend),
           par.settings=list(superpose.symbol=list(pch=1:5)),
           panel = function(x, y, ...) {
             panel.abline(h = 0, v = 0, col = "gray", lty = 2)
             if (PCAobj$type == "par") {
               panel.text(x, y, names.df$Participant)
             } else {
               panel.xyplot(x, y, ...)
             }
           })
  }  
}

