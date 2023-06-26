## Various possibilities for visualizing the raw data. Variable datf
## should be the raw data object.
showRawData <- function(datf,
                        type = c("aa", "participant", "sequential"),
                        relevantAAs, 
                        what = c("all", "aas", "essentials", "totals"),
                        scale = list(y = "free"),
                        ...) {
  type <- match.arg(type)
  what <- match.arg(what)

  Intervention <- PP <- NULL # to avoid R CMD check NOTEs

  if (missing(relevantAAs) )
    relevantAAs <- switch(what,
                          aas = attr(datf, "aanames"),
                          totals = attr(datf, "totalnames"),
                          essentials = intersect(attr(datf, "aanames"),
                                                 aaessentials()),
                          c(attr(datf, "aanames"), attr(datf, "totalnames")))
  
  idvars <- c("Participant", "Period", "Time", "Intervention")
  datf.df <- melt(datf[, c(idvars, relevantAAs)], variable.name = "AA",
                  id.vars = idvars)
  datf.df$AA <- factor(datf.df$AA)
  datf.df$Period <- factor(datf.df$Period)

  datf.df2 <- datf.df[datf.df$Time == min(datf.df$Time),]
  datf.df2$Time <- min(datf.df$Time) - 1
  datf.df2$value <- NA
  datf.df3 <- rbind(datf.df2, datf.df)
  datf.df3 <- datf.df3[order(datf.df3$Participant,
                             datf.df3$Period,
                             datf.df3$AA,
                             datf.df3$Time),]
  
  maxTime <- max(datf$Time)
  nProt <- nlevels(datf$Intervention)

  switch(
    type,
    aa = {
      xyplot(value ~ Time | AA, data = datf.df3, groups = Intervention,
             auto.key = list(points = FALSE, lines = TRUE, columns = nProt),
             xlab = "Time (mins)", ylab = "AA level", ...,
             as.table = TRUE, type = "l", scale = scale)
    },
    participant =
      useOuterStrips(
        combineLimits(
          xyplot(value ~ Time | Participant + AA, data = datf.df3,
                 groups = Intervention, as.table = TRUE, ...,
                 xlab = "Time (mins)", ylab = "AA level",
                 auto.key = list(points = FALSE, lines = TRUE,
                                 columns = nProt),
                 type = "l", scale = scale)
        )),
    sequential =
      useOuterStrips(
        combineLimits(
          xyplot(value ~ Time | Participant + AA, groups = Intervention,
                 auto.key = list(points = FALSE, lines = TRUE, columns = nProt),
                 data = datf.df3, type = "l", as.table = TRUE,
                 xlab = "Time (mins)", ylab = "AA level",
                 scale = scale, ...,
                 prepanel = function(...) {
                   list(xlim = c(0, nProt * maxTime))
                 },
                 panel = function(x, y, ..., groups, subscripts) {
                   ##                   x[is.na(y)] <- NA
                   x2 <- x +
                     (as.integer(datf.df3$Period[subscripts]) - 1) *
                     maxTime
                   ## browser()
                   panel.xyplot(x2, y, groups = groups,
                                subscripts = subscripts, ...)
                 })
  ))
  )
}


#############################################################################

## Need to insert a check that the relevant what stuff is also present
## in pardf, probably more important than comparing to datf
## Add check if relvars is empty
showDataFits <- function(datf, pardf, relevantAAs,
                         what = c("all", "aas", "essentials", "totals"),
                         xlab = "Time (min.)", ylab = "AA level",
                         scale = list(y = "free"), baseLineCorr = FALSE,
                         points = TRUE, noLegend = FALSE, ...) {
  what <- match.arg(what)

  Intervention <- NULL # to avoid R CMD check NOTEs

  if (missing(relevantAAs))
    relevantAAs <- switch(what,
                          aas = attr(datf, "aanames"),
                          totals = attr(datf, "totalnames"),
                          essentials = intersect(attr(datf, "aanames"),
                                                 aaessentials()),
                          levels(pardf$AA))
  relevantAAs <- intersect(relevantAAs, levels(pardf$AA))
  idvars <- c("Participant", "Period", "Intervention", "Time")

  ## Allow for subsets of participants - only include participant
  ## levels from pardf in datf
  idx <- datf$Participant %in% levels(pardf$Participant)  
  datf.df <- melt(datf[idx, c(idvars, relevantAAs)], variable.name = "AA",
                  id.vars = idvars)
  datf.df$AA <- factor(datf.df$AA)
#  datf.df$Participant <- factor(datf.df$Participant)

  mytype <- ifelse(points, "p", "n")
  parameters = c("a", "m", "c", "d")

  if (noLegend) {
    auto.k  <- FALSE
  } else {
    auto.k  <- list(column = nlevels(datf.df$Intervention),
                    points = points, lines = !points)
  }

  if (baseLineCorr) {
    for (ii in 1:nrow(pardf)) {
      idx <- which(datf.df$Participant == pardf$Participant[ii] &
                   datf.df$Period == pardf$Period[ii] &
                   datf.df$AA == pardf$AA[ii])
      datf.df$value[idx] <- datf.df$value[idx] - pardf$d[ii]
    }
    pardf$d <- 0.0
  }

  ## added selection on period in params
  mypanel <- function(x, y, ..., groups, subscripts) {
    xx <- seq(min(x), max(x), length = 50)/15
    panel.xyplot(x, y, groups = groups, type = mytype,
                 subscripts = subscripts, ...)
    aa <- as.character(datf.df[subscripts[1], "AA"])
    participant <- datf.df[subscripts[1], "Participant"]
    periods <- unique(datf.df[subscripts, "Period"])
    mycols <- trellis.par.get("superpose.line")$col
    
    params <-
      pardf[pardf$Participant == participant &
            pardf$AA == aa &
            pardf$Period %in% periods,]
    for (pp in 1:nrow(params)) {
      if (!any(is.na(params[pp, parameters])))
        panel.lines(xx*15, woodFun(params[pp, parameters], xx),
                    col = mycols[as.integer(params[pp, "Intervention"])])
    }
  }

  useOuterStrips(
    combineLimits(
      xyplot(value ~ Time | Participant + AA, groups = Intervention,
             xlab = xlab, ylab = ylab, scale = scale, ...,
             data = datf.df, as.table = TRUE,
             par.settings =
               list(superpose.symbol =
                      list(pch=1:nlevels(datf.df$Intervention))),
             auto.key = auto.k, panel = mypanel)
    )
  )
}

showParameters <- function(params,
                           relevantAAs,
                           what = c("all", "aas", "essentials", "totals"),
                           ...) {
  what <- match.arg(what)

  params <- params[names(params) != "RMS"]

  if (missing(relevantAAs))
    relevantAAs <- switch(what,
                          aas = attr(params, "aanames"),
                          totals = attr(params, "totalnames"),
                          essentials = intersect(attr(params, "aanames"),
                                                 aaessentials()),
                          levels(params$AA))
  relevantAAs <- intersect(relevantAAs, levels(params$AA))

  Intervention <- NULL # to avoid R CMD check NOTEs

  params.df <- melt(params[params$AA %in% relevantAAs,],
                    id.vars = c("Participant", "AA", "Intervention", "Period"))
  params.df$AA <- factor(params.df$AA)
  
  useOuterStrips(
    combineLimits(
      xyplot(Participant ~ value | variable + AA, groups = Intervention,
             par.settings =
               list(superpose.symbol =
                      list(pch=1:nlevels(params.df$Intervention))),
             between = list(x = .2),
             scales = list(x = list(relation = "free"),
                           y = list(relation="free", draw = FALSE)),
             data = params.df, scale = "free", as.table = TRUE, ...,
             auto.key = list(place = "top",
                             column = nlevels(params$Intervention)))
    )
  )
}

prepanel.ci <- function(x, y, lx, ux, subscripts, standard = NULL) {
  lx <- as.numeric(lx[subscripts])
  ux <- as.numeric(ux[subscripts])

  x <- as.numeric(x)
  if (is.null(standard)) standard <- NA
  
  list(xlim = range(standard, x, ux, lx, na.rm = TRUE))
}

panel.ci <- function(x, y, lx, ux, subscripts, pch = 16, 
                     standard = 0, truth = NULL, ...) {
  if (length(subscripts) > 0) {
    lx <- lx[subscripts]
    ux <- ux[subscripts]
    if (length(standard) > 1)
      standard <- unique(standard[subscripts])
    
    mycols <- 2 - (standard > lx & standard < ux)
    panel.dotplot(x, y, pch = pch, col = mycols, ...)
    if (!is.null(standard))
      panel.abline(v = standard, col = "grey", lty = 2)
    if (!is.null(truth))
      panel.abline(v = unique(truth[subscripts]), col = "blue", lty = 3)
    panel.arrows(lx, y, ux, y, col = mycols,
                 length = .02, unit = "npc",
                 angle = 90, code = 3, lwd = 2)
  }
}

panel.ci2 <- function(x, y, lx, ux, subscripts, standard = 0, ...) {
  if (length(subscripts) > 0) {
    lx <- lx[subscripts]
    ux <- ux[subscripts]
    if (length(standard) > 1)
      standard <- unique(standard[subscripts])

    args <- list(...)
    panel.abline(h = 1:length(x), col = "gray")
    if (!is.null(standard))
      panel.abline(v = standard, col = "grey", lty = 2)
    panel.arrows(lx, y, ux, y, lwd = 2,
                 length = .02, unit = "npc",
                 angle = 90, code = 3, col = args$col.line)
    panel.xyplot(x, y, ...)
  }
}

showPoIs <- function(AUC = NULL, Height = NULL,
                     Time2Max = NULL, ...) {
  Type <- NULL # to avoid errors during R CMD check
  pois <- list(AUC = AUC, Height = Height, Time2Max = Time2Max)
  pois <- pois[!sapply(pois, is.null)]
  for (ii in seq(along = pois)) {
    pois[[ii]]$Type <- names(pois)[ii]
    colnames(pois[[ii]])[colnames(pois[[ii]]) == names(pois)[ii]] <- "value"
  }
  allpois <- do.call(rbind, pois)
  allpois$Intervention <- factor(allpois$Intervention)

  allcols <- c("orange", "purple", "darkgreen")
  names(allcols) <- c("AUC", "Height", "Time2Max")
  allcols <- allcols[names(pois)]


  if (length(pois) > 1) {
  pl <- xyplot(Intervention ~ value | Type + AA, data = allpois,
               panel = panel.superpose, pch = 16,
               groups = Type, prepanel = prepanel.ci, panel.groups = panel.ci2, 
               par.settings = list(superpose.symbol = list(col = allcols),
                                   superpose.line = list(col = allcols)),
               lx = allpois$low, ux = allpois$up,
               scales = list(x = "free"),
               ...)
    combineLimits(useOuterStrips(pl))
  } else {
    xyplot(Intervention ~ value | AA, data = allpois,
           panel = panel.superpose, pch = 16, main = names(pois)[1],
           groups = Type, prepanel = prepanel.ci, panel.groups = panel.ci2, 
           par.settings = list(superpose.symbol = list(col = allcols),
                               superpose.line = list(col = allcols)),
           lx = allpois$low, ux = allpois$up,
           scales = list(x = "free"),
           ...)
  }
}

showCIs <- function(resultsTable, ...) {
  standard <- ifelse(regexpr("/", resultsTable$contrast[1])[1] < 0, 0, 1)
  xlab <- ifelse(standard == 1, "Intervention / Ref", "Intervention - Ref")
  xyplot(AA ~ Estimate | contrast, data = resultsTable,
         main = resultsTable$Target[1],
         ylab = "Amino acid", xlab = xlab, standard = standard,
         panel = panel.ci, prepanel = prepanel.ci, as.table = TRUE,
         ux = resultsTable$upper.CL, lx = resultsTable$lower.CL, ...)
}

showCombinedCIs <- function(lAUC = NULL, AUC = NULL,
                            Height = NULL, Time2Max = NULL,
                            relevantAAs,
                            what = c("all", "aas", "essentials", "totals"),
                            between, subset, ...) {
  what <- match.arg(what)

  pois <- list(lAUC = lAUC, AUC = AUC, Height = Height, Time2Max = Time2Max)
  pois <- pois[!sapply(pois, is.null)]
  pois.df <- do.call(rbind, pois)

  if (missing(subset)) {
    if (missing(relevantAAs) )
      relevantAAs <- switch(what,
                            aas = aminoacids(),
                            totals = aatotals(),
                            essentials = aaessentials(),
                            levels(pois.df$AA))
    relevantAAs <- intersect(relevantAAs, levels(pois.df$AA))
    if (length(relevantAAs) == 0)
      stop("No overlap between relevantAAs/what arguments and PoI tables")

    subset <- pois.df$AA %in% relevantAAs
  }
  
  pois.df$standard <- 0
  pois.df$standard[sapply(pois.df$contrast,
                          function(x) regexpr("/", x)[1] > 0)] <- 1
  pois.df$Protein <-
    factor(sapply(strsplit(as.character(pois.df$contrast), "[/-]"), "[[", 1))

  if (missing(between)) {
    between <- rep(0, length(pois))
    if (names(pois)[1] == "lAUC") between[1] <- 0.5
    between <- list(x = between)
  }
  
  useOuterStrips(
    combineLimits(
      xyplot(AA ~ Estimate | Target + Protein, data = pois.df,
             scale = list(x = "free"), as.table = TRUE,
             xlab = "Intervention, reference comparison", ylab = "Amino acid",
             panel = panel.ci, , prepanel = prepanel.ci,
             standard = pois.df$standard, between = between,
             subset = subset, ...,
             lx = pois.df$lower.CL, ux = pois.df$upper.CL)
    )
  )
}

