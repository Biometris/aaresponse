aminoacids <- function() {
  c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Gln", "Gln/Arg",
    "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro",
    "Ser", "Thr", "Trp", "Tryp", "Tyr", "Val")
}

aatotals <- function() {
  c("TAA", "TEAA", "EAA")
}

aaessentials <- function() {
    c("His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Tryp", "Val")
}

checkAAdata <- function(df, aanames, totalnames, quiet = FALSE) {  
  cnames <- colnames(df)
  rnames <- c("Participant", "Period", "Intervention", "Time")
  missingVars <- setdiff(rnames, cnames)
  ## are all required names present?
  if (length(missingVars) > 0)
    stop("Required variables absent:", paste(missingVars, collapse = ", "))

  if (any(is.na(df[,rnames])))
    stop("Missing values present in required variables.")
  
  if (missing(aanames)) aanames <- aminoacids()
  if (missing(totalnames)) totalnames <- aatotals()
  ## sometimes people specify totals within aanames
  aanames <- setdiff(aanames, totalnames)
  allnames <- c(aanames, totalnames)

  foundAAs <- intersect(aanames, cnames)
  foundTotals <- intersect(totalnames, cnames)
  
  allAAs <- c(foundAAs, foundTotals)
  NAA <- length(allAAs)
  if (NAA == 0)
    stop("No amino-acid information found.")

  if (max(df$Time) < 30)
    stop("Time should be recorded in minutes")

  maxNAA <- 10
  nlines <- ifelse(NAA %% maxNAA == 0, NAA %/% maxNAA, NAA %/% maxNAA + 1)
  AAstring <-
    lapply(1:nlines,
           function(ii)
             paste(allAAs[((ii-1)*10 + 1):min(NAA, ii*10)], collapse = ", "))

  if (!quiet) {
    cat("Experimental setup (interventions and periods):")
    xod <- summary.XOdesign(df, verbose = FALSE)
    print(xod)
    if (max(xod) - min(xod) > 1)
      cat("Warning: interventions not evenly distributed.\n")
  }

  for (var in c("Participant", "Period", "Intervention")) 
    if (!is.factor(df[[var]]))
      df[[var]] <- factor(df[[var]])

  if (!quiet) {
    cat("\nNumber of participants:", nlevels(df$Participant))
    cat("\nInterventions:", paste(levels(df$Intervention), sep = ", "))
    cat("\nTime range: from", min(df$Time), "to", max(df$Time), "minutes.\n")
    cat("\nInformation on", NAA, "amino acids and aa totals:\n\t",
        do.call(paste, args = list(AAstring, collapse = "\n\t ")), "\n\n")
  }
  
  attr(df, "aanames") <- foundAAs
  attr(df, "totalnames") <- foundTotals

  class(df) <- c("aar", "data.frame")
  df
}

## use getAnywhere("[.data.frame") to see source code...
'[.aar' <- function(x, ...) {
  r <- NextMethod("[")
  attr(r, "aanames") <- intersect(attr(x, "aanames"), colnames(r))
  attr(r, "totalnames") <- intersect(attr(x, "totalnames"), colnames(r))
  r
}

## '[[.aar' <- function(x, ...) {
##   r <- NextMethod("[[")
##   attr(r, "aanames") <- intersect(attr(x, "aanames"), colnames(r))
##   attr(r, "totalnames") <- intersect(attr(x, "totalnames"), colnames(r))
##   r
## }

## next two functions are not necessary - attributes are kept upon assignment
## '[<-.aar' <- function(x, i, j, value) {
##   r <- NextMethod("[<-", i=i, j=j, value = value)
##   mostattributes(r) <- attributes(x)
##   r
## }

## '[[<-.aar' <- function(x, ...) {
##   r <- NextMethod("[[<-", ...)
##   mostattributes(r) <- attributes(x)
##   r
## }

## subset.aar <- function(x, ...) {
##   r <- NextMethod("subset")
##   attr(r, "aanames") <- intersect(attr(x, "aanames"), colnames(r))
##   attr(r, "totalnames") <- intersect(attr(x, "totalnames"), colnames(r))
##   r
## }
