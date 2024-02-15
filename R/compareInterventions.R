fitMixedModels <-
  function(dt, refIntervention,
           model.formula, lm.alternative,
           mainFun = lmer, respondersOnly = TRUE)
{
  imputedValues <- "AUC.orig" %in% colnames(dt)

  formulaText <- as.character(model.formula)
  target <- strsplit(formulaText, " ")[1]
  
  if (!respondersOnly) {
    if (target == "Time2Max") {
      warning("No imputation possible for Time2Max - presenting results for respondersOnly")
      respondersOnly <- TRUE
    }
    if (!imputedValues) {
      warning("No imputed values found - presenting results for respondersOnly")
      respondersOnly <- TRUE
    }
  }

  if (respondersOnly) {
    if (imputedValues) {
      dt$AUC <- dt$AUC.orig
      dt$Height <- dt$Height.orig
    }
  }
  
  if (!missing(refIntervention)) {
    dt$Intervention <- relevel(dt$Intervention, ref = refIntervention)
  } else {
    refIntervention <- levels(dt$Intervention)[1]
  }

  doFit <- function(dtf) {
    modelFit <- mainFun(model.formula, data = dtf)
    
    if (isSingular(modelFit)) { # lm model
      lm(lm.alternative, data = dtf)
    } else {
      modelFit
    }
  }

  by(dt, dt["AA"], doFit)
}

testPeriodEffect <- function(fitModels) {
  doPeriodTests <- function(mFit) {
    mFit <- emmeans(mFit, c("Intervention", "Period"))
    joint_tests(mFit)
  }
  
  emTests <- lapply(fitModels, doPeriodTests)
  pvals <- lapply(emTests, "[[", "p.value")
  pvalsPeriod <- sapply(pvals, function(x) if (length(x) == 2) x[2] else NA)
  if (any(naP <- which(is.na(pvalsPeriod))))
    warning(paste("No period p values could be calculated for",
                  names(emTests)[naP], "\n"))
  if (any(smallP <- which(p.adjust(pvalsPeriod, "BH") < .05)))
    warning(paste("Period significant in these cases:",
                  names(emTests)[smallP]))

  pvalsPeriod
}  

doComparisons <- function(fitModels, logTransform = FALSE, ...) {
  doEMmeans <- function(mFit, logT) {
    if (logT) {
      regrid(emmeans(mFit, "Intervention"), transform = "log")
    } else {
      emmeans(mFit, "Intervention")
    }
  }

  emMods <- lapply(fitModels, doEMmeans, logT = logTransform)
  results <- lapply(emMods, contrast, method = "trt.vs.ctrl",
                    type = "response")
  ## sometimes p values are NaNs leading to errors, so we have to
  ## check for that
  myconfint <- function(xxx, ...) {
    result <- try(confint(xxx, ...), silent = TRUE)
    if (inherits(result, "try-error")) 
      result <- NA

    result
  }
  
  allConfInts <- lapply(results, myconfint, ...)
  allConfInts <- allConfInts[!sapply(allConfInts, function(x) all(is.na(x)))]
  result.df <- as.data.frame(do.call(rbind, allConfInts))
  names(result.df)[2] <- "Estimate"
  result.df <-
    cbind(
      data.frame(AA = factor(sapply(strsplit(rownames(result.df), "\\."),
                                    "[[", 1)),
                 Target = as.character(terms(fitModels[[1]]))[2]),
      result.df)
  rownames(result.df) <- NULL

  result.df
}

## Would be nice to be able to add additional arguments for the
## mixed-model functions, especially when using functions from other
## packages than lme4... Can we use the ellipses for that?
compareInterventions <-
  function(dt, refIntervention,
           target = c("AUC", "Height", "Time2Max"),
           logTransform = FALSE,
           model.formula = "~ Period + Intervention + (1 | Participant)",
           lm.alternative = "~ Period + Intervention",
           mainFun = lmer, singularFun = lm,
           respondersOnly = TRUE, ...)
{
  model.terms <- trimws(strsplit(model.formula, "\\+")[[1]])

  model.formula <- formula(paste(target, model.formula, collapse = " "))
  lm.alternative <- formula(paste(target, lm.alternative, collapse = " "))

  fitModels <-
    fitMixedModels(dt = dt,
                   refIntervention = refIntervention,
                   model.formula = model.formula,
                   lm.alternative = lm.alternative,
                   mainFun = mainFun,
                   respondersOnly = respondersOnly)

  if ("Period" %in% model.terms)
    periodPvals <- testPeriodEffect(fitModels)
  
  doComparisons(fitModels, logTransform, ...)
}
