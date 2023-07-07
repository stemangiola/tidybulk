
#' Get matrix from tibble
#' 
#' @importFrom purrr map2_dfc
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom glue glue
#' @importFrom dplyr join_by
#'
#' @keywords internal
#' @noRd
#' 
lmer_to_confidence_intervals_random_effects = function(fit){
  
  
  ster = parameters::standard_error(fit, effects = "random")
  ster = map2_dfc(
    ster, names(ster),
    ~ .x |> 
      as_tibble(rownames = "group_id") |>
      setNames( 
        "group_id" |>
          c(glue("{.y}__{colnames(.x)}"))
      )) |> 
    pivot_longer(-group_id, names_to = "parameter", values_to = "CI")    
  
  mod = lme4::ranef(fit, condVar=T) 
  mod = map2_dfc(
    mod, names(mod),
    ~ .x |> 
      as_tibble(rownames = "group_id") |>
      setNames( 
        "group_id" |>
          c(glue("{.y}__{colnames(.x)}"))
      )) |> 
    pivot_longer(-group_id, names_to = "parameter", values_to = "mode")   
  
  mod |>
    left_join(ster, join_by(group_id, parameter)) |> 
    mutate(lower = mode - CI, upper = mode + CI) |> 
    select(-CI) |> 
    pivot_wider(names_from = parameter, values_from = c(lower, mode, upper), names_glue = "{parameter}__{.value}") |> 
    pivot_wider(names_from = group_id, values_from = -group_id, names_glue = "{group_id}__{.value}")  
}


glmerCore = function (geneList, fullFormula, reduced, data, control, offset, 
          modelData, designMatrix, hyp.matrix, ...) 
{
  data[, "count"] <- geneList$y
  disp <- geneList$dispersion
  fit <- try(suppressMessages(suppressWarnings(lme4::glmer(fullFormula, 
                                                           data = data, control = control, offset = offset, family = MASS::negative.binomial(theta = 1/disp), 
                                                           ...))), silent = TRUE)
  
  
  # If errors return
  if (inherits(fit, "try-error"))
    return(list(stats = NA, coef = NA, stdErr = NA, chisq = NA, 
                df = NA, predict = NA, optinfo = NA, tryErrors = fit[1]))
  
  # Otherwise keep going
  if (length(attr(fit@pp$X, "msgRankdrop")) > 0) {
    return(list(stats = NA, predict = NA, optinfo = NA, 
                tryErrors = attr(fit@pp$X, "msgRankdrop")))
  }
  stdErr <- suppressWarnings(coef(summary(fit))[, 2])
  singular <- as.numeric(lme4::isSingular(fit))
  conv <- length(slot(fit, "optinfo")$conv$lme4$messages)
  vcov. <- suppressWarnings(as.matrix(vcov(fit, complete = FALSE)))
  fixedEffects <- fixef(fit)
  stats <- setNames(c(disp, AIC(fit), as.numeric(logLik(fit))), 
                    c("Dispersion", "AIC", "logLik"))
  if (is.null(reduced)) {
    test <- glmmSeq:::lmer_wald(fixedEffects, hyp.matrix, vcov.)
  }
  else {
    fit2 <- try(suppressMessages(suppressWarnings(lme4::glmer(reduced, 
                                                              data = data, control = control, offset = offset, 
                                                              family = MASS::negative.binomial(theta = 1/disp), 
                                                              ...))), silent = TRUE)
    if (!inherits(fit2, "try-error")) {
      lrt <- anova(fit, fit2)
      test <- list(chisq = setNames(lrt$Chisq[2], "LRT"), 
                   df = lrt$Df[2])
    }
    else {
      test <- list(chisq = NA, df = NA)
    }
  }
  newY <- predict(fit, newdata = modelData, re.form = NA)
  a <- designMatrix %*% vcov.
  b <- as.matrix(a %*% t(designMatrix))
  predVar <- diag(b)
  newSE <- sqrt(predVar)
  newLCI <- exp(newY - newSE * 1.96)
  newUCI <- exp(newY + newSE * 1.96)
  predictdf <- c(exp(newY), newLCI, newUCI)
  return(list(
    stats = stats, 
    coef = fixedEffects, 
    stdErr = stdErr, 
    chisq = test$chisq,
    df = test$df,
    predict = predictdf, 
    optinfo = c(singular, conv), 
    fit = fit,
    tryErrors = ""
  ))
  

}

glmmSeq = function (modelFormula, countdata, metadata, id = NULL, dispersion = NA, 
                    sizeFactors = NULL, reduced = NULL, modelData = NULL, designMatrix = NULL, 
                    method = c("lme4", "glmmTMB"), control = NULL, family = nbinom2, 
                    cores = 1, removeSingles = FALSE, zeroCount = 0.125, verbose = TRUE, 
                    returnList = FALSE, progress = FALSE, ...) 
{
  glmmcall <- match.call(expand.dots = TRUE)
  method <- match.arg(method)
  if (is.null(control)) {
    control <- switch(method, lme4 = lme4::glmerControl(optimizer = "bobyqa"), 
                      glmmTMB = glmmTMBControl())
  }
  countdata <- as.matrix(countdata)
  if (length(lme4::findbars(modelFormula)) == 0) {
    stop("No random effects terms specified in formula")
  }
  if (ncol(countdata) != nrow(metadata)) {
    stop("countdata columns different size to metadata rows")
  }
  if (!is.null(sizeFactors) & ncol(countdata) != length(sizeFactors)) {
    stop("Different sizeFactors length")
  }
  if (!is.numeric(zeroCount)) 
    stop("zeroCount must be numeric")
  if (zeroCount < 0) 
    stop("zeroCount must be >= 0")
  if (zeroCount > 0) 
    countdata[countdata == 0] <- zeroCount
  fullFormula <- update.formula(modelFormula, count ~ ., simplify = FALSE)
  subFormula <- lme4::subbars(modelFormula)
  variables <- rownames(attr(terms(subFormula), "factors"))
  subsetMetadata <- metadata[, variables]
  if (is.null(id)) {
    fb <- lme4::findbars(modelFormula)
    id <- sub(".*[|]", "", fb)
    id <- gsub(" ", "", id)
  }
  ids <- as.character(metadata[, id])
  if (removeSingles) {
    nonSingles <- names(table(ids))[table(ids) > 1]
    nonSingleIDs <- ids %in% nonSingles
    countdata <- countdata[, nonSingleIDs]
    sizeFactors <- sizeFactors[nonSingleIDs]
    subsetMetadata <- subsetMetadata[nonSingleIDs, ]
    ids <- ids[nonSingleIDs]
  }
  if (!is.null(sizeFactors)) 
    offset <- log(sizeFactors)
  else offset <- NULL
  if (verbose) 
    cat(paste0("\nn = ", length(ids), " samples, ", length(unique(ids)), 
               " individuals\n"))
  FEformula <- nobars(modelFormula)
  if (is.null(modelData)) {
    reducedVars <- rownames(attr(terms(FEformula), "factors"))
    varLevels <- lapply(reducedVars, function(x) {
      if (is.factor(metadata[, x])) {
        return(levels(subsetMetadata[, x]))
      }
      else {
        sort(unique(subsetMetadata[, x]))
      }
    })
    modelData <- expand.grid(varLevels)
    colnames(modelData) <- reducedVars
    if (method == "glmmTMB") {
      modelData <- cbind(modelData, .id = NA)
      colnames(modelData)[which(colnames(modelData) == 
                                  ".id")] <- id
    }
  }
  if (is.null(designMatrix)) {
    designMatrix <- model.matrix(FEformula, modelData)
  }
  if (is.null(reduced)) {
    test.stat <- "Wald"
    hyp.matrix <- glmmSeq:::hyp_matrix(fullFormula, metadata, "count")
  }
  else {
    if (length(lme4::findbars(reduced)) == 0) {
      stop("No random effects terms specified in reduced formula")
    }
    subReduced <- lme4::subbars(reduced)
    redvars <- rownames(attr(terms(subReduced), "factors"))
    if (any(!redvars %in% variables)) {
      stop("Extra terms in reduced formula not found full formula")
    }
    reduced <- update.formula(reduced, count ~ ., simplify = FALSE)
    test.stat <- "LRT"
    hyp.matrix <- NULL
  }
  start <- Sys.time()
  if (method == "lme4") {
    if (!all(rownames(countdata) %in% names(dispersion))) {
      stop("Some dispersion values are missing")
    }
    fullList <- lapply(rownames(countdata), function(i) {
      list(y = countdata[i, ], dispersion = dispersion[i])
    })
    if (Sys.info()["sysname"] == "Windows" & cores > 1) {
      cl <- paralle::makeCluster(cores)
      on.exit(stopCluster(cl))
      dots <- list(...)
      varlist <- c("glmerCore", "fullList", "fullFormula", 
                   "reduced", "subsetMetadata", "control", "modelData", 
                   "offset", "designMatrix", "hyp.matrix", "dots")
      clusterExport(cl, varlist = varlist, envir = environment())
      if (progress) {
        resultList <- pbapply::pblapply(fullList, function(geneList) {
          args <- c(list(geneList = geneList, fullFormula = fullFormula, 
                         reduced = reduced, data = subsetMetadata, 
                         control = control, offset = offset, modelData = modelData, 
                         designMatrix = designMatrix, hyp.matrix = hyp.matrix), 
                    dots)
          do.call(glmerCore, args)
        }, cl = cl)
        
        names(resultList) <- rownames(countdata)
        noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)
        
        ci_random_effect_df = 
          do.call(
            rbind,
            pbapply::pblapply(resultList[noErr], function(x) lmer_to_confidence_intervals_random_effects(x$fit), cl = cl)
          )
        
      }
      else {
        resultList <- parLapply(cl = cl, fullList, function(geneList) {
          args <- c(list(geneList = geneList, fullFormula = fullFormula, 
                         reduced = reduced, data = subsetMetadata, 
                         control = control, offset = offset, modelData = modelData, 
                         designMatrix = designMatrix, hyp.matrix = hyp.matrix), 
                    dots)
          do.call(glmerCore, args)
        })
        
        names(resultList) <- rownames(countdata)
        noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)
        
        ci_random_effect_df = 
          do.call(
            rbind,
            parLapply(cl = cl, resultList[noErr], function(x) lmer_to_confidence_intervals_random_effects(x$fit))
          )
        
      }
    }
    else {
      if (progress) {
        resultList <- pbmcapply::pbmclapply(fullList, function(geneList) {
          glmerCore(geneList, fullFormula, reduced, 
                    subsetMetadata, control, offset, modelData, 
                    designMatrix, hyp.matrix, ...)
        }, mc.cores = cores)
        if ("value" %in% names(resultList)) resultList <- resultList$value
        
        names(resultList) <- rownames(countdata)
        noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)
        
        ci_random_effect_df = 
          do.call(
            rbind,
            pbmcapply::pbmclapply(resultList[noErr], function(x) lmer_to_confidence_intervals_random_effects(x$fit), mc.cores = cores)
          )
        
      }
      else {
        resultList <- mclapply(fullList, function(geneList) {
          glmerCore(geneList, fullFormula, reduced, 
                    subsetMetadata, control, offset, modelData, 
                    designMatrix, hyp.matrix, ...)
        }, mc.cores = cores)
        
        names(resultList) <- rownames(countdata)
        noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)
        
        ci_random_effect_df = 
          do.call(
            rbind,
            mclapply(resultList[noErr], function(x) lmer_to_confidence_intervals_random_effects(x$fit), mc.cores = cores)
          )
        
      }
    }
  }
  else {
    fullList <- lapply(rownames(countdata), function(i) {
      countdata[i, ]
    })
    if (Sys.info()["sysname"] == "Windows" & cores > 1) {
      cl <- makeCluster(cores)
      on.exit(stopCluster(cl))
      dots <- list(...)
      varlist <- c("glmmTMBcore", "fullList", "fullFormula", 
                   "reduced", "subsetMetadata", "family", "control", 
                   "modelData", "offset", "designMatrix", "hyp.matrix", 
                   "dots")
      parallel::clusterExport(cl, varlist = varlist, envir = environment())
      if (progress) {
        resultList <- pbapply::pblapply(fullList, function(geneList) {
          args <- c(list(geneList = geneList, fullFormula = fullFormula, 
                         reduced = reduced, data = subsetMetadata, 
                         family = family, control = control, offset = offset, 
                         modelData = modelData, designMatrix = designMatrix, 
                         hyp.matrix = hyp.matrix), dots)
          do.call(glmmSeq:::glmmTMBcore, args)
        }, cl = cl)
        
        names(resultList) <- rownames(countdata)
        noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)
        
        ci_random_effect_df = 
          do.call(
            rbind,
            pbapply::pblapply(resultList[noErr], function(x) lmer_to_confidence_intervals_random_effects(x$fit), cl = cl)
          )
        
      }
      else {
        resultList <- parLapply(cl = cl, fullList, function(geneList) {
          args <- c(list(geneList = geneList, fullFormula = fullFormula, 
                         reduced = reduced, data = subsetMetadata, 
                         family = family, control = control, offset = offset, 
                         modelData = modelData, designMatrix = designMatrix, 
                         hyp.matrix = hyp.matrix), dots)
          do.call(glmmSeq:::glmmTMBcore, args)
        })
        
        names(resultList) <- rownames(countdata)
        noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)
        
        ci_random_effect_df = 
          do.call(
            rbind,
            parLapply(cl = cl, resultList[noErr], function(x) lmer_to_confidence_intervals_random_effects(x$fit))
          )
        
      }
    }
    else {
      if (progress) {
        resultList <- pbmcapply::pbmclapply(fullList, function(geneList) {
          glmmSeq:::glmmTMBcore(geneList, fullFormula, reduced, 
                      subsetMetadata, family, control, offset, 
                      modelData, designMatrix, hyp.matrix, ...)
        }, mc.cores = cores)
        if ("value" %in% names(resultList)) resultList <- resultList$value
        
        names(resultList) <- rownames(countdata)
        noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)
        
        ci_random_effect_df = 
          do.call(
            rbind,
            pbmcapply::pbmclapply(resultList[noErr], function(x) lmer_to_confidence_intervals_random_effects(x$fit), mc.cores = cores)
          )
        
      }
      else {
        resultList <- mclapply(fullList, function(geneList) {
          glmmSeq:::glmmTMBcore(geneList, fullFormula, reduced, 
                      subsetMetadata, family, control, offset, 
                      modelData, designMatrix, hyp.matrix, ...)
        }, mc.cores = cores)
        
        names(resultList) <- rownames(countdata)
        noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)
        
        ci_random_effect_df = 
          do.call(
            rbind,
            mclapply(resultList[noErr], function(x) lmer_to_confidence_intervals_random_effects(x$fit), mc.cores = cores)
          )
        
      }
    }
  }
  if (returnList) 
    return(resultList)

  
  if (sum(noErr) == 0) {
    message("All genes returned an error. Check call. Check sufficient data in each group")
    outputErrors <- vapply(resultList[!noErr], function(x) {
      x$tryErrors
    }, FUN.VALUE = c("test"))
    print(outputErrors[1])
    return(outputErrors)
  }
  end <- Sys.time()
  if (verbose) 
    print(end - start)
  predList <- lapply(resultList[noErr], "[[", "predict")
  outputPredict <- do.call(rbind, predList)
  outLabels <- apply(modelData, 1, function(x) paste(x, collapse = "_"))
  colnames(outputPredict) <- c(paste0("y_", outLabels), paste0("LCI_", 
                                                               outLabels), paste0("UCI_", outLabels))
  optInfo <- t(vapply(resultList[noErr], function(x) {
    setNames(x$optinfo, c("Singular", "Conv"))
  }, FUN.VALUE = c(1, 1)))
  s <- glmmSeq:::organiseStats(resultList[noErr], "Wald")
  meanExp <- rowMeans(log2(countdata[noErr, , drop = FALSE] + 
                             1))
  
  # Add coeff for random effects
  s$coef = s$coef |> cbind(ci_random_effect_df )
  
  s$res <- cbind(s$res, meanExp)
  if (method == "lme4") {
    if (sum(!noErr) != 0) {
      if (verbose) 
        cat(paste("Errors in", sum(!noErr), "gene(s):", 
                  paste(names(noErr)[!noErr], collapse = ", ")))
      outputErrors <- vapply(resultList[!noErr], function(x) {
        x$tryErrors
      }, FUN.VALUE = c("test"))
    }
    else {
      outputErrors <- c("No errors")
    }
  }
  else {
    err <- is.na(s$res[, "AIC"])
    if (any(err)) {
      if (verbose) 
        cat(paste("Errors in", sum(err), "gene(s):", 
                  paste(names(err)[err], collapse = ", ")))
      outputErrors <- vapply(resultList[err], function(x) {
        x$message
      }, FUN.VALUE = character(1))
    }
    else outputErrors <- c("No errors")
  }
  new("GlmmSeq", info = list(call = glmmcall, offset = offset, 
                             designMatrix = designMatrix, method = method, control = control, 
                             family = substitute(family), test.stat = test.stat, 
                             dispersion = dispersion), formula = fullFormula, stats = s, 
      predict = outputPredict, reduced = reduced, countdata = countdata, 
      metadata = subsetMetadata, modelData = modelData, optInfo = optInfo, 
      errors = outputErrors, vars = list(id = id, removeSingles = removeSingles))
}
