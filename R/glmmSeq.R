
# From https://github.com/easystats/parameters
lme4_standard_error = function  (model)
{
  rand.se <- lme4::ranef(model, condVar = TRUE)
  n.groupings <- length(rand.se)
  for (m in 1:n.groupings) {
    vars.m <- attr(rand.se[[m]], "postVar")
    K <- dim(vars.m)[1]
    J <- dim(vars.m)[3]
    names.full <- dimnames(rand.se[[m]])
    rand.se[[m]] <- array(NA, c(J, K))
    for (j in 1:J) {
      rand.se[[m]][j, ] <- sqrt(diag(as.matrix(vars.m[,
                                                      , j])))
    }
    dimnames(rand.se[[m]]) <- list(names.full[[1]], names.full[[2]])
  }
  rand.se
}

glmmTMB_standard_error = function (model){
  rand.se <- glmmTMB::ranef(model, condVar = TRUE)$cond
  n.groupings <- length(rand.se)
  for (m in 1:n.groupings) {
    vars.m <- attr(rand.se[[m]], "condVar")
    K <- dim(vars.m)[1]
    J <- dim(vars.m)[3]
    names.full <- dimnames(rand.se[[m]])
    rand.se[[m]] <- array(NA, c(J, K))
    for (j in 1:J) {
      rand.se[[m]][j, ] <- sqrt(diag(as.matrix(vars.m[,
                                                      , j])))
    }
    dimnames(rand.se[[m]]) <- list(names.full[[1]], names.full[[2]])
  }
  rand.se
}

#' Get matrix from tibble
#'
#' @importFrom purrr map2_dfc
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr join_by
#'
#' @keywords internal
#' @noRd
#'
glmmTMB_to_confidence_intervals_random_effects = function(fit){


  ster = glmmTMB_standard_error(fit)
  ster = map2_dfr(
    ster, names(ster),
    ~ .x |>
      as_tibble(rownames = "group_id") |>
      setNames(
        "group_id" |>
          c(sprintf("%s__%s", .y, colnames(.x)))
      ) |>
      pivot_longer(-group_id, names_to = "parameter", values_to = "CI")
  )

  mod = glmmTMB::ranef(fit, condVar=T)$cond
  mod = map2_dfr(
    mod, names(mod),
    ~ .x |>
      as_tibble(rownames = "group_id") |>
      setNames(
        "group_id" |>
          c(sprintf("%s__%s", .y, colnames(.x)))
      ) |>
      pivot_longer(-group_id, names_to = "parameter", values_to = "mode")
  )

  mod |>
    left_join(ster, join_by(group_id, parameter)) |>
    mutate(lower = mode - CI, upper = mode + CI) |>
    select(-CI) |>
    tidyr::unite("parameter", c(group_id, parameter), sep="_") |>
    pivot_wider(names_from = parameter, values_from = c(lower, mode, upper), names_glue = "{parameter}__{.value}")
}


#' Get matrix from tibble
#'
#' @importFrom purrr map2_dfc
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr join_by
#'
#' @keywords internal
#' @noRd
#'
lmer_to_confidence_intervals_random_effects = function(fit){


  ster = lme4_standard_error(fit)
  ster = map2_dfr(
    ster, names(ster),
    ~ .x |>
      as_tibble(rownames = "group_id") |>
      setNames(
        "group_id" |>
          c(sprintf("%s__%s", .y, colnames(.x)))
      ) |>
      pivot_longer(-group_id, names_to = "parameter", values_to = "CI")
  )

  mod = lme4::ranef(fit, condVar=T)
  mod = map2_dfr(
    mod, names(mod),
    ~ .x |>
      as_tibble(rownames = "group_id") |>
      setNames(
        "group_id" |>
          c(sprintf("%s__%s", .y, colnames(.x)))
      ) |>
      pivot_longer(-group_id, names_to = "parameter", values_to = "mode")
  )

  mod |>
    left_join(ster, join_by(group_id, parameter)) |>
    mutate(lower = mode - CI, upper = mode + CI) |>
    select(-CI) |>
    tidyr::unite("parameter", c(group_id, parameter), sep="_") |>
    pivot_wider(names_from = parameter, values_from = c(lower, mode, upper), names_glue = "{parameter}__{.value}")
}

glmmTMBcore = function (geneList, fullFormula, reduced, data, family, control,
          offset, modelData, designMatrix, hyp.matrix, ...)
{
  data[, "count"] <- geneList
  fit <- try(suppressMessages(suppressWarnings(glmmTMB::glmmTMB(fullFormula,
                                                       data, family, control = control, offset = offset, ...))),
             silent = TRUE)
  if (!inherits(fit, "try-error")) {
    singular <- conv <- NA
    stdErr <- suppressWarnings(coef(summary(fit))$cond[, 2])
    vcov. <- vcov(fit)$cond
    fixedEffects <- glmmTMB::fixef(fit)$cond
    disp <- glmmTMB::sigma(fit)
    msg <- fit$fit$message
    stats <- setNames(c(disp, AIC(fit), as.numeric(logLik(fit))),
                      c("Dispersion", "AIC", "logLik"))
    if (is.null(reduced)) {
      test <- glmmSeq:::lmer_wald(fixedEffects, hyp.matrix, vcov.)
    }
    else {
      fit2 <- try(suppressMessages(suppressWarnings(glmmTMB::glmmTMB(reduced,
                                                            data, family, control = control, offset = offset,
                                                            ...))), silent = TRUE)
      if (!inherits(fit2, "try-error")) {
        lrt <- anova(fit, fit2)
        test <- list(chisq = setNames(lrt$Chisq[2], "LRT"),
                     df = lrt[2, "Chi Df"])
      }
      else {
        test <- list(chisq = NA, df = NA)
      }
    }
    newY <- predict(fit, newdata = modelData, re.form = NA)
    a <- designMatrix %*% vcov.
    b <- as.matrix(a %*% t(designMatrix))
    predVar <- diag(b)
    newSE <- suppressWarnings(sqrt(predVar))
    newLCI <- exp(newY - newSE * 1.96)
    newUCI <- exp(newY + newSE * 1.96)
    predictdf <- c(exp(newY), newLCI, newUCI)

    ci_random_effect_df = glmmTMB_to_confidence_intervals_random_effects(fit)


    return(list(
      stats = stats,
      coef = fixedEffects,
      stdErr = stdErr,
      chisq = test$chisq,
      df = test$df,
      predict = predictdf,
      optinfo = c(singular, conv),
      ci_random_effect_df = ci_random_effect_df,
      message = msg,
      tryErrors = ""))
  }
  else {
    return(list(stats = NA, coef = NA, stdErr = NA, chisq = NA,
                df = NA, predict = NA, optinfo = NA, ci_random_effect_df = NA, tryErrors = fit[1]))
  }
}

glmerCore = function (geneList, fullFormula, reduced, data, control, offset,
          modelData, designMatrix, hyp.matrix, max_rows_for_matrix_multiplication = Inf, return_fit =FALSE, ...)
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
  fixedEffects <- lme4::fixef(fit)
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


  ci_random_effect_df = lmer_to_confidence_intervals_random_effects(fit)

  return_list =
    list(
      stats = stats,
      coef = fixedEffects,
      stdErr = stdErr,
      chisq = test$chisq,
      df = test$df,
      predict = predictdf,
      optinfo = c(singular, conv),
      ci_random_effect_df = ci_random_effect_df,
      tryErrors = ""
    )

  # If fit not needed do not return
  if(!return_fit) {
    rm(fit)
    gc()
  } else {
    return_list$fit = fit

  }

  return(return_list)


}

glmmSeq = function (modelFormula, countdata, metadata, id = NULL, dispersion = NA,
                    sizeFactors = NULL, reduced = NULL, modelData = NULL, designMatrix = NULL,
                    method = c("lme4", "glmmTMB"), control = NULL, family = glmmTMB::nbinom2,
                    cores = 1, removeSingles = FALSE, zeroCount = 0.125, verbose = TRUE,
                    returnList = FALSE, progress = FALSE, max_rows_for_matrix_multiplication = Inf, ...)
{
  glmmcall <- match.call(expand.dots = TRUE)
  method <- match.arg(method)
  if (is.null(control)) {
    control <- switch(method, lme4 = lme4::glmerControl(optimizer = "bobyqa"),
                      glmmTMB = glmmTMB::glmmTMBControl())
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
  subsetMetadata <- metadata[, variables, drop=FALSE]
  if (is.null(id)) {
    fb <- lme4::findbars(modelFormula)
    id <- sub(".*[|]", "", fb)
    id <- gsub(" ", "", id)
  }
  ids <- as.character(metadata[, id])
  if (removeSingles) {
    nonSingles <- names(table(ids))[table(ids) > 1]
    nonSingleIDs <- ids %in% nonSingles
    countdata <- countdata[, nonSingleIDs, drop=FALSE]
    sizeFactors <- sizeFactors[nonSingleIDs]
    subsetMetadata <- subsetMetadata[nonSingleIDs, , drop=FALSE]
    ids <- ids[nonSingleIDs]
  }
  if (!is.null(sizeFactors))
    offset <- log(sizeFactors)
  else offset <- NULL
  if (verbose)
    message(paste0("\nn = ", length(ids), " samples, ", length(unique(ids)),
               " individuals\n"))
  FEformula <- lme4::nobars(modelFormula)
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
      na_matrix = rep(NA, nrow(modelData) * length(id)) |> matrix(ncol = length(id))
      colnames(na_matrix) = id
      modelData <- modelData |> cbind(na_matrix)
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

  #-----------------------------------------------#
  # If matrix is too big because model is too big
  #-----------------------------------------------#
  rows_to_sample = sample(seq_len(nrow(designMatrix)), min(max_rows_for_matrix_multiplication, nrow(designMatrix)))
  if(length(rows_to_sample) < nrow(designMatrix)) warning(sprintf("tidybulk says: for calculating p-value the combination of covariates has been limited to 10000 rather than %d othwerwise a matrix multiplication in the glmmSeq would overflow the momery available for arrays.", nrow(designMatrix)))
  designMatrix = designMatrix[rows_to_sample,,drop=FALSE]
  modelData = modelData[rows_to_sample,,drop=FALSE]
  #----------------------------------------------#


  if (method == "lme4") {

    # FASTER - Stefano
    control$calc.derivs = FALSE

    if (!all(rownames(countdata) %in% names(dispersion))) {
      stop("Some dispersion values are missing")
    }
    fullList <- lapply(rownames(countdata), function(i) {
      list(y = countdata[i, ], dispersion = dispersion[i])
    })
    if(cores == 1){
      resultList <- lapply(fullList, function(geneList) {
        glmerCore(geneList, fullFormula, reduced,
                  subsetMetadata, control, offset, modelData,
                  designMatrix, hyp.matrix,, max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication, ...)
      })
    }
    else if (Sys.info()["sysname"] == "Windows" & cores > 1) {
      cl <- parallel::makeCluster(cores)
      on.exit(stopCluster(cl))
      dots <- list(...)
      varlist <- c("glmerCore", "fullList", "fullFormula",
                   "reduced", "subsetMetadata", "control", "modelData",
                   "offset", "designMatrix", "hyp.matrix", "dots")
      clusterExport(cl, varlist = varlist, envir = environment())
      if (progress) {

        # Check if package is installed, otherwise install
        if (find.package("pblapply", quiet = TRUE) %>% length %>% equals(0)) {
          message("tidybulk says: Installing pblapply needed for differential transcript abundance analyses")
          if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos = "https://cloud.r-project.org")
          BiocManager::install("pblapply", ask = FALSE)
        }

        resultList <- pbapply::pblapply(fullList, function(geneList) {
          args <- c(list(geneList = geneList, fullFormula = fullFormula,
                         reduced = reduced, data = subsetMetadata,
                         control = control, offset = offset, modelData = modelData,
                         designMatrix = designMatrix, hyp.matrix = hyp.matrix, max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication),
                    dots)
          do.call(glmerCore, args)
        }, cl = cl)

        names(resultList) <- rownames(countdata)
        noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)

      }
      else {
        resultList <- parLapply(cl = cl, fullList, function(geneList) {
          args <- c(list(geneList = geneList, fullFormula = fullFormula,
                         reduced = reduced, data = subsetMetadata,
                         control = control, offset = offset, modelData = modelData,
                         designMatrix = designMatrix, hyp.matrix = hyp.matrix, max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication),
                    dots)
          do.call(glmerCore, args)
        })

        names(resultList) <- rownames(countdata)
        noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)

      }
    }
    else {
      if (progress) {

        # Check if package is installed, otherwise install
        if (find.package("pbmcapply", quiet = TRUE) %>% length %>% equals(0)) {
          message("tidybulk says: Installing pbmcapply needed for differential transcript abundance analyses")
          if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos = "https://cloud.r-project.org")
          BiocManager::install("pbmcapply", ask = FALSE)
        }

        resultList <- pbmcapply::pbmclapply(fullList, function(geneList) {
          glmerCore(geneList, fullFormula, reduced,
                    subsetMetadata, control, offset, modelData,
                    designMatrix, hyp.matrix, , max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication, ...)
        }, mc.cores = cores)
        if ("value" %in% names(resultList)) resultList <- resultList$value

      }
      else {
        resultList <- mclapply(fullList, function(geneList) {
          glmerCore(geneList, fullFormula, reduced,
                    subsetMetadata, control, offset, modelData,
                    designMatrix, hyp.matrix,, max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication, ...)
        }, mc.cores = cores)

      }
    }
  }
  else {

    # FASTER - Stefano
    control$profile=TRUE

    fullList <- lapply(rownames(countdata), function(i) {
      countdata[i, ]
    })

    if(cores == 1){
      resultList <- lapply(fullList, function(geneList) {
        glmmTMBcore(geneList, fullFormula, reduced,
                              subsetMetadata, family, control, offset,
                              modelData, designMatrix, hyp.matrix, ...)
      })
    }
    else if (Sys.info()["sysname"] == "Windows" & cores > 1) {
      cl <- makeCluster(cores)
      on.exit(stopCluster(cl))
      dots <- list(...)
      varlist <- c("glmmTMBcore", "fullList", "fullFormula",
                   "reduced", "subsetMetadata", "family", "control",
                   "modelData", "offset", "designMatrix", "hyp.matrix",
                   "dots")
      parallel::clusterExport(cl, varlist = varlist, envir = environment())
      if (progress) {

        # Check if package is installed, otherwise install
        if (find.package("pblapply", quiet = TRUE) %>% length %>% equals(0)) {
          message("tidybulk says: Installing pblapply needed for differential transcript abundance analyses")
          if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos = "https://cloud.r-project.org")
          BiocManager::install("pblapply", ask = FALSE)
        }

        resultList <- pbapply::pblapply(fullList, function(geneList) {
          args <- c(list(geneList = geneList, fullFormula = fullFormula,
                         reduced = reduced, data = subsetMetadata,
                         family = family, control = control, offset = offset,
                         modelData = modelData, designMatrix = designMatrix,
                         hyp.matrix = hyp.matrix), dots)
          do.call(glmmTMBcore, args)
        }, cl = cl)

      }
      else {
        resultList <- parLapply(cl = cl, fullList, function(geneList) {
          args <- c(list(geneList = geneList, fullFormula = fullFormula,
                         reduced = reduced, data = subsetMetadata,
                         family = family, control = control, offset = offset,
                         modelData = modelData, designMatrix = designMatrix,
                         hyp.matrix = hyp.matrix), dots)
          do.call(glmmTMBcore, args)
        })

      }
    }
    else {
      if (progress) {

        # Check if package is installed, otherwise install
        if (find.package("pbmcapply", quiet = TRUE) %>% length %>% equals(0)) {
          message("tidybulk says: Installing pbmcapply needed for differential transcript abundance analyses")
          if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos = "https://cloud.r-project.org")
          BiocManager::install("pbmcapply", ask = FALSE)
        }

        resultList <- pbmcapply::pbmclapply(fullList, function(geneList) {
          glmmTMBcore(geneList, fullFormula, reduced,
                      subsetMetadata, family, control, offset,
                      modelData, designMatrix, hyp.matrix, ...)
        }, mc.cores = cores)
        if ("value" %in% names(resultList)) resultList <- resultList$value

      }
      else {
        resultList <- mclapply(fullList, function(geneList) {
          glmmTMBcore(geneList, fullFormula, reduced,
                      subsetMetadata, family, control, offset,
                      modelData, designMatrix, hyp.matrix, ...)
        }, mc.cores = cores)

      }
    }
  }

  # Premature return
  if (returnList) return(resultList)

  names(resultList) <- rownames(countdata)
  noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)

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
  meanExp <- rowMeans(log2(countdata[noErr, , drop = FALSE] + 1))

  # Add coeff for random effects
  ci_random_effect_df =
    do.call(
      rbind,
     lapply(resultList[noErr], function(x) x$ci_random_effect_df)
    )
  s$coef = s$coef |> cbind(ci_random_effect_df )

  s$res <- cbind(s$res, meanExp)
  if (method == "lme4") {
    if (sum(!noErr) != 0) {
      if (verbose)
        message(paste("Errors in", sum(!noErr), "gene(s):",
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
        message(paste("Errors in", sum(err), "gene(s):",
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
