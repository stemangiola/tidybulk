## Source car:::ConjComp
car_ConjComp <- function (X, Z = diag(nrow(X)), ip = diag(nrow(X))){
  xq <- qr(t(Z) %*% ip %*% X)
  if (xq$rank == 0)
    return(Z)
  Z %*% qr.Q(xq, complete = TRUE)[, -(1:xq$rank)]
}

## Source car:::relatives
car_relatives <- function (term, names, factors)
{
  is.relative <- function(term1, term2) {
    all(!(factors[, term1] & (!factors[, term2])))
  }
  if (length(names) == 1)
    return(NULL)
  which.term <- which(term == names)
  (1:length(names))[-which.term][sapply(names[-which.term],
                                        function(term2) is.relative(term, term2))]
}

#' @importFrom stats pchisq
organiseStats = function (resultList, test.stat)
{
  statsList <- lapply(resultList, "[[", "stats")
  s <- do.call(rbind, statsList)
  coefList <- lapply(resultList, "[[", "coef")
  cf <- do.call(rbind, coefList)
  SEList <- lapply(resultList, "[[", "stdErr")
  stdErr <- do.call(rbind, SEList)
  if (test.stat == "Wald") {
    chisqList <- lapply(resultList, "[[", "chisq")
    chisq <- do.call(rbind, chisqList)
    dfList <- lapply(resultList, "[[", "df")
    df <- do.call(rbind, dfList)
    pvals <- pchisq(chisq, df = df, lower.tail = FALSE)
    colnames(df) <- colnames(chisq)
    colnames(pvals) <- colnames(chisq)
    s <- list(res = s, coef = cf, stdErr = stdErr, Chisq = chisq,
              Df = df, pvals = pvals)
  }
  else if (test.stat == "F") {
    NumDF <- lapply(resultList, function(x) x$test[, 1])
    NumDF <- do.call(rbind, NumDF)
    DenDF <- lapply(resultList, function(x) x$test[, 2])
    DenDF <- do.call(rbind, DenDF)
    Fval <- lapply(resultList, function(x) x$test[, 3])
    Fval <- do.call(rbind, Fval)
    pvals <- lapply(resultList, function(x) x$test[, 4])
    pvals <- do.call(rbind, pvals)
    if (ncol(NumDF) == 1) {
      colnames(NumDF) <- colnames(DenDF) <- colnames(Fval) <- colnames(pvals) <- rownames(resultList[[1]]$test)
    }
    s <- list(res = s, coef = cf, stdErr = stdErr, NumDF = NumDF,
              DenDF = DenDF, Fval = Fval, pvals = pvals)
  }
  else {
    LRT <- lapply(resultList, "[[", "test")
    LRT <- do.call(rbind, LRT)
    chisq <- LRT[, 1, drop = FALSE]
    df <- LRT[, 2, drop = FALSE]
    pvals <- LRT[, 3, drop = FALSE]
    colnames(chisq) <- colnames(df) <- colnames(pvals) <- "LRT"
    s <- list(res = s, coef = cf, stdErr = stdErr, Chisq = chisq,
              Df = df, pvals = pvals)
  }
}

hyp_matrix = function (fullFormula, metadata, LHS)
{
  reduced2 <- lme4::nobars(fullFormula)
  fac <- attr(terms(reduced2), "factors")
  data2 <- metadata
  data2[, LHS] <- rep(0, nrow(data2))
  dm2 <- model.matrix(reduced2, data2)
  assign <- attr(dm2, "assign")
  term.labels <- attr(terms(reduced2), "term.labels")
  p <- length(assign)
  I.p <- diag(p)
  n.terms <- length(term.labels)
  hyp.matrix.1 <- hyp.matrix.2 <- list()
  for (i in seq_len(n.terms)) {
    which.term <- i
    subs.term <- which(assign == which.term)
    relatives <- car_relatives(term.labels[i], term.labels,
                               fac)
    subs.relatives <- NULL
    for (relative in relatives) subs.relatives <- c(subs.relatives,
                                                    which(assign == relative))
    hyp.matrix.1[[i]] <- I.p[subs.relatives, , drop = FALSE]
    hyp.matrix.2[[i]] <- I.p[c(subs.relatives, subs.term),
                             , drop = FALSE]
  }
  names(hyp.matrix.1) <- term.labels
  return(list(hyp.matrix.1 = hyp.matrix.1, hyp.matrix.2 = hyp.matrix.2))
}

lmer_wald  = function (fixef, hyp.matrix, vcov.)
{
  hyp.matrix.1 <- hyp.matrix[[1]]
  hyp.matrix.2 <- hyp.matrix[[2]]
  hyp.list <- lapply(seq_along(hyp.matrix.1), function(i) {
    hyp.matrix.term <- if (nrow(hyp.matrix.1[[i]]) == 0) {
      hyp.matrix.2[[i]]
    }
    else t(car_ConjComp(t(hyp.matrix.1[[i]]), t(hyp.matrix.2[[i]]),
                        vcov.))
    hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term,
                                              1, function(x) all(x == 0)), , drop = FALSE]
    hyp.matrix.term
  })
  b <- fixef
  V <- vcov.
  chi_val <- lapply(hyp.list, function(L) {
    as.vector(t(L %*% b) %*% solve(L %*% V %*% t(L)) %*%
                (L %*% b))
  })
  df <- unlist(lapply(hyp.list, NROW))
  list(chisq = setNames(unlist(chi_val), names(hyp.matrix.1)),
       df = df)
}


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

    if(vcov. |> as.numeric() |> is.nan() |> all())
      return(list(stats = NA, coef = NA, stdErr = NA, chisq = NA,
                  df = NA, predict = NA, optinfo = NA, ci_random_effect_df = NA,
                  message = "vcov. failed to calculate", tryErrors = fit[1]))

    fixedEffects <- glmmTMB::fixef(fit)$cond
    disp <- glmmTMB::sigma(fit)
    msg <- fit$fit$message
    stats <- setNames(c(disp, AIC(fit), as.numeric(logLik(fit))),
                      c("Dispersion", "AIC", "logLik"))
    if (is.null(reduced)) {
      test <- lmer_wald(fixedEffects, hyp.matrix, vcov.)
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
                df = NA, predict = NA, optinfo = NA, ci_random_effect_df = NA, message = "COMPLETE FIT ERROR", tryErrors = fit[1]))
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
    test <- lmer_wald(fixedEffects, hyp.matrix, vcov.)
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

# From glmmSeq

# Class definitions

setClassUnion("character_or_list", c("character", "list"))
setClassUnion("df_or_matrix", c("data.frame", "matrix"))
setClassUnion("list_or_matrix", c("list", "matrix"))
setClassUnion("formulaOrNULL", c("formula", "NULL"))

#' An S4 class to define the glmmSeq output
#'
#' @keywords internal
#' @noRd
#'
#' @slot info List including the matched call, dispersions, offset, designMatrix
#' @slot formula The model formula
#' @slot stats Statistics from fitted models
#' @slot predict Predicted values
#' @slot reduced Optional reduced formula for LRT
#' @slot countdata The input expression data with count data in rows
#' @slot metadata The input metadata
#' @slot modelData Model data for predictions
#' @slot optInfo Information on whether the model was singular or converged
#' @slot errors Any errors
#' @slot vars List of variables stored from the original call, including the
#'   `id` variable (by default automatically identified from the random effect
#'   term in the model) and `removeSingles` argument

setClass("GlmmSeq", slots = list(
  info = "list",
  formula = "formula",
  stats = "list_or_matrix",
  predict = "df_or_matrix",
  reduced = "formulaOrNULL",
  countdata = "df_or_matrix",
  metadata = "df_or_matrix",
  modelData = "df_or_matrix",
  optInfo = "matrix",
  errors = "character_or_list",
  vars = "list"
))

glmmSeq = function (modelFormula, countdata, metadata, id = NULL, dispersion = NA,
                    sizeFactors = NULL, reduced = NULL, modelData = NULL, designMatrix = NULL,
                    method = c("lme4", "glmmTMB"), control = NULL, family = glmmTMB::nbinom2,
                    cores = 1, removeSingles = FALSE, zeroCount = 0.125, verbose = TRUE,
                    returnList = FALSE, progress = FALSE, max_rows_for_matrix_multiplication = Inf, avoid_forking = FALSE, ...)
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
    hyp.matrix <- hyp_matrix(fullFormula, metadata, "count")
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
      if(avoid_forking){
        library(parallel)
        cl = parallel::makeCluster(cores, type = "PSOCK")
        #parallel::clusterEvalQ(cl,c(library(dplyr),library(glmmSeq)))
        #clusterExport(cl, list("varname1", "varname2"),envir=environment())
        resultList <- parallel::clusterApply(
          cl,
          fullList,
          function(geneList) {
            glmerCore(geneList, fullFormula, reduced,
                      subsetMetadata, control, offset, modelData,
                      designMatrix, hyp.matrix,, max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication, ...)
          }
        )
      } 
      else if (progress) {
        
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
  s <- organiseStats(resultList[noErr], "Wald")
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

# From GlmmSeq

#' Summarise a 'glmmSeq'/'lmmSeq' object
#'
#' Summarise results from [glmmSeq] or [lmmSeq] analysis
#'
#' @keywords internal
#' @noRd
#'
#' @param object an object of class `"GlmmSeq"` or `"lmmSeq"`
#' @param gene an optional character value specifying a single gene whose
#'   results are summarised
#' @param digits integer, used for number formatting
#' @param ... arguments to be passed to other methods
#' @return
#' If `gene=NULL` a dataframe of results for all genes is returned. Otherwise
#' the output of GLMM or LMM model results for a single gene including
#' coefficients, test statistics, p-values is printed and the dataframe for all
#' genes is returned invisibly.
#' @seealso [glmmSeq()], [lmmSeq()]
summary_lmmSeq <- function(object,
                           gene = NULL,
                           digits = max(3L, getOption("digits") - 3L), ...) {
  if (is.null(gene)) {
    statSet <- names(object@stats)
    gp <- lapply(statSet, function(i) {
      out <- object@stats[[i]]
      if (i %in% c("Chisq", "Fval", "Df", "NumDF", "DenDF")) colnames(out) <- paste(i, colnames(out), sep = "_")
      if (i == "pvals") colnames(out) <- paste0("P_", colnames(out))
      if (i == "stdErr") colnames(out) <- paste0("se_", colnames(out))
      out
    })
    do.call(cbind, gp)
  } else {
    out <- lapply(object@stats, function(i) i[gene, ])
    if (is(object, "GlmmSeq")) {
      cat("Generalised linear mixed model\n")
      cat(paste0("Method: ", object@info$method, "\n"))
      if (object@info$method == "lme4") {
        cat("Family: Negative Binomial\n")
      } else {
        cat(paste0("Family: ", object@info$family, "\n"))
      }
    } else {
      cat("Linear mixed model\n")
    }
    cat("Formula: ")
    print(object@formula)
    print(out$res)
    cat("\nFixed effects:\n")
    cfdf <- data.frame(Estimate = out$coef,
                       `Std. Error` = out$stdErr, check.names = FALSE)
    print(cfdf, digits = digits)
    if (object@info$test.stat == "Wald") {
      cat("\nAnalysis of Deviance Table (Type II Wald chisquare tests)\n")
      testdf <- data.frame(Chisq = out$Chisq, Df = out$Df,
                           `Pr(>Chisq)` = out$pvals,
                           row.names = colnames(object@stats$Chisq),
                           check.names = FALSE)
    } else if (object@info$test.stat == "LRT") {
      cat("\nLikelihood ratio test\nReduced formula: ")
      print(object@reduced)
      testdf <- data.frame(Chisq = out$Chisq, Df = out$Df,
                           `Pr(>Chisq)` = out$pvals,
                           row.names = " ",
                           check.names = FALSE)
    } else {
      cat("\nType III Analysis of Variance Table with Satterthwaite's method\n")
      testdf <- data.frame(NumDF = out$NumDF, DenDF = out$DenDF,
                           `F value` = out$Fval, `Pr(>F)` = out$pvals,
                           row.names = colnames(object@stats$Fval),
                           check.names = FALSE)
    }
    print(testdf, digits = digits)
    invisible(out)
  }
}
