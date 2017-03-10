# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# surrosurv() ##################################################################
surrosurv <- function(
  data, 
  models        = c('Clayton', 'Plackett', 'Hougaard',
                    'Poisson I', 'Poisson T', 'Poisson TI', 'Poisson TIa'),
  intWidth      = NULL,
  nInts         = 8,
  cop.OPTIMIZER = 'bobyqa',
  poi.OPTIMIZER = 'bobyqa',
  verbose       = TRUE,
  twoStage       = FALSE,
  keep.data     = TRUE) {
  # ************************************************************************** #
  data$trialref <- factor(data$trialref)
  models <- tolower(noSpP(models))
  if ('poisson' %in% models) {
    models <- setdiff(models, 'poisson')
    models <- unique(c(models, paste0('poisson', c('i', 't', 'ti', 'tia'))))
  }
  models <- match.arg(models, several.ok = TRUE, choices = c(
    'clayton', 'plackett', 'hougaard', paste0('poisson', c('i', 't', 'ti', 'tia'))))
  Poissons <- models[grepl('poisson', models)]
  begin <- Sys.time()
  W <- options()$warn
  options(warn = -1)
  
  ####  Parameter estimation ####
  if (any(!grepl('poisson', models))) {
    # library('SurvCorr')
    INIrho <- survcorr(Surv(timeS, statusS) ~ 1, 
                       Surv(timeT, statusT) ~ 1, data = data)$rho
  } else INIkTau <- NULL
  
  # Copula approach
  copulas <- function(cop) {
    f <- function(){
      if (verbose) message(paste0(
        '- Estimating model: ', toupper(substr(cop, 1, 1)), substr(cop, 2, 100)),
        appendLF = FALSE)
      res <- try(copuSurr(data = data, family = cop, 
                          optimx.method = cop.OPTIMIZER, 
                          varcor1 = TRUE,
                          #INIkTau = INIkTau
                          INIrho = INIrho,
                          twoStage = twoStage
      ), silent = FALSE) #TRUE)
      if (class(res) == 'try-error') {
        res <- list(kTau = NA, alpha = NA, beta = NA,
                    R2 = NA, ranef = NA, 
                    VarCor2 = NA, VarCor1 = NA)
        res <- list(unadj = res, adj = res)
      }
      if (verbose)
        message(paste0(' (', format(res$adj$runTime, digits = 2), ')'))
      return(res)
    }
    return(f)
  }
  clayton = copulas('Clayton')
  plackett = copulas('Plackett')
  hougaard = copulas('Hougaard')
  
  # Poisson approach
  poisson = function(poissons = Poissons){
    res <- try(poisSurr(data = data, 
                        intWidth = intWidth, nInts = nInts,
                        OPTIMIZER = poi.OPTIMIZER, models = poissons,
                        verbose = verbose), 
               silent = TRUE)
    if (class(res) == 'try-error') 
      res <- list(modelT = list(R2 = NA, kTau = NA), 
                  modelI = list(R2 = NA, kTau = NA), 
                  modelTI = list(R2 = NA, kTau = NA), 
                  modelTIa = list(R2 = NA, kTau = NA))[
                    sub('poisson ', 'model', Poissons)]
    return(res)
  }
  # -------------------------------------------------------------------------- #
  models <- c(models[!grepl('poisson', models)],
              if (any(grepl('poisson', models))) 'poisson')
  if (verbose) message('Computation may take very long. Please wait...')
  fitRES <- lapply(models, function(x) {
    eval(call(paste(x)))
  })
  names(fitRES) <- models
  
  if ('poisson' %in% models)
    fitRES <- c(fitRES[!(names(fitRES) == 'poisson')], 
                fitRES$poisson)
  names(fitRES) <- sapply(names(fitRES), function(x) paste0(
    toupper(substr(x, 1, 1)), substr(x, 2, 100)))
  
  if (all(c('kTau', 'R2') %in% names(attributes(data)))) {
    RES <- list('True Values' = attributes(data)[c('kTau', 'R2')])
  } else 
    RES <- list()
  
  RES <- c(RES, fitRES)
  class(RES) <- c('surrosurv', class(RES))
  attributes(RES) <- c(
    attributes(RES),
    list(intWidth   = intWidth,
         runTime    = Sys.time() - begin,
         trialSizes = table(data$trialref)))
  if (keep.data) attr(RES, 'data') <- data
  options(warn = W)
  return(RES)
}

# print.surrosurv ##############################################################
print.surrosurv <- function(x, silent=FALSE, digits=2, na.print='-.--', ...) {
  models <- c('Clayton', 'Plackett', 'Hougaard', 'Poisson')
  models <- names(which(sapply(models, function(y) any(grepl(y, names(x))))))
  
  # Copula approach
  copulas <- function(cop) {
    f <- function() {
      res <- rbind(x[[cop]][['unadj']][c('kTau', 'R2')],
                   x[[cop]][['adj']][c('kTau', 'R2')])
      rownames(res) <- paste(cop, c('unadj', 'adj'))
      res[sapply(res, is.null)] <- NA
      return(res)
    }
    return(f)
  }
  Clayton = copulas('Clayton')
  Plackett = copulas('Plackett')
  Hougaard = copulas('Hougaard')
  # Poisson approach
  Poisson = function(){
    res <- t(sapply(x[grepl('Poisson', names(x ))], function(y) y[c('kTau', 'R2')]))
    res[sapply(res, is.null)] <- NA
    return(res)
  }
  
  fitRES <- do.call(rbind, lapply(models, function(x) eval(call(paste(x)))))
  
  RES <- fitRES
  if (('True Values' %in% names(x)) &
      (!all(sapply(x[['True Values']], is.null))))
    RES <- rbind('True Values' = x[['True Values']][c('kTau', 'R2')], RES)
  
  if (silent) return(RES) else print(
    RES, na.print = na.print, digits = digits, ...)
}


# confint.surrosurv ############################################################
confint.surrosurv <- function(
  object,
  parm,
  level = 0.95,
  method = c('boot'),
  nsim = 100,
  parallel = TRUE, 
  nCores, 
  models = names(object),
  intWidth,
  keep.allRes = FALSE,
  ...) {
  models <- tolower(noSpP(models))
  if ('poisson' %in% models) {
    models <- setdiff(models, 'poisson')
    models <- unique(c(models, paste0('poisson', c('i', 't', 'ti', 'tia'))))
  }
  models <- match.arg(models, several.ok = TRUE, choices = c(
    'clayton', 'plackett', 'hougaard', paste0('poisson', c('i', 't', 'ti', 'tia'))))
  intWidth <- attr(object, 'intWidth')
  
  # library('parallel')
  
  if (parallel) {
    totCores <- detectCores()
    
    if (missing(nCores)) {
      nCores <- min(nsim, totCores)
      message(paste0(
        'Parallel computing on ', nCores,
        ' cores (the total number of ', 
        ifelse(nsim < totCores, 'trials', 'cores detected'),
        ')'))
    } else {
      if (nCores > min(nsim, totCores))
        message(paste0(
          'The number of cores (nCores=', nCores, ') is greater than',
          'the number of ', 
          ifelse(nsim < totCores, 'trials', 'cores detected'),
          ')'))
      
      nCores <- min(nCores, nsim, totCores)
      message(paste('Parallel computing on', nCores,'cores'))
    }
  } else nCores <- 1
  
  if ('data' %in% names(attributes(object))) {
    data <- attr(object, 'data')
  } else {
    stop(paste(
      "The fitted models 'obect' must have an attribute 'data'",
      "containing the original data.",
      "See the option 'keep.data' for the surrosurv() function."))
  }
  
  if (Sys.info()[1] == "Windows") {
    cl <- makeCluster(nCores, type = 'PSOCK')
    clusterExport(cl, c('data', 'models', 'intWidth'), envir = environment())
  } else {
    cl <- makeCluster(nCores, type = 'FORK')
  }
  clusterEvalQ(cl, library('survival'))
  clusterEvalQ(cl, library('surrosurv'))
  ciRES <- clusterApplyLB(cl, 1:nsim, function(x) {
    return(surrosurv(data, models = models, intWidth = intWidth))
  }, ...) 
  stopCluster(cl)
  rm(cl)
  
  ciRES <- lapply(ciRES, print, silent = TRUE)
  arrayRES <- array(unlist(ciRES), dim = c(nrow(ciRES[[1]]), 
                                           ncol(ciRES[[1]]), 
                                           length(ciRES)))
  LCL <- apply(arrayRES, 1:2, quantile, p = .025)
  UCL <- apply(arrayRES, 1:2, quantile, p = .975)
  dimnames(LCL) <- dimnames(UCL) <- dimnames(ciRES[[1]])
  
  res <- list(LCL = LCL, UCL = UCL)
  if (keep.allRes) res <- c(res, allRes = ciRES)
  
  class(res) <- c('ciSurrosurv', class(res))
  return(res)
}

print.ciSurrosurv <- function(x, ...) {
  res <- list(kTau = cbind(' 2.5 %' = x$LCL[, 1],
                           '97.5 %' = x$UCL[, 1]),
              R2   = cbind(' 2.5 %' = x$LCL[, 2],
                           '97.5 %' = x$UCL[, 2]))
  lapply(names(res), function(p) {
    cat(paste('\n', p, '\n'))
    print(res[[p]])
    })
}
