#####################################################################################
#####################################################################################
#####################################################################################
surrosurv <- function(
  data, 
  models = c('Clayton', 'Plackett', 'Hougaard', 'Poisson'),
  intWidth = NULL,  nInts = 8,
  cop.OPTIMIZER = 'bobyqa',
  poi.OPTIMIZER = 'bobyqa',
  verbose = FALSE) {
  # ******************************************************************************* #
  models <- match.arg(models, several.ok = TRUE)
  begin <- Sys.time()
  W <- options()$warn
  options(warn=-1)
  
  ### *** Parameter estimation *** ###############################################
  #   if (missing(intWidth)) {
  #     times <- c(data$timeT, data$timeS)
  #     status <- c(data$statusT, data$statusS)
  #     sfit <- survfit(Surv(times, status) ~ 1)
  #     intWidth <- #quantile(c(data$timeT, data$timeS), .4)
  #       sfit$time[which(sfit$surv <= .6)[1]]
  #     rm(times, status, sfit)
  #   }  
  
  if (length(setdiff(models, 'Poisson'))) {
    # library('SurvCorr')
    INIrho <- survcorr(Surv(timeS, statusS) ~ 1, 
                        Surv(timeT, statusT) ~ 1, data = data)$rho
    # library('NADA')
      #INIkTau <- cenken(-data$timeS, (data$statusS == 0),
      #                       -data$timeT, (data$statusT == 0))$tau
  } else INIkTau <- NULL
  
  # Copula approach
  copulas <- function(cop) {
    f <- function(){
      res <- try(copuSurr(data = data, family = cop, 
                          optimx.method = cop.OPTIMIZER, 
                          varcor1 = TRUE,
                          #INIkTau = INIkTau
                          INIrho = INIrho
                          ), silent = TRUE)
      if (class(res) == 'try-error'){
        res <- list(kTau = NA, alpha = NA, beta = NA,
                    R2 = NA, ranef = NA, 
                    VarCor2 = NA, VarCor1 = NA)
        res <- list(unadj = res, adj = res)
      }
      return(res)
    }
    return(f)
  }
  Clayton = copulas('Clayton')
  Plackett = copulas('Plackett')
  Hougaard = copulas('Hougaard')
  
  # Poisson approach
  Poisson = function(){
    res <- try(poisSurr(data = data, 
                        intWidth = intWidth, nInts = nInts,
                        OPTIMIZER = poi.OPTIMIZER), silent = TRUE)
    if (class(res) == 'try-error') 
      res <- list(modelT = list(R2 = NA, kTau = NA), 
                  modelI = list(R2 = NA, kTau = NA), 
                  modelTI = list(R2 = NA, kTau = NA), 
                  modelTIa = list(R2 = NA, kTau = NA))
    return(res)
  }
  ################################################################################
  fitRES <- lapply(models, function(x, verb=verbose) {
    if (verbose) message(paste('Estimating model:', x))
    eval(call(paste(x)))
  })
  names(fitRES) <- models
  if ('Poisson' %in% models)
    fitRES <- c(fitRES[!(names(fitRES) == 'Poisson')], list(
      PoissonT = fitRES$Poisson$modelT,
      PoissonI = fitRES$Poisson$modelI,
      PoissonTI  = fitRES$Poisson$modelTI,
      PoissonTIa = fitRES$Poisson$modelTIa))
  
  if (all(c('kTau', 'R2') %in% names(attributes(data)))) {
    RES <- list('True Values' = attributes(data)[c('kTau', 'R2')])
  } else 
    RES <- list()
  
  RES <- c(RES, fitRES, list(intWidth = intWidth,
                             runTime = Sys.time() - begin))
  class(RES) <- c('surrosurv', class(RES))
  options(warn = W)
  return(RES)
}

#####################################################################################
#####################################################################################
#####################################################################################
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
    res <- rbind(PoissonT = x[['PoissonT']][c('kTau', 'R2')],
                 PoissonI = x[['PoissonI']][c('kTau', 'R2')],
                 PoissonTI  = x[['PoissonTI']][c('kTau', 'R2')],
                 PoissonTIa  = x[['PoissonTIa']][c('kTau', 'R2')])
    res[sapply(res, is.null)] <- NA
    return(res)
  }
  
  fitRES <- do.call(rbind, lapply(models, function(x) eval(call(paste(x)))))
  
  RES <- fitRES
  if (('True Values' %in% names(x)) &
        (!all(sapply(x[['True Values']], is.null))))
    RES <- rbind('True Values' = x[['True Values']][c('kTau', 'R2')], RES)
  
  #   RES <- lapply(RES, function(x) {x[is.null(x)] <- NA; return(x)})
  if(silent) return(RES) else print(RES, na.print=na.print, digits=digits, ...)
}
