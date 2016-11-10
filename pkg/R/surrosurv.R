#####################################################################################
#####################################################################################
#####################################################################################
surrosurv <- function(
  data, 
  models = c('Clayton', 'Plackett', 'Hougaard',
             'Poisson I', 'Poisson T', 'Poisson TI', 'Poisson TIa'),
  intWidth = NULL,  nInts = 8,
  cop.OPTIMIZER = 'bobyqa',
  poi.OPTIMIZER = 'bobyqa',
  verbose = FALSE) {
  # ************************************************************************** #
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
  options(warn=-1)
  
  ### *** Parameter estimation *** #############################################
  #   if (missing(intWidth)) {
  #     times <- c(data$timeT, data$timeS)
  #     status <- c(data$statusT, data$statusS)
  #     sfit <- survfit(Surv(times, status) ~ 1)
  #     intWidth <- #quantile(c(data$timeT, data$timeS), .4)
  #       sfit$time[which(sfit$surv <= .6)[1]]
  #     rm(times, status, sfit)
  #   }  
  
  if (any(!grepl('poisson', models))) {
    # library('SurvCorr')
    INIrho <- survcorr(Surv(timeS, statusS) ~ 1, 
                        Surv(timeT, statusT) ~ 1, data = data)$rho
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
      if (class(res) == 'try-error') {
        res <- list(kTau = NA, alpha = NA, beta = NA,
                    R2 = NA, ranef = NA, 
                    VarCor2 = NA, VarCor1 = NA)
        res <- list(unadj = res, adj = res)
      }
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
                        OPTIMIZER = poi.OPTIMIZER, models = poissons), 
               silent = TRUE)
    if (class(res) == 'try-error') 
      res <- list(modelT = list(R2 = NA, kTau = NA), 
                  modelI = list(R2 = NA, kTau = NA), 
                  modelTI = list(R2 = NA, kTau = NA), 
                  modelTIa = list(R2 = NA, kTau = NA))[
                    sub('poisson ', 'model', Poissons)]
    return(res)
  }
  ##############################################################################
  models <- c(models[!grepl('poisson', models)],
              ifelse(any(grepl('poisson', models)), 'poisson', ''))
  fitRES <- lapply(models, function(x, verb = verbose) {
    if (verbose) message(paste('Estimating model:', x))
    eval(call(paste(x)))
  })
  names(fitRES) <- models
  
  if ('poisson' %in% models)
    fitRES <- c(fitRES[!(names(fitRES) == 'poisson')], 
                # list(
                #   PoissonT = fitRES$Poisson$modelT,
                #   PoissonI = fitRES$Poisson$modelI,
                #   PoissonTI  = fitRES$Poisson$modelTI,
                #   PoissonTIa = fitRES$Poisson$modelTIa)
                fitRES$poisson)
  names(fitRES) <- sapply(names(fitRES), function(x) paste0(
    toupper(substr(x, 1, 1)), substr(x, 2, 100)))
  
  if (all(c('kTau', 'R2') %in% names(attributes(data)))) {
    RES <- list('True Values' = attributes(data)[c('kTau', 'R2')])
  } else 
    RES <- list()
  
  RES <- c(RES, fitRES, list(intWidth = intWidth,
                             runTime = Sys.time() - begin))
  class(RES) <- c('surrosurv', class(RES))
  attr(RES, 'trialSizes') <- table(data$trialref)
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
    res <- 
      # rbind(PoissonT = x[['PoissonT']][c('kTau', 'R2')],
      #            PoissonI = x[['PoissonI']][c('kTau', 'R2')],
      #            PoissonTI  = x[['PoissonTI']][c('kTau', 'R2')],
      #            PoissonTIa  = x[['PoissonTIa']][c('kTau', 'R2')])
      t(sapply(x[grepl('Poisson', names(x ))], function(y) y[c('kTau', 'R2')]))
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
