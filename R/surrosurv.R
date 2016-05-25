#####################################################################################
#####################################################################################
#####################################################################################
surrosurv <- function(
  data, 
  models = c('Clayton', 'Plackett', 'Hougaard', 'Poisson'),
  intWidth, 
  cop.OPTIMIZER='bobyqa',
  poi.OPTIMIZER='bobyqa',
  verbose=FALSE) {
  # ******************************************************************************* #
  models <- match.arg(models, several.ok = TRUE)
  begin <- Sys.time()
  ### *** Parameter estimation *** ###############################################
  if (missing(intWidth))
    intWidth <- quantile(c(data$timeT, data$timeS), .4)
  
  # Copula approach
  copulas <- function(cop) {
    f <- function(){
      res <- try(copuSurr(data = data, family = cop, 
                          optimx.method = cop.OPTIMIZER, 
                          varcor1 = TRUE))
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
                        intWidth = intWidth, 
                        OPTIMIZER = poi.OPTIMIZER))
    if (class(res) == 'try-error') 
      res <- list(model2a = list(R2 = NA, kTau = NA), 
                  model2b = list(R2 = NA, kTau = NA), 
                  model3 = list(R2 = NA, kTau = NA), 
                  model5b = list(R2 = NA, kTau = NA))
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
      Poisson2a = fitRES$Poisson$model2a,
      Poisson2b = fitRES$Poisson$model2b,
      Poisson3  = fitRES$Poisson$model3,
      Poisson5b = fitRES$Poisson$model5b))
  
  if (all(c('R2', 'kTau') %in% names(attributes(data)))) {
    RES <- list('True Values' = attributes(data)[c('R2', 'kTau')])
  } else 
    RES <- list()
  
  RES <- c(RES, fitRES, list(intWidth = intWidth,
                             runTime = Sys.time() - begin))
  class(RES) <- c('surrosurv', class(RES))
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
      res <- rbind(x[[cop]][['unadj']][c('R2', 'kTau')],
                   x[[cop]][['adj']][c('R2', 'kTau')])
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
    res <- rbind(Poisson2a = x[['Poisson2a']][c('R2', 'kTau')],
                 Poisson2b = x[['Poisson2b']][c('R2', 'kTau')],
                 Poisson3  = x[['Poisson3']][c('R2', 'kTau')],
                 Poisson5b  = x[['Poisson5b']][c('R2', 'kTau')])
    res[sapply(res, is.null)] <- NA
    return(res)
  }
  
  fitRES <- do.call(rbind, lapply(models, function(x) eval(call(paste(x)))))
  
  RES <- fitRES
  if (('True Values' %in% names(x)) &
        (!all(sapply(x[['True Values']], is.null))))
    RES <- rbind('True Values' = x[['True Values']], RES)
  
  #   RES <- lapply(RES, function(x) {x[is.null(x)] <- NA; return(x)})
  if(silent) return(RES) else print(RES, na.print=na.print, digits=digits, ...)
}


#####################################################################################
#####################################################################################
#####################################################################################
convergence <- function(x, kkttol = 1e-8) {
  if (!'surrosurv' %in% class(x))
    stop('x must be of class surrosurv')
  
  models <- c('Clayton', 'Plackett', 'Hougaard', 'Poisson')
  models <- models[sapply(models, function(y) any(grepl(y, names(x))))]
  
  kktCtrl <- list(kkttolA = 1/kkttol,
                  kkt2tolA = kkttol)
  
  copulas <- function(cop) {
    f <- function() {
      res <- t(sapply(x[[cop]], function(y) {
        #                 library('optextras')
        kkt <- try(kktc(
          par = as.numeric(y$optimxRES[, 1:attr(y$optimxRES, 'npar')]),
          fval = y$optimxRES$value,
          ngr = attr(y$optimxRES, 'grad'), 
          nHes = attr(y$optimxRES, 'hessian'), 
          nbm=0, control = kktCtrl))
        if (class(kkt) == 'try-error') kkt <- c(kkt1 = NA, kkt2 = NA)
        return(rbind(c(
          unlist(kkt[c('kkt1', 'kkt2')]),
          VCpd = tryCatch({
            #                         library('corpcor')
            is.positive.definite(y$VarCor2, tol=kkttol)
          }, error=function(x) return(NA)))))
      }))
      rownames(res) <- paste(cop, rownames(res))
      return(res)
    }
    return(f)
  }
  Clayton = copulas('Clayton')
  Plackett = copulas('Plackett')
  Hougaard = copulas('Hougaard')
  Poisson = function() {
    #         library('optextras')
    which_models <- x[grepl('Poisson', names(x))]
    convres <- t(sapply(
      which_models,
      function (y) tryCatch({
        kkt <- kktc(par = y$optinfo$val,
                    fval = y$optinfo$control$maxfun,
                    ngr = y$optinfo$derivs$gradient, 
                    nHes = y$optinfo$derivs$Hessian, 
                    nbm=0, control = kktCtrl)
        #                 library('matrixcalc')
        return(c(unlist(kkt[c('kkt1', 'kkt2')]),
                 VCpd = tryCatch(
                   is.positive.definite(y$VarCor$trialref, tol=kkttol),
                   error = function(x) return(NA))))
      }, error=function(x) return(rep(NA, 3))))
    )
    names(convres) <- names(which_models)
    return(convres)
  }
  
  fitRES <- do.call(rbind, lapply(models, function(x) eval(call(paste(x)))))
  return(fitRES)
}