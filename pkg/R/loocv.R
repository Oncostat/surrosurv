################################################################################
loocv <- function(data, nCores, parallel = TRUE, 
                  models = c('Clayton', 'Plackett', 'Hougaard',
                             'Poisson I', 'Poisson T', 'Poisson TI', 'Poisson TIa'),
                  ...) {
  # ************************************************************************** #
  models <- tolower(noSpP(models))
  if ('poisson' %in% models) {
    models <- setdiff(models, 'poisson')
    models <- unique(c(models, paste0('poisson', c('i', 't', 'ti', 'tia'))))
  }
  models <- match.arg(models, several.ok = TRUE, choices = c(
    'clayton', 'plackett', 'hougaard', paste0('poisson', c('i', 't', 'ti', 'tia'))))
  data$trialref <- factor(data$trialref)
  
  # library('parallel')
  
  if (parallel) {
    totCores <- detectCores()
    nTrials <- nlevels(data$trialref)
    
    if (missing(nCores)) {
      nCores <- min(nTrials, totCores)
      message(paste0(
        'Parallel computing on ', nCores,
        ' cores (the total number of ', 
        ifelse(nTrials < totCores, 'trials', 'cores detected'),
        ')'))
    } else {
      if (nCores > min(nTrials, totCores))
        message(paste0(
          'The number of cores (nCores=', nCores, ') is greater than',
          'the number of ', 
          ifelse(nTrials < totCores, 'trials', 'cores detected'),
          ')'))
      
      nCores <- min(nCores, nTrials, totCores)
      message(paste('Parallel computing on', nCores,'cores'))
    }
  } else nCores <- 1
  
  loof <- function(TRIAL, models2predict = models, ...) {
    alpha <- coef(coxph(Surv(timeS, statusS) ~ trt, 
                        data = data[data$trialref == TRIAL, ]))[['trt']]
    beta <- coef(coxph(Surv(timeT, statusT) ~ trt, 
                       data = data[data$trialref == TRIAL, ]))[['trt']]
    reddata <- data[data$trialref != TRIAL, ]
    reddata$trialref <- factor(reddata$trialref)
    predint <- tryCatch(
      expr = {
        surrofit <- surrosurv(data = reddata, models = models2predict, ...)
        lapply(attr(predict(surrofit), 'predf'), function(x) x(alpha))
      },
      error = function(e) {
        for (cop in c('Clayton', 'Plackett', 'Hougaard')) {
          coppos <- which(models2predict == cop)
          if (length(coppos) == 1) {
            models2predict <- c(
              models2predict[1:coppos - 1],
              paste(cop, c('unadj', 'adj'), sep='.'),
              models2predict[(coppos + 1):length(models2predict)])
          }
        }
        poipos <- which(models2predict == 'Poisson')
        if (length(poipos) == 1) {
          models2predict <- c(
            models2predict[1:poipos - 1],
            paste0('Poisson', c('T', 'TI', 'TIa')))
        }
        models2predict
        return(
          mapply(function(i) return(c(fit = NA, lwr = NA, upr = NA)), 
                 models2predict, SIMPLIFY = FALSE))
      })
    RES <- c(list(margPars = c(alpha = alpha, beta = beta)),
             predint)
    return(RES)
  }
  
  cl <- makeCluster(nCores)
  clusterExport(cl, 'data')
  clusterEvalQ(cl, library('survival'))
  clusterEvalQ(cl, library('surrosurv'))
  loocvRES <- clusterApplyLB(cl, levels(data$trialref), loof, ...) 
  stopCluster(cl)
  rm(cl)
  
  names(loocvRES) <- levels(data$trialref)
  class(loocvRES) <- c('loocvSurrosurv', class(loocvRES))
  return(loocvRES)
}
################################################################################


################################################################################
print.loocvSurrosurv <- function(x, n = 6, silent = FALSE, ...) {
  # ************************************************************************** #
  models <- setdiff(names(x[[1]]), 'margPars')
  RES <- lapply(models, function(y){
    preds <- sapply(x, function(trial)
      c(obsBeta = trial$margPars['beta'], trial[[y]][-1])
    )
    rownames(preds) <- c('obsBeta', 'lwr', 'upr')
    return(preds)
  }
  )
  names(RES) <- models
  
  if (silent) return(RES)
  
  for (i in 1:length(RES)) {
    method <- format.methodNames(RES)[i]
    cat('\n  ', method, '\n')
    res2print <- format(RES[[i]][, 1:n], digits = 1, na.encode = FALSE)
    if (nrow(RES[[i]] > n))
      res2print <- cbind(res2print, '  ' = c('...', '...'))
    # rownames(res2print) <- paste0(
    #   sub('trt', '    Treatment effects on ', rownames(res2print)), ':')
    print(res2print, quote = FALSE, ...)
  }
}
################################################################################


################################################################################
plot.loocvSurrosurv <- function(x, 
                                models, 
                                exact.models, ...) {
  # ************************************************************************** #
  x <- print(x, silent = TRUE)
  if (missing(models)) models <- setdiff(names(x), 'margPars')
  
  if (missing(exact.models))
    exact.models <- any(tolower(noSpP(names(x))) %in% tolower(noSpP(models)))
  
  if (exact.models) {
    ind <- which(tolower(noSpP(names(x))) %in% tolower(noSpP(models)))
  } else {
    ind <- which(sapply(tolower(noSpP(names(x))), function(mod)
      !all(is.na(pmatch(tolower(noSpP(models)), mod)))))
  }
  
  if (length(ind)) {
    par(mfrow = n2mfrow(length(ind)))
    for (i in ind) {
      plot(NA, main = format.methodNames(x)[i],
           xlim = c(0, ncol(x[[i]])) + .5, 
           xlab = 'Trials', xaxt = 'n', 
           ylim = range(x[[i]], na.rm = TRUE) + 
             c(-1, 1) * diff(range(x[[i]], na.rm = TRUE)) / 20,
           ylab = "Hazard ratio on the true endpoint", yaxt = 'n',
           panel.first = abline(h = 0, col = 'grey'))
      axis(1, 1:ncol(x[[i]]), labels = colnames(x[[i]]))
      axis(2, axTicks(2), format(round(exp(axTicks(2)), 2)), las = 1)
      segments(1:ncol(x[[i]]), x[[i]][2, ], y1 = x[[i]][3, ], 
               col = rgb(.8, .8, .8, .8), lwd = 5)
      COLs <- 2 - (colSums(apply(x[[i]], 2, function(x)
        order(x) == c(2, 1, 3))) == 3)
      NC <- is.na(x[[i]][2, ]) | is.na(x[[i]][3, ])
      COLs[NC] <- 0
      points(1:ncol(x[[i]]), x[[i]][1, ], pch = 16, cex = 1.4, col = COLs)
      if (any(NC))
        mtext('x', 1, -1, at = which(NC), col=2, font=2)
    }
  }
}
################################################################################
