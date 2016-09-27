################################################################################
loocv <- function(data, nCores, parallel = TRUE, ...) {
  # ************************************************************************** #
  # models <- match.arg(models, several.ok = TRUE)
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
  
  loof <- function(TRIAL, ...) {
    alpha <- coef(coxph(Surv(timeS, statusS) ~ trt, 
                        data = data[data$trialref == TRIAL, ]))[['trt']]
    beta <- coef(coxph(Surv(timeT, statusT) ~ trt, 
                        data = data[data$trialref == TRIAL, ]))[['trt']]
    reddata <- data[data$trialref != TRIAL, ]
    reddata$trialref <- factor(reddata$trialref)
    predint <- tryCatch({
      surrofit <- surrosurv(data = reddata, ...)
      lapply(attr(predict(RES), 'predf'), function(x) x(alpha))
    }, error = function(e) return(c(fit = NA, lwr = NA, upr = NA)))
    RES <- c(list(margPars = c(alpha = alpha, beta = beta)),
             predint)
    return(RES)
  }
  
  cl <- makeCluster(nCores)
  clusterExport(cl, 'data')
  clusterEvalQ(cl, library('survival'))
  clusterEvalQ(cl, library('surrosurv'))
  loocvRES <- clusterApplyLB(cl, levels(data$trialref), loof)  
  stopCluster(cl)
  rm(cl)
  
  class(loocvRES) <- c('loocvSurrosurv', class(loocvRES))
  return(loocvRES)
}
################################################################################


################################################################################
print.loocvSurrosurv <- function(
  x, 
  silent = FALSE, 
  digits = 2, 
  na.print = '-.--', 
  ...) {
  # ************************************************************************** #
  RES <- x
  
  if(silent) return(RES) else print(RES, na.print=na.print, digits=digits, ...)
}
################################################################################


################################################################################
plot.loocvSurrosurv <- function(
  x,
  ...) {
  # ************************************************************************** #
}
################################################################################
