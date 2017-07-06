################################################################################
loocv <- function(object, ...)
  UseMethod('loocv')
################################################################################
loocv.surrosurv <- function(object,
                            models,
                            nCores,
                            parallel = TRUE,
                            ...) {
  if (missing(models)) {
    allmodels <- c(
      'Clayton',
      'Plackett',
      'Hougaard',
      'Poisson I',
      'Poisson T',
      'Poisson TI',
      'Poisson TIa'
    )
    models <- allmodels[allmodels %in% names(object)]
  }
  loocv.data.frame(attr(object, 'data'), models = models, ...)
}
################################################################################
loocv.data.frame <- function(object,
                             models = c('Clayton', 'Poisson TI'),
                             nCores,
                             parallel = TRUE,
                             ...) {
  # ************************************************************************** #
  models <- tolower(noSpP(models))
  if ('poisson' %in% models) {
    models <- setdiff(models, 'poisson')
    models <-
      unique(c(models, paste0('poisson', c(
        'i', 't', 'ti', 'tia'
      ))))
  }
  models <- match.arg(
    models,
    several.ok = TRUE,
    choices = c('clayton', 'plackett', 'hougaard', paste0('poisson', c(
      'i', 't', 'ti', 'tia'
    )))
  )
  object$trialref <- factor(object$trialref)
  
  nTrials <- nlevels(object$trialref)
  if (nTrials < 3)
    stop('At least three trials are needed for cross-validation.')
  # library('parallel')
  
  if (parallel) {
    totCores <- detectCores()
    
    if (missing(nCores))
      nCores <- NA
    if (is.na(as.integer(nCores))) {
      nCores <- min(nTrials, totCores)
      message(
        paste0(
          'Parallel computing on ',
          nCores,
          ' cores (the total number of ',
          ifelse(nTrials < totCores, 'trials', 'cores detected'),
          ')'
        )
      )
    } else {
      if (nCores > min(nTrials, totCores))
        message(
          paste0(
            'The number of cores (nCores=',
            nCores,
            ') is greater than',
            'the number of ',
            ifelse(nTrials < totCores, 'trials', 'cores detected'),
            ' (',
            min(nTrials, totCores),
            ')'
          )
        )
      
      nCores <- min(nCores, nTrials, totCores)
      message(paste('Parallel computing on', nCores, 'cores'))
    }
  } else
    nCores <- 1
  
  loof <- function(TRIAL, models2predict = models, ...) {
    alpha <- coef(coxph(Surv(timeS, statusS) ~ trt,
                        data = object[object$trialref == TRIAL,]))[['trt']]
    beta <- coef(coxph(Surv(timeT, statusT) ~ trt,
                       data = object[object$trialref == TRIAL,]))[['trt']]
    redobject <- object[object$trialref != TRIAL,]
    redobject$trialref <- factor(redobject$trialref)
    predint <- tryCatch(
      expr = {
        surrofit <- surrosurv(data = redobject, models = models2predict, ...)
        res <- lapply(attr(predict(surrofit), 'predf'), function(x) x(alpha))
        names <- names(res)
        res <- lapply(1:length(res), function(i)
          rbind(res[[i]], t(print(
            surrofit, silent = TRUE
          ))[, i, drop = FALSE]))
        names(res) <- names
        attr(res, 'convergence') <- convergence(surrofit)
        attr(res, 'convals') <- convals(surrofit)
        res
      },
      error = function(e) {
        for (cop in c('clayton', 'plackett', 'hougaard')) {
          coppos <- which(models2predict == cop)
          if (length(coppos) == 1) {
            models2predict <- c(models2predict[1:coppos - 1],
                                paste(cop, c('unadj', 'adj'), sep = '.'),
                                as.character(na.omit(models2predict[
                                  (coppos + 1):(length(models2predict) + 1)])))
          }
        }
        poipos <- which(models2predict == 'Poisson')
        if (length(poipos) == 1) {
          models2predict <- c(models2predict[0:(poipos - 1)],
                              paste0('Poisson', c('T', 'TI', 'TIa')))
        }
        res <- mapply(function(i)
          return(t(t(c(fit = NA,
                       lwr = NA,
                       upr = NA,
                       kTau = NA,
                       R2 = NA)
          ))),
          models2predict, SIMPLIFY = FALSE)
        attr(res, 'convergence') <- attr(res, 'convals') <- 
          matrix(NA, length(models2predict), 3)
        return(res)
      }
    )
    RES <- c(list(margPars = c(alpha = alpha, beta = beta)),
             predint)
    attr(RES, 'convergence') <- attr(predint, 'convergence')
    attr(RES, 'convals') <- attr(predint, 'convals')
    return(RES)
  }
  
  if (Sys.info()[1] == "Windows") {
    cl <- makeCluster(nCores, type = 'PSOCK')
    clusterExport(cl, 'object', environment())
  } else {
    cl <- makeCluster(nCores, type = 'FORK')
  }
  clusterEvalQ(cl, library('survival'))
  clusterEvalQ(cl, library('surrosurv'))
  loocvRES <- clusterApplyLB(cl, levels(object$trialref), loof, ...)
  stopCluster(cl)
  rm(cl)
  
  names(loocvRES) <- levels(object$trialref)
  class(loocvRES) <- c('loocvSurrosurv', class(loocvRES))
  return(loocvRES)
}
################################################################################


################################################################################
print.loocvSurrosurv <- function(x,
                                 n = min(length(x), 6),
                                 silent = FALSE,
                                 ...) {
  # ************************************************************************** #
  models <- setdiff(names(x[[1]]), 'margPars')
  RES <- lapply(models, function(y) {
    preds <- sapply(x, function(trial) {
      trialRes <- if (is.null(trial[[y]])) {
        rep(NA, 5)
      } else {
        trial[[y]]
      }
      c(obsAlpha = trial$margPars['alpha'], 
        obsBeta = trial$margPars['beta'], 
        trialRes)
    })
    preds <- matrix(unlist(preds),
                    nrow = 7,
                    dimnames = list(
                      c('obsAlpha', 'obsBeta', 'predict', 'lwr', 'upr', 'kTau', 'R2'),
                      names(x)
                    ))
    return(preds)
  })
  names(RES) <- models
  
  if (silent)
    return(RES)
  
  n <- min(n, ncol(RES[[1]]))
  for (i in 1:length(RES)) {
    method <- format.methodNames(RES)[i]
    cat('\n  ', method, '\n')
    res2print <-
      format(
        as.data.frame(RES[[i]][, 1:n, drop = FALSE]),
        digits = 1,
        nsmall = 2,
        na.encode = FALSE
      )
    if (ncol(RES[[i]]) > n)
      res2print <-
      cbind(res2print, '  ' = rep('...', nrow(res2print)))
    # rownames(res2print) <- paste0(
    #   sub('trt', '    Treatment effects on ', rownames(res2print)), ':')
    print(res2print, quote = FALSE, ...)
  }
}
################################################################################


################################################################################
plot.loocvSurrosurv <- function(x,
                                models,
                                exact.models, 
                                plot.type = c('classic', 'regression'),
                                main, ylab, xlab, ...) {
  # ************************************************************************** #
  x <- print(x, silent = TRUE)
  
  if (missing(models))
    models <- setdiff(names(x), 'margPars')
  
  if (missing(exact.models))
    exact.models <-
      any(tolower(noSpP(names(x))) %in% tolower(noSpP(models)))
  
  if (exact.models) {
    ind <- which(tolower(noSpP(names(x))) %in% tolower(noSpP(models)))
  } else {
    ind <- which(sapply(tolower(noSpP(names(
      x
    ))), function(mod)
      !all(is.na(
        pmatch(tolower(noSpP(models)), mod)
      ))))
  }
  
  plot.type <- match.arg(plot.type)
  
  mypalette <-
    c(
      orange = rgb(234, 104, 14, maxColorValue = 255, alpha = 180),
      blue = rgb(28, 170, 155, maxColorValue = 255, alpha = 180),
      magenta = rgb(191,  22, 120, maxColorValue = 255, alpha = 180),
      green = rgb(
        184 * .7,
        201,
        27 * .5,
        maxColorValue = 255,
        alpha = 180
      ),
      grey = rgb(79, 71, 68, maxColorValue = 255, alpha = 180)
    )
  palette(mypalette[c(2:4, 1, 5)])
  
  if (length(ind)) {
    par(mfrow = n2mfrow(length(ind)))
    if (!missing(main)) {
      if (length(main) == 1) mains <- rep(main, length(ind))
      else mains <- main
    } else mains <- format.methodNames(x)
    
    if (plot.type == 'classic') {
      XLAB <- ifelse(missing(xlab), 'Trials', xlab)
    } else {
      XLAB <- ifelse(missing(xlab), 'Hazard ratio on the surrogate endpoint', xlab)
    }
    
    for (i in ind) {
      if (plot.type == 'classic') {
        xlims <- c(0, ncol(x[[i]])) + .5
        Xs <- 1:ncol(x[[i]])
      } else {
        xlims <- range(x[[i]]['obsAlpha', ])
        Xs <- x[[i]]['obsAlpha', ]
      }
      plot(
        NA,
        main = mains[i],
        xlim = xlims,
        xlab = XLAB,
        xaxt = 'n',
        ylim = range(x[[i]], na.rm = TRUE) +
          c(-1, 1) * diff(range(x[[i]], na.rm = TRUE)) / 20,
        ylab = ifelse(missing(ylab), 'Hazard ratio on the true endpoint', ylab),
        yaxt = 'n',
        panel.first = abline(h = 0, col = 'grey')
      )
      if (plot.type == 'classic') {
        axis(1, 1:ncol(x[[i]]), labels = colnames(x[[i]]))
      } else axis(1, axTicks(1), format(round(exp(axTicks(1)), 2)), las = 1)
      axis(2, axTicks(2), format(round(exp(axTicks(2)), 2)), las = 1)
      segments(
        Xs,
        x[[i]]['lwr',],
        y1 = x[[i]]['upr',],
        col = rgb(.6, .6, .6, .8),
        lwd = 5
      )
      # segments(1:ncol(x[[i]]) - .01, x[[i]]['predict', ], 1:ncol(x[[i]]) + .01,
      #          col = rgb(.2, .2, .2), lwd = 3, lend = 2)
      points(
        Xs,
        x[[i]]['predict',],
        pch = 15,
        col = rgb(.2, .2, .2),
        lwd = 3,
        lend = 2
      )
      COLs <-
        2 - (colSums(apply(x[[i]][c('obsBeta', 'lwr', 'upr'),], 2,
                           function(x)
                             order(x) == c(2, 1, 3))) == 3)
      NC <- is.na(x[[i]]['lwr',]) | is.na(x[[i]]['upr',])
      COLs[NC] <- 0
      points(Xs,
             x[[i]]['obsBeta',],
             pch = 16,
             cex = 1.4,
             col = COLs)
      if (any(NC))
        mtext('x',
              1,
              -1,
              at = which(NC),
              col = 2,
              font = 2)
    }
  }
}
################################################################################