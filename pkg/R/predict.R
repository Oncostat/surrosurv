predict.surrosurv <- function(object,
                              models = names(object),
                              exact.models,
                              ...) {
  if (missing(exact.models))
    exact.models <- all(sapply(tolower(noSpP(models)), function(m)
      any(tolower(noSpP(names(object))) %in% m)))
  
  if (exact.models) {
    ind <- which(tolower(noSpP(names(object))) %in% tolower(noSpP(models)))
  } else {
    ind <- which(sapply(tolower(noSpP(names(object))), function(mod)
      !all(is.na(pmatch(tolower(noSpP(models)), mod)))))
  }
  models <- names(object)[ind]
  
  copulas <- intersect(c('Clayton', 'Plackett', 'Hougaard'), models)
  poissons <- grep('PoissonT', models, value = TRUE)
  
  allRES <- c(
    # Copula models
    sapply(copulas, function(cop) {
      list(
        unadj = as.data.frame(object[[cop]]$unadj$step1[c('alpha', 'beta')]),
        adj = as.data.frame(object[[cop]]$adj$ranef))
    }),
    # Poisson models
    lapply(poissons, function(poi) {
      # object[[poi]]$ranef$trialref[, c('trtS', 'trtT')]
      as.data.frame(t(
        t(object[[poi]]$ranef$trialref[, c('trtS', 'trtT')]) + 
          unlist(object[[poi]][c('alpha', 'beta')])))
    })
  )
  
  allRES <- lapply(allRES, function(object) {
    colnames(object) <- c('trtS', 'trtT')
    return(object)
  })
  names(allRES) <- c(
    paste(rep(copulas, each = 2), rep(c('unadj', 'adj'), length(copulas)), sep='.'),
    poissons)
  
  class(allRES)  <- c('predictSurrosurv', class(allRES))
  attr(allRES, 'trialSizes') <- attr(object, 'trialSizes')
  
  attr(allRES, 'predf') <- c(
    # Copula models
    sapply(copulas, function(cop) {
      list(
        unadj = Vectorize(function(x) {
          # est <- coef(object[[cop]]$unadj$step2) %*% c(1, x)
          # n <- 2 + object[[cop]]$unadj$step2$df.residual
          # t <- qt(c(.025, .975), n - 2)
          # S <- summary(object[[cop]]$unadj$step2)$sigma
          # sqrt(
          #     1 + 1 / n + n * (
          #       x - mean(object[[cop]]$unadj$step1$alpha)
          #     )^2 / (
          #       n * sum(object[[cop]]$unadj$step1$alpha^2) -
          #         sum(object[[cop]]$unadj$step1$alpha)^2
          #     )
          # ) * S * t + est
          res <- as.vector(predict(object[[cop]]$unadj$step2,
                                   data.frame(alpha = x),
                                   interval = 'prediction'))
          names(res) <- c('fit', 'lwr', 'upr')
          return(res)
        }),
        adj = Vectorize(function(x) {
          cMEAN <- object[[cop]]$adj$step2$gamma %*% c(1, x) 
          # The following variance is given by equation (18.27), page 332
          # Burzykowski and Buyse, An alternative measure for surrogate 
          # endpoint validation, in Burzykowski, Molenberghs, Buyse, 2005
          # The evaluation of surrogate endpoints, New York Springer
          # library(Matrix)
          pdVCOV <- nearPD(object[[cop]]$adj$step2$gamma.vcov)$mat
          pVAR <- pdVCOV[1, 1] + 2 * x * pdVCOV[1, 2] + x^2 * pdVCOV[2, 2] +
            object[[cop]]$adj$step2$sigma
          res <- rep(cMEAN, 3) +
            qnorm(c(fit = .5, lwr = .025, upr = .975)) * sqrt(pVAR)
          names(res) <- c('fit', 'lwr', 'upr')
          return(res)
        }))
    }),
    # Poisson models
    lapply(poissons, function(poi) {
      Vectorize(function(x) {
        # VCOV <- as.matrix(
        #   object[[poi]]$VarCor$trialref[paste0('trt', c('S', 'T')),
        #                                 paste0('trt', c('S', 'T'))] +
        #     object[[poi]]$fixVarCor[paste0('trt', c('S', 'T')),
        #                             paste0('trt', c('S', 'T'))])
        # CORR <- VCOV[1, 2] / sqrt(prod(diag(VCOV)))
        # AVG <- unlist(object[[poi]][c('alpha', 'beta')])
        # cMEAN <- AVG['beta.trtT'] + (x - AVG['alpha.trtS']) * CORR *
        #   # sqrt(VCOV['trtT', 'trtT'] / VCOV['trtS', 'trtS'])
        #   sqrt(exp(diff(log(diag(VCOV)))))
        # res <- as.vector(
        #   rep(cMEAN, 3) + qnorm(c(fit=.5, lwr=.025, upr=.975)) * 
        #     sqrt(VCOV['trtS', 'trtS'] * (1 - VCOV[1, 2] / 
        #                                    prod(sqrt(diag(VCOV))))))
        cMEAN <- object[[poi]]$regr$gamma %*% c(1, x) 
        # The following varianve is given by equation (18.27), page 332
        # Burzykowski and Buyse, An alternative measure for surrogate 
        # endpoint validation, in Burzykowski, Molenberghs, Buyse, 2005
        # The evaluation of surrogate endpoints, New York Springer
        # library(Matrix)
        pdVCOV <- nearPD(object[[poi]]$regr$gamma.vcov)$mat
        pVAR <- pdVCOV[1, 1] + 2 * x *  pdVCOV[1, 2] + x^2 * pdVCOV[2, 2] +
          object[[poi]]$regr$sigma
        res <- rep(cMEAN, 3) +
          qnorm(c(fit = .5, lwr = .025, upr = .975)) * sqrt(pVAR)
        names(res) <- c('fit', 'lwr', 'upr')
        return(res)
      })
    })
  )
  attr(allRES, 'surro.stats') <- print(object, silent = TRUE)
  attr(allRES, 'surro.stats') <- attr(allRES, 'surro.stats')[
    noSpP(rownames(attr(allRES, 'surro.stats'))) %in%
      noSpP(names(allRES)), , drop=FALSE]
  noSpP(names(allRES))
  names(attr(allRES, 'predf')) <- names(allRES)
  return(allRES)
}

format.methodNames <- function(x) {
  sub('.unadj', ' copula (Unadjusted)', sub(
    '.adj', ' copula (Adjusted)', sub(
      'Poisson', 'Poisson ', names(x)),
    fixed = TRUE), fixed = TRUE)
}

noSpP <- function(x) gsub('[\\. ]', '', x)

print.predictSurrosurv <- function(x, n = 6, ...) {
  cat('Treatment effect prediction for surrosurv object\n')
  for (i in 1:length(x)) {
    method <- format.methodNames(x)[i]
    cat('\n  ', method, '\n')
    res2print <- format(t(head(x[[i]], n)), digits=1)
    if (nrow(x[[i]] > n))
      res2print <- cbind(res2print, '  ' = c('...', '...'))
    rownames(res2print) <- paste0(
      sub('trt', '    Treatment effects on ', rownames(res2print)), ':')
    print(res2print, quote = FALSE, ...)
  }
}

##############################################################################

ste <- function(x, models = names(x), exact.models) {
  if (class(x)[1] == 'surrosurv') x <- predict(x)
  if (missing(exact.models))
    exact.models <- all(sapply(tolower(noSpP(models)), function(m)
      any(tolower(noSpP(names(x))) %in% m)))
  
  if (exact.models) {
    ind <- which(tolower(noSpP(names(x))) %in% tolower(noSpP(models)))
  } else {
    ind <- which(sapply(tolower(noSpP(names(x))), function(mod)
      !all(is.na(pmatch(tolower(noSpP(models)), mod)))))
  }
  
  res <- sapply(ind, function(i) {
    f <- function(y) attr(x, 'predf')[[i]](y)[3,]^2
    ste <- optimize(f, c(-1e8, 1e8))$minimum
    return(ste)
  })
  names(res) <- names(x)[ind]
  class(res) <- 'steSurrosurv'
  return(res)
}

print.steSurrosurv <- function(x, digits = 2, ...) {
  res2print <- cbind(
    beta = round(x, digits),
    HR = round(exp(x), digits))
  print(res2print, quote = FALSE, ...)
}


##############################################################################

plot.surrosurv <- function(x, ...) {
  plot(predict(x), ...)
}


plot.predictSurrosurv <- function(
  x, 
  models = names(x), 
  exact.models,
  pred.ints = TRUE,
  show.ste = TRUE,
  surro.stats = TRUE,
  xlab, 
  ylab,
  xlim,
  ylim,
  mfrow,
  main,
  ...) {
  # ************************************************************************** #
  if (missing(xlab)) xlab <- 'Treatment effect (HR) on S'
  if (missing(ylab)) ylab <- 'Treatment effect (HR) on T'
  if (missing(exact.models))
    exact.models <- all(sapply(tolower(noSpP(models)), function(m)
      any(tolower(noSpP(names(x))) %in% m)))
  w <- attr(x, 'trialSizes')
  if (var(w)) {
    w <- .8 + 3 * (w - min(w)) / (max(w) - min(w))
  } else w <- 1
  
  if (exact.models) {
    ind <- which(tolower(noSpP(names(x))) %in% tolower(noSpP(models)))
  } else {
    ind <- which(sapply(tolower(noSpP(names(x))), function(mod)
      !all(is.na(pmatch(tolower(noSpP(models)), mod)))))
  }
  PREDF <- attr(x, 'predf')[ind]
  if (show.ste)
    STE <- ste(x, names(x)[ind], exact.models = TRUE)
  if (surro.stats){
    SURRO.STATS <- round(matrix(as.numeric(attr(x, 'surro.stats')), ncol = 2), 
                         2)[ind, , drop = FALSE]
    SURRO.STATS[is.na(SURRO.STATS)] <- '\u2014'
  }
  x <- x[ind]
  
  if(!missing('xlim')) xlims <- log(xlim)
  if(!missing('ylim')) ylims <- log(ylim)
  
  if (length(x)) {
    if (missing(mfrow)) {
      par(mfrow = n2mfrow(length(x)))
    } else {
      par(mfrow = mfrow)
    }
    for (i in 1:length(x)) {
      # abcoeff <- try(coef(lm(x[[i]][, 2:1])), silent = TRUE)
      predf.fit <- Vectorize(function(y) {PREDF[[i]](y)[1]})
      set0in <- function(xint) {
        if (xint[1] > 0) xint[1] <- 0
        if (xint[2] < 0) xint[2] <- 0
        five <- diff(range(xint)) / 20
        return(xint + c(-1, 1) * five)
      }
      if(missing('xlim'))
        xlims <- set0in(range(x[[i]][, 1]))
      if(missing('ylim'))
        ylims <- set0in(range(x[[i]][, 2]))
      plot(x[[i]], asp = 1, xlim = xlims, ylim = ylims,
           cex = w, pch = 21, bg = rgb(.5, .5, .5, .5),
           panel.first = {
             if (pred.ints) {
               Xs <- seq(axTicks(1)[1] - 1, rev(axTicks(1))[1] + 1, 
                         length.out = 1e2)
               polygon(c(Xs, rev(Xs)),
                       c(PREDF[[i]](Xs)[2, ], rev(PREDF[[i]](Xs)[3, ])),
                       col = rgb(.8, .8, .8, .75), border = NA)
             }
             abline(h=0, v=0, col='grey')
             # if (all(is.finite(abcoeff)))
             # abline(lm(x[[i]][, 2:1]), col=rgb(.2, .2, .2, .8), lwd = 2)
             curve(predf.fit, from = axTicks(1) - 1,
                   to = rev(axTicks(1))[1] + 1,
                   col=rgb(.2, .2, .2, .8), lwd = 2, add=TRUE)
           },
           xaxt = 'n', yaxt = 'n',
           main =  ifelse(missing(main), format.methodNames(x)[i], main),
           xlab = xlab, ylab = ylab, ...)
      if (show.ste) {
        points(STE[i], 0, col = 2, pch = '|', font = 2)
        text(STE[i], 0, col = 2, font = 2, adj = c(1.1, -.5),
             labels = paste('STE =', round(exp(STE[i]), 2)))
      }
      if (surro.stats) {
        mtext(bquote(paste(R^2, ' = ', .(SURRO.STATS[i, 2]))),
              side = 1, line = -1.2, adj = .9, font = 2)
        mtext(bquote(paste(tau, ' = ', .(SURRO.STATS[i, 1]))),
              side = 1, line = -2.8, adj = .9, font = 2)
      }
      axis(1, at = axTicks(1), 
           labels = sapply(
             exp(axTicks(1)), function(x)
               ifelse(x < .1, format(x, nsmall=3, digits=2),
                      ifelse(x < 1, format(x, nsmall=2, digits=2),
                             ifelse(x < 10, format(x, nsmall=1, digits=2),
                                    format(x, digits=2, drop0trailing=FALSE)
                             )
                      )
               )
           )
      )
      axis(2, at = axTicks(2), 
           labels = sapply(
             exp(axTicks(2)), function(x)
               ifelse(x < .1, format(x, nsmall=3, digits=2),
                      ifelse(x < 1, format(x, nsmall=2, digits=2),
                             ifelse(x < 10, format(x, nsmall=1, digits=2),
                                    format(x, digits=2, drop0trailing=FALSE)
                             )
                      )
               )
           )
      )
    }
  } else warning('No model selected', immediate. = TRUE)
}