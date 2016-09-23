predict.surrosurv <- function(object, ...) {
  copulas <- intersect(c('Clayton', 'Plackett', 'Hougaard'), names(object))
  poissons <- grep('PoissonT', names(object), value = TRUE)
  
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
    paste(rep(copulas, each = 2), c('unadj', 'adj'), sep='.'),
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
          predict(object[[cop]]$unadj$step2,
                  data.frame(alpha = x),
                  interval = 'prediction')
        }),
        adj = Vectorize(function(x) {
          VCOV <- object[[cop]]$adj$step2$Psi + object[[cop]]$adj$step2$vcov
          CORR <- VCOV[1, 2] / sqrt(prod(diag(VCOV)))
          AVG <- object[[cop]]$adj$step2$coefficients
          cMEAN <- AVG[, 'beta'] + (x - AVG[, 'alpha']) * CORR * 
            # sqrt(VCOV['beta', 'beta'] / VCOV['alpha', 'alpha'])
            sqrt(exp(diff(log(diag(VCOV)))))
          return(rep(cMEAN, 3) + qnorm(c(fit=.5, lwr=.025, upr=.975)) * 
                   sqrt(VCOV['beta', 'beta'] * (
                     1 - VCOV['alpha', 'beta'] / prod(sqrt(diag(VCOV))))))
        }))
    }),
    # Poisson models
    lapply(poissons, function(poi) {
      Vectorize(function(x) {
        VCOV <- as.matrix(
          object[[poi]]$VarCor$trialref[paste0('trt', c('S', 'T')),
                                        paste0('trt', c('S', 'T'))] +
            object[[poi]]$fixVarCor[paste0('trt', c('S', 'T')),
                                    paste0('trt', c('S', 'T'))])
        CORR <- VCOV[1, 2] / sqrt(prod(diag(VCOV)))
        AVG <- unlist(object[[poi]][c('alpha', 'beta')])
        cMEAN <- AVG['beta.trtT'] + (x - AVG['alpha.trtS']) * CORR *
          # sqrt(VCOV['trtT', 'trtT'] / VCOV['trtS', 'trtS'])
          sqrt(exp(diff(log(diag(VCOV)))))
        res <- rep(cMEAN, 3) + qnorm(c(fit=.5, lwr=.025, upr=.975)) * 
          sqrt(VCOV['trtS', 'trtS'] * (1 - VCOV[1, 2] / prod(sqrt(diag(VCOV)))))
        names(res) <- c('fit', 'lwr', 'upr')
        return(cbind(res))
      })
    })
  )
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

plot.surrosurv <- function(x, ...)
  plot(predict(x), ...)


plot.predictSurrosurv <- function(
  x, 
  models = names(x), 
  exact.models,
  pred.ints = TRUE,
  xlab, 
  ylab, ...) {
  if (missing(xlab)) xlab <- 'Treatment effect (HR) on S'
  if (missing(ylab)) ylab <- 'Treatment effect (HR) on T'
  if (missing(exact.models))
    exact.models <- any(tolower(noSpP(names(x))) %in% tolower(noSpP(models)))
  
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
  x <- x[ind]
  
  if (length(x)) {
    par(mfrow = n2mfrow(length(x)))
    for (i in 1:length(x)) {
      # abcoeff <- try(coef(lm(x[[i]][, 2:1])), silent = TRUE)
      predf.fit <- Vectorize(function(y) {PREDF[[i]](y)[1]})
      set0in <- function(xint) {
        if (xint[1] > 0) xint[1] <- 0
        if (xint[2] < 0) xint[2] <- 0
        return(xint)
      }
      xlims <- set0in(range(x[[i]][, 1]))
      ylims <- set0in(range(x[[i]][, 2]))
      plot(x[[i]], asp = 1, xlim = xlims, ylim = ylims,
           cex = w, pch = 21, bg = rgb(.5, .5, .5, .5),
           panel.first = {
             if (pred.ints) {
               Xs <- seq(axTicks(1)[1] - 1, rev(axTicks(1))[1] + 1, length.out = 1e2)
               polygon(c(Xs, rev(Xs)),
                       c(PREDF[[i]](Xs)[2, ], rev(PREDF[[i]](Xs)[3, ])),
                       col = rgb(.8, .8, .8, .5), border = NA)
             }
             abline(h=0, v=0, col='grey')
             # if (all(is.finite(abcoeff)))
             # abline(lm(x[[i]][, 2:1]), col=rgb(.2, .2, .2, .8), lwd = 2)
             curve(predf.fit, from = axTicks(1) - 1,
                   to = rev(axTicks(1))[1] + 1,
                   col=rgb(.2, .2, .2, .8), lwd = 2, add=TRUE)
           },
           xaxt = 'n', yaxt = 'n',
           main =  format.methodNames(x)[i], 
           xlab = xlab, ylab = ylab, ...)
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