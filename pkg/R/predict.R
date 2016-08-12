predict.surrosurv <- function(x, ...) {
  copulas <- intersect(c('Clayton', 'Plackett', 'Hougaard'), names(x))
  poissons <- grep('PoissonT', names(x), value = TRUE)
  
  allRES <- c(
    # Copula models
    sapply(copulas, function(cop) {
      list(
        unadj = as.data.frame(x[[cop]]$unadj$step1[c('alpha', 'beta')]),
        adj = as.data.frame(x[[cop]]$adj$ranef))
    }),
    # Poisson models
    lapply(poissons, function(poi) {
      x[[poi]]$ranef$trialref[, c('trtS', 'trtT')]
    })
  )
  
  allRES <- lapply(allRES, function(x) {
    colnames(x) <- c('trtS', 'trtT')
    return(x)
  })
  names(allRES) <- c(
    paste(rep(copulas, each = 2), c('unadj', 'adj'), sep='.'),
    poissons)
  
  class(allRES)  <- c('predictSurrosurv', class(allRES))
  attr(allRES, 'trialSizes') <- attr(x, 'trialSizes')
  return(allRES)
}

format.methodNames <- function(x) {
  sub('.unadj', ' copula (Unadjusted)', sub(
    '.adj', ' copula (Adjusted)', sub(
      'Poisson', 'Poisson ', names(x)),
    fixed = TRUE), fixed = TRUE)
}

print.predictSurrosurv <- function(x, n = 6, ...) {
  cat('Treatment effect prediction for surrosurv object\n' )
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
  x, models = names(x), ...) {
  if (missing(xlab)) xlab <- 'Treatment effect (HR) on S'
  if (missing(ylab)) ylab <- 'Treatment effect (HR) on T'
  
  w <- attr(x, 'trialSizes')
  if (var(w)) {
    w <- .8 + 3 * (w - min(w)) / (max(w) - min(w))
  } else w <- 1
  
  x <- x[which(sapply(tolower(names(x)), function(mod)
    !all(is.na(pmatch(tolower(models), mod)))))]
  
  if (length(x)) {
    par(mfrow = n2mfrow(length(x)))
    for (i in 1:length(x)) {
      plot(x[[i]], asp = 1,
           cex = w, pch = 21, bg = rgb(.5, .5, .5, .5),
           panel.first = {
             abline(h=0, v=0, col='grey')
             abline(lm(x[[i]][, 2:1]), col=rgb(.2, .2, .2, .8), lwd = 2)
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
    }
  } else warning('No model selected', immediate. = TRUE)
}