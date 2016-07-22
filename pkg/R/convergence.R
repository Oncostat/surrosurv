#####################################################################################
#####################################################################################
#####################################################################################
convals <- function(x, restricted = FALSE) {
  if (!'surrosurv' %in% class(x))
    stop('x must be of class surrosurv')
  
  models <- c('Clayton', 'Plackett', 'Hougaard', 'Poisson')
  models <- models[sapply(models, function(y) any(grepl(y, names(x))))]
  
  copulas <- function(cop) {
    f <- function() {
      res <- t(sapply(x[[cop]], function(y) {
        vals <- try({
          if (restricted) {
            c(maxAgrad = abs(attr(y$optimxRES, 'grad')[1]),
              minev1 = tryCatch(
                min(eigen(solve(
                  attr(y$optimxRES, 'hessian')[1, 1, drop=FALSE])
                )$values) > 0,
              error = function(x) return(NA)))
          } else {
            c(maxAgrad = max(abs(attr(y$optimxRES, 'grad'))),
              minev1 = tryCatch(
                min(eigen(solve(attr(y$optimxRES, 'hessian')))$values),
                error = function(x) return(NA)))
          }
        }, silent = TRUE)
        if (class(vals) == 'try-error') vals <- c(maxAgrad = NA, minev1 = NA)
        return(rbind(c(vals,
                       minev2 = tryCatch(
                         min(eigen(y$VarCor2)$values),
                         error=function(x) return(NA)))))
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
    which_models <- x[grepl('Poisson', names(x))]
    convres <- t(sapply(
      which_models,
      function (y) tryCatch({
        vals <- try({
          c(maxAgrad = max(abs(y$optinfo$derivs$gradient)), 
            minev1 = tryCatch(
              min(eigen(solve(y$optinfo$derivs$Hessian))$values),
              error = function(x) return(NA)))
        })
        if (class(vals) == 'try-error') vals <- c(maxAgrad = NA, minev1 = NA)
        return(c(vals,
                 minev2 = tryCatch(min(eigen(y$VarCor$trialref)$values),
                                   error = function(x) return(NA))))
      }, error=function(x) return(rep(NA, 3))))
    )
    names(convres) <- names(which_models)
    if (restricted)
      convres[, 1:2] <- NA
    return(convres)
  }
  
  fitRES <- do.call(rbind, lapply(models, function(x) eval(call(paste(x)))))
  return(fitRES)
}

#####################################################################################
#####################################################################################
#####################################################################################
convergence <- function(x, kkttol = 1e-3, kkt2tol = 1e-6, restricted = FALSE) {
  if (!'surrosurv' %in% class(x))
    stop('x must be of class surrosurv')
    
  ConvRES <- convals(x, restricted)
  checkConvRES <- cbind(
    ConvRES[, 1] <= kkttol,
    ConvRES[, 2:3] >= kkt2tol)
  
  return(checkConvRES)
}