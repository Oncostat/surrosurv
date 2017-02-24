#####################################################################################
#####################################################################################
#####################################################################################
convals <- function(x) {
  if (!'surrosurv' %in% class(x))
    stop('x must be of class surrosurv')
  
  models <- c('Clayton', 'Plackett', 'Hougaard', 'Poisson')
  models <- models[sapply(models, function(y) any(grepl(y, names(x))))]
  
  copulas <- function(cop) {
    f <- function() {
      res <- t(sapply(x[[cop]], function(y) {
        vals <- c(
          maxSgrad = tryCatch(
            max(abs(with(attributes(y$optimxRES), solve(hessian, grad)))),
            error = function(x) return(Inf)),
          minHev = tryCatch(
            min(eigen(attr(y$optimxRES, 'hessian'))$values),
            error = function(x) return(-1)),
          minREev = tryCatch(
            min(eigen(y$VarCor2)$values),
            error = function(x) return(-1))
        )
        return(vals)
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
      function(y) tryCatch({
        vals <- c(
          maxSgrad = tryCatch(
            max(abs(with(y$optinfo$derivs,solve(Hessian, gradient)))),
            error = function(x) return(Inf)),
          minHev = tryCatch(
            min(eigen(y$optinfo$derivs$Hessian)$values),
            error = function(x) return(-1)),
          minREev = tryCatch(
            min(eigen(y$VarCor$trialref)$values),
            error = function(x) return(-1))
        )
        return(vals)
      }, error = function(x) return(c(Inf, -1, -1))))
    )
    names(convres) <- names(which_models)
    return(convres)
  }
  
  fitRES <- do.call(rbind, lapply(models, function(x) eval(call(paste(x)))))
  fitRES[grepl('unadj|PoissonI', rownames(fitRES)), 'minREev'] <- NA
  class(fitRES) <- 'conv'
  return(fitRES)
}

#####################################################################################
#####################################################################################
#####################################################################################
convergence <- function(x, kkttol = 1e-2, kkt2tol = 1e-8) {
  if (!'surrosurv' %in% class(x))
    stop('x must be of class surrosurv')
  
  ConvRES <- convals(x)
  checkConvRES <- cbind(
    ConvRES[, 1, drop = FALSE] <= kkttol,
    ConvRES[, 2:3] > kkt2tol)
  
  class(checkConvRES) <- 'conv'
  return(checkConvRES)
}

#####################################################################################
#####################################################################################
#####################################################################################

print.conv <- function(x, na.print = '---', ...){
  class(x) <- 'matrix'
  print(x, na.print = na.print, ... )
}