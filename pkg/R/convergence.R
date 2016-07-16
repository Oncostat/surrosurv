#####################################################################################
#####################################################################################
#####################################################################################
convergence <- function(x, kkttol = 1e-3, kkt2tol = 1e-6, restricted = FALSE) {
  if (!'surrosurv' %in% class(x))
    stop('x must be of class surrosurv')
  
  models <- c('Clayton', 'Plackett', 'Hougaard', 'Poisson')
  models <- models[sapply(models, function(y) any(grepl(y, names(x))))]
  
  kktCtrl <- list(kkttolA = kkttol,
                  kkt2tolA = kkt2tol)
  
  copulas <- function(cop) {
    f <- function() {
      res <- t(sapply(x[[cop]], function(y) {
        # library('optextras')
        kkt <- try({
          if (restricted) {
            kktc(
              par = as.numeric(y$optimxRES[, 1]),
              fval = y$optimxRES$value,
              ngr = attr(y$optimxRES, 'grad')[1], 
              nHes = attr(y$optimxRES, 'hessian')[1, 1, drop=FALSE], 
              nbm=0, control = kktCtrl)
          } else {
            kktc(
              par = as.numeric(y$optimxRES[, 1:attr(y$optimxRES, 'npar')]),
              fval = y$optimxRES$value,
              ngr = attr(y$optimxRES, 'grad'), 
              nHes = attr(y$optimxRES, 'hessian'), 
              nbm=0, control = kktCtrl)
          }
        }, silent = TRUE)
        if (class(kkt) == 'try-error') kkt <- c(kkt1 = NA, kkt2 = NA)
        return(rbind(c(
          unlist(kkt[c('kkt1', 'kkt2')]),
          VCpd = tryCatch({
            # library('corpcor')
            is.positive.definite(y$VarCor2, tol = kkttol)
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
                   is.positive.definite(y$VarCor$trialref, tol = kkttol),
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