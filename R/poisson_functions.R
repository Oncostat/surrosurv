poisSurr <- function(data,
                     OPTIMIZER,
                     # = 'bobyqa',
                     MAXFUN = 1e8,
                     intWidth = NULL,
                     nInts = NULL,
                     models = paste('Poisson', c('I', 'T', 'TI', 'TIa')),
                     verbose = TRUE) {
  models <- tolower(noSpP(models))
  models <- match.arg(models,
                      several.ok = TRUE,
                      choices = paste0('poisson', c('i', 't', 'ti', 'tia')))
  #     library('parfm')
  #     library('survival')
  #     library('msm')
  #     library('lme4')
  
  GLMERCONTROL <- glmerControl(
    optimizer = 'optimx',
    # calc.derivs = FALSE,
    boundary.tol = 0,
    optCtrl = list(
      maxfun = MAXFUN,
      dowarn = FALSE,
      method = OPTIMIZER,
      starttests = FALSE,
      kkt = FALSE
    )
  )
  
  ############################################################################
  ### *** Data poissonization *** ############################################
  ############################################################################
  if (verbose) {
    message('- Data poissonization', appendLF = FALSE)
    startTime <- Sys.time()
  }
  
  T <- S <- JOINT <- list()
  
  if ((min(data$trt) == 0) & (max(data$trt) == 1)) {
    message('Recoding treatment from 0/1 to -0.5/0.5')
    data$trt <- data$trt - .5
  }
  
  if (is.null(intWidth)) {
    times <- unlist(data[, grep('time', names(data))])
    stata <- unlist(data[, grep('status', names(data))])
    # library('survival')
    smod <- survfit(Surv(times, stata) ~ 1)
    smin <- min(smod$surv)
    all.breaks <- c(0, sapply(1 - 1:nInts / nInts, function(p)
      smod$time[(smod$surv - smin) / (1 - smin) <= p][1]))
    rm(times, stata, smod, smin)
  } else
    all.breaks <- NULL
  
  colnames(data)[colnames(data) == "timeS"] <- 'time'
  colnames(data)[colnames(data) == "statusS"] <- 'status'
  S$poidata <- poissonize(
    data,
    all.breaks = all.breaks,
    interval.width = intWidth,
    nInts = nInts,
    factors = c('trialref', 'trt', 'id'),
    compress = TRUE
  )
  colnames(data)[colnames(data) == "time"] <- 'timeS'
  colnames(data)[colnames(data) == "status"] <- 'statusS'
  
  colnames(data)[colnames(data) == "timeT"] <- 'time'
  colnames(data)[colnames(data) == "statusT"] <- 'status'
  T$poidata <- poissonize(
    data,
    all.breaks = all.breaks,
    interval.width = intWidth,
    nInts = nInts,
    factors = c('trialref', 'trt', 'id'),
    compress = TRUE
  )
  colnames(data)[colnames(data) == "time"] <- 'timeT'
  colnames(data)[colnames(data) == "status"] <- 'statusT'
  
  if (verbose) {
    execTime <- Sys.time() - startTime
    message(paste0(' (', format(execTime, digits = 2), ')'))
  }
  ############################################################################
  
  
  ############################################################################
  ### *** Mixed Poisson model estimation *** #################################
  ############################################################################
  JOINT$poidata <- rbind(cbind(S$poidata, ep = 'S'),
                         cbind(T$poidata, ep = 'T'))
  
  JOINT$poidata$trtT <-
    JOINT$poidata$trt * (JOINT$poidata$ep == "T")
  JOINT$poidata$trtS <-
    JOINT$poidata$trt * (JOINT$poidata$ep == "S")
  
  JOINT$poidata$trialrefT <- JOINT$poidata$trialrefS <-
    JOINT$poidata$trialref
  JOINT$poidata$trialrefT[JOINT$poidata$ep != "T"]  <-
    levels(JOINT$poidata$trialref)[1]
  JOINT$poidata$trialrefS[JOINT$poidata$ep != "S"]  <-
    levels(JOINT$poidata$trialref)[1]
  
  # W <- options()$warn
  # options(warn = -1)
  
  if (nlevels(JOINT$poidata$interval) > 1) {
    baseform <- m ~ -1 + interval * ep - ep + trtT + trtS + offset(log(Rt))
  } else  {
    baseform <- m ~ -1 + ep + trtT + trtS + offset(log(Rt))
  }
  
  RES <- list()
  
  ##############################################################################
  ##############################################################################
  #     library('optimx')
  # * Model T: Poisson model with
  #             - random treatment-trial interaction
  if ('poissont' %in% models) {
    if (verbose)
      message('- Estimating model: Poisson T', appendLF = FALSE)
    
    system.time({
      JOINT$poifitT <- glmer(
        update(baseform, . ~ . +
                 (-1 + trtT + trtS | trialref)),
        data = JOINT$poidata,
        family = poisson,
        control = GLMERCONTROL
      )
    }) -> attr(JOINT$poifitT, 'exec.time')
    R2 <- VarCorr(JOINT$poifitT)$trialref[1, 2] ^ 2 /
      prod(diag(VarCorr(JOINT$poifitT)$trialref))
    
    # ------------------------------------------------------------------------ #
    # Regression parameters ####
    PARS <- c(
      alpha = fixef(JOINT$poifitT)[['trtS']],
      beta  = fixef(JOINT$poifitT)[['trtT']],
      Psi.aa = VarCorr(JOINT$poifitT)$trialref[['trtS', 'trtS']],
      Psi.ab = VarCorr(JOINT$poifitT)$trialref[['trtS', 'trtT']],
      Psi.bb = VarCorr(JOINT$poifitT)$trialref[['trtT', 'trtT']]
    )
    # library(msm)
    Psi.vcov <- deltamethod(
      g = list( ~ x1 ^ 2, ~ x1 * x2, ~ x2 ^ 2 + x3 ^ 2),
      # for this tranformation, see the Remark on page 612 of
      # van Houwelingen, Arends, Stijnen, Stat Med 2002; 21:589-624
      mean = JOINT$poifitT@theta,
      cov = solve(-JOINT$poifitT@optinfo$derivs$Hessian[1:3, 1:3]),
      ses = FALSE
    )
    VCOV <- diag(5)
    VCOV[1:2, 1:2] <-
      vcov(JOINT$poifitT)[c('trtS', 'trtT'), c('trtS', 'trtT')]@x
    VCOV[3:5, 3:5] <- Psi.vcov
    rownames(VCOV) <- colnames(VCOV) <- names(PARS)
    
    # The following parameter estimates are given by equations (18.23), page 331
    # Burzykowski and Buyse, An alternative measure for surrogate
    # endpoint validation, in Burzykowski, Molenberghs, Buyse, 2005
    # The evaluation of surrogate endpoints, New York Springer
    regr <- list(
      gamma = c(
        gamma0 = as.numeric(PARS['beta'] - PARS['alpha'] *
                              PARS['Psi.ab'] / PARS['Psi.aa']),
        gamma1 = as.numeric(PARS['Psi.ab'] / PARS['Psi.aa'])
      ),
      gamma.vcov = deltamethod(
        g = list( ~ x2 - x1 * x4 / x3, ~ x4 / x3),
        mean = PARS,
        cov = VCOV,
        ses = FALSE
      )
    )
    regr$sigma <- PARS['Psi.bb'] * (1 - R2)
    # ------------------------------------------------------------------------ #
    
    RES$PoissonT <- list(
      kTau = NULL,
      alpha = fixef(JOINT$poifitT)['trtS'],
      beta = fixef(JOINT$poifitT)['trtT'],
      fixVarCor = vcov(JOINT$poifitT)[c('trtS', 'trtT'), c('trtS', 'trtT')],
      R2 = R2,
      ranef = ranef(JOINT$poifitT),
      VarCor = VarCorr(JOINT$poifitT),
      optinfo = JOINT$poifitT@optinfo,
      regr = regr,
      runTime = as.difftime(round(
        attr(JOINT$poifitT, 'exec.time')[3] / 60, 1
      ), units = 'mins')
    )
    
    if (verbose)
      message(paste0(' (', format(RES$PoissonT$runTime, digits = 2), ')'))
  }#############################################################################
  
  ##############################################################################
  # * Model I: Poisson model with
  #             - individual random intercept
  if ('poissoni' %in% models) {
    if (verbose)
      message('- Estimating model: Poisson I', appendLF = FALSE)
    
    system.time({
      JOINT$poifitI <-  glmer(
        update(baseform, . ~ . +
                 (1 | id)) ,
        data = JOINT$poidata,
        family = poisson,
        control = GLMERCONTROL
      )
    }) -> attr(JOINT$poifitI, 'exec.time')
    
    attr(JOINT$poifitI, 'kTau') <- fr.lognormal(what = 'tau',
                                                        sigma2 = as.double(summary(JOINT$poifitI)$varcor$id))
    
    RES$PoissonI <- list(
      kTau = attr(JOINT$poifitI, 'kTau'),
      alpha = fixef(JOINT$poifitI)['trtS'],
      beta = fixef(JOINT$poifitI)['trtT'],
      fixVarCor = vcov(JOINT$poifitI)[c('trtS', 'trtT'), c('trtS', 'trtT')],
      R2 = NULL,
      ranef = ranef(JOINT$poifitI),
      VarCor = VarCorr(JOINT$poifitI),
      optinfo = JOINT$poifitI@optinfo,
      runTime = as.difftime(round(
        attr(JOINT$poifitI, 'exec.time')[3] / 60, 1
      ), units = 'mins')
    )
    
    if (verbose)
      message(paste0(' (', format(RES$PoissonI$runTime, digits = 2), ')'))
  }#############################################################################
  
  # marpars <- surrosurv:::margFits(data)
  # vcovmar <- matrix(c(var(marpars$alpha), cov(marpars$alpha, marpars$beta),
  #                     cov(marpars$alpha, marpars$beta), var(marpars$beta)), 2)
  # theta <- as.numeric((chol(vcovmar)))[-2] /
  #   var(glm(baseform, data = JOINT$poidata, family = poisson)$residual)
             
  
  ##############################################################################
  # * Model TI: Poisson model with
  #            - random treatment-trial interaction
  #            - individual random intercept
  if ('poissonti' %in% models) {
    if (verbose)
      message('- Estimating model: Poisson TI', appendLF = FALSE)
    
    system.time({
      JOINT$poifitTI <- glmer(
        update(baseform, . ~ . +
                 (-1 + trtT + trtS | trialref) +
                 (1 | id)),
        data = JOINT$poidata,
        family = poisson,
        control = GLMERCONTROL
        # , start = theta
      )
    }) -> attr(JOINT$poifitTI, 'exec.time')
    
    
    attr(JOINT$poifitTI, 'kTau') <- fr.lognormal(what = 'tau',
                                                         sigma2 = as.double(summary(JOINT$poifitTI)$varcor$id))
    R2 <- VarCorr(JOINT$poifitTI)$trialref[1, 2] ^ 2 /
      prod(diag(VarCorr(JOINT$poifitTI)$trialref))
    
    # ------------------------------------------------------------------------ #
    # Regression parameters ####
    PARS <- c(
      alpha = fixef(JOINT$poifitTI)[['trtS']],
      beta  = fixef(JOINT$poifitTI)[['trtT']],
      Psi.aa = VarCorr(JOINT$poifitTI)$trialref[['trtS', 'trtS']],
      Psi.ab = VarCorr(JOINT$poifitTI)$trialref[['trtS', 'trtT']],
      Psi.bb = VarCorr(JOINT$poifitTI)$trialref[['trtT', 'trtT']]
    )
    # library(msm)
    Psi.vcov <- deltamethod(
      g = list( ~ x1 ^ 2, ~ x1 * x2, ~ x2 ^ 2 + x3 ^ 2),
      # for this tranformation, see the Remark on page 612 of
      # van Houwelingen, Arends, Stijnen, Stat Med 2002; 21:589-624
      mean = JOINT$poifitTI@theta[2:4],
      cov = solve(-JOINT$poifitTI@optinfo$derivs$Hessian[2:4, 2:4]),
      ses = FALSE
    )
    VCOV <- diag(5)
    VCOV[1:2, 1:2] <-
      vcov(JOINT$poifitTI)[c('trtS', 'trtT'), c('trtS', 'trtT')]@x
    VCOV[3:5, 3:5] <- Psi.vcov
    rownames(VCOV) <- colnames(VCOV) <- names(PARS)
    
    # The following parameter estimates are given by equations (18.23), page 331
    # Burzykowski and Buyse, An alternative measure for surrogate
    # endpoint validation, in Burzykowski, Molenberghs, Buyse, 2005
    # The evaluation of surrogate endpoints, New York Springer
    regr <- list(
      gamma = c(
        gamma0 = as.numeric(PARS['beta'] - PARS['alpha'] *
                              PARS['Psi.ab'] / PARS['Psi.aa']),
        gamma1 = as.numeric(PARS['Psi.ab'] / PARS['Psi.aa'])
      ),
      gamma.vcov = deltamethod(
        g = list( ~ x2 - x1 * x4 / x3, ~ x4 / x3),
        mean = PARS,
        cov = VCOV,
        ses = FALSE
      )
    )
    regr$sigma <- PARS['Psi.bb'] * (1 - R2)
    # ------------------------------------------------------------------------ #
    
    RES$PoissonTI <- list(
      kTau = attr(JOINT$poifitTI, 'kTau'),
      alpha = fixef(JOINT$poifitTI)['trtS'],
      beta = fixef(JOINT$poifitTI)['trtT'],
      fixVarCor = vcov(JOINT$poifitTI)[c('trtS', 'trtT'), c('trtS', 'trtT')],
      R2 = R2,
      ranef = ranef(JOINT$poifitTI),
      VarCor = VarCorr(JOINT$poifitTI),
      optinfo = JOINT$poifitTI@optinfo,
      regr = regr,
      runTime = as.difftime(round(
        attr(JOINT$poifitTI, 'exec.time')[3] / 60, 1
      ), units = 'mins')
    )
    
    if (verbose)
      message(paste0(' (', format(RES$PoissonTI$runTime, digits = 2), ')'))
  }#############################################################################
  
  ##############################################################################
  # * Model TIa: Poisson model with
  #             - random treatment-trial interaction
  #             - individual random intercept
  #             - random trial intercept (shared by the two EndPoints)
  if ('poissontia' %in% models) {
    if (verbose)
      message('- Estimating model: Poisson TIa', appendLF = FALSE)
    
    system.time({
      JOINT$poidata$trialrefSH <- JOINT$poidata$trialref
      JOINT$poifitTIa <-  glmer(
        update(
          baseform,
          . ~ . +
            (-1 + trtT + trtS | trialref) +
            (1 | id) +
            (1 | trialrefSH)
        ),
        data = JOINT$poidata,
        family = poisson,
        control = GLMERCONTROL
      )
    }) -> attr(JOINT$poifitTIa, 'exec.time')
    
    attr(JOINT$poifitTIa, 'kTau') <- fr.lognormal(what = 'tau',
                                                          sigma2 = as.double(summary(JOINT$poifitTIa)$varcor$id))
    R2 <- VarCorr(JOINT$poifitTIa)$trialref[1, 2] ^ 2 /
      prod(diag(VarCorr(JOINT$poifitTIa)$trialref))
    
    # ------------------------------------------------------------------------ #
    # Regression parameters ####
    PARS <- c(
      alpha = fixef(JOINT$poifitTIa)[['trtS']],
      beta  = fixef(JOINT$poifitTIa)[['trtT']],
      Psi.aa = VarCorr(JOINT$poifitTIa)$trialref[['trtS', 'trtS']],
      Psi.ab = VarCorr(JOINT$poifitTIa)$trialref[['trtS', 'trtT']],
      Psi.bb = VarCorr(JOINT$poifitTIa)$trialref[['trtT', 'trtT']]
    )
    # library(msm)
    Psi.vcov <- deltamethod(
      g = list( ~ x1 ^ 2, ~ x1 * x2, ~ x2 ^ 2 + x3 ^ 2),
      # for this tranformation, see the Remark on page 612 of
      # van Houwelingen, Arends, Stijnen, Stat Med 2002; 21:589-624
      mean = JOINT$poifitTIa@theta[3:5],
      cov = solve(-JOINT$poifitTIa@optinfo$derivs$Hessian[3:5, 3:5]),
      ses = FALSE
    )
    VCOV <- diag(5)
    VCOV[1:2, 1:2] <-
      vcov(JOINT$poifitTIa)[c('trtS', 'trtT'), c('trtS', 'trtT')]@x
    VCOV[3:5, 3:5] <- Psi.vcov
    rownames(VCOV) <- colnames(VCOV) <- names(PARS)
    
    # The following parameter estimates are given by equations (18.23), page 331
    # Burzykowski and Buyse, An alternative measure for surrogate
    # endpoint validation, in Burzykowski, Molenberghs, Buyse, 2005
    # The evaluation of surrogate endpoints, New York Springer
    regr <- list(
      gamma = c(
        gamma0 = as.numeric(PARS['beta'] - PARS['alpha'] *
                              PARS['Psi.ab'] / PARS['Psi.aa']),
        gamma1 = as.numeric(PARS['Psi.ab'] / PARS['Psi.aa'])
      ),
      gamma.vcov = deltamethod(
        g = list( ~ x2 - x1 * x4 / x3, ~ x4 / x3),
        mean = PARS,
        cov = VCOV,
        ses = FALSE
      )
    )
    regr$sigma <- PARS['Psi.bb'] * (1 - R2)
    # ------------------------------------------------------------------------ #
    
    RES$PoissonTIa <- list(
      kTau = attr(JOINT$poifitTIa, 'kTau'),
      alpha = fixef(JOINT$poifitTIa)['trtS'],
      beta = fixef(JOINT$poifitTIa)['trtT'],
      fixVarCor = vcov(JOINT$poifitTIa)[c('trtS', 'trtT'), c('trtS', 'trtT')],
      R2 = R2,
      ranef = ranef(JOINT$poifitTIa),
      VarCor = VarCorr(JOINT$poifitTIa),
      optinfo = JOINT$poifitTIa@optinfo,
      regr = regr,
      runTime = as.difftime(attr(JOINT$poifitTIa, 'exec.time')[3] / 60,
                            units = 'mins')
    )
    
    if (verbose)
      message(paste0(' (', format(RES$PoissonTIa$runTime, digits = 2), ')'))
  }#############################################################################
  ##############################################################################
  
  # options(warn = W)
  
  attr(RES, 'intWidth') <- intWidth
  attr(RES, 'nInts') <- attr(RES, 'nInts')
  return(RES)
}
