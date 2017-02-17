poisSurr <- function(data,
                     OPTIMIZER, # = 'bobyqa',
                     MAXFUN = 1e8,
                     intWidth = NULL,
                     nInts = NULL,
                     models = paste('Poisson', c('I', 'T', 'TI', 'TIa')),
                     verbose = TRUE) {
  models <- tolower(noSpP(models))
  models <- match.arg(models, several.ok = TRUE, 
                      choices = paste0('poisson', c('i', 't', 'ti', 'tia')))
  #     library('parfm')
  #     library('survival')
  #     library('msm')
  #     library('lme4')
  
  GLMERCONTROL <- glmerControl(optimizer = 'optimx', 
                               # calc.derivs = FALSE,
                               boundary.tol = 0,
                               optCtrl = list(maxfun = MAXFUN,
                                              dowarn = FALSE,
                                              method = OPTIMIZER,
                                              starttests = FALSE,
                                              kkt = FALSE))
  
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
  } else all.breaks <- NULL
  
  colnames(data)[colnames(data) == "timeS"] <- 'time'
  colnames(data)[colnames(data) == "statusS"] <- 'status'
  S$poidata <- poissonize(data, all.breaks = all.breaks, 
                          interval.width = intWidth, nInts = nInts,
                          factors = c('trialref', 'trt', 'id'), compress = TRUE)
  colnames(data)[colnames(data) == "time"] <- 'timeS'
  colnames(data)[colnames(data) == "status"] <- 'statusS'
  
  colnames(data)[colnames(data) == "timeT"] <- 'time'
  colnames(data)[colnames(data) == "statusT"] <- 'status'
  T$poidata <- poissonize(data, all.breaks = all.breaks, 
                          interval.width = intWidth, nInts = nInts, 
                          factors = c('trialref', 'trt', 'id'), compress = TRUE)
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
  
  JOINT$poidata$trtT <- JOINT$poidata$trt * (JOINT$poidata$ep == "T")
  JOINT$poidata$trtS <- JOINT$poidata$trt * (JOINT$poidata$ep == "S")
  
  JOINT$poidata$trialrefT <- JOINT$poidata$trialrefS <-
    JOINT$poidata$trialref
  JOINT$poidata$trialrefT[JOINT$poidata$ep != "T"]  <- 
    levels(JOINT$poidata$trialref)[1]
  JOINT$poidata$trialrefS[JOINT$poidata$ep != "S"]  <- 
    levels(JOINT$poidata$trialref)[1]
  
  W <- options()$warn
  options(warn = -1)
  
  if (nlevels(JOINT$poidata$interval) > 1) {
    baseform <- m~ -1 +interval*ep-ep + trtT + trtS + offset(log(Rt))
  } else  {
    baseform <- m~ -1 + ep + trtT + trtS + offset(log(Rt))
  }
  
  RES <- list()
  
  ##############################################################################
  ##############################################################################
  #     library('optimx')
  # * Model T: Poisson model with
  #             - random treatment-trial interaction
  if ('poissont' %in% models) {
    if (verbose) message('- Estimating model: Poisson T', appendLF = FALSE)
    
    system.time({
      JOINT$poifitT <- glmer(
        update(baseform, .~. +
                 (-1 + trtT + trtS | trialref)),
        data = JOINT$poidata, family = poisson,
        control = GLMERCONTROL)
    }) -> attr(JOINT$poifitT, 'exec.time')
    
    RES$PoissonT <- list(
      kTau = NULL,
      alpha = fixef(JOINT$poifitT)['trtS'],
      beta = fixef(JOINT$poifitT)['trtT'],
      fixVarCor = vcov(JOINT$poifitT)[c('trtS', 'trtT'), c('trtS', 'trtT')],
      R2 = VarCorr(JOINT$poifitT)$trialref[1, 2]^2 /
        prod(diag(VarCorr(JOINT$poifitT)$trialref)),
      ranef = ranef(JOINT$poifitT),
      VarCor = VarCorr(JOINT$poifitT),
      optinfo = JOINT$poifitT@optinfo,
      runTime = as.difftime(
        round(attr(JOINT$poifitT, 'exec.time')[3] / 60, 1), units = 'mins'))
    
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
        update(baseform, .~. +
                 (1 | id)) ,
        data = JOINT$poidata, family = poisson,
        control = GLMERCONTROL)
    }) -> attr(JOINT$poifitI, 'exec.time')
    
    attr(JOINT$poifitI, 'kTau') <- parfm:::fr.lognormal(
      what = 'tau',
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
      runTime = as.difftime(
        round(attr(JOINT$poifitI, 'exec.time')[3] / 60, 1), units = 'mins'))
    
    if (verbose) 
      message(paste0(' (', format(RES$PoissonI$runTime, digits = 2), ')'))
  }#############################################################################
  
  ##############################################################################
  # * Model TI: Poisson model with
  #            - random treatment-trial interaction
  #            - individual random intercept 
  if ('poissonti' %in% models) {
    if (verbose) 
      message('- Estimating model: Poisson TI', appendLF = FALSE)
    
    system.time({
      JOINT$poifitTI <- glmer(
        update(baseform, .~. +
                 (-1 + trtT + trtS | trialref) +
                 (1 | id)),
        data = JOINT$poidata, family = poisson,
        control = GLMERCONTROL)
    }) -> attr(JOINT$poifitTI, 'exec.time')
    
    attr(JOINT$poifitTI, 'kTau') <- parfm:::fr.lognormal(
      what = 'tau',
      sigma2 = as.double(summary(JOINT$poifitTI)$varcor$id))
    
    RES$PoissonTI <- list(
      kTau = attr(JOINT$poifitTI, 'kTau'),
      alpha = fixef(JOINT$poifitTI)['trtS'],
      beta = fixef(JOINT$poifitTI)['trtT'],
      fixVarCor = vcov(JOINT$poifitTI)[c('trtS', 'trtT'), c('trtS', 'trtT')],
      R2 = VarCorr(JOINT$poifitTI)$trialref[1, 2]^2 /
        prod(diag(VarCorr(JOINT$poifitTI)$trialref)),
      ranef = ranef(JOINT$poifitTI),
      VarCor = VarCorr(JOINT$poifitTI),
      optinfo = JOINT$poifitTI@optinfo,
      runTime = as.difftime(
        round(attr(JOINT$poifitTI, 'exec.time')[3] / 60, 1), units = 'mins'))
    
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
        update(baseform, .~. +
                 (-1 + trtT + trtS | trialref) +
                 (1 | id) +
                 (1 | trialrefSH)),
        data = JOINT$poidata, family = poisson,
        control = GLMERCONTROL)
    }) -> attr(JOINT$poifitTIa, 'exec.time')
    
    attr(JOINT$poifitTIa, 'kTau') <- parfm:::fr.lognormal(
      what = 'tau',
      sigma2 = as.double(summary(JOINT$poifitTIa)$varcor$id))
    
    RES$PoissonTIa <- list(
      kTau = attr(JOINT$poifitTIa, 'kTau'),
      alpha = fixef(JOINT$poifitTIa)['trtS'],
      beta = fixef(JOINT$poifitTIa)['trtT'],
      fixVarCor = vcov(JOINT$poifitTIa)[c('trtS', 'trtT'), c('trtS', 'trtT')],
      R2 = VarCorr(JOINT$poifitTIa)$trialref[1, 2]^2 /
        prod(diag(VarCorr(JOINT$poifitTIa)$trialref)),
      ranef = ranef(JOINT$poifitTIa),
      VarCor = VarCorr(JOINT$poifitTIa),
      optinfo = JOINT$poifitTIa@optinfo,
      runTime = as.difftime(attr(JOINT$poifitTIa, 'exec.time')[3] / 60,
                            units = 'mins'))
    
    if (verbose) 
      message(paste0(' (', format(RES$PoissonTIa$runTime, digits = 2), ')'))
  }#############################################################################
  ##############################################################################
  
  options(warn = W)
  
  attr(RES, 'intWidth') <- intWidth
  attr(RES, 'nInts') <- attr(RES, 'nInts')
  return(RES)
}
