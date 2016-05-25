poisSurr <- function(data,
                     OPTIMIZER = 'bobyqa',
                     MAXFUN = 1e8,
                     intWidth = 365.25/4) {
  #     library('phmm')
  #     library('parfm')
  #     library('survival')
  #     library('msm')
  #     library('lme4')
  
  GLMERCONTROL <- glmerControl(optimizer='optimx', 
                               boundary.tol = 0,
                               optCtrl = list(maxfun = MAXFUN,
                                              dowarn = FALSE,
                                              method = OPTIMIZER))
  
  #   GLMERCONTROL <- glmerControl(optimizer='bobyqa',
  #                                optCtrl = list(maxfun=MAXFUN, dowarn=FALSE),
  #                                boundary.tol=0)
  
  ############################################################################
  ### *** Data poissonization *** ############################################
  ############################################################################
  T <- S <- JOINT <- list()
  
  if ((min(data$trt) == 0) & (max(data$trt) == 1)) {
    message('Recoding treatment from 0/1 to -0.5/0.5')
    data$trt <- data$trt - .5
  }
  
  colnames(data)[colnames(data) == "timeS"] <- 'time'
  colnames(data)[colnames(data) == "statusS"] <- 'status'
  S$poidata <- poissonize(data, interval.width = intWidth, 
                          factors = c('trialref', 'trt', 'id'), compress=TRUE)
  colnames(data)[colnames(data) == "time"] <- 'timeS'
  colnames(data)[colnames(data) == "status"] <- 'statusS'
  
  colnames(data)[colnames(data) == "timeT"] <- 'time'
  colnames(data)[colnames(data) == "statusT"] <- 'status'
  T$poidata <- poissonize(data, interval.width = intWidth, 
                          factors = c('trialref', 'trt', 'id'), compress=TRUE)
  colnames(data)[colnames(data) == "time"] <- 'timeT'
  colnames(data)[colnames(data) == "status"] <- 'statusT'
  
  ############################################################################
  
  
  ############################################################################
  ### *** Mixed Poisson model estimation *** #################################
  ############################################################################
  JOINT$poidata <- rbind(cbind(S$poidata, ep='S'),
                         cbind(T$poidata, ep='T'))
  
  JOINT$poidata$trtT <- JOINT$poidata$trt * (JOINT$poidata$ep=="T")
  JOINT$poidata$trtS <- JOINT$poidata$trt * (JOINT$poidata$ep=="S")
  
  JOINT$poidata$trialrefT <- JOINT$poidata$trialrefS <-
    JOINT$poidata$trialref
  JOINT$poidata$trialrefT[JOINT$poidata$ep != "T"]  <- 
    levels(JOINT$poidata$trialref)[1]
  JOINT$poidata$trialrefS[JOINT$poidata$ep != "S"]  <- 
    levels(JOINT$poidata$trialref)[1]
  
  W <- options()$warn
  options(warn=-1)
  
  if (nlevels(JOINT$poidata$interval) > 1) {
    baseform <- m~-1+interval*ep-ep + trtT + trtS + offset(log(Rt))
  } else  {
    baseform <- m~-1+ ep + trtT + trtS + offset(log(Rt))
  }
  
  #     library('optimx')
  # * Model 2a: Poisson model with
  #             - random treatment-trial interaction
  system.time({
    JOINT$poifit2a <- glmer(
      update(baseform, .~. + (-1+trtT+trtS|trialref)),
      data=JOINT$poidata, family=poisson,
      control=GLMERCONTROL)
  }) -> attr(JOINT$poifit2a, 'exec.time')
  
  # * Model 2b: Poisson model with
  #             - individual random intercept
  system.time({
    JOINT$poifit2b <-  glmer(
      update(baseform, .~. + (1|id)) ,
      data=JOINT$poidata, family=poisson,
      control=GLMERCONTROL)
  }) -> attr(JOINT$poifit2b, 'exec.time')
  
  attr(JOINT$poifit2b, 'kTau') <- parfm:::fr.lognormal(
    what = 'tau',
    sigma2 = as.double(summary(JOINT$poifit2b)$varcor$id))
  
  # * Model 3: Poisson model with
  #            - random treatment-trial interaction
  #            - individual random intercept 
  system.time({
    JOINT$poifit3 <- update(JOINT$poifit2a, .~. + (1|id))
  }) -> attr(JOINT$poifit3, 'exec.time')
  
  attr(JOINT$poifit3, 'kTau') <- parfm:::fr.lognormal(
    what = 'tau',
    sigma2 = as.double(summary(JOINT$poifit3)$varcor$id))
  
  # * Model 5b: Poisson model with
  #             - random treatment-trial interaction
  #             - individual random intercept 
  #             - random trial intercept (shared by the two EndPoints)
  system.time({
    JOINT$poidata$trialrefSH <- JOINT$poidata$trialref
    JOINT$poifit5b <- update(JOINT$poifit3, .~. + (1|trialrefSH))
  }) -> attr(JOINT$poifit5b, 'exec.time')
  
  attr(JOINT$poifit5b, 'kTau') <- parfm:::fr.lognormal(
    what = 'tau',
    sigma2 = as.double(summary(JOINT$poifit5b)$varcor$id))
  ############################################################################
  
  options(warn=W)
  
  RES <- list(
    model2a = list(
      kTau = NULL,
      alpha = fixef(JOINT$poifit2a)['trtS'],
      beta = fixef(JOINT$poifit2a)['trtT'],
      R2 = VarCorr(JOINT$poifit2a)$trialref[1, 2]^2 /
        prod(diag(VarCorr(JOINT$poifit2a)$trialref)),
      ranef = ranef(JOINT$poifit2a),
      VarCor = VarCorr(JOINT$poifit2a),
      optinfo = JOINT$poifit2a@optinfo
    ),
    model2b = list(
      kTau = attr(JOINT$poifit2b, 'kTau'),
      alpha = fixef(JOINT$poifit2b)['trtS'],
      beta = fixef(JOINT$poifit2b)['trtT'],
      R2 = NULL,
      ranef = ranef(JOINT$poifit2b),
      VarCor = VarCorr(JOINT$poifit2b),
      optinfo = JOINT$poifit2b@optinfo
    ),
    model3 = list(
      kTau = attr(JOINT$poifit3, 'kTau'),
      alpha = fixef(JOINT$poifit3)['trtS'],
      beta = fixef(JOINT$poifit3)['trtT'],
      R2 = VarCorr(JOINT$poifit3)$trialref[1, 2]^2 /
        prod(diag(VarCorr(JOINT$poifit3)$trialref)),
      ranef = ranef(JOINT$poifit3),
      VarCor = VarCorr(JOINT$poifit3),
      optinfo = JOINT$poifit3@optinfo
    ),
    model5b = list(
      kTau = attr(JOINT$poifit5b, 'kTau'),
      alpha = fixef(JOINT$poifit5b)['trtS'],
      beta = fixef(JOINT$poifit5b)['trtT'],
      R2 = VarCorr(JOINT$poifit5b)$trialref[1, 2]^2 /
        prod(diag(VarCorr(JOINT$poifit5b)$trialref)),
      ranef = ranef(JOINT$poifit5b),
      VarCor = VarCorr(JOINT$poifit5b),
      optinfo = JOINT$poifit5@optinfo
    )
  )
  return(RES)
}
