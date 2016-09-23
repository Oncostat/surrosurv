###############################################################################
# Functions for 1st step
############################################################################### 
fitClayton <- function(data, varcor=FALSE, 
                       optimx.method,
                       MAXFUN = 1e8, 
                       ini_rho = 0.5,
                       ini_kTau = 0.5) {
  OPTIMX_CONTROL <- list(maxit = MAXFUN,
                         dowarn = FALSE)
  Mpars <- margFits(data)
  
  INIpars <- c(log(copula::iRho(claytonCopula(), ini_rho)),
               #log(copula::iTau(claytonCopula(), ini_kTau)),
               log(Mpars$lambdaS), 
               log(Mpars$rhoS),
               Mpars$alpha, 
               log(Mpars$lambdaT),
               log(Mpars$rhoT),
               Mpars$beta)
  nTrials <- (length(INIpars) - 1) / 6
  
  #     library('optimx')
  RES <- optimx(INIpars, mloglik, data=data, family='clayton',
                method = optimx.method, #itnmax=5,
                control= OPTIMX_CONTROL)
  
  if (varcor) {
    #         library('optextras')
    attr(RES, 'grad') <- grnd(as.numeric(RES[, 1:attr(RES, 'npar')]),
                              mloglik, data=data, family='clayton')
    attr(RES, 'hessian') <- optimHess(RES[1:length(INIpars)],
                                      mloglik, data=data, family='clayton',
                                      control= OPTIMX_CONTROL)
    
    #         library('msm')
    tranform <- eval(parse(text = paste(
      # theta
      'list(~ exp(x1),\n', 
      # lambdaS and rhoS
      paste0('~ exp(x', 1 + 1:(2*nTrials), '),', collapse='\n'),
      # alpha
      paste0('~ x', 1 + 2*nTrials + 1:nTrials, ',', collapse='\n'), 
      # lambdaT and rhoT
      paste0('~ exp(x', 1 + 3*nTrials + 1:(2*nTrials), '),', collapse='\n'),
      # beta
      paste0('~ x', 1 + 5*nTrials + 1:nTrials, collapse=',\n'),
      ')'
    )))
    VarCor <- try(deltamethod(tranform,
                              mean = as.numeric(RES[1:length(INIpars)]),
                              cov = solve(attr(RES, 'hessian')),
                              ses = FALSE))
  }
  
  RES <- list(
    theta   = as.numeric(exp(RES[,1])),
    lambdaS = as.numeric(exp(RES[,1 + 1:nTrials])),
    rhoS    = as.numeric(exp(RES[,1 + nTrials + 1:nTrials])),
    alpha   = as.numeric(RES[,1 + 2 * nTrials + 1:nTrials]),
    lambdaT = as.numeric(exp(RES[,1 + 3 * nTrials + 1:nTrials])),
    rhoT    = as.numeric(exp(RES[,1 + 4 * nTrials + 1:nTrials])),
    beta    = as.numeric(RES[,1 + 5 * nTrials + 1:nTrials]),
    optimxRES = RES
  )
  if (varcor) RES <- c(RES, list(VarCor=VarCor))
  return(RES)
}

kTau.clay <- function(x) {
  #     library('copula')
  copula::tau(claytonCopula(x$theta))
}

############################################################################### 
fitPlackett <- function(data, varcor=FALSE, 
                        optimx.method,
                        MAXFUN = 1e8, 
                        ini_kTau = 0.5, 
                        ini_rho = 0.5) {
  OPTIMX_CONTROL <- list(maxit = MAXFUN,
                         dowarn = FALSE)

  Mpars <- margFits(data)
  
  INIpars <- c(log(copula::iRho(plackettCopula(), ini_rho)),
               #log(copula::iTau(plackettCopula(), ini_kTau)),
               log(Mpars$lambdaS), 
               log(Mpars$rhoS),
               Mpars$alpha, 
               log(Mpars$lambdaT),
               log(Mpars$rhoT),
               Mpars$beta)
  nTrials <- (length(INIpars) - 1) / 6
  
  #     library('optimx')
  RES <- optimx(INIpars, mloglik, data=data, family='plackett',
                method = optimx.method,
                control= OPTIMX_CONTROL)
  
  if (varcor) {
    #         library('optextras')
    attr(RES, 'grad') <- grnd(as.numeric(RES[, 1:attr(RES, 'npar')]),
                              mloglik, data=data, family='plackett')
    attr(RES, 'hessian') <- optimHess(RES[1:length(INIpars)], 
                                      mloglik, data=data, family='plackett',
                                      control= OPTIMX_CONTROL)
    #             library('msm')
    tranform <- eval(parse(text = paste(
      # theta
      'list(~ exp(x1),\n', 
      # lambdaS and rhoS
      paste0('~ exp(x', 1 + 1:(2*nTrials), '),', collapse='\n'),
      # alpha
      paste0('~ x', 1 + 2*nTrials + 1:nTrials, ',', collapse='\n'), 
      # lambdaT and rhoT
      paste0('~ exp(x', 1 + 3*nTrials + 1:(2*nTrials), '),', collapse='\n'),
      # beta
      paste0('~ x', 1 + 5*nTrials + 1:nTrials, collapse=',\n'),
      ')'
    )))
    VarCor <- try(deltamethod(tranform,
                              mean = as.numeric(RES[1:length(INIpars)]),
                              cov = solve(attr(RES, 'hessian')),
                              ses=FALSE))
  }
  
  RES <- list(
    theta   = as.numeric(exp(RES[,1])),
    lambdaS = as.numeric(exp(RES[,1 + 1:nTrials])),
    rhoS    = as.numeric(exp(RES[,1 + nTrials + 1:nTrials])),
    alpha   = as.numeric(RES[,1 + 2 * nTrials + 1:nTrials]),
    lambdaT = as.numeric(exp(RES[,1 + 3 * nTrials + 1:nTrials])),
    rhoT    = as.numeric(exp(RES[,1 + 4 * nTrials + 1:nTrials])),
    beta    = as.numeric(RES[,1 + 5 * nTrials + 1:nTrials]),
    optimxRES = RES
  )
  if (varcor) RES <- c(RES, list(VarCor=VarCor))
  return(RES)
}

kTau.plack <- function(x) {
  #     library('copula')
  copula::tau(plackettCopula(x$theta))
}

############################################################################### 
fitHougaard <- function(data, varcor=FALSE, 
                        optimx.method,
                        MAXFUN = 1e8, 
                        ini_rho = 0.5, 
                        ini_kTau = 0.5) {
  OPTIMX_CONTROL <- list(maxit = MAXFUN,
                         dowarn = FALSE)
  
  Mpars <- margFits(data)
  
  INIpars <- c(log(1 / copula::iRho(gumbelCopula(), ini_rho)),
               #log(1 - ini_kTau),
               log(Mpars$lambdaS), 
               log(Mpars$rhoS),
               Mpars$alpha, 
               log(Mpars$lambdaT),
               log(Mpars$rhoT),
               Mpars$beta)
  nTrials <- (length(INIpars) - 1) / 6
  
  #     library('optextras')
  RES <- optimx(INIpars, mloglik, data=data, family='hougaard',
                method = optimx.method,
                control= OPTIMX_CONTROL)
  
  if (varcor) {
    #         library('optimx')
    attr(RES, 'grad') <- grnd(as.numeric(RES[, 1:attr(RES, 'npar')]),
                              mloglik, data=data, family='hougaard')
    attr(RES, 'hessian') <- optimHess(RES[1:length(INIpars)], 
                                      mloglik, data=data, family='hougaard',
                                      control= OPTIMX_CONTROL)
    #             library('msm')
    tranform <- eval(parse(text = paste(
      # theta
      'list(~ exp(-exp(x1)),\n', 
      # lambdaS and rhoS
      paste0('~ exp(x', 1 + 1:(2*nTrials), '),', collapse='\n'),
      # alpha
      paste0('~ x', 1 + 2*nTrials + 1:nTrials, ',', collapse='\n'), 
      # lambdaT and rhoT
      paste0('~ exp(x', 1 + 3*nTrials + 1:(2*nTrials), '),', collapse='\n'),
      # beta
      paste0('~ x', 1 + 5*nTrials + 1:nTrials, collapse=',\n'),
      ')'
    )))
    VarCor <- try(deltamethod(tranform,
                              mean = as.numeric(RES[1:length(INIpars)]),
                              cov = solve(attr(RES, 'hessian')),
                              ses=FALSE))
  }
  
  RES <- list(
    theta   = as.numeric(exp(-exp(RES[,1]))),
    lambdaS = as.numeric(exp(RES[,1 + 1:nTrials])),
    rhoS    = as.numeric(exp(RES[,1 + nTrials + 1:nTrials])),
    alpha   = as.numeric(RES[,1 + 2 * nTrials + 1:nTrials]),
    lambdaT = as.numeric(exp(RES[,1 + 3 * nTrials + 1:nTrials])),
    rhoT    = as.numeric(exp(RES[,1 + 4 * nTrials + 1:nTrials])),
    beta    = as.numeric(RES[,1 + 5 * nTrials + 1:nTrials]),
    optimxRES = RES
  )
  if (varcor) RES <- c(RES, list(VarCor=VarCor))
  return(RES)
}

kTau.houg <- function(x) {
  1 - x$theta
}
############################################################################### 


margFits <- function(data) {
  tS <- data$timeS
  tT <- data$timeT
  cS <- data$statusS
  cT <- data$statusT
  
  #     library('eha')
  # Surrogate S
  sreg <- survreg(Surv(timeS, statusS) ~ 
                    trt:trialref + strata(trialref) + trialref, data=data)
  lS <- exp(-(
    sreg$coeff['(Intercept)'] + 
      c(0, sreg$coeff[(names(sreg$coeff) != '(Intercept)') &
                        !grepl('trt', names(sreg$coeff))])
  ) / sreg$scale)
  rS <- 1/sreg$scale
  bS <- -sreg$coeff[grep('trt', names(sreg$coef))] / sreg$scale
  
  #     sS <- exp(-lS[data$trialref] * tS^(rS[data$trialref]) * 
  #                   exp(lS[data$trialref] * data$trt))
  #     l_fS <- log(lS[data$trialref] * exp(lS[data$trialref] * data$trt) * 
  #                     rS[data$trialref] * 
  #                     tS^(rS[data$trialref]-1) *
  #                     sS)
  
  # True endpoint T
  treg <- survreg(Surv(timeT, statusT) ~ 
                    trt:trialref + strata(trialref) + trialref, data=data)
  lT <- exp(-(
    treg$coeff['(Intercept)'] + 
      c(0, treg$coeff[(names(treg$coeff) != '(Intercept)') &
                        !grepl('trt', names(treg$coeff))])
  ) / treg$scale)
  rT <- 1/treg$scale
  bT <- -treg$coeff[grep('trt', names(treg$coef))] / treg$scale
  
  #     sT <- exp(-lT[data$trialref] * tT^(rT[data$trialref]) * 
  #                   exp(lT[data$trialref] * data$trt))
  #     l_fT <- log(lT[data$trialref] * exp(lT[data$trialref] * data$trt) *
  #                     rT[data$trialref] * 
  #                     tT^(rT[data$trialref]-1) * 
  #                     sT)
  
  RES <- list(
    lambdaS = as.numeric(lS), 
    rhoS = as.numeric(rS), 
    alpha = as.numeric(bS),
    lambdaT = as.numeric(lT), 
    rhoT = as.numeric(rT), 
    beta = as.numeric(bT))
  return(RES)
}


mloglik <- function(pars, data, family=c('clayton', 'plackett', 'hougaard')) {
  family <- match.arg(family)
  # pars contains
  #  - one value for log(theta)
  #  - nTrials values for log(lambdaS)
  #  - nTrials values for log(rhoS)
  #  - nTrials values for alphaS
  #  - nTrials values for log(lambdaT)
  #  - nTrials values for log(rhoT)
  #  - nTrials values for betaT
  theta <- exp(pars[1])
  nTrials <- (length(pars) - 1) / 6
  lambdaS <- exp(pars[1 + 1:nTrials])
  rhoS <- exp(pars[1 + nTrials + 1:nTrials])
  alpha <- pars[1 + 2 * nTrials + 1:nTrials]
  lambdaT <- exp(pars[1 + 3 * nTrials + 1:nTrials])
  rhoT <- exp(pars[1 + 4 * nTrials + 1:nTrials])
  beta <- pars[1 + 5 * nTrials + 1:nTrials]
  
  ll <- 0
  
  tS <- data$timeS
  tT <- data$timeT
  cS <- data$statusS
  cT <- data$statusT
  
  lS <- lambdaS[data$trialref] * exp(alpha[data$trialref] * data$trt)
  rS <- rhoS[data$trialref]
  sS <- exp(-lS * tS^rS)
  l_fS <- log(lS * rS * tS^(rS-1) * sS)
  
  lT <- lambdaT[data$trialref] * exp(beta[data$trialref] * data$trt)
  rT <- rhoT[data$trialref]
  sT <- exp(-lT * tT^rT)
  l_fT <- log(lT * rT * tT^(rT-1) * sT)
  
  if (family == 'clayton') {
    ### S_ST(s, t) = (S_S(s)^(-theta) + S_T(t)^(-theta) - 1)^(-1/theta)
    l_Sst <- - log(sS^(-theta) + sT^(-theta) - 1) / theta
    
    # Both censored
    sel <- !cS & !cT
    if (any(sel)) {
      ll <- ll + sum(l_Sst[sel])
    }
    # S censored only
    sel <- !cS & cT
    if (any (sel)) {
      ll <- ll + sum(l_fT[sel]) + 
        (1 + theta) * sum(l_Sst[sel] - log(sT[sel]))
    }
    # T censored only
    sel <- cS & !cT
    if (any (sel)) {
      ll <- ll + sum(l_fS[sel]) + 
        (1 + theta) * sum(l_Sst[sel] - log(sS[sel]))
    }
    # both measured
    sel <- cS & cT
    if (any (sel)) {
      ll <- ll + sum(sel) * log(1+theta)  +
        sum(l_fS[sel] + l_fT[sel]) +
        (1+2*theta) *  sum(l_Sst[sel]) +
        - (1+theta) * sum(log(sS[sel] * sT[sel]))
    }
  } else if (family == 'plackett') {
    if (theta == 1) {
      l_Sst <- log(sS) + log(sT)
      
      if (any(!cS & !cT)) {
        ll <- ll + sum(l_Sst[!cS & !cT])
      }
      if (any (!cS & cT)) {
        ll <- ll + sum(log(sT[!cS & cT]) + l_fS[!cS & cT])
      }
      if (any (cS & !cT)) {
        ll <- ll + sum(log(sS[cS & !cT]) + l_fT[cS & !cT])
      }
      if (any (cS & cT)) {
        ll <- ll + sum(l_fS[cS & cT] + l_fT[cS & cT])
      }
    } else {
      ### S_ST(s, t) = (Q - sqrt(R))/(2 *(theta-1))
      ### Q = 1 + (theta-1) * (S_S(s)+S_T(t))
      ### R = Q^2 - 4 * S_S(s) * S_T(t) * theta * (theta-1)
      Q <- 1 + (theta-1) * (sS + sT)
      R <- Q^2 - 4 * theta * (theta-1) * sS * sT
      
      l_Sst <- log((Q - R^(1/2)) / (2 * (theta-1)))
      
      # Both censored
      sel <- !cS & !cT
      if (any(sel)) {
        ll <- ll + sum(l_Sst[sel])
      }
      # S censored only
      sel <- !cS & cT
      if (any (sel)) {
        ll <- ll - sum(sel) * log(2) + 
          sum(log(1 - (
            Q[sel] - 2 * theta * sS[sel]) * R[sel]^(-1/2)
          ) + l_fT[sel])
      }
      # T censored only
      sel <- cS & !cT
      if (any(sel)) {
        ll <- ll - sum(sel) * log(2) + 
          sum(log(1 - (
            Q[sel] - 2 * theta * sT[sel]) * R[sel]^(-1/2)
          ) + l_fS[sel])
      }
      # both measured
      sel <- cS & cT
      if (any (sel)) {
        ll <- ll + sum(sel) * log(theta) + sum(
          -3/2 * log(R[sel]) + 
            log(Q[sel] - 2 * (theta-1) * sS[sel] * sT[sel]
            ) + l_fS[sel] + l_fT[sel])
      }
    }
  } else if (family == 'hougaard') {
    ### S_ST(s, t) = exp(-Q^theta)
    ### Q = (-log u)^1/theta + (-log v)^1/theta
    Q = (-log(sS))^(1/theta) + (-log(sT))^(1/theta)
    l_Sst <- - Q^theta
    
    # Both censored
    sel <- !cS & !cT
    if (any(sel)) {
      ll <- ll + sum(l_Sst[sel])
    }
    # S censored only
    sel <- !cS & cT
    if (any (sel)) {
      ll <- ll + sum(l_fT[sel]) + 
        (1/theta-1) * sum(log(-log(sT[sel]))) +
        -sum(log(sT[sel])) + sum(l_Sst[sel]) +
        (theta - 1) * sum(log(Q[sel]))
    }
    # T censored only
    sel <- cS & !cT
    if (any (sel)) {
      ll <- ll + sum(l_fS[sel]) + 
        (1/theta-1) * sum(log(-log(sS[sel]))) +
        -sum(log(sS[sel])) + sum(l_Sst[sel]) +
        (theta - 1) * sum(log(Q[sel]))
    }
    # both measured
    sel <- cS & cT
    if (any (sel)) {
      ll <- ll + sum(l_fS[sel] + l_fT[sel]) +
        (1/theta-1) * sum(log(-log(sS[sel])) + log(-log(sT[sel]))) +
        -sum(log(sS[sel]) + log(sT[sel])) + sum(l_Sst[sel]) +
        (theta - 2) * sum(log(Q[sel])) +
        sum(log(1/theta - 1 + Q[sel]^theta))
    }
  }
  return(-ll)
}

###############################################################################
# Functions for 2nd step
############################################################################### 
copuSurr <- function(data, family=c('clayton', 'plackett', 'hougaard'), 
                     varcor1=FALSE, 
                     optimx.method,
                     INIrho = 0.5,
                     INIkTau = 0.5) {
  family <- tolower(family)
  family <- match.arg(family)
  margPars <- margFits(data)
  
  if (family == 'clayton'){
    step1 <- fitClayton(data, varcor=varcor1,
                        optimx.method = optimx.method,
                        ini_rho = INIrho
                        #ini_kTau = INIkTau
                        )
    kTau <- kTau.clay(step1)
  } else if (family == 'plackett'){
    step1 <- fitPlackett(data, varcor=varcor1,
                         optimx.method = optimx.method,
                         ini_rho = INIrho
                         #ini_kTau = INIkTau
                         )
    kTau <- kTau.plack(step1)
  } else if (family == 'hougaard'){
    step1 <- fitHougaard(data, varcor=varcor1,
                         optimx.method = optimx.method,
                         ini_rho = INIrho
                         #ini_kTau = INIkTau
                         )
    kTau <- kTau.houg(step1)
  }
  
  step1est <- data.frame(trteff = c(step1$alpha, step1$beta),
                         ep = rep(c('S', 'T'), each=length(step1$alpha)),
                         trial = rep(1:length(step1$alpha), 2), 
                         weights = as.vector(table(data$trialref)))
  
  ############################################################################
  # 2nd step using the mvmeta package                                        #
  #                                                                          #
  # Gasparrini A, Armstrong B, Kenward MG. Multivariate meta-analysis for    #
  #   non-linear abd other multi-parameter associations. Stat Med 2012; 31.  #
  #                                                                          #
  # the same as performed in SAS by                                          #
  # Burzykowski T,Molenberghs G, Buyse M, Geys H, Renard D. Validation of    #
  #   surrogate and points in multiple randomized clinical trials with fail- #
  #   ure time end points. Appl Staist 2001; 50.                             #
  # and                                                                      #
  # van Houwelingen mstrong B, Kenward MG. Multivariate meta-analysis for    #
  #   non-linear abd other multi-parameter associations. Stat Med 2012; 31.  #
  ############################################################################
  poss <- cumsum(sapply(step1, length))
  ind_AlphaBeta <- which(names(poss) %in% c('alpha', 'beta'))
  ind_AlphaBeta <- mapply(':', poss[ind_AlphaBeta-1]+1, poss[ind_AlphaBeta])
  myS <- with(step1est, cbind(
    alpha_Var = diag(step1$VarCor[ind_AlphaBeta[, 1], ind_AlphaBeta[, 1]]),
    albe_Covar = diag(step1$VarCor[ind_AlphaBeta[, 1], ind_AlphaBeta[, 2]]),
    beta_Var = diag(step1$VarCor[ind_AlphaBeta[, 2], ind_AlphaBeta[, 2]])))
  #     library('mvmeta')
  bvmodel <- mvmeta(cbind(alpha, beta), S = myS, 
                    data = as.data.frame(step1[c('alpha', 'beta')]))
  alpha <- bvmodel$coefficients[1]
  beta <- bvmodel$coefficients[2]
  R2 <- bvmodel$Psi[1, 2]^2 / prod(diag(bvmodel$Psi))    
  
  RES <- list(
    kTau = kTau,
    alpha = alpha,
    beta = beta,
    step1 = step1,
    step2 = bvmodel,
    R2 = R2,
    ranef = blup(bvmodel),
    optimxRES = step1$optimxRES,
    VarCor2 = bvmodel$Psi
  )
  if (varcor1) RES <- c(RES, list(VarCor1 = step1$VarCor))
  
  step2 <- lm(beta ~ alpha,
              data = as.data.frame(step1[c('alpha', 'beta')]))
  alpha <- mean(step1$alpha)
  beta <- mean(step1$beta)
  R2 <- summary(step2)$r.squared
  
  RES <- list(
    unadj = list(
      kTau = kTau,
      alpha = alpha,
      beta = beta,
      step1 = step1,
      step2 = step2,
      R2 = R2,
      ranef = NULL,
      optimxRES = step1$optimxRES,
      VarCor2 = NULL),
    adj = RES)
  if (varcor1) RES$unadj <- c(RES$unadj, list(VarCor1 = step1$VarCor))
  
  return(RES)
}

