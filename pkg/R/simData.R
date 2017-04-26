# ################################################################################
# ################################################################################
# ################################################################################
# # Simualtion procedure based on a mixture of exponential and truncated normal  #
# #   random variables, described by Shi et al (2011) and Renfro et al (2012)    #
# #                                                                              #
# #   Shi Q, Renfro LA, Bot BM, Burzykowski T, Buyse M, Sargent DJ (2011).       #
# #   Comparative assessment of trial-level surrogacy measures for candidate     #
# #   time-to-event surrogate endpoints in clinical trials. Computational        #
# #   statistics & data analysis, 55(9), 2748-2757.                              #
# #                                                                              #
# #   Renfro LA, Shi Q, Sargent DJ, Carlin BP (2012). Bayesian adjusted R2 for   #
# #   the meta-analytic evaluation of surrogate time-to-event endpoint           #
# #   clinical trials. Statistics in medicine, 31(8), 743-761.                   #
# ################################################################################
# ################################################################################
# ################################################################################
# simData.mx <- function(
#   ############################################################################
#   ############################ *** PARAMETERS *** ############################
#   R2 = 0.6,           # adjusted trial-level R^2
#   N = 30,             # number of trials
#   ni = 2000,          # number of patients per trial
#   nifix = TRUE,       # is ni fix? (or average)
#   gammaWei = c(1, 1), # shape parameters of the Weibull distributions
#   censorT,            # censoring rate for the true endpoint T
#   censorA,            # administrative censoring at time censorA
#   indCorr= TRUE,      # S and T correlated?
#   baseCorr = 0.5,     # correlation between baseline hazards (rho_basehaz)
#   baseVars = c(.5,.5),# variances of baseline random effects (S and T)
#   alpha = 0,          # average treatment effect on S
#   beta = 0,           # average treatment effect on T
#   alphaVar = 1,       # variance of a_i (theta_aa^2)
#   betaVar = 1,        # variance of b_i (theta_bb^2)
#   mstS = 4*365.25,    # median survival time for S in the control arm
#   mstT = 8*365.25     # median survival time for T in the control arm
#   ############################################################################
# ) {
#   if (length(gammaWei) == 1)
#     gammaWei <- rep(gammaWei, 2)
#   
#   if (nifix) {
#     trialref <- rep(1:N, each = ni)
#   } else {
#     trialref <- sort(sample(1:N, N * ni, replace = TRUE, prob = rpois(N, 5)))
#   }
#   
#   data <- data.frame(
#     trialref = factor(trialref),
#     trt = rbinom(n = ni * N, size = 1, prob = 0.5) - .5
#   )
#   data$id <- factor(mapply(paste, data$trialref, 
#                            unlist(lapply(table(data$trialref), function(x) 1:x)),
#                            #rep(1:ni, N), 
#                            sep='.'))
#   
#   # First stage: regression coefficients *********************************** #
#   muS <- log(mstS / log(2))
#   muT <- log(mstT / log(2))
#   d_ab <- sqrt(R2 * alphaVar * betaVar)
#   d_ST <- sqrt(baseCorr) * sqrt(prod(baseVars))
#   Sigma_trial <- matrix(c(baseVars[1],        d_ST,        0,       0,
#                           d_ST,        baseVars[2],        0,       0,
#                           0,                     0, alphaVar,    d_ab,
#                           0,                     0,     d_ab, betaVar ), 4)
#   #   library('MASS')
#   pars <- mvrnorm(n= N, mu = c(muS, muT, -alpha, -beta), 
#                   Sigma = Sigma_trial)
#   
#   ATTRs <- list(
#     'N'=N,
#     'ni'=ni,
#     'gammaWei'=gammaWei,
#     'censorT'= ifelse(missing(censorT), NA, censorT),
#     'censorA'= ifelse(missing(censorA), NA, censorA),
#     'indCorr'=indCorr,
#     'baseCorr'=baseCorr,
#     'baseVars'=baseVars,
#     'alphaVar'= alphaVar,
#     'betaVar'=betaVar,
#     'mstS'=mstS,
#     'mstT'=mstT,
#     'alpha'=alpha,
#     'beta'=beta,
#     'pars'=pars,
#     'R2'=R2
#   )
#   pars <- pars[data$trialref, ]
#   # ************************************************************************ #
#   
#   
#   # Second stage: survival times ******************************************* #
#   # Y, the truncated normal random variables
#   #   library('msm')
#   Y <- rtnorm(n = ni * N, lower=0)
#   if (indCorr) {
#     Y <- cbind(Y, Y)
#   } else {
#     Y <- cbind(Y, rtnorm(n = ni * N, lower=0))
#   }
#   
#   # lambda, the exponential random variables
#   lambdaS <- rexp(n = ni * N, rate = 1)
#   lambdaT <- rexp(n = ni * N, rate = 1)
#   
#   # S and T, the times
#   deltaS <- exp(pars[, 1] + pars[, 3] * data$trt)
#   deltaT <- exp(pars[, 2] + pars[, 4] * data$trt)
#   data$S <- deltaS * (Y[, 1] * sqrt(2 * lambdaS))^(1 / gammaWei[1])
#   data$T <- deltaT * (Y[, 2] * sqrt(2 * lambdaT))^(1 / gammaWei[2])
#   # ************************************************************************ #
#   
#   
#   # *** Censoring ********************************************************** #
#   data$C <- rep(Inf, nrow(data))
#   # Random censoring
#   if (!missing(censorT)) {
#     findK <- Vectorize(function(logK) {
#       (mean(runif(10 * length(data$T), 0, exp(logK)) < data$T) - censorT)^2
#     })
#     suppressWarnings({k <- exp(optim(log(max(data$T)), findK)$par)})
#     data$C <- runif(length(data$T), 0, k)
#   }
#   # Administrative censoring
#   if (!missing(censorA)) {
#     data$C <- pmin(data$C, censorA)
#   }
#   
#   data$timeT <- pmin(data$C, data$T)
#   data$statusT <- mapply('>=', data$C, data$T) * 1
#   data$timeS <- pmin(data$C, data$S)
#   data$statusS <- mapply('>=', data$C, data$S) * 1
#   # kTau <- cor(data$S, data$T, method='kendall')
#   kTau <- mean(
#     cor(data$S[data$trt == -.5], data$T[data$trt == -.5], method='kendall'),
#     cor(data$S[data$trt == 0.5], data$T[data$trt == 0.5], method='kendall'))
#   data <- data[, !(names(data) %in% c('C', 'T', 'S'))]
#   # ************************************************************************ #
#   
#   attributes(data) <- c(
#     attributes(data), ATTRs, list(
#       'kTau' = kTau
#     ))
#   return(data)
# }
# ################################################################################
# # summary(DATA <- simData())
# # c(mean(DATA$C<DATA$S), mean(DATA$C<DATA$T))
# 
# 
# ################################################################################
# ################################################################################
# ################################################################################
# # Simualtion procedure based on proportional hazard models with random effects #
# #   both at individual and at trial level                                      #
# ################################################################################
# ################################################################################
# ################################################################################
logsigma2 <- seq(-2, 4.5, by = .25)
kTaus <- Vectorize(function(lnSigma2)
  parfm:::fr.lognormal(what = 'tau', sigma2 = exp(lnSigma2)))
# plot(exp(logsigma2), kTaus_ln(logsigma2))

find.sigma2 <- function(tau) {
  loss <- function(x)
    (kTaus(x) - tau) ^ 2
  RES <- optimize(loss, c(-10, 4))
  return(exp(RES$minimum))
}
# 
# simData.re <- function(
#   ############################################################################
#   ############################ *** PARAMETERS *** ############################
#   R2 = 0.6,           # adjusted trial-level R^2
#   N = 30,             # number of trials
#   ni = 2000,          # number of patients per trial
#   nifix = TRUE,       # is ni fix? (or average)
#   gammaWei = c(1, 1), # shape parameters of the Weibull distributions
#   censorT,            # censoring rate for the true endpoint T
#   censorA,            # administrative censoring at time censorA
#   kTau= 0.6,          # individual dependence between S and T (Kendall's tau)
#   baseCorr = 0.5,     # correlation between baseline hazards (rho_basehaz)
#   baseVars = c(.5,.5),# variances of baseline random effects (S and T)
#   alpha = 0,          # average treatment effect on S
#   beta = 0,           # average treatment effect on T
#   alphaVar = 1,       # variance of a_i (theta_aa^2)
#   betaVar = 1,        # variance of b_i (theta_bb^2)
#   mstS = 4*365.25,    # median survival time for S in the control arm
#   mstT = 8*365.25     # median survival time for T in the control arm
#   ############################################################################
# ) {
#   if (length(gammaWei) == 1)
#     gammaWei <- rep(gammaWei, 2)
#   
#   if (nifix) {
#     trialref <- rep(1:N, each = ni)
#   } else {
#     trialref <- sort(sample(1:N, N * ni, replace = TRUE, prob = rpois(N, 5)))
#   }
#   
#   data <- data.frame(
#     trialref = factor(trialref),
#     trt = rbinom(n = ni * N, size = 1, prob = 0.5) - .5
#   )
#   data$id <- factor(mapply(paste, data$trialref, 
#                            unlist(lapply(table(data$trialref), function(x) 1:x)),
#                            #rep(1:ni, N), 
#                            sep='.'))
#   
#   # First stage: random effects ******************************************** #
#   # Trial-level
#   muS <- log(log(2) / mstS)
#   muT <- log(log(2) / mstT)
#   d_ab <- sqrt(R2 * alphaVar * betaVar)
#   d_ST <- sqrt(baseCorr) * sqrt(prod(baseVars))
#   Sigma_trial <- matrix(c(baseVars[1],        d_ST,        0,       0,
#                           d_ST,        baseVars[2],        0,       0,
#                           0,                     0, alphaVar,    d_ab,
#                           0,                     0,     d_ab, betaVar ), 4)
#   #   library('MASS')
#   pars <- mvrnorm(n= N, mu = c(muS, muT, alpha, beta), 
#                   Sigma = Sigma_trial)
#   
#   ATTRs <- list(
#     'N'=N,
#     'ni'=ni,
#     'gammaWei'=gammaWei,
#     'censorT'= ifelse(missing(censorT), NA, censorT),
#     'censorA'= ifelse(missing(censorA), NA, censorA),
#     'kTau'=kTau,
#     'baseCorr'=baseCorr,
#     'baseVars'=baseVars,
#     'alphaVar'= alphaVar,
#     'betaVar'=betaVar,
#     'mstS'=mstS,
#     'mstT'=mstT,
#     'alpha'=alpha,
#     'beta'=beta,
#     'pars'=pars,
#     'R2'=R2
#   )
#   pars <- pars[data$trialref, ]
#   
#   # Individ2ual-level
#   sigma2 <- find.sigma2(kTau)
#   indFrailty <- rnorm(nrow(data), mean = 0, sd = sqrt(sigma2))
#   # ************************************************************************ #
#   
#   # Second stage: survival times ******************************************* #
#   # Weibull hazard:
#   # h(t) = lambda gamma x^(gamma-1)
#   lambda.S <- exp(pars[, 1] + pars[, 3] * data$trt + indFrailty)
#   lambda.T <- exp(pars[, 2] + pars[, 4] * data$trt + indFrailty)
#   
#   # Weibull hazard, parametrization as in rweibull:
#   # h(t) = shape scale^(-shape) x^(shape-1)
#   # shape = gamma
#   # scale = lambda^(-1/gamma)
#   shape.S <- gammaWei[1]
#   scale.S <- lambda.S^-(1/shape.S)
#   shape.T <- gammaWei[2]    
#   scale.T <- lambda.T^-(1/shape.T)
#   
#   # S and T, the times
#   data$S <- rweibull(nrow(data), shape = shape.S, scale = scale.S)
#   data$T <- rweibull(nrow(data), shape = shape.T, scale = scale.T)
#   # ************************************************************************ #
#   
#   
#   # Censoring ************************************************************** #
#   data$C <- rep(Inf, nrow(data))
#   # Random censoring
#   if (!missing(censorT)) {
#     findK <- Vectorize(function(logK) {
#       (mean(runif(10 * length(data$T), 0, exp(logK)) < data$T) - censorT)^2
#     })
#     suppressWarnings({k <- exp(optim(log(max(data$T)), findK)$par)})
#     data$C <- runif(length(data$T), 0, k)
#   }
#   # Administrative censoring
#   if (!missing(censorA)) {
#     data$C <- pmin(data$C, censorA)
#   }
#   
#   data$timeT <- pmin(data$C, data$T)
#   data$statusT <- mapply('>=', data$C, data$T) * 1
#   data$timeS <- pmin(data$C, data$S)
#   data$statusS <- mapply('>=', data$C, data$S) * 1
#   # kTau <- cor(data$S, data$T, method='kendall')
#   kTau <- mean(
#     cor(data$S[data$trt == -.5], data$T[data$trt == -.5], method='kendall'),
#     cor(data$S[data$trt == 0.5], data$T[data$trt == 0.5], method='kendall'))
#   data <- data[, !(names(data) %in% c('C', 'T', 'S'))]
#   # ************************************************************************ #
#   
#   attributes(data) <- c(
#     attributes(data), ATTRs, list(
#       'kTau.empir' = kTau
#     ))
#   return(data)
# }
# 
# 
# ################################################################################
# ################################################################################
# ################################################################################
# # Simualtion procedure based on proportional hazard models                     #
# #   with Clayton copulas                                                       #
# #                                                                              #
# #   Burzykowski T, Abrahantes JC (2005). Validation in the case of two         #
# #   failure-time endpoints. In The Evaluation of Surrogate Endpoints           #
# #   (pp. 163-194). Springer New York.                                          #
# ################################################################################
# ################################################################################
# ################################################################################
# simData.cc <- function(
#   ############################################################################
#   ############################ *** PARAMETERS *** ############################
#   R2 = 0.6,           # adjusted trial-level R^2
#   N = 30,             # number of trials
#   ni = 2000,          # number of patients per trial
#   nifix = TRUE,       # is ni fix? (or average)
#   gammaWei = c(1, 1), # shape parameters of the Weibull distributions
#   censorT,            # censoring rate for the true endpoint T
#   censorA,            # administrative censoring at time censorA
#   kTau= 0.6,          # individual dependence between S and T (Kendall's tau)
#   baseCorr = 0.5,     # correlation between baseline hazards (rho_basehaz)
#   baseVars = c(.5,.5),# variances of baseline random effects (S and T)
#   alpha = 0,          # average treatment effect on S
#   beta = 0,           # average treatment effect on T
#   alphaVar = 1,       # variance of a_i (theta_aa^2)
#   betaVar = 1,        # variance of b_i (theta_bb^2)
#   mstS = 4*365.25,    # median survival time for S in the control arm
#   mstT = 8*365.25     # median survival time for T in the control arm
#   ############################################################################
# ) {
#   if (length(gammaWei) == 1)
#     gammaWei <- rep(gammaWei, 2)
#   
#   if (nifix) {
#     trialref <- rep(1:N, each = ni)
#   } else {
#     trialref <- sort(sample(1:N, N * ni, replace = TRUE, prob = rpois(N, 5)))
#   }
#   
#   data <- data.frame(
#     trialref = factor(trialref),
#     trt = rbinom(n = ni * N, size = 1, prob = 0.5) - .5
#   )
#   data$id <- factor(mapply(paste, data$trialref, 
#                            unlist(lapply(table(data$trialref), function(x) 1:x)),
#                            #rep(1:ni, N), 
#                            sep='.'))
#   
#   # First stage: random effects ******************************************** #
#   # Trial-level
#   muS <- log(log(2) / mstS)
#   muT <- log(log(2) / mstT)
#   d_ab <- sqrt(R2 * alphaVar * betaVar)
#   d_ST <- sqrt(baseCorr) * sqrt(prod(baseVars))
#   Sigma_trial <- matrix(c(baseVars[1],        d_ST,        0,       0,
#                           d_ST,        baseVars[2],        0,       0,
#                           0,                     0, alphaVar,    d_ab,
#                           0,                     0,     d_ab, betaVar ), 4)
#   #   library('MASS')
#   pars <- mvrnorm(n= N, mu = c(muS, muT, alpha, beta), 
#                   Sigma = Sigma_trial)
#   
#   ATTRs <- list(
#     'N'=N,
#     'ni'=ni,
#     'gammaWei'=gammaWei,
#     'censorT'= ifelse(missing(censorT), NA, censorT),
#     'censorA'= ifelse(missing(censorA), NA, censorA),
#     'kTau'=kTau,
#     'baseCorr'=baseCorr,
#     'baseVars'=baseVars,
#     'alphaVar'= alphaVar,
#     'betaVar'=betaVar,
#     'mstS'=mstS,
#     'mstT'=mstT,
#     'alpha'=alpha,
#     'beta'=beta,
#     'pars'=pars,
#     'R2'=R2
#   )
#   pars <- pars[data$trialref, ]
#   
#   # Individ2ual-level
#   theta <- 2 * kTau / (1 - kTau)
#   # ************************************************************************ #
#   
#   # Second stage: survival times ******************************************* #
#   # Weibull hazard:
#   # h(t) = lambda gamma x^(gamma-1)
#   lambda.S <- exp(pars[, 1] + pars[, 3] * data$trt)
#   lambda.T <- exp(pars[, 2] + pars[, 4] * data$trt)
#   
#   # S times
#   US <- runif(nrow(data), 0, 1)
#   data$S <- (-log(US) / lambda.S)^(1/gammaWei[1])
#   
#   # T times | S times
#   UT <- runif(nrow(data), 0, 1)
#   UT_prime <- ((UT^(-theta / (1+theta)) - 1) * US^(-theta) + 1)^(-1/theta)
#   data$T <- (-log(UT_prime) / lambda.T)^(1/gammaWei[2])
#   # ************************************************************************ #
#   #     plot(data$S, data$T, log='xy', col=data$trt+1.5)
#   
#   # Censoring ************************************************************** #
#   data$C <- rep(Inf, nrow(data))
#   # Random censoring
#   if (!missing(censorT)) {
#     findK <- Vectorize(function(logK) {
#       (mean(runif(10 * length(data$T), 0, exp(logK)) < data$T) - censorT)^2
#     })
#     suppressWarnings({k <- exp(optim(log(max(data$T)), findK)$par)})
#     data$C <- runif(length(data$T), 0, k)
#   }
#   # Administrative censoring
#   if (!missing(censorA)) {
#     data$C <- pmin(data$C, censorA)
#   }
#   
#   data$timeT <- pmin(data$C, data$T)
#   data$statusT <- mapply('>=', data$C, data$T) * 1
#   data$timeS <- pmin(data$C, data$S)
#   data$statusS <- mapply('>=', data$C, data$S) * 1
#   # kTau <- cor(data$S, data$T, method='kendall')
#   kTau <- mean(
#     cor(data$S[data$trt == -.5], data$T[data$trt == -.5], method='kendall'),
#     cor(data$S[data$trt == 0.5], data$T[data$trt == 0.5], method='kendall'))
#   data <- data[, !(names(data) %in% c('C', 'T', 'S'))]
#   # ************************************************************************ #
#   
#   attributes(data) <- c(
#     attributes(data), ATTRs, list(
#       'kTau.empir' = kTau
#     ))
#   return(data)
# }



simData <- function(method = c('re', 'cc', 'mx')) {
  method <- match.arg(method)
  
  simfun <- function(
    ############################################################################
    ############################ *** PARAMETERS *** ############################
    R2 = 0.6,           # adjusted trial-level R^2
    N = 30,             # number of trials
    ni = 200,           # number of patients per trial
    nifix = TRUE,       # is ni fix? (or average)
    gammaWei = c(1, 1), # shape parameters of the Weibull distributions
    censorT,            # censoring rate for the true endpoint T
    censorA,            # administrative censoring at time censorA
    kTau = 0.6,         # individual dependence between S and T (Kendall's tau)
    indCorr = TRUE,     # S and T correlated?
    baseCorr = 0.5,     # correlation between baseline hazards (rho_basehaz)
    baseVars = c(.2,.2),# variances of baseline random effects (S and T)
    alpha = 0,          # average treatment effect on S
    beta = 0,           # average treatment effect on T
    alphaVar = 0.1,     # variance of a_i (theta_aa^2)
    betaVar = 0.1,      # variance of b_i (theta_bb^2)
    mstS = 4*365.25,    # median survival time for S in the control arm
    mstT = 8*365.25     # median survival time for T in the control arm
    ############################################################################
  ) {
    
    if (length(gammaWei) == 1)
      gammaWei <- rep(gammaWei, 2)
    
    if (nifix) {
      trialref <- rep(1:N, each = ni)
    } else {
      trialref <-
        sort(sample(1:N, N * ni, replace = TRUE, prob = rpois(N, 5)))
    }
    
    data <- data.frame(trialref = factor(trialref),
                       trt = rbinom(n = ni * N, size = 1, prob = 0.5) - .5)
    data$id <- factor(mapply(paste, data$trialref,
                             unlist(lapply(table(data$trialref), function(x)
                               1:x)),
                             #rep(1:ni, N),
                             sep = '.'))
    
    # First stage: random effects ******************************************** #
    # Trial-level
    muS <- log(log(2) / mstS) * (-1) ^ (method == 'mx')
    muT <- log(log(2) / mstT) * (-1) ^ (method == 'mx')
    d_ab <- sqrt(R2 * alphaVar * betaVar)
    d_ST <- sqrt(baseCorr) * sqrt(prod(baseVars))
    Sigma_trial <-
      matrix(c(
        baseVars[1],
        d_ST,
        0,
        0,
        d_ST,
        baseVars[2],
        0,
        0,
        0,
        0,
        alphaVar,
        d_ab,
        0,
        0,
        d_ab,
        betaVar
      ),
      4)
    #   library('MASS')
    pars <- mvrnorm(n = N,
                    mu = c(muS, muT, alpha, beta),
                    Sigma = Sigma_trial)
    pars <- pars[data$trialref,]
    
    ATTRs <- list(
      'N' = N,
      'ni' = ni,
      'gammaWei' = gammaWei,
      'censorT' = ifelse(missing(censorT), NA, censorT),
      'censorA' = ifelse(missing(censorA), NA, censorA),
      'baseCorr' = baseCorr,
      'baseVars' = baseVars,
      'alphaVar' = alphaVar,
      'betaVar' = betaVar,
      'mstS' = mstS,
      'mstT' = mstT,
      'alpha' = alpha,
      'beta' = beta,
      'pars' = pars,
      'R2' = R2
    )
    
    if (method == 'mx') {
      ATTRs$indCorr <- indCorr
      
      # Second stage: survival times ***************************************** #
      # Y, the truncated normal random variables
      #   library('msm')
      Y <- rtnorm(n = ni * N, lower = 0)
      if (indCorr) {
        Y <- cbind(Y, Y)
      } else {
        Y <- cbind(Y, rtnorm(n = ni * N, lower = 0))
      }
      
      # lambda, the exponential random variables
      lambdaS <- rexp(n = ni * N, rate = 1)
      lambdaT <- rexp(n = ni * N, rate = 1)
      
      # S and T, the times
      deltaS <- exp(pars[, 1] + pars[, 3] * data$trt)
      deltaT <- exp(pars[, 2] + pars[, 4] * data$trt)
      data$S <-
        deltaS * (Y[, 1] * sqrt(2 * lambdaS)) ^ (1 / gammaWei[1])
      data$T <-
        deltaT * (Y[, 2] * sqrt(2 * lambdaT)) ^ (1 / gammaWei[2])
      # ********************************************************************** #
      
      ATTRs$kTau <- mean(cor(data$S[data$trt == -.5], data$T[data$trt == -.5], method =
                               'kendall'),
                         cor(data$S[data$trt == 0.5], data$T[data$trt == 0.5], method = 'kendall'))
    } else {
      indCorr <- NULL
      ATTRs$kTau <- kTau
      
      if (method == 're') {
        # Individ2ual-level
        sigma2 <- find.sigma2(kTau)
        indFrailty <- rnorm(nrow(data), mean = 0, sd = sqrt(sigma2))
        # ******************************************************************** #
        
        # Second stage: survival times *************************************** #
        # Weibull hazard:
        # h(t) = lambda gamma x^(gamma-1)
        lambda.S <-
          exp(pars[, 1] + pars[, 3] * data$trt + indFrailty)
        lambda.T <-
          exp(pars[, 2] + pars[, 4] * data$trt + indFrailty)
        
        # Weibull hazard, parametrization as in rweibull:
        # h(t) = shape scale^(-shape) x^(shape-1)
        # shape = gamma
        # scale = lambda^(-1/gamma)
        shape.S <- gammaWei[1]
        scale.S <- lambda.S ^ -(1 / shape.S)
        shape.T <- gammaWei[2]
        scale.T <- lambda.T ^ -(1 / shape.T)
        
        # S and T, the times
        data$S <-
          rweibull(nrow(data), shape = shape.S, scale = scale.S)
        data$T <-
          rweibull(nrow(data), shape = shape.T, scale = scale.T)
        # ******************************************************************** #
      } else {
        # Individ2ual-level
        theta <- 2 * kTau / (1 - kTau)
        # ******************************************************************** #
        
        # Second stage: survival times *************************************** #
        # Weibull hazard:
        # h(t) = lambda gamma x^(gamma-1)
        lambda.S <- exp(pars[, 1] + pars[, 3] * data$trt)
        lambda.T <- exp(pars[, 2] + pars[, 4] * data$trt)
        
        # S times
        US <- runif(nrow(data), 0, 1)
        data$S <- (-log(US) / lambda.S) ^ (1 / gammaWei[1])
        
        # T times | S times
        UT <- runif(nrow(data), 0, 1)
        UT_prime <-
          ((UT ^ (-theta / (1 + theta)) - 1) * US ^ (-theta) + 1) ^ (-1 / theta)
        data$T <- (-log(UT_prime) / lambda.T) ^ (1 / gammaWei[2])
        # ******************************************************************** #
        #     plot(data$S, data$T, log='xy', col=data$trt+1.5)
      }
    }
    
    # Censoring ************************************************************** #
    data$C <- rep(Inf, nrow(data))
    # Random censoring
    if (!missing(censorT)) {
      findK <- Vectorize(function(logK) {
        (mean(runif(10 * length(data$T), 0, exp(logK)) < data$T) - censorT) ^ 2
      })
      suppressWarnings({
        k <- exp(optim(log(max(data$T)), findK)$par)
      })
      data$C <- runif(length(data$T), 0, k)
    }
    # Administrative censoring
    if (!missing(censorA)) {
      data$C <- pmin(data$C, censorA)
    }
    
    data$timeT <- pmin(data$C, data$T)
    data$statusT <- mapply('>=', data$C, data$T) * 1
    data$timeS <- pmin(data$C, data$S)
    data$statusS <- mapply('>=', data$C, data$S) * 1
    # kTau <- cor(data$S, data$T, method='kendall')
    kTau <- mean(cor(data$S[data$trt == -.5], data$T[data$trt == -.5], method =
                       'kendall'),
                 cor(data$S[data$trt == 0.5], data$T[data$trt == 0.5], method = 'kendall'))
    data <- data[, !(names(data) %in% c('C', 'T', 'S'))]
    # ************************************************************************ #
    attributes(data) <- c(attributes(data), ATTRs)
    return(data)
  }
  
  if (method == 'mx') {
    formals(simfun)$kTau <- NULL
  } else {
    formals(simfun)$indCorr <- NULL
  }
  return(simfun)
}

simData.cc <- simData('cc')
simData.re <- simData('re')
simData.mx <- simData('mx')
