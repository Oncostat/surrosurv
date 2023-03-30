# code from Georg Heinze used for starting values, originally in survcorr package . Ref: Schemper M, Kaider A, Wakounig S, Heinze G. Estimating the correlation of bivariate failure times under censoring. Stat Med. 2013 Nov 30;32(27):4781-90. doi: 10.1002/sim.5874. PMID: 23775542.


## Helper function for bivariate pearson correlation.
pearson = function(data) {
  cormatrix = cor(data, method="pearson")
  cormatrix[1, 2]
}

survcorr = function
(
  formula1,  # Formula for defining time-to-event 1, e.g. Surv(TIME1, STATUS1) ~ 1.
  formula2,  # Formula for defining time-to-event 2, =Surv(TIME2, STATUS2) ~ 1.
  data,    # Data set 
  methods="imi",  # Vector of method names. Allowed are "imi" and "copula", but copula not yet implemented.
  alpha=0.05,     # Confidence level; interval is [alpha, 1-alpha].
  intra=FALSE,    # TRUE if the survival times are symmetric and interchangeable.
  M=10,           # Imputation size.
  MCMCSteps=10,   # Number of MCMC steps for double-censored observations.
  epsilon= 0.001, # Precision epsilon.
  maxiter=100     # Maximum number of iterations.
)
### IMI (Iterative Multiple Imputation)
### MP, 2013-12
{
  ## Extract times ('t') and events ('delta') of data set.
  ## For <intra=T> mirror the data. 
  tmp1 = as.matrix(model.frame(formula1, data))
  if (intra) {
    tmp = as.matrix(model.frame(formula2, data))
    tmp2 = rbind(tmp, tmp1)
    tmp1 = rbind(tmp1, tmp)
  }
  else
    tmp2 = as.matrix(model.frame(formula2, data))
  
  colnames(tmp1) = colnames(tmp2) = c("time", "status")
  tmp1 = as.data.frame(tmp1)
  tmp2 = as.data.frame(tmp2)
  
  t1 = as.vector(tmp1$time)
  delta1 = as.vector(tmp1$status)
  order1 = order(-t1)
  
  t2 = as.vector(tmp2$time)
  delta2 = as.vector(tmp2$status)
  order2 = order(-t2)
  
  n = length(t1)
  minTime = 0.00000001
  
  ## Fill data set. 
  data1 = data.frame(index=1:n, t1=t1, delta1=delta1, t1unif=rep(NA, n), z1=rep(NA, n))[order1, ]
  obj = survfit(coxph(Surv(time, status) ~ 1, data=tmp1, model=TRUE), type="aalen", se.fit=FALSE, conf.type="none")
  # Get one row per object.
  data1$t1unif = rev(rep(obj$surv, obj$n.event + obj$n.censor))
  data1$t1unif = pmin(pmax(data1$t1unif, minTime), 1 - minTime)
  data1$z1 = qnorm(data1$t1unif)
  
  data2 = data.frame(index=1:n, t2=t2, delta2=delta2, t2unif=rep(NA, n), z2=rep(NA, n))[order2, ]
  obj = survfit(coxph(Surv(time, status) ~ 1, data=tmp2, model=TRUE), type="aalen", se.fit=FALSE, conf.type="none")
  # Get one row per object.
  data2$t2unif = rev(rep(obj$surv, obj$n.event + obj$n.censor))
  data2$t2unif = pmin(pmax(data2$t2unif, minTime), 1 - minTime)
  data2$z2 = qnorm(data2$t2unif)
  
  ##d = cbind(data1[order(data1$index), -1], data2[order(data2$index), -1])
  ##print(d)
  
  z1orig = data1$z1[order(data1$index)]
  z2orig = data2$z2[order(data2$index)]
  
  ## 1) Get initial correlation coefficient estimate.
  r0 = pearson(cbind(z1orig, z2orig)[delta1 & delta2, , drop=F])
  rj.all = rep(NA, M)
  
  ## Calc. indicators for censoring.
  ind1 = !delta1 & delta2
  ind2 = delta1 & !delta2
  indBoth = !delta1 & !delta2
  nBoth = sum(indBoth)
  
  z1M = z2M = matrix(NA, n, M)
  for (iImp in 1:M) {
    rj = r0
    
    ## Generate random data for 2).
    runif1 = runif(n)
    runif2 = runif(n)
    runif1mcmc = matrix(runif(nBoth*MCMCSteps), nBoth, MCMCSteps)
    runif2mcmc = matrix(runif(nBoth*MCMCSteps), nBoth, MCMCSteps)
    
    for (iter in seq(length=maxiter)) {
      ## Reset the z's
      z1 = z1orig
      z2 = z2orig
      
      ## a) Only z1i is censored.
      p1 = 1 - pnorm(z1, mean=rj * z2, sd=sqrt(1 - rj^2))     # n
      unif1 = 1 - pmin(pmax(p1 + runif1 * (1 - p1), minTime), 1 - minTime)     # n
      z1[ind1] = qnorm(unif1, mean=rj * z2, sd=sqrt(1 - rj^2))[ind1]
      
      ## b) Only z2i is censored.
      p2 = 1 - pnorm(z2, mean=rj * z1, sd=sqrt(1 - rj^2))     # n
      unif2 = 1 - pmin(pmax(p2 + runif2 * (1 - p2), minTime), 1 - minTime)     # n
      z2[ind2] = qnorm(unif2, mean=rj * z1, sd=sqrt(1 - rj^2))[ind2]
      
      ## c) Both are censored.
      z1new = z1[indBoth]
      z2new = z2[indBoth]
      for (MCMCStep in seq(length=MCMCSteps)) {
        p1 = 1 - pnorm(z1[indBoth], mean=rj * z2new, sd=sqrt(1 - rj^2))     # n
        unif1 = 1 - pmin(pmax(p1 + runif1mcmc[, MCMCStep] * (1 - p1), minTime), 1 - minTime)     # n
        z1new = qnorm(unif1, mean=rj * z2new, sd=sqrt(1 - rj^2))
        
        p2 = 1 - pnorm(z2[indBoth], mean=rj * z1new, sd=sqrt(1 - rj^2))     # n
        unif2 = 1 - pmin(pmax(p2 + runif2mcmc[, MCMCStep] * (1 - p2), minTime), 1 - minTime)     # n
        z2new = qnorm(unif2, mean=rj * z1new, sd=sqrt(1 - rj^2))
      }
      z1[indBoth] = z1new
      z2[indBoth] = z2new
      
      ## Save old correlations.
      rj.old = rj
      
      ## 3) Calculate the new correlation.
      rj = pearson(cbind(z1, z2))
      rj.all[iImp] = rj
      
      ## 5) Stop condition: all correlations converge.
      if (abs(rj - rj.old) < epsilon) {
        z1M[, iImp] = z1
        z2M[, iImp] = z2
        break
      }
    }
  }
  
  ## 6) Calculate point and interval estimates.
  n = nrow(data)  # Above n was twice as large for intra=TRUE.
  if (intra)
    df0 = n - 3/2
  else
    df0 = n - 3
  rj.trans = atanh(rj.all)
  rj.t.mean = mean(rj.trans)
  rj.t.var = 1 / df0
  BetwVar = var(rj.trans)
  TotalVar = rj.t.var + BetwVar * (M + 1) / M
  r.W.hat = tanh(rj.t.mean)
  
  df = (M - 1) * (1 + M / (BetwVar * (M + 1) * df0))^2 
  limits = tanh(rj.t.mean + c(-1, 1) * qt(1 - alpha/2, df) * sqrt(TotalVar))
  
  ## Create an object of class <survcorr>
  simData = list(M=M, z1M=z1M, z2M=z2M, delta1=delta1, delta2=delta2, t1=t1, t2=t2)
  obj = list(rho=r.W.hat, ci.lower=limits[1], ci.upper=limits[2], simData=simData, M=M, MCMCSteps=MCMCSteps, rj.trans=rj.trans, rj.t.mean=rj.t.mean, 
             var=c(within=rj.t.var, between=BetwVar, total=TotalVar), df=df, intra=intra, alpha=alpha, call=match.call())
  class(obj) = "survcorr"
  
  obj
}
