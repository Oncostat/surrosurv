

test_that("Data", {
  data(gastadv)
  expect_snapshot({
    nrow(gastadv)
    names(gastadv)
  })
})

#takes around 6 min to fit
test_that("surrosurv - lite", {
  skip_on_cran()
  data(gastadv)
  # dplyr::count(gastadv, trialref, trt, sort=TRUE)
  nrow(df) #400
  
  tictoc::tic('one')
  set.seed(42)
  df = dplyr::slice_head(gastadv, n=10, by=c(trialref, trt))
  allSurroRes <- surrosurv(df, verbose=TRUE) #150s
  loocvRes <- loocv(df, models=c('Clayton', 'PoissonTI'), nCores=5) #205s
  tictoc::toc() #205
  
  saveRDS(allSurroRes, test_path("allSurroRes_lite.rds"))
  saveRDS(loocvRes, test_path("loocvRes_lite.rds"))
  if(FALSE){
    allSurroRes = readRDS(test_path("allSurroRes_lite.rds"))
    loocvRes = readRDS(test_path("loocvRes_lite.rds"))
  }
  
  expect_snapshot({
    allSurroRes
    convergence(allSurroRes)
    convals(allSurroRes)
    predict(allSurroRes, models='PoissonTI')
    # plot(allSurroRes, c('Clayton adj', 'PoissonTI'))
    ste(allSurroRes)
    
    loocvRes
    # plot(loocvRes)
  })
})



test_that("Data poissonization", {
  skip("TODO")
  data(gastadv)
  gastadv.poi <- gastadv
  gastadv.poi$time <- gastadv.poi$timeT / 365.25
  gastadv.poi$status <- gastadv.poi$statusT
  
  fitcox <- coxph(Surv(time, status) ~ trt, data=gastadv.poi)
  cox.base <- basehaz(fitcox, centered = FALSE)
  
  
  plot(
    stepfun(cox.base$time [-nrow (cox.base)], 
            exp(-cox.base$hazard)),
    ylim = 0:1 ,
    yaxs = 'i',
    xlim = c (0 , 5) ,
    xaxs = 'i',
    col = 1 ,
    lwd = 2 ,
    bty = 'l',
    do.points = FALSE ,
    verticals = FALSE ,
    main = 'Overall Survival\nAdvanced GASTRIC meta - analysis ',
    xlab = 'Years ', ylab = 'Survival probability '
  )
  
  lines(
    stepfun(cox.base$time[-nrow(cox.base)], 
            exp(-cox.base$hazard*exp(coef(fitcox)['trt']))) ,
    col = 2 ,
    pch = '',
    lwd = 2
  )
  
})
