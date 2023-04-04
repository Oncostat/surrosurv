skip_on_cran()

#Fitting the whole dataset takes 1 hour
#It should be made
#snapshot file "surrosurv_total.md" is 


# 3.1. Fitting the surrogacy models
test_that("surrosurv - total", {
  cache = test_path("allSurroRes.rds")
  if(FALSE){
    tictoc::tic("surrosurv")
    set.seed(42)
    allSurroRes <- surrosurv(gastadv, verbose=TRUE) #1500s = 25min
    tictoc::toc()
    beepr::beep()
    saveRDS(allSurroRes, cache)
  } else {
    allSurroRes = readRDS(cache)
  }
  
  expect_snapshot({
    allSurroRes
    convergence(allSurroRes)
    convals(allSurroRes)
    predict(allSurroRes, models='PoissonTI')
    # plot(allSurroRes, c('Clayton adj', 'PoissonTI'))
    ste(allSurroRes)
  })
  
})


# 3.2.1. Leave-one-trial-out cross-validation
test_that("loocv - total", {
  cache = test_path("loocvRes.rds")
  if(FALSE){
    tictoc::tic("loocv")
    set.seed(42)
    loocvRes <- loocv(gastadv, models=c('Clayton', 'PoissonTI'), nCores=3) #2900s = 50min
    saveRDS(loocvRes, cache)
    tictoc::toc()
    beepr::beep()
  } else {
    loocvRes = readRDS(cache)
  }
  
  
  expect_snapshot({
    loocvRes
    # plot(loocvRes)
  })
  
})