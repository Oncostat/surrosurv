library(surrosurv)
data('gastadv')
cvRes <- loocv(gastadv, 'Clayton', nCores = 2)
cvRes
plot(cvRes)
save(cvRes, file = './cvRes.RData')
