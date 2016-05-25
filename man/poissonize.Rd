\name{poissonize}
\alias{poissonize}
\alias{plotsson}
\title{
  Tranform survival data for fitting a Poisson model
}
\description{
  This function  tranform survival data into a format compatible with
  the \code{\link{glm}()} function for fitting a Poisson model.
}
\usage{
poissonize(data, interval.width = max(data[, grep("time", names(data))]),
           factors = NULL, compress = TRUE)
plotsson(x, type = "haz", add = FALSE, xscale = 1, ...)
}
\arguments{
  \item{data}{
    a data frame with columns:
    \itemize{
     \item id       : the patient identifyier
     \item time     : the event/censoring time
     \item status   : the event(1) or censoring(0) indicator
     \item ...      : other factors such like the covariables needed
      in the regression model   
    }
}
  \item{interval.width}{
    the width of the time intervals on wich the risks will be assumed consant
}
  \item{factors}{
    a vector of characters, containing the names of the factors
    to be kept in the transformed data set
}
  \item{compress}{
    a logical, indicating whether the record with the same
    factor profile should be summarized into one record,
    i.e. whether the data should be expressed in a short form
}

  \item{x}{
  The fitted Poisson model on the poissonized data
}
  \item{type}{
  the type of plot, either 'haz' for the hazard function
  or 'Surv', for the survival curve
}
  \item{add}{
  should the plot added to the active device?
}
  \item{xscale}{
  scaling factor for the time (x) axis
}
  \item{\dots}{
  other \link[=par]{graphical parameters}
}
}
\references{
  Crowther MJ, Riley RD, Staessen JA, Wang J, Gueyffier F, Lambert PC.
  Individual patient data meta-analysis of survival data
  using Poisson regression models.
  \emph{BMC Medical Research Methodology} 2012; \bold{12}:34.
  \doi{10.1186/1471-2288-12-34}.
}
\author{
  Federico Rotolo
}
\note{
  This code is hugely inspired by original code made publicly available by Stephanie Kovalchik
  \href{http://r.789695.n4.nabble.com/exponential-proportional-hazard-model-td805536.html}{[web link]}
}

\examples{
################################################################################
# Example 1 - KIDNEY data                                                      #
################################################################################
library(survival)
data(kidney)
kidney <- kidney[1:(nrow(kidney)/2)*2,]
head(kidney) 

par(mfrow=c(1, 3))
for (int in c(50, 20, 10)) {
    head(wdata1 <- poissonize(kidney, int, factors=c('disease'), compress=FALSE))
    head(wdata2 <- poissonize(kidney, int, factors=c('disease')))
    
    fitcox <- (coxph(Surv(time, status)~disease, data=kidney)) 
    fitpoi1 <- glm(event~-1+interval+disease+offset(log(time)), data=wdata1, fam=poisson)
    fitpoi2 <- glm(m~-1+interval+offset(log(Rt))+disease, data=wdata2, family=poisson)
    cox.base <- basehaz(fitcox, centered = FALSE)
    plot(stepfun(cox.base$time[-nrow(cox.base)], exp(-cox.base$hazard)),
         ylim=0:1, xlim=c(0, max(cox.base$time)),
         do.points=FALSE, verticals=FALSE, xaxs='i',
         main=paste0('KIDNEY data set\nInterval width=', int),
         xlab='Time', ylab='Survival probability')
    plotsson(fitpoi1, 'Surv', add=TRUE, col=2, lty=2)
    plotsson(fitpoi2, 'Surv', add=TRUE, col=3, lty=3) 
    legend('topright', col=1:3, lty=1:3,
           legend = c('Breslow (Cox)', 'Poisson', 'Poisson (compressed dataset)'))
}
print(cbind(Cox=coef(fitcox),
            Poisson_=rev(rev(coef(fitpoi1))[1:3]),
            Poisson_Compressed=rev(rev(coef(fitpoi2))[1:3])), digits=2)



#################################################################################
#  Example 2 - COLON data                                                       #
#################################################################################
library(survival)
data(colon)
head(wdata1 <- poissonize(subset(colon, etype==1), 365.25, 
                          factors=c('surg', 'sex', 'age'), compress=FALSE))
head(wdata2 <- poissonize(subset(colon, etype==1), 365.25, 
                          factors=c('surg', 'sex', 'age')))

fitcox <- coxph(Surv(time,status)~surg+sex+age, 
                data=subset(colon, etype==1))

system.time({
    fitpoi1 <- glm(event~-1+interval+surg+sex+age+offset(log(time)),
                   data=wdata1, fam="poisson")
})
system.time({
    fitpoi2 <- glm(m~-1+interval+offset(log(Rt))+surg+sex+age,
                   data=wdata2, family=poisson)
})
{
    cox.base <- basehaz(fitcox, centered = FALSE)
    plot(stepfun(cox.base$time[-nrow(cox.base)], exp(-cox.base$hazard)),
         ylim=0:1, xlim=c(0, max(cox.base$time)),
         do.points=FALSE, verticals=FALSE, xaxs='i',
         main='COLON data set', xlab='Time', ylab='Survival probability')
    plotsson(fitpoi1, 'Surv', add=TRUE, col=2, lty=2) 
    plotsson(fitpoi2, 'Surv', add=TRUE, col=3, lty=3) 
    legend('topright', col=1:3, lty=1:3,
           legend = c('Cox', 'Poisson', 'Poisson (compressed dataset)'))
}
print(cbind(Cox=coef(fitcox),
            Poisson_=rev(rev(coef(fitpoi1))[1:3]),
            Poisson_Compressed=rev(rev(coef(fitpoi2))[1:3])), digits=2)



#################################################################################
#  Example 3 - LUNG data                                                       #
#################################################################################
library(survival)
data(lung)
lung$status <- lung$status-1
lung$id <- 1:nrow(lung)
head(wdata1 <- poissonize(lung, 365.25/12, 
                          factors=c('pat.karno', 'sex', 'age'), compress=FALSE))
head(wdata2 <- poissonize(lung, 365.25/12, factors=c('pat.karno', 'sex', 'age')))

fitcox <- coxph(Surv(time,status)~pat.karno+sex+age, data=lung)

system.time({
    fitpoi1 <- glm(event~-1+interval+pat.karno+sex+age+offset(log(time)),
                   data=wdata1, fam="poisson")
})
system.time({
    fitpoi2 <- glm(m~-1+interval+offset(log(Rt))+pat.karno+sex+age,
                   data=wdata2, family=poisson)
})
{
    cox.base <- basehaz(fitcox, centered = FALSE)
    plot(stepfun(cox.base$time[-nrow(cox.base)], exp(-cox.base$hazard)),
         ylim=0:1, xlim=c(0, max(cox.base$time)),
         do.points=FALSE, verticals=FALSE, xaxs='i',
         main='LUNG data set', xlab='Time', ylab='Survival probability')
    plotsson(fitpoi1, 'Surv', add=TRUE, col=2, lty=2) 
    plotsson(fitpoi2, 'Surv', add=TRUE, col=3, lty=3) 
    legend('topright', col=1:3, lty=1:3,
           legend = c('Cox', 'Poisson', 'Poisson (compressed dataset)'))
}
print(cbind(Cox=coef(fitcox),
            Poisson_=rev(rev(coef(fitpoi1))[1:3]),
            Poisson_Compressed=rev(rev(coef(fitpoi2))[1:3])), digits=2)

}
\keyword{ Poisson }
\keyword{ Survival data }
