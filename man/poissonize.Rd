\name{poissonize}
\alias{poissonize}
\alias{plotsson}
\title{
  Transform survival data for fitting a Poisson model
}
\description{
  This function transform survival data into a format compatible with
  the \code{\link{glm}()} function for fitting an auxiliary Poisson model,
  providing the parameter estimates of the associated proportional hazard model.
}
\usage{
poissonize(data,
           all.breaks = NULL, interval.width = NULL, nInts = 8,
           factors = NULL, compress = TRUE)
plotsson(x, type = c('survival', 'hazard'),
         add = FALSE, xscale = 1, by, col, ...)
}
\arguments{
  \item{data}{
    a data frame with columns:
    \itemize{
     \item id       : the patient identifier
     \item time     : the event/censoring time
     \item status   : the event(1) or censoring(0) indicator
     \item ...      : other factors such like the covariables needed
      in the regression model   
    }
}
  \item{all.breaks}{
    the breakpoints between time intervals
}
  \item{interval.width}{
    the width of the time intervals on which the risks will be assumed constant,
    in case of intervals of the same length.
    This parameter is ignored if \code{all.breaks} is specified
}
  \item{nInts}{
    the number of intervals containing the same expected number of events
    (used only if \code{is.null(interval.width)}, see Details).
    This parameter is ignored if either \code{all.breaks}
    or \code{interval.width} is specified
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
  \item{by}{
  covariate for which a different curve per level has to be plotted
}
  \item{col, \dots}{
  other \link[=par]{graphical parameters}
}
}

\details{
  If \code{interval.width} is not null, the study period is divided into
    equal-length intervals of length \code{interval.width}.
  Otherwise, \code{nInts} intervals are used, and the location of their bounds
    is computed based on the empirical quantiles of the survival function.
}
\references{
  Whitehead, J.
    Fitting Cox's regression model to survival
    data using GLIM.
    \emph{J Roy Stat Soc C Appl Stat} 1980; \bold{29}(3):268-275.
    \url{https://www.jstor.org/stable/2346901}.
  
  Crowther MJ, Riley RD, Staessen JA, Wang J, Gueyffier F, Lambert PC.
    Individual patient data meta-analysis of survival data
    using Poisson regression models.
    \emph{BMC Medical Research Methodology} 2012; \bold{12}:34.
    \doi{10.1186/1471-2288-12-34}.
}
\author{
  \packageAuthor{surrosurv}
}
\note{
  This code is hugely inspired by original code made publicly available by Stephanie Kovalchik.
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
    head(wdata1 <- poissonize(kidney, interval.width = int, 
                              factors = c('disease'), compress = FALSE))
    head(wdata2 <- poissonize(kidney, interval.width = int, 
                              factors = c('disease'), compress = TRUE))
    
    fitcox <- (coxph(Surv(time, status) ~ disease, data = kidney))
    fitpoi1 <- glm(event ~ -1 + interval + disease + offset(log(time)),
                   data = wdata1, family = 'poisson')
    fitpoi2 <- glm(m ~ -1 + interval + offset(log(Rt)) + disease,
                   data = wdata2, family = 'poisson')
    cox.base <- basehaz(fitcox, centered = FALSE)
    plot(stepfun(cox.base$time[-nrow(cox.base)], exp(-cox.base$hazard)),
         ylim = 0:1, xlim = c(0, max(cox.base$time)),
         do.points = FALSE, verticals = FALSE, xaxs = 'i',
         main = paste0('KIDNEY data set\nInterval width = ', int),
         xlab = 'Time', ylab = 'Survival probability')
    plotsson(fitpoi1, 'Surv', add = TRUE, col = 2, lty = 2)
    plotsson(fitpoi2, 'Surv', add = TRUE, col = 3, lty = 3) 
    legend('topright', col = 1:3, lty = 1:3,
           legend = c('Breslow (Cox)', 'Poisson',
                      'Poisson (compressed dataset)'))
}
print(cbind(Cox                = coef(fitcox),
            Poisson            = rev(rev(coef(fitpoi1))[1:3]),
            Poisson_Compressed = rev(rev(coef(fitpoi2))[1:3])), digits = 2)



################################################################################
#  Example 2 - COLON data                                                      #
################################################################################
library(survival)
data(colon)
head(wdata1 <- poissonize(subset(colon, etype == 1), interval.width = 365.25, 
                          factors=c('surg', 'sex', 'age'), compress = FALSE))
head(wdata2 <- poissonize(subset(colon, etype == 1), interval.width = 365.25, 
                          factors=c('surg', 'sex', 'age'), compress = TRUE))

fitcox <- coxph(Surv(time, status) ~ surg + sex + age, 
                data = subset(colon, etype == 1))

system.time({
    fitpoi1 <- glm(event ~ -1 + interval + surg + sex + age + offset(log(time)),
                   data = wdata1, fam = 'poisson')
})
system.time({
    fitpoi2 <- glm(m ~ -1 + interval + offset(log(Rt)) + surg + sex + age,
                   data = wdata2, family = 'poisson')
})
{
    cox.base <- basehaz(fitcox, centered = FALSE)
    par(mfrow = c(1, 1))
    plot(stepfun(cox.base$time[-nrow(cox.base)], exp(-cox.base$hazard)),
         ylim = 0:1, xlim = c(0, max(cox.base$time)),
         do.points = FALSE, verticals = FALSE, xaxs = 'i',
         main = 'COLON data set', xlab = 'Time', ylab = 'Survival probability')
    plotsson(fitpoi1, 'Surv', add = TRUE, col = 2, lty = 2) 
    plotsson(fitpoi2, 'Surv', add = TRUE, col = 3, lty = 3) 
    legend('topright', col = 1:3, lty = 1:3,
           legend = c('Cox', 'Poisson', 'Poisson (compressed dataset)'))
}
print(cbind(Cox                = coef(fitcox),
            Poisson            = rev(rev(coef(fitpoi1))[1:3]),
            Poisson_Compressed = rev(rev(coef(fitpoi2))[1:3])), digits = 2)



################################################################################
#  Example 3 - LUNG data                                                       #
################################################################################
library(survival)
data(lung)
lung$status <- lung$status - 1
lung$id <- 1:nrow(lung)
head(wdata1 <- poissonize(lung, interval.width = 365.25/12, 
                          factors = c('pat.karno', 'sex', 'age'), 
                          compress = FALSE))
head(wdata2 <- poissonize(lung, interval.width = 365.25/12,
                          factors = c('pat.karno', 'sex', 'age'),
                          compress = TRUE))

fitcox <- coxph(Surv(time, status) ~ pat.karno + sex + age, data = lung)

system.time({
    fitpoi1 <- glm(event ~ -1 + interval + pat.karno + sex + age + 
                              offset(log(time)),
                   data = wdata1, family = 'poisson')
})
system.time({
    fitpoi2 <- glm(m ~ -1 + interval + pat.karno + sex + age + offset(log(Rt)),
                   data = wdata2, family = 'poisson')
})
{
    cox.base <- basehaz(fitcox, centered = FALSE)
    plot(stepfun(cox.base$time[-nrow(cox.base)], exp(-cox.base$hazard)),
         ylim = 0:1, xlim = c(0, max(cox.base$time)),
         do.points = FALSE, verticals = FALSE, xaxs = 'i',
         main = 'LUNG data set', xlab = 'Time', ylab = 'Survival probability')
    plotsson(fitpoi1, 'Surv', add = TRUE, col = 2, lty = 2) 
    plotsson(fitpoi2, 'Surv', add = TRUE, col = 3, lty = 3) 
    legend('topright', col = 1:3, lty = 1:3,
           legend = c('Cox', 'Poisson', 'Poisson (compressed dataset)'))
}
print(cbind(Cox                = coef(fitcox),
            Poisson            = rev(rev(coef(fitpoi1))[1:3]),
            Poisson_Compressed = rev(rev(coef(fitpoi2))[1:3])), digits = 2)

}
\keyword{ Poisson }
\keyword{ Survival data }
