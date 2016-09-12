\name{surrosurv}
\alias{surrosurv}
\alias{print.surrosurv}
\alias{plot.surrosurv}
\alias{predict.surrosurv}
\alias{plot.predictSurrosurv}
\alias{print.predictSurrosurv}
\title{
  Fit and print the models for evaluating the surrogacy strength 
  of a candidate surrogate endpoint
}
\description{
  The function \code{surrosurv} fits (all or a subset of) statistical
  models to evaluate a surrogate endpoint S for a given true endpoint T,
  using individual data from a meta-analysis of randomized controlled trials.
}
\usage{
surrosurv(data,
          models = c("Clayton", "Plackett", "Hougaard", "Poisson"), 
          intWidth,  nInts = 8,
          cop.OPTIMIZER = "bobyqa",
          poi.OPTIMIZER = "bobyqa",
          verbose = FALSE)

\method{predict}{surrosurv}(object, ...)

\method{print}{surrosurv}(x, silent = FALSE, 
      digits = 2, na.print = "-.--", ...)
          
\method{print}{predictSurrosurv}(x, n = 6, ...)
          
\method{plot}{surrosurv}(x, ...)

\method{plot}{predictSurrosurv}(x, models = names(x), exact.models = TRUE,
                      xlab, ylab, ...)
}

\arguments{
  \item{data}{
  A \link{data.frame} with columns
  \itemize{
    \item \code{trialref}, the trial reference
    \item \code{trt}, the treatment arm (-0.5 or 0.5)
    \item \code{id}, the patient id
    \item \code{timeT}, the value of the true endpoint T
    \item \code{statusT}, the censoring/event (0/1) indicator of the true endpoint T
    \item \code{timeS}, the value of the surrogate endpoint S
    \item \code{statusS}, the censoring/event (0/1) indicator of the surrogate endpoint S
  }
}
  \item{models}{
  For \code{surrosurv()}, which models should be fitted: 
    Clayton copula (unadjusted and adjusted),
    Plackett (unadjusted and adjusted),
    Hougaard (unadjusted and adjusted),
    Poisson (with individual-level heterogeneity only,
              with trial-level heterogeneity only,
              with both individual- and trial-level heterogeneity,
              with both individual- and trial-level heterogeneity 
              and with random per-trial intercept)
              
  For \code{plot.surrosurv()}, 
    the models for which plots must be drawn.
}
  \item{intWidth}{
    the width of time intervals for data Poissonization
    (see \link{poissonize})
}
  \item{nInts}{
    the number of time intervals for data Poissonization
    (see \link{poissonize})
}
  \item{cop.OPTIMIZER}{
    the optimizer for copula models (see \link{optimx})
}
  \item{poi.OPTIMIZER}{
    the optimizer for Poisson models (see \link{optimx})
}
  \item{verbose}{
    should the function print out the model being fitted
}

  \item{x, object}{
    The fitted models, an object of class \code{\link{surrosurv}}
}
  \item{silent}{
    Should the results be return for storing without printing them?
}
  \item{digits, na.print, xlab, ylab, ...}{
    other parameters for \link{print} or \link{plot}
}
  \item{n}{
    the number of row to print
}
  \item{exact.models}{
    If \code{TRUE}, only the elements of \code{x} which match exactly
      any of \code{models} are plotted.
    If \code{exact.models = TRUE}, partial matching is used.
    
}
}

\details{
  Three copula models can be fit: Clayton (1978), Plackett (1965), and
  Hougaard (1986).
  For all of them the linear regression at the second step
  is computed both via simple LS regression
  and via a linear model adjusted for measurement error
  of the log-hazard ratios estimated at the first step.
  This adjusted model is the one described by
  Burzykowski et al. (2001), which relies on the results
  by van Houwelingen et al. (2002).
  
  The moxed Poisson models that can be fit
  are used to estimate parameters of mixed proportional hazard models,
  as described for instance by Crowther et al (2014).
  The statistical details are provided in Rotolo et al (WP).
}
\value{
  The fitted models, an object of class \code{\link{surrosurv}}, 
  with possible elements
}
  
\references{
  Burzykowski T, Molenberghs G, Buyse M et al.
  Validation of surrogate end points in multiple randomized clinical trials
  with failure time end points.
  \emph{Journal of the Royal Statistical Society C} 2001;
  \bold{50}:405--422.
  \doi{10.1111/1467-9876.00244}
  
  Clayton DG.
  A model for association in bivariate life tables
  and its application in epidemiological studies of familial tendency 
  in chronic disease incidence.
  \emph{Biometrika} 1978; \bold{65}:141--151.
  \doi{10.1093/biomet/65.1.141}
  
  Crowther MJ, Riley RD, Staessen JA, Wang J, Gueyffier F, Lambert PC.
  Individual patient data meta-analysis of survival data
  using Poisson regression models.
  \emph{BMC Medical Research Methodology} 2012; \bold{12}:34.
  \doi{10.1186/1471-2288-12-34}.
  
  Gasparrini A, Armstrong B, Kenward MG.
  Multivariate meta-analysis for non-linear and other multi-parameter associations.
  \emph{Statistics in Medicine} 2012; \bold{31}:3821--39.
  \doi{10.1002/sim.5471}
  
  Hougaard P.
  A class of multivariate failure time distributions.
  \emph{Biometrika} 1986; \bold{73}:671--678.
  \doi{10.1093/biomet/73.3.671}
  
  Plackett RL.
  A class of bivariate distributions. 
  \emph{Journal of the America Statistical Association} 1965; \bold{60}:516--522.
  \doi{10.1080/01621459.1965.10480807}
  
  Rotolo F, Paoletti X, Burzykowski T, Buyse M, Michiels S.
  Evaluation of failure time surrogate endpoints in individual patient data
  meta-analyses of randomized clinical trials. A Poisson approach.
  \emph{Working paper}.
  
  van Houwelingen HC, Arends LR, Stijnen T.
  Advanced methods in meta-analysis: multivariate approach and meta-regression.
  \emph{Statistics in Medicine} 2002; \bold{21}:589--624.
  \doi{10.1002/sim.1040}

}

\author{
  \packageAuthor{surrosurv}
}

\examples{
  set.seed(150)
  data <- simData.re(J = 20, Nj = 250,
                     R2 = 0.8, kTau = 0.4,
                     alpha = log(0.95), beta = log(0.85),
                     censorA = 15 * 365.25)
  library(survival)
  par(mfrow = 1:2)
  plot(survfit(Surv(timeS, statusS) ~ trt, data = data), lty = 1:2, 
       xscale = 365.25, main = 'Progression-Free Survival\n(S)', col = 2)
  plot(survfit(Surv(timeT, statusT) ~ trt, data = data), lty = 1:2,
       xscale = 365.25, main = 'Overall Survival\n(T)')
       
  \dontrun{
    surrores <- surrosurv(data, verbose = TRUE)
    surrores
    convergence(surrores)
  }
}

\keyword{ surrogate endpoint }
\keyword{ survival }
\keyword{ meta-analysis }
\keyword{ randomized controlled trial }
\keyword{ copula }
\keyword{ Clayton }
\keyword{ Plackett }
\keyword{ Hougaard }
\keyword{ proportional hazard model }
\keyword{ generalized linear mixed model }
\keyword{ Poisson }
