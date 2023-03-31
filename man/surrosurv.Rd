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
          models = c('Clayton', 'Plackett', 'Hougaard',
                     'Poisson I', 'Poisson T', 'Poisson TI', 'Poisson TIa'),
          intWidth = NULL,  nInts = 8,
          cop.OPTIMIZER = "bobyqa",
          poi.OPTIMIZER = "bobyqa",
          verbose = TRUE,
          twoStage = FALSE,
          keep.data = TRUE)

\method{predict}{surrosurv}(object, models = names(object), exact.models, ...)

\method{print}{surrosurv}(x, silent = FALSE, 
      digits = 2, na.print = "-.--", ...)
          
\method{print}{predictSurrosurv}(x, n = 6, ...)
          
\method{plot}{surrosurv}(x, ...)

\method{plot}{predictSurrosurv}(x, models = names(x), exact.models,
                      pred.ints = TRUE,
                      show.ste = TRUE,
                      surro.stats = TRUE, 
                      xlab, ylab, 
                      xlim, ylim, mfrow, main, ...)
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
  For \code{surrosurv()}, the models should be fitted/plotted/predicted.
  Possible models are:
    Clayton copula (unadjusted and adjusted),
    Plackett copula (unadjusted and adjusted),
    Hougaard copula (unadjusted and adjusted),
    Poisson (with individual-level heterogeneity only,
              with trial-level heterogeneity only,
              with both individual- and trial-level heterogeneity,
              with both individual- and trial-level heterogeneity 
              and with random per-trial intercept).
  }
  \item{exact.models}{
    If \code{TRUE}, plots or predictions are generated only
      for the elements of \code{x} which match exactly any of \code{models}.
    If \code{exact.models = TRUE}, partial matching is used.
    By default, \code{exact.models = TRUE} if all the \code{models}
      match exactly any of the \code{names(x)} (or  \code{names(object)})
      and \code{exact.models = FALSE} otherwise.
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
    the optimizer for copula models (see \code{\link{optimx}})
  }
  \item{poi.OPTIMIZER}{
    the optimizer for Poisson models (see \code{\link{optimx}})
  }
  \item{verbose}{
    should the function print out the model being fitted
  }
  \item{twoStage}{
    should the parameters of the baseline hazard functions
    fixed to their marginal estimates (Shih and Louis, 1995)
  }
  \item{keep.data}{
    should the data object be kept as attribute of the returned results?
    (this is needed for \code{confint.surrosurv()})
  }

  \item{x, object}{
    The fitted models, an object of class \code{\link{surrosurv}}
  }
  \item{silent}{
    Should the results be return for storing without printing them?
  }
  \item{digits, na.print, xlab, ylab, xlim, ylim, main, ...}{
    other parameters for \code{\link[base]{print}} or \code{\link[graphics]{plot}}
  }
  \item{mfrow}{
    the number of rows and columns for displaying the plots
    (see \code{\link[graphics]{par}}).
    If missing, the default is computed using the function
      \code{\link[grDevices]{n2mfrow}}
  }
  \item{n}{
    the number of rows to print
  }
  \item{pred.ints}{
    Should the prediction intervals be plotted?
  }
  \item{show.ste}{
    Should the surrogate threshold effect be showed?
  }
  \item{surro.stats}{
    Should the surrogacy statistics be showed?
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
  
  The mixed Poisson models that can be fit
  are used to estimate parameters of mixed proportional hazard models,
  as described for instance by Crowther et al (2014).
  The statistical details are provided in Rotolo et al (WP).
  
  The function \code{predict()} returns the estimated values
    of the log-hazard ratios on the true and the surrogate endpoints.
  The list of the prediction functions (for all the \code{models})
    is available as \code{attr(predict.surrosurv(...), 'predf')}.
}

\value{
  The fitted models, an object of class \code{\link{surrosurv}}.
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
  A Poisson approach for the validation of failure time surrogate endpoints
  in individual patient data meta-analyses.
  \emph{Statistical Methods in Medical Research} 2017;
  \bold{In Press}.
  \doi{10.1177/0962280217718582}

  Shih JH, Louis TA.
  Inferences on the Association Parameter in Copula Models
  for Bivariate Survival Data.
  \emph{Biometrics} 1995; \bold{51}:1384--1399.
  \doi{10.2307/2533269}
  
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
  data <- simData.re(N = 20, ni = 250,
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
    # Long computation time!
    surrores <- surrosurv(data, verbose = TRUE)
    convergence(surrores)
    surrores
  }
  
  # Advanced GASTRIC data
  \dontrun{
    # Long computation time!
    data('gastadv')
    allSurroRes <- surrosurv(gastadv, c('Clayton', 'Poisson'), verbose = TRUE)
    convergence(allSurroRes)
    allSurroRes
    predict(allSurroRes)
    plot(allSurroRes)
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
