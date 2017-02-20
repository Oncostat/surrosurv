\name{loocv}
\alias{loocv}
\alias{loocv.surrosurv}
\alias{print.loocvSurrosurv}
\alias{plot.loocvSurrosurv}
\title{
  Leave-one-trial-out cross-validation for treatment effect prediction
}
\description{
  The function \code{loocv()} computed leave-one-out predcition of the treatment
    effect on the true endpoint for each trial,
    based on the observed effect on the surrogate endpoint in the trial itself
    and based on the meta-analytic model fitted on the remaining trials
    (Michiels et al, 2009).
}
\usage{
\method{loocv}{surrosurv}(object, nCores, parallel = TRUE, 
                             models = c('Clayton', 'Plackett', 'Hougaard',
                                        'Poisson I', 'Poisson T',
                                        'Poisson TI', 'Poisson TIa'), ...)

\method{print}{loocvSurrosurv}(x, n = 6, silent = FALSE, ...)
          
\method{plot}{loocvSurrosurv}(x, models, exact.models, ...)
}

\arguments{
  \item{object}{
    Either an object of class \code{surrosurv} with an attribute \code{data}
      of class \link{data.frame} or a \link{data.frame} with columns
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
  \item{nCores}{
    The number of cores for parallel computing
  }
  \item{parallel}{
    Should results be computed using parallelization?
  }
  \item{models, exact.models}{
    Which models should be fitted (see \code{\link{surrosurv}()})
  }
  \item{x}{
    The fitted models, an object of class \code{\link{surrosurv}}
  }
  \item{n}{
    the number of rows to print
  }
  \item{silent}{
    Should the results be return for storing without printing them?
  }
  \item{...}{
    Further parameters to be passed to \code{surrosurv}
     (for \code{loocv()})
     or to the generics \code{print()} and \code{plot()}
  }
}

\value{
  An object of class \code{loocvSurrosurv} 
  containing, for each trial:
  \item{margPars}{
    the observed treatment effects 
    on the surrogate ednpoint (\code{alpha}) and
    on the true endpoint (\code{beta})
  }
  \item{...}{
    for each method in \code{models}
    the predicted value and prediction interval for \code{beta}.
  }
}
  
\references{
  Michiels S, Le Maitre A, Buyse M, et al.
    Surrogate endpoints for overall survival
    in locally advanced head and neck cancer:
    meta-analyses of individual patient data.
    \emph{Lancet Oncol}. 2009;10(4):341-50.
    doi:\href{http://dx.doi.org/10.1016/S1470-2045(09)70023-3}{10.1016/S1470-2045(09)70023-3}
}

\author{
  \packageAuthor{surrosurv}
}

\examples{
  \donttest{
  # Possibly long computation time!
  data('gastadv')
  cvRes <- loocv(gastadv)
  cvRes
  plot(cvRes)
  }
}

\keyword{ cross-validation }
\keyword{ leave-one-out }
