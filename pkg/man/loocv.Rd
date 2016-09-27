\name{loocv}
\alias{loocv}
\alias{loocv.surrosurv}
\alias{print.loocvSurrosurv}
\alias{plot.loocvSurrosurv}
\title{
  Performs leave-one-trial-out cross-validation
}
\description{
  The function \code{loocv} ...
}
\usage{
loocv(data, nCores, parallel = TRUE, ...)

\method{print}{loocvSurrosurv}(x, silent = FALSE, 
      digits = 2, na.print = "-.--", ...)
          
\method{plot}{loocvSurrosurv}(x, ...)
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
  \item{nCores}{
  The number of cores for parallel computing
}
  \item{parallel}{
  Should results be computed using parallelization?
}
  \item{...}{
  Further parameters to be passed to \code{surrosurv}
   (for \code{loocv()}, in particular \code{models})
   or to the generics \code{print()} and \code{plot()}.
              
  For \code{plot.surrosurv()}, 
    the models for which plots must be drawn.
}
  \item{x}{
    The fitted models, an object of class \code{\link{surrosurv}}
}
  \item{silent}{
    Should the results be return for storing without printing them?
}
  \item{digits, na.print}{
    parameters to be passed to the generic \link{print}
}
}

\details{
  Three
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
  Michiels S, Le Maître A, Buyse M, et al.
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
  \dontrun{
  data('gastadv')
  cvRes <- loocv(gastadv)
  cvRes
  plot(cvRes)
  }
}

\keyword{ cross-validation }
\keyword{ leave-one-out }
