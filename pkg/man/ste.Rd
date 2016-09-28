\name{ste}
\alias{ste}
\alias{print.steSurrosurv}
\title{
  Surrogate threshold effect
}
\description{
  The function \code{ste()} computes the surrogate threshold effect (STE)
    of a .
}
\usage{
ste(x, models = names(x), exact.models)

\method{print}{steSurrosurv}(x, digits = 2, ...)
}

\arguments{
  \item{x}{
    The fitted models, an object of class \code{\link{surrosurv}}
  }
  \item{models, exact.models}{
    Which models should be fitted (see \code{\link{surrosurv}()})
  }
  \item{digits}{
    the number of digits
  }
  \item{...}{
    Further parameters to be passed to 
    the generic \code{print()} function
  }
}

\value{
  An object of class \code{steSurrosurv} 
  containing, for each trial:
  \item{}{
    
  }
  \item{}{
    
  }
}
  
\references{
  Burzykowski T, Buyse M.
    Surrogate threshold effect: an alternative measure
    for meta-analytic surrogate endpoint validation.
    \emph{Pharm Stat}. 2006;5(3):173-86.
    doi:\href{http://dx.doi.org/10.1002/pst.207}{10.1002/pst.207}
}

\author{
  \packageAuthor{surrosurv}
}

\examples{
  \donttest{
  # Possibly long computation time!
  data('gastadv')
  mod <- surrosurv(gastadv, 'Clayton')
  ste(mod)
  }
}

\keyword{ surrogate threshold effect }
\keyword{ ste }
