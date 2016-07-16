\name{convergence}
\alias{convergence}
\title{
  Assesses the convergence of fitted models for surrogacy evaluation
}
\description{
  This function evaluates whether the fitted models for evaluating 
  the surrogacy of a candidate endpoint have converged.
  Convergence is assessed by checking whether the
  two Kuhn--Karush--Tucker conditions
  (small gradient + positive definite Hessian) are fulfilled
  and whether the estiamted variance-covariance matrix of trial-level random effects
  is definie positive.
}
\usage{
  convergence(x, kkttol = 1e-3, kkt2tol = 1e-6, restricted = FALSE)
}
\arguments{
  \item{x}{
  The fitted models, an object of class \code{\link{surrosurv}}
}
  \item{kkttol, kkt2tol}{
  The tolerance thresholds for the KKT conditions
    (see \code{\link[optextras]{kktc}})
    and for positive definitiveness of the VCOV matrix of random effects
    (see \code{\link[matrixcalc]{is.positive.definite}}).
  The default values \code{kkttol = 1e-3} and \code{kkt2tol = 1e-6}
    are the same values used by \code{\link[optextras]{kktc}}.
}
  \item{restricted}{
  If \code{FALSE} (default), the convergence is checked on all the parameters.
  If \code{TRUE}, the convergence is checked only on parameters
    which contribute to the computation of surrogacy measures.
}
}
\author{
  \packageAuthor{surrosurv}
}

\seealso{\code{\link[optextras]{kktc}}, 
         \code{\link[matrixcalc]{is.positive.definite}}}

\keyword{ convergence }
\keyword{ KKT }
\keyword{ Kuhn--Karush--Tucker conditions }
