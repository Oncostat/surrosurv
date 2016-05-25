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
  convergence(x, kkttol = 1e-8)
}
\arguments{
  \item{x}{
  The fitted models, an object of class \code{\link{surrosurv}}
}
  \item{kkttol}{
  The tolerance for the KKT conditions (see \code{\link[optextras]{kktc}})
    and for positive definitiveness of the VCOV matrix of random effects
    (see \code{\link[matrixcalc]{is.positive.definite}})
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
