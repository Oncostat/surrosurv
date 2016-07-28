\name{convergence.surrosurv}
\alias{convergence.surrosurv}
\alias{convals}
\title{
  Assesses the convergence of fitted models for surrogacy evaluation
}
\description{
  This function evaluates whether the fitted models for evaluating 
    the surrogacy of a candidate endpoint have converged.
  Convergence is assessed by checking whether the maximum absolute gradient 
    is small enough, and whether the Hessian matrix
    and the variance-covariance matrix of random treatment effects
    are definie positive.
}
\usage{
  convals(x)
  \method{convergence}{surrosurv}(x, kkttol = 1e-3, kkt2tol = 1e-6)
}
\arguments{
  \item{x}{
    The fitted models, an object of class \code{\link{surrosurv}}.
  }
  \item{kkttol}{
    The tolerance thresholds for the assessing whether the maximum absolute
      gradient is small enough.
  }
  \item{kkt2tol}{
    The tolerance thresholds for checking whether
      the Hessian matrix and
      the variance-covariance matrix of random treatment effects
      are positive definite.
    The threshold is for the minimum of the eigenvalues.
  }
}

\value{
  The function \code{convals()} returns a matrix with one row per model and three columans,
    reporting the values of the maximum absolute gradient, 
    of the minimum eigenvalue of the Hessian matrix, and
    of the minimum eigenvalue of the estimated variance-covariance matrix of random treatment effects.
  The function \code{convergence()} returns a matrix with the same structure as \code{convals()},
    with \code{TRUE}/\code{FALSE} values for the test of the results of \code{convals()}
    against the given thresholds \code{kkttol} and \code{kkt2tol}.
}
\author{
  \packageAuthor{surrosurv}
}


\keyword{ convergence }
\keyword{ KKT }
\keyword{ Kuhn--Karush--Tucker conditions }
