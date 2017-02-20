\name{convergence}
\alias{convergence}
\alias{convals}
\title{
  Assesses the convergence of fitted models for surrogacy evaluation
}
\description{
  This function evaluates whether the fitted models for evaluating 
    the surrogacy of a candidate endpoint have converged.
  Convergence is assessed by checking whether the maximum gradient is small enough,
    and whether the Hessian matrix
    and the variance-covariance matrix of random treatment effects
    are positive definite.
}
\usage{
  convals(x)
  convergence(x, kkttol = 1e-2, kkt2tol = 0)
}
\arguments{
  \item{x}{
    The fitted models, an object of class \code{\link{surrosurv}}.
  }
  \item{kkttol}{
    The tolerance thresholds for the assessing whether
      the maximum scaled gradient is small enough.
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
    reporting the values of the maximum scaled gradient (\code{maxSgrad}), 
    of the minimum eigenvalue of the Hessian matrix (\code{minHev}), and
    of the minimum eigenvalue of the estimated variance-covariance matrix
    of random treatment effects (\code{minREev}).
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
