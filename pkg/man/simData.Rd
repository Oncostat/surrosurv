\name{simData}
\alias{simData.re}
\alias{simData.cc}
\alias{simData.mx}
\title{
  Generate survival times for two endpoints in a meta-analysis of randomized trials
}
\description{
  Data are generated from a mixed proportional hazard model,
  a Clayton copula model (Burzykowski and Cortinas Abrahantes, 2005),
  or a mixture of half-normal and exponential random variables
  (Shi et al., 2011).
}
\usage{
simData.re(R2 = 0.6, N = 30, ni = 200,
           nifix = TRUE, gammaWei = c(1, 1), censorT, censorA, 
           kTau= 0.6, baseCorr = 0.5, baseVars = c(0.2, 0.2),
           alpha = 0, beta = 0,
           alphaVar = 0.1, betaVar = 0.1,
           mstS = 4 * 365.25, mstT = 8 * 365.25)
           
simData.cc(R2 = 0.6, N = 30, ni = 200,
           nifix = TRUE, gammaWei = c(1, 1), censorT, censorA,
           kTau= 0.6, baseCorr = 0.5, baseVars = c(0.2, 0.2),
           alpha = 0, beta = 0,
           alphaVar = 0.1, betaVar = 0.1,
           mstS = 4 * 365.25, mstT = 8 * 365.25)
           
simData.mx(R2 = 0.6, N = 30, ni = 200,
           nifix = TRUE, gammaWei = c(1, 1), censorT, censorA,
           indCorr = TRUE, baseCorr = 0.5, baseVars = c(0.2, 0.2),
           alpha = 0, beta = 0,
           alphaVar = 0.1, betaVar = 0.1,
           mstS = 4 * 365.25, mstT = 8 * 365.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R2}{
  The desired trial-level surrogacy \eqn{R^2}
}
  \item{N}{
  The number of trials
}
  \item{ni}{
  The (fixed or average) number of patients per trial
}
  \item{nifix}{
  Should all trials have the same size (if \code{nifix = TRUE})
    of should the \code{N * ni} patients be randomly assigned to trials
    with random probabilities (if \code{nifix = FALSE})?
}
  \item{gammaWei}{
    The shape parameter(s) of the Weibull distributions.
    Either one or two values. If one value is provided, it is used for both endpoints
}
  \item{censorT}{
    censoring rate for the true endpoint T
    (before adding administrative censoring)
}
  \item{censorA}{
    administrative censoring at time censorA
}
  \item{kTau}{
    The desired individual-level dependence between S and T (Kendall's tau)
}
  \item{indCorr}{
    Should S and T be correlated or not? (for \code{.mx} method)
}
  \item{baseCorr}{
    correlation between baseline hazards (\eqn{\rho_{basehaz}}{rho_basehaz})
}
  \item{baseVars}{
    variances of baseline random effects (S and T)
}
  \item{alpha}{
    average treatment effect on S
}
  \item{beta}{
    average treatment effect on T
}
  \item{alphaVar}{
    variance of \eqn{a_i} (\eqn{\theta_a^2}{theta_a^2})
}
  \item{betaVar}{
    variance of \eqn{b_i} (\eqn{\theta_b^2}{theta_b^2})
}
  \item{mstS}{
    median survival time for S in the control arm
}
  \item{mstT}{
    median survival time for T in the control arm
}
}
\details{
  
  The function \code{simData.re} generates data from a proportional hazard model
  with random effects at individual level and 
  random effects and random treatment effects at trial level.
  Individual dependence can be tuned in terms of Kendall's \eqn{tau}
  (\code{kTau}).
  
  The function \code{simData.cc} generates data from a Copula function
  as shown by Burzykowski and Cortinas Abrahantes (2005).
  Individual dependence can be tuned in terms of Kendall's \eqn{tau}
  (\code{kTau}).
  
  The function \code{simData.mx} implements the simulation method by Shi et al. (2011).
  This model is based on a mixture of half-normal and exponential random variables.
  Under this model, individual dependence can be induced by using the same
  half-normal random variable for S and T.
  This is obtained by setting \code{indCorr = TRUE},
  but the amount of correlation is not dependent on a single parameter.
}
\value{
  A \link{data.frame} with columns
  \item{trialref}{the trial reference}
  \item{trt}{the treatment arm (-0.5 or 0.5)}
  \item{id}{the patient id}
  \item{timeT}{the value of the true endpoint T}
  \item{statusT}{the censoring/event (0/1) indicator of the true endpoint T}
  \item{timeS}{the value of the surrogate endpoint S}
  \item{statusS}{the censoring/event (0/1) indicator of the surrogate endpoint S}
}
\references{
Burzykowski T, Cortinas Abrahantes J (2005). Validation in the case of two failure-time endpoints.
In \emph{The Evaluation of Surrogate Endpoints} (pp. 163-194). Springer, New York.

Shi Q, Renfro LA, Bot BM, Burzykowski T, Buyse M, Sargent DJ.
Comparative assessment of trial-level surrogacy measures for candidate
time-to-event surrogate endpoints in clinical trials. 
\emph{Computational Statistics & Data Analysis} 2011; 55: 2748--2757.
}

\author{
  \packageAuthor{surrosurv}
}

\examples{
  set.seed(1)
  simData.re(N = 2, ni = 5)
  simData.cc(N = 2, ni = 5)
  simData.mx(N = 2, ni = 5)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
