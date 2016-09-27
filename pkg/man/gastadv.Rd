\name{gastadv}
\docType{data}
\alias{gastadv}
\title{Individual data from the advanced GASTRIC meta-analysis}
\description{
  
}
\usage{data(gastadv)}
\format{A dataframe containing data of 4069 patients from 20 randomized
  trials of chemotherapy.
  \describe{
    \item{timeT:}{Overall survival time (days).}
    \item{statusT:}{Overall survival indicator (0=censored, 1=death).}
    \item{timeS:}{Progression-free survival time (days).}
    \item{statusS:}{Progression-free survival indicator (0=censored, 1=prgoression on death).}
    \item{trialref:}{Trial indicator}
    \item{trt:}{Treatment arm (-0.5 = control, 0.5=chemomtherapy).}
    \item{id:}{Patient identifier.}
    }
  }
\references{
The GASTRIC group.
  Role of chemotherapy for advanced/recurrent gastric cancer:
  An individual-patient-data meta-analysis. 
  \emph{Eur J Cancer}, 49(7):1565-77, 2013.
  doi:\href{http://dx.doi.org/10.1016/j.ejca.2012.12.016}{10.1016/j.ejca.2012.12.016}.
  
Paoletti X, Oba K, Bang Y-J, et al.
  Progression-free survival as a surrogate for overall survival
  in advanced/recurrent gastric cancer trials: a meta-analysis.
  \emph{J Ntl Cancer Inst}, 105(21):1667-70, 2013.
  doi:\href{http://dx.doi.org/10.1093/jnci/djt269}{10.1093/jnci/djt269}.
}
\examples{
  \dontrun{
    data('gastadv')
    allSurroRes <- surrosurv(gastadv, c('Clayton', 'Poisson'), verbose = TRUE)
    convergence(allSurroRes)
    allSurroRes
    predict(allSurroRes)
    plot(allSurroRes)
  }
}

\keyword{datasets}
\keyword{gastadv}
\keyword{advanced}
\keyword{gastric}
\keyword{cancer}
\keyword{surrogate}

      
      