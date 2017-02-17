\name{gastadv}
\docType{data}
\alias{gastadv}
\title{Individual data from the advanced GASTRIC meta-analysis}
\description{
  The \code{gastadv} dataset contains individual data
    (overall and progression-free survival)
    of 4069 patients with advanced/recurrent gastric cancer
    from 20 randomized trials of chemotherapy.
}

\usage{data(gastadv)}
\format{A dataframe with variables:
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
  
\source{
The authors thank the GASTRIC
  (Global Advanced/Adjuvant Stomach Tumor Research International Collaboration)
  Group for permission to use their data.
The investigators who contributed to GASTRIC are listed in 
  Paoletti et al (2013) and GASTRIC (2013).
The GASTRIC Group data are available within the \inline{surrosurv} package
  for research purposes, under the conditions that
  (1) the research be scientifically appropriate,
  (2) the confidentiality of individual patient data be protected,
  (3) the results of the analyses be shared with the GASTRIC Group prior to public communication,
  (4) the source of data be fully acknowledged as above, and
  (5) resulting data and results be further shared with the research community.
}
\references{
  Paoletti X, Oba K, Bang Y-J, et al.
    Progression-free survival as a surrogate for overall survival
    in advanced/recurrent gastric cancer trials: a meta-analysis.
    \emph{J Ntl Cancer Inst}, 105(21):1667-70, 2013.
    doi:\href{http://dx.doi.org/10.1093/jnci/djt269}{10.1093/jnci/djt269}.
    
  The GASTRIC group.
    Role of chemotherapy for advanced/recurrent gastric cancer:
    An individual-patient-data meta-analysis. 
    \emph{Eur J Cancer}, 49(7):1565-77, 2013.
    doi:\href{http://dx.doi.org/10.1016/j.ejca.2012.12.016}{10.1016/j.ejca.2012.12.016}.
    
  Buyse M, Molenberghs G, Paoletti Xavier et al.
    Statistical evaluation of surrogate endpoints
    with examples from cancer clinical trials.
    \emph{Biom J}, 58(1):104-32, 2016.
    doi:\href{http://dx.doi.org/10.1002/bimj.201400049}{10.1002/bimj.201400049}
}
\examples{
  \dontrun{
    data('gastadv')
    allSurroRes <- surrosurv(gastadv, c('Clayton', 'PoissonTIa'), verbose = TRUE)
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

      
      