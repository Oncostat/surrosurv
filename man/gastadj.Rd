\name{gastadj}
\docType{data}
\alias{gastadj}
\title{Individual data from the adjuvant GASTRIC meta-analysis}
\description{
  The \code{gastadj} dataset contains individual data
    (overall and disease-free survival)
    of 3288 patients with resectable gastric cancer
    from 14 randomized trials of adjuvant chemotherapy.
}

\usage{data(gastadj)}
\format{A dataframe with variables:
  \describe{
    \item{timeT:}{Overall survival time (days).}
    \item{statusT:}{Overall survival indicator (0=censored, 1=death).}
    \item{timeS:}{Disease-free survival time (days).}
    \item{statusS:}{Disease-free survival indicator (0=censored, 1=progression on death).}
    \item{trialref:}{Trial indicator}
    \item{trt:}{Treatment arm (-0.5 = control, 0.5=chemotherapy).}
    \item{id:}{Patient identifier.}
    }
  }
  
\source{
The authors thank the GASTRIC
  (Global Advanced/Adjuvant Stomach Tumor Research International Collaboration)
  Group for permission to use their data.
The investigators who contributed to GASTRIC are listed in 
  Oba et al (2013) and GASTRIC (2010).
The GASTRIC Group data are available within the \code{surrosurv} package
  for research purposes, under the conditions that
  (1) the research be scientifically appropriate,
  (2) the confidentiality of individual patient data be protected,
  (3) the results of the analyses be shared with the GASTRIC Group prior to public communication,
  (4) the source of data be fully acknowledged as above, and
  (5) resulting data and results be further shared with the research community.
}
\references{
  Paoletti X, Oba K, Bang Y-J, et al.
    Disease-free survival as a surrogate for overall survival in adjuvant trials
    of gastric cancer: a meta-analysis.
    \emph{J Ntl Cancer Inst}, 105(21):1600-7, 2013.
    \doi{10.1093/jnci/djt270}.
    
  The GASTRIC group.
    Benefit of adjuvant chemotherapy for resectable gastric cancer:
    a meta-analysis. 
    \emph{JAMA}, 303(17):1729-37, 2010.
    \doi{10.1001/jama.2010.534}.
    
  Buyse M, Molenberghs G, Paoletti Xavier et al.
    Statistical evaluation of surrogate endpoints
    with examples from cancer clinical trials.
    \emph{Biom J}, 58(1):104-32, 2016.
    \doi{10.1002/bimj.201400049}
}
\examples{
  \dontrun{
    data('gastadj')
    allSurroRes <- surrosurv(gastadj, c('Clayton', 'PoissonTIa'), verbose = TRUE)
    convergence(allSurroRes)
    allSurroRes
    predict(allSurroRes)
    plot(allSurroRes)
  }
}

\keyword{datasets}
\keyword{gastadj}
\keyword{adjuvant}
\keyword{gastric}
\keyword{cancer}
\keyword{surrogate}

      
      