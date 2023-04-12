[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/surrosurv)](https://cran.r-project.org/package=surrosurv)

# surrosurv
Evaluation of Failure Time Surrogate Endpoints in Individual Patient Data Meta-Analyses

## Description
Provides functions for the evaluation of
    surrogate endpoints when both the surrogate and the true endpoint are failure
    time variables. The approaches implemented are: (1) the two-step approach
    (Burzykowski et al, 2001) [<DOI:10.1111/1467-9876.00244>](http://dx.doi.org/10.1111/1467-9876.00244)
    with a copula model (Clayton, Plackett, Hougaard) at
    the first step and either a linear regression of log-hazard ratios at the second
    step (either adjusted or not for measurement error); (2) mixed proportional
    hazard models estimated via mixed Poisson GLM.