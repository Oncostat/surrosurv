# surrosurv - total

    Code
      allSurroRes
    Output
                     kTau R2  
      Clayton unadj  0.61 0.45
      Clayton adj    0.61 0.41
      Plackett unadj 0.62 0.45
      Plackett adj   0.62 0.4 
      Hougaard unadj 0.32 0.45
      Hougaard adj   0.32 0.38
      PoissonT       -.-- 1   
      PoissonI       0.51 -.--
      PoissonTI      0.51 0.63
      PoissonTIa     0.51 0.83
    Code
      convergence(allSurroRes)
    Output
                     maxSgrad minHev minREev
      Clayton unadj     FALSE  FALSE     ---
      Clayton adj       FALSE  FALSE    TRUE
      Plackett unadj    FALSE  FALSE     ---
      Plackett adj      FALSE  FALSE    TRUE
      Hougaard unadj    FALSE   TRUE     ---
      Hougaard adj      FALSE   TRUE    TRUE
      PoissonT           TRUE   TRUE   FALSE
      PoissonI           TRUE   TRUE     ---
      PoissonTI          TRUE   TRUE    TRUE
      PoissonTIa         TRUE  FALSE    TRUE
    Code
      convals(allSurroRes)
    Output
                         maxSgrad        minHev      minREev
      Clayton unadj  1.500047e+00 -6.109032e+00          ---
      Clayton adj    1.500047e+00 -6.109032e+00 1.010648e-02
      Plackett unadj 3.247321e+02 -5.225147e+00          ---
      Plackett adj   3.247321e+02 -5.225147e+00 8.886121e-03
      Hougaard unadj 1.412709e+01  7.740882e-01          ---
      Hougaard adj   1.412709e+01  7.740882e-01 8.015472e-03
      PoissonT       2.519186e-05  1.299472e+02 6.688100e-12
      PoissonI       2.949779e-05  6.798776e+01          ---
      PoissonTI      3.214986e-06  6.700705e+01 2.041672e-02
      PoissonTIa     1.995411e-04 -1.456454e+07 1.024372e-01
    Code
      predict(allSurroRes, models = "PoissonTI")
    Output
      Treatment effect prediction for surrosurv object
      
         Poisson TI 
                                  1     2     3     4     5     6        
          Treatment effects on S: -0.52 -0.42 -0.38 -0.08 -0.51 -0.38 ...
          Treatment effects on T: -0.26 -0.08 -0.27  0.41 -0.41 -0.15 ...
    Code
      ste(allSurroRes)
    Output
                      beta   HR
      Clayton.unadj  -0.61 0.54
      Clayton.adj    -0.44 0.65
      Plackett.unadj -0.61 0.54
      Plackett.adj   -4.36 0.01
      Hougaard.unadj -0.61 0.54
      Hougaard.adj   -1.30 0.27
      PoissonT       -0.17 0.84
      PoissonTI      -0.65 0.52
      PoissonTIa     -1.16 0.31

# loocv - total

    Code
      loocvRes
    Output
      
         Clayton copula (Unadjusted) 
                   1     2     3     4     5     6    
      obsAlpha -0.65 -0.52 -0.12 -0.28 -0.24 -0.45 ...
      obsBeta  -0.31 -0.21 -0.09 -0.02 -0.22 -0.34 ...
      predict  -0.40 -0.31 -0.07 -0.17 -0.14 -0.27 ...
      lwr      -0.76 -0.65 -0.42 -0.51 -0.48 -0.62 ...
      upr      -0.05  0.02  0.28  0.17  0.21  0.09 ...
      kTau      0.60  0.60  0.61  0.60  0.60  0.60 ...
      R2        0.49  0.49  0.45  0.46  0.46  0.44 ...
      
         Clayton copula (Adjusted) 
                   1     2     3     4     5      6    
      obsAlpha -0.65 -0.52 -0.12 -0.28 -0.24 -0.453 ...
      obsBeta  -0.31 -0.21 -0.09 -0.02 -0.22 -0.342 ...
      predict  -0.39 -0.31 -0.09 -0.18 -0.14 -0.261 ...
      lwr      -0.69 -0.57 -0.35 -0.41 -0.39 -0.517 ...
      upr      -0.09 -0.04  0.17  0.06  0.10 -0.004 ...
      kTau      0.60  0.60  0.61  0.60  0.60  0.604 ...
      R2        0.46  0.45  0.42  0.46  0.43  0.411 ...
      
         Poisson TI 
                   1     2     3     4     5     6    
      obsAlpha -0.65 -0.52 -0.12 -0.28 -0.24 -0.45 ...
      obsBeta  -0.31 -0.21 -0.09 -0.02 -0.22 -0.34 ...
      predict  -0.69 -0.44  0.09 -0.20 -0.02 -0.38 ...
      lwr      -1.22 -1.00 -0.74 -0.83 -0.60 -0.91 ...
      upr      -0.15  0.11  0.91  0.43  0.56  0.14 ...
      kTau      0.51  0.52  0.51  0.52  0.51  0.52 ...
      R2        0.70  0.60  0.48  0.39  0.62  0.74 ...

