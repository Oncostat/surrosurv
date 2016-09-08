<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

<head>
  <title>surrosurv: Evaluation of Failure Time Surrogate Endpoints in Individual Patient Data
    Meta-Analyses in R</title>

  <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
  <meta name="author" content="F. Rotolo" />
  <meta name="description" content="Evaluation of Failure Time Surrogate Endpoints in Individual Patient Data
    Meta-Analyses in R" />
  <meta name="keywords" content="surrogate, endpoint, evaluation, surrogacy, copula, buyse, burzykowski, poisson, poissonize, survival, R, package" /> 

  <link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  <link rel="stylesheet" type="text/css" href="https://r-forge.r-project.org/themes/css/gforge.css" />
  <link rel="stylesheet" type="text/css" href="https://r-forge.r-project.org/themes/rforge/css/theme.css" />
</head>
<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- project title  -->

<h2>Evaluation of Failure Time Surrogate Endpoints in Individual Patient Data
    Meta-Analyses in R</h2>
<em>F.Rotolo</em>
<p><a href="http://cran.r-project.org/web/packages/surrosurv/index.html">Information about the current plublic release on CRAN can be found here.</a>
<br><a href="http://r-forge.r-project.org/projects/surrosurv/">Information about the project development can be found here.</a>


...


<!--References--><hr>
<h3>References</h3>
<p> 
  Burzykowski T, Molenberghs G, Buyse M et al.
  Validation of surrogate end points in multiple randomized clinical trials
  with failure time end points.
  <em>Journal of the Royal Statistical Society C</em> 2001;
  <b>50</b>:405--422.
  <a href="http://dx.doi.org/10.1111/1467-9876.00244">10.1111/1467-9876.00244</a>
</p>
  
<p>
  Clayton DG.
  A model for association in bivariate life tables
  and its application in epidemiological studies of familial tendency 
  in chronic disease incidence.
  <em>Biometrika}</em>1978; 
  <b>65</b>:141--151.
  <a href="http://dx.doi.org/10.1093/biomet/65.1.141">10.1093/biomet/65.1.141</a>
</p>
  
<p>
  Crowther MJ, Riley RD, Staessen JA, Wang J, Gueyffier F, Lambert PC.
  Individual patient data meta-analysis of survival data
  using Poisson regression models.
  <em>BMC Medical Research Methodology</em> 2012;
  <b>12</b>:34.
  <a href="http://dx.doi.org/10.1080/01621459.1965.10480807">10.1080/01621459.1965.10480807</a>
  \doi{10.1186/1471-2288-12-34}.
</p>
  
<p>
  Gasparrini A, Armstrong B, Kenward MG.
  Multivariate meta-analysis for non-linear and other multi-parameter associations.
  <em>Statistics in Medicine</em> 2012; 
  <b>31</b>:3821--39.
  <a href="http://dx.doi.org/10.1080/01621459.1965.10480807">10.1080/01621459.1965.10480807</a>
  \doi{10.1002/sim.5471}
</p>
  
<p>
  Hougaard P.
  A class of multivariate failure time distributions.
  <em>Biometrika</em> 1986; 
  <b>73</b>:671--678.
  <a href="http://dx.doi.org/10.1080/01621459.1965.10480807">10.1080/01621459.1965.10480807</a>
  \doi{10.1093/biomet/73.3.671}
</p>
  
<p>
  Plackett RL.
  A class of bivariate distributions. 
  <em>Journal of the America Statistical Association</em> 1965;
  <b>60</b>:516--522.
  <a href="http://dx.doi.org/10.1080/01621459.1965.10480807">10.1080/01621459.1965.10480807</a>
</p>
  
<p>
  Rotolo F, Paoletti X, Burzykowski T, Buyse M, Michiels S.
  Evaluation of failure time surrogate endpoints in individual patient data
  meta-analyses of randomized clinical trials. A Poisson approach.
  <em>Working paper</em>.
</p>
  
<p>
  van Houwelingen HC, Arends LR, Stijnen T.
  Advanced methods in meta-analysis: multivariate approach and meta-regression.
  <em>Statistics in Medicine</em> 2002;
  <b>21</b>:589--624.
  <a href="http://dx.doi.org/10.1002/sim.1040">10.1002/sim.1040</a>
</p>
  
</body>
</html>
