#
--
A reference page for the **LatticeKrig** package DOI 
--

 DOI version of source code Version 8.4 :
[LatticeKrig_5.5.tar.gz](www.image.ucar.edu/LatticeKrig/LatticeKrig_5.5.tar.gz) 

 **doi:XXXXX** 
 
 For the most recent version of **LatticeKrig** please use 
  the [R studio CRAN mirror](http://cran.rstudio.com) to download and install this package.
 
 [LatticeKrig home page](http://www.image.ucar.edu/~nychka/LatticeKrig)
 
 
###Description:
This is the DESCRIPTION file distributed with the package.

```
Package: LatticeKrig
Version: 5.5
Date: 2016-04-20
Title: Multiresolution Kriging Based on Markov Random Fields
Author: Douglas Nychka [aut, cre], Dorit Hammerling [aut], Stephan Sain [aut], Nathan Lenssen [aut]
Authors@R: c( 
            person("Douglas", "Nychka", role = c("aut", "cre"),
               email = "nychka@ucar.edu"),
            person("Dorit", "Hammerling", role = c("aut"),
               email = "hammerling@samsi.info"),
            person("Stephan", "Sain", role = "aut",
               email = "ssain@ucar.edu"),
            person("Nathan", "Lenssen", role = "aut",
               email = "lenssen@ucar.edu")) 
Maintainer: Douglas Nychka <nychka@ucar.edu>
```
```Description: ```
Functions for the interpolation of large spatial
  datasets. This package follows a "fixed rank Kriging" approach using
  a large number of basis functions and provides spatial estimates
  that are comparable to standard families of covariance functions.
  Using a large number of basis functions allows for estimates that
  can come close to interpolating the observations (a spatial model
  with a small nugget variance.)  Moreover, the covariance model for this method
  can approximate the Matern covariance family but also allows for a
  multi-resolution model and supports efficient computation of the
  profile likelihood for estimating covariance parameters. This is
  accomplished through compactly supported basis functions and a
  Markov random field model for the basis coefficients. These features
  lead to sparse matrices for the computations. An extension of this 
  version over previous ones ( < 5.4 ) is the support for different 
  geometries besides a rectangular domain. 
  One benefit of the LatticeKrig model/approach 
  is the facility to do unconditional and conditional
  simulation of the field for large numbers of arbitrary points. There
  is also the flexibility for estimating non-stationary covariances. Included are
  generic methods for prediction, standard errors for prediction,
  plotting of the estimated surface and conditional and unconditional
  simulation.

```

License: GPL (>= 2)
URL: http://www.r-project.org
Depends: R (>= 3.0.1), methods, spam, fields (>= 6.9.1)
Packaged: 2016-04-20; nychka
NeedsCompilation: yes
```

*Dates:*	2016 [Copyrighted]

*Rights:* 	This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version. 

--

<div> 
 <p>&copy; 2016 UCAR | 
<a href="http://www.ucar.edu/legal/privacy_policy.shtml">Privacy Policy</a> | 
<a href="http://www.ucar.edu/legal/terms_of_use.shtml">Terms of Use</a> | 
<a href="http://www.nsf.gov">Sponsored by NSF</a> | 
<a href="http://www.ucar.edu">Managed by UCAR</a> |  
<a href="mailto:kconrad@ucar.edu">Webmaster/Feedback</a> | 
<em>Postal Address: P.O. Box 3000, Boulder, CO 80307-3000  &bull;  Shipping Address: 1850 Table Mesa Drive, Boulder, CO 80305 &bull; <a href="http://www.ncar.ucar.edu/organization/about/">Contact</a></em>
</p>
</div>

<div>
<p>The National Center for Atmospheric Research is sponsored by the National Science Foundation.
Any opinions, findings and conclusions or recommendations expressed in this publication are those 
of the author(s) and do not necessarily reflect the views of the National Science Foundation.
</p>
</div>
