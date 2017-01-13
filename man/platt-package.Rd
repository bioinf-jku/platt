\name{platt-package}
\alias{platt-package}
\alias{platt}
\docType{package}
\title{
This package performs Platt-scaling.
}
\description{
This package implements a fast and robust function to estimate the
parameters of a sigmoid function that can map prediction values from 
arbitrary machine learning methods to probabilistic predictions.
}
\details{
\tabular{ll}{
Package: \tab platt\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2015-01-07\cr
License: \tab LGPL (>= 2.0)\cr
}
The package contains a testing data set MMP. The core function ist
"plattScaling" that fits the sigmoid function. Additionally the 
functions "predictProb" is provided to map new values to the sigmoid and
the function "ensemble" combines multiple probabilstic predicitons into 
one predcition.
}
\author{
Guenter Klambauer and Andreas Mayr.

Maintainer:  <klambauer@bioinf.jku.at>
}
\keyword{Platt scaling}
