\name{PACVB-package}
\alias{PACVB-package}
\alias{PACVB}
\docType{package}
\title{
Variational Bayes (VB) Approximation of Gibbs Posteriors with Hinge Losses
}
\description{
The package computes a Gaussian approximation of a Gibbs posterior for convexified risks.
We use a hinge loss, and a hinge variation of the AUC loss. The implementation follows the lines of Alquier et al. [2015]. The authors obtain PAC-Bayesian bounds for the variational approximation of an exponential weighted average estimator. The optimal bound is obtain through a gradient descent on a convex objective.  
}
\details{
\tabular{ll}{
Package: \tab PACVB\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-07-29\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Author: James Ridgway\cr
Maintainer: James Ridgway <james.ridgway@bristol.ac.uk>
}
\references{
Alquier, P., Ridgway, J., and  Chopin, N. On the properties of variational approximations of Gibbs posteriors. arXiv preprint, 2015
}
\keyword{PAC-Bayesian, classification, ranking}

