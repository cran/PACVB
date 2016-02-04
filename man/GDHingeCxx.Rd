\name{GDHingeCxx}
\alias{GDHingeCxx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
C++ internal function to compute the VB approximation  with convexified classification loss (use GDHinge instead).
%%  ~~function to do ... ~~
}
\description{
The function is the C++ internal function used by GDHinge}
\usage{
GDHingeCxx(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{
See arguments of GDHinge
%%     ~~Describe \code{X} here~~
}
}
\value{
\item{m}{ Mean of the Gaussian approximation}
\item{V}{ Variance matrix of the Gaussian approximation}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Alquier, P., Ridgway, J., and  Chopin, N. On the properties of variational approximations of Gibbs posteriors. arXiv preprint, 2015
}
\author{
James Ridgway
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ VB }
\keyword{ Classification }
\keyword{ PAC-Bayesian Bounds }% __ONLY ONE__ keyword per line
