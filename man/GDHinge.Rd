\name{GDHinge}
\alias{GDHinge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the VB (variational Bayes) approximation of a Gibbs measure with convexified classification loss
using a gradient descent.
}
\description{
The function computes a VB approximation of a Gibbs measure with a convexified 01-loss using a gradient descent.
The user must specify the design matrix, the response vector and the inverse temperature parameter. The number of iteration of the convex solver is fixed a priori (the default value is 10), an informed choice can be made using theorem 6.3 in 
Alquier et al. [2015].  
%The function computes the EP approximation of a logit regression with Gaussian prior. The user must specify the design matrix, the response vector
%The function computes the EP approximation of a logit regression with Gaussian prior. The user must specify the design matrix, the response vector
%and the prior variance. For more information on the default prior variance see Chopin and Ridgway [2015].
}
\usage{
GDHinge(X,Y,lambda,theta=0,K=100,v=10,ls=FALSE, 
                            B=0,family="F1",eps=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Design matrix. The matrix should include a constant column if a bias is to be considered. In addition
the Gradient descent has been calibrated considering a centered and scale design matrix.
}
\item{Y}{
Response vector. The vector should take values in \{-1,1\}. 
}
\item{lambda}{
Inverse temperature of the Gibbs posterior (See Alquier et al. (2015) for guide lines)
}
\item{theta}{
Initial value of the gradient descent. In the case of the "F1" family the last value of the vector is the initial log-variance. The vector should be of size p+1, 
where p is the number of columns of the design matrix.
If no initial value are chosen the algorithm is initialized to a Gaussian random vector. 
}

\item{K}{
Number of iteration of the gradient descent. Default value K=10. An informed choice can be made using theorem 6.3 in 
Alquier et al. [2015].
}
  \item{v}{
Prior variance. The prior is taken to be spherical Gaussian, the variance must therefore be specified in the form of a scalar.
For default choices see Alquier et al. [2015].  
The default is arbitrarily set to 10.
}
\item{ls}{
Logical value. Indicates if a linesearch should be used to find an optimal step length. Default value is FALSE. The 
option is not available for stochastic gradient descent.
}
\item{B}{
 Batch sizes when considering a stochastic gradient descent. B=0 corresponds to standard gradient descent.
}
\item{family}{
    Approximate family to consider when implementing VB. Possible values are: "F0" variance is fixed to 1/(sample size) times identity; "F1" spherical Gaussian. (see Alquier et al. [2015] for details)
}
\item{eps}{
    Probability of the empirical bound of the theoretical risk to be considered. Default is 0.05
}
}
\details{
The implementation is based on theorem 6.3 of Alquier et al. [2015] using convex solver presented in Nesterov [2004] (section 3.2.3).
The calibration depends on an upper bound on the l2 distance between the solution and the initial value. We use an arbitrary value of sqrt (p+1) 
We also give the possibility to use a linesearch algorithm satisfying the Wolfe conditions.
}
\value{
\item{m}{ Mean of the Gaussian approximation}
\item{s}{ Variance of the Gaussian approximation}
\item{bound}{Empirical bound on the aggregated risk. A negative value indicates that temperaure was taken outside of the admissible interval. The bound assumes that each element of the design is bounded by 1. c_x=1 in Alquier et al. [2015]}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\section{Warning}{
The columns of the design matrix should be centered and scale fo proper behaviour of the algorithm.
}
\references{
Alquier, P., Ridgway, J., and  Chopin, N. On the properties of variational approximations of Gibbs posteriors. arXiv preprint, 2015.

Nesterov, Y. Introductory lectures on convex optimization, volume 87. Springer Science and Buisness Media, 2004. 
}
\author{
James Ridgway
}
% \section{Warning }{
%The current implementation does not include damping or the possibility to use fractional EP (hopefully it will in a future version).
%This might results in poor performance for large datasets.  
%}


\examples{
data(Pima.tr)
Y<-2*as.matrix(as.numeric(Pima.tr[,8]))-3% formating the response
X<-data.matrix(Pima.tr[,1:7])
m<-apply(X,2,mean)
v<-apply(X,2,sd)
X<-t(apply(t(apply(X,1,"-",m)),1,"/",v))%centered normalized
X<-cbind(1,X)% adding a constant to the model
l<-45
Sol<-GDHinge(X,Y,l)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ VB }
\keyword{ classification }
\keyword{ PAC-Bayesian Bounds }% __ONLY ONE__ keyword per line
