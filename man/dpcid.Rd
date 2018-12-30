\name{dpcid}
\alias{dpcid}
\title{
Differential partial correlation identification with the ridge and the fusion penalties}
\usage{
dpcid(A,B,lambda1,lambda2,niter=1000,tol=1e-6,scaling=FALSE)
}
\arguments{
  \item{A}{An observed dataset from the first condition.}

  \item{B}{An observed dataset from the second condition.}

  \item{lambda1}{A tuning parameter for the ridge penalty.}

  \item{lambda2}{A a tuning parameter for the fusion penalty between two precision matrices.}

  \item{niter}{A total number of iterations in the block-wise coordinate descent.}

  \item{tol}{A tolerance for the convergence.}

  \item{scaling}{a logical flag for scaling variable to have unit variance. Default is FALSE.}
}

\value{
  \item{rho1}{An estimated partial correlatioin matrix of the first condition.}
  \item{rho2}{An estimated partial correlatioin matrix of the second condition.}

  \item{wd1}{A vector of estimated diagonal elements of the first precision matrices.}
  \item{wd2}{A vector of estimated diagonal elements of the second precision matrices.}
  \item{diff_edge}{An index matrix of different edges between two conditions.}
  \item{n_diff}{The number of different edges between two conditions.}
}
\description{
DPCID is a procedure for the differential partial correlation identification
with the ridge and the fusion penalties.
This function conducts the two stage procedure (diagonal and partial correlation steps).
}

\details{
In the first step (lshr.cov), each precision matrix is estimated
from the optimal linear shrinkage covariance matrix.
In the second step (dpcid_core), two partial correlation matrices
are jointly estimated with a given tuning parameters lambda1 and lambda2
and fixed diagonal elements of two precision matrices.
}

\references{
Yu, D., Lee, S. H., Lim, J., Xiao, G., Craddock, R. C., and Biswal, B. B. (2018). Fused Lasso
Regression for Identifying Differential Correlations in Brain Connectome Graphs. \emph{Statistical Analysis and Data Mining}, 11, 203--226.
}

\examples{
library(MASS)

## True precision matrix
omega1 <- matrix(0,5,5)
omega1[1,2] <- omega1[1,3] <- omega1[1,4] <- 1
omega1[2,3] <- omega1[3,4] <- 1.5
omega1 <- t(omega1) + omega1
diag(omega1) <- 3

omega2 <- matrix(0,5,5)
omega2[1,3] <- omega2[1,5] <- 1.5
omega2[2,3] <- omega2[2,4] <- 1.5
omega2 <- t(omega2) + omega2
diag(omega2) <- 3


Sig1 = solve(omega1)
Sig2 = solve(omega2)
X1 = mvrnorm(50,rep(0,5),Sig1)
X2 = mvrnorm(50,rep(0,5),Sig2)

lambda1 = 0.2
lambda2 = 0.2

res = dpcid(X1,X2,lambda1,lambda2,niter=1000,tol=1e-6)

}

