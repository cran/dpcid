\name{dpcid_l1_core}
\alias{dpcid_l1_core}
\title{
Identification of two partial correlation matrices having sparse differences}
\usage{
dpcid_l1_core(A,B,lambda1,lambda2,wd1,wd2,rho1_init,rho2_init,niter=1000,tol=1e-6)

}
\arguments{
  \item{A}{An observed dataset from the first condition.}

  \item{B}{An observed dataset from the second condition.}

  \item{lambda1}{A tuning parameter for the lasso penalty.}

  \item{lambda2}{A a tuning parameter for the fusion penalty between two precision matrices.}

  \item{wd1}{The estimate of diagonal elements of the precision matrix of the first condition.}
  \item{wd2}{The estimate of diagonal elements of the precision matrix of the second condition.}
  \item{rho1_init}{An initial value for the partial correlation matrix of the first condition.}
  \item{rho2_init}{An initial value for the partial correlation matrix of the second condition.}

  \item{niter}{A total number of iterations in the block-wise coordinate descent.}

  \item{tol}{A tolerance for the convergence.}
}

\value{
  \item{rho1}{An estimated partial correlatioin matrix of the first condition.}
  \item{rho2}{An estimated partial correlatioin matrix of the second condition.}
  \item{resid1}{Residuals of the first condtion.}
  \item{resid2}{Residuals of the second condtion.}
}
\description{
dpcid_core estimates two partial correlation matrices by
applying the regression approach with the lasso penalty and the fusion penalty.}

\details{
Dpcid_l1_core is the partial correlation step of the differential partial correlation identification
method with the lasso penalty instead of the ridge penalty by Yu et al.(2018).
The dpcid_l1_core estimates two partial correlation matrices
with the estimated diagonal elements of two precision matrices from
the optimal linear shrinkage estimates.
The estimated precision matrices by the optimal linear shrinkage estimates
can simply be used as the initial values of PM1 and PM2 in the dpcid_core.
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

A = scale(X1,center=TRUE,scale=FALSE)
B = scale(X2,center=TRUE,scale=FALSE)

lambda1 = 1
lambda2 = 1

shr_res = lshr.cov(A)
PM1 = shr_res$shr_inv

shr_res = lshr.cov(B)
PM2 = shr_res$shr_inv

wd1 = diag(PM1)
wd2 = diag(PM2)

rho1_init = -(1/sqrt(wd1))*PM1
rho1_init = t( 1/sqrt(wd1)*t(rho1_init))
diag(rho1_init) = 1

rho2_init = -(1/sqrt(wd2))*PM2
rho2_init = t( 1/sqrt(wd2)*t(rho2_init))
diag(rho2_init) = 1

res = dpcid_l1_core(A, B, lambda1, lambda2, wd1,wd2, rho1_init, rho2_init)
    
}

