\name{crit.dpcid}
\alias{crit.dpcid}
\title{Akaike information criterion (AIC) and Bayesian information criterion (BIC) for the DPCID}
\usage{
crit.dpcid(A,B,l1,seq_l2,wd1,wd2,rho1_init,rho2_init,niter=1000,tol=1e-6,scaling=FALSE)
}
\arguments{
  \item{A}{An observed dataset from the first condition.}

  \item{B}{An observed dataset from the second condition.}

  \item{l1}{The selected lambda1 in cv.lambda1.}

  \item{seq_l2}{A sequence of tuning parameter lambda2 for the fusion penalty}
  \item{wd1}{The estimate of diagonal elements of the precision matrix of the first condition.}
  \item{wd2}{The estimate of diagonal elements of the precision matrix of the second condition.}
  \item{rho1_init}{An initial value for the partial correlation matrix of the first condition.}
  \item{rho2_init}{An initial value for the partial correlation matrix of the second condition.}

  \item{niter}{A total number of iterations in the block-wise coordinate descent.}

  \item{tol}{A tolerance for the convergence.}

  \item{scaling}{a logical flag for scaling variable to have unit variance. Default is FALSE.}
}

\value{
  \item{aic}{A vector of aic values corresponding to a given sequence of tuning paramters.}
  \item{bic}{A vector of bic values corresponding to a given sequence of tuning paramters.}
}
\description{
aic.dpcid returns the AIC values corresponding to the given lambda1 and lambda2 values
for the DPCID.
}

\details{
crit.dpcid needs the estimates of the diagonal elements of two precision matrices.
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

l1 = 0.3
seq_l2 = seq(0.1,1,by=0.2)

crit =crit.dpcid(A,B,l1,seq_l2,wd1,wd2,rho1_init,rho2_init)
crit$aic
crit$bic
}

