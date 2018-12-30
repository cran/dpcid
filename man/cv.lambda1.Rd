\name{cv.lambda1}
\alias{cv.lambda1}
\title{K-fold Crossvalidation for the ridge step}
\usage{
cv.lambda1(A,B, nfold,seq_lambda1,niter=1000,tol=1e-6,scaling=FALSE)
}
\arguments{
  \item{A}{An observed dataset from the first condition.}

  \item{B}{An observed dataset from the second condition.}

  \item{nfold}{the number of folds in the crossvalidation(i.e., K in K-fold cross validation)}

  \item{seq_lambda1}{A sequence of tuning parameters for the ridge penalty}

  \item{niter}{A total number of iterations in the block-wise coordinate descent.}

  \item{tol}{A tolerance for the convergence.}

  \item{scaling}{a logical flag for scaling variable to have unit variance. Default is FALSE.}
}

\value{
  \item{cv}{A vector of crossvalidated errors corresponding to a given sequence of tuning paramters.}
  \item{pm1}{A matrix of the inverse of the linear shrinkage covariance estimates for the first condition.}
  \item{pm2}{A matrix of the inverse of the linear shrinkage covariance estimates for the second condition.}
}
\description{
K-fold crossvalidation for the choice of the lambda1.
}

\details{
cv.lambda1 returns a vector of the K-fold crossvalidated errors and matrices of the initial estimates of the precision matrices.
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

nfold = 5
seq_lam1 = seq(0.5,3,by=0.5)
cv_vec = cv.lambda1(X1,X2,nfold,seq_lam1,niter=1000,tol=1e-6)$cv

}

