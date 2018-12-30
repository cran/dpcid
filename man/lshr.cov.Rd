\name{lshr.cov}
\alias{lshr.cov}
\title{Linear shrinkage estimates of covariance and inverse covariance matrix}
\usage{
lshr.cov(X,scaling=FALSE)
}
\arguments{
  \item{X}{An observed dataset from a specific condition.}
  \item{scaling}{a logical flag for scaling variable to have unit variance. Default is FALSE.}
}

\value{
  \item{shr_cov}{Linear shrinkage estimate of the covariance matrix.}
  \item{shr_inv}{Linear shrinkage estimate of the inverse covariance matrix.}
}
\description{
Linear shrinkage estimates of covariance and inverse covariance matrix.
}

\details{
shr_covp returns the optimal linear shrinkage parameter, the linear
shrinkage estimates of the  covariance and the precision matrix.
}

\references{
Ledoit, O. and M.~Wolf, M. (2004). A well-conditioned estimator for large-dimensional covariance matrices, \emph{Journal of Multivariate Analysis}, 88, 365--411.

Yu, D., Lee, S. H., Lim, J., Xiao, G., Craddock, R. C., and Biswal, B. B. (2018). Fused Lasso
Regression for Identifying Differential Correlations in Brain Connectome Graphs. \emph{Statistical Analysis and Data Mining}, 11, 203--226.
}

\examples{
library(MASS)

## True precision matrix
omega <- matrix(0,5,5)
omega[1,2] <- omega[1,3] <- omega[1,4] <- 1
omega[2,3] <- omega[3,4] <- 1.5
omega <- t(omega) + omega
diag(omega) <- 3

Sig = solve(omega)
X = mvrnorm(50,rep(0,5),Sig)
lshr.cov(X)
}

