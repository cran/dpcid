# Linear shrinkage estimates of covariance and inverse covariance matrix
# Arguments
# X: observed dataset from a specific condition
# scaling: logical value for scaling columns to have unit variance estimates
# Values
# shr_cov: Linear shrinkage estimate of the covariance matrix
# shr_inv: Linear shrinkage estimate of the inverse covariance matrix
lshr.cov <- function(X,scaling=FALSE)
{
  A = scale(X,center=T,scale=scaling)
  n = nrow(A)
  p = ncol(A)
  
  SA = var(A)*((n-1)/n)
  mk = sum(diag(SA))/p
  temp = 2*sum(SA[upper.tri(SA)]^2)
  dk2 = (sum((diag(SA) - mk)^2) + temp)/p
  
  temp2 = apply(A,1,function(x) sum(x^2))
  bbk2 = ( sum(temp2^2) - n*(sum(diag(SA)^2)+temp))/(n^2*p)
  bk2 = min(bbk2,dk2)
  ak2 = dk2 - bk2
  if(ak2==0)
     shr_cov = mk*diag(p)
  else
     shr_cov = (bk2/dk2)*mk*diag(p)+(ak2/dk2)*SA
     
  inv_shr = solve(shr_cov)

  return(list(shr_cov=shr_cov,shr_inv=inv_shr))
}


# DPCID for two precision matrices
# Arguments
# A: observed dataset from the first condition
# B: observed dataset from the second condition
# lambda1: a tuning parameter for the ridge penalty
# lambda2: a tuning parameter for the fusion penalty between two precision matrices
# wd1: the estimate of diagonal elements of the precision matrix of the first condition
# wd2: the estimate of diagonal elements of the precision matrix of the second condition
# rho1_init: an initial value for the partial correlation matrix of the first condition
# rho2_init: an initial value for the partial correlation matrix of the second condition
# niter: a total number of iterations in the block-wise coordinate descent
# tol: a tolerance for the convergence
# Values
# rho1: an estimated partial correlatioin matrix of the first condition
# rho2: an estimated partial correlatioin matrix of the first condition
dpcid_core <- function(A,B,lambda1,lambda2,wd1,wd2,rho1_init,rho2_init,niter=1000,tol=1e-6)
{
   n1 = nrow(A)
   n2 = nrow(B)
   p = ncol(A)

   out = .C('dpcid_c', as.integer(n1), as.integer(n2), as.integer(p),
   as.double(A), as.double(B), as.double(lambda1), as.double(lambda2),
   as.double(wd1),as.double(wd2),
   as.double(rho1_init), as.double(rho2_init), as.integer(niter), as.double(tol),
   rho1 = as.double(matrix(0,p,p)), rho2 = as.double(matrix(0,p,p)),
   resid1 = as.double(rep(0,n1*p)),resid2=as.double(rep(0,n2*p)))

   out1 = matrix(out$rho1,p,p)
   out2 = matrix(out$rho2,p,p)
   resid1 = matrix(out$resid1,n1,p)
   resid2 = matrix(out$resid2,n2,p)

   return(list(rho1 = out1,rho2=out2,resid1=resid1,resid2=resid2))
}

# DPCID for two precision matrices  (with l1 penalty instead of l2 norm)
# Arguments
# A: observed dataset from the first condition
# B: observed dataset from the second condition
# lambda1: a tuning parameter for the lasso penalty
# lambda2: a tuning parameter for the fusion penalty between two precision matrices
# wd1: the estimate of diagonal elements of the precision matrix of the first condition
# wd2: the estimate of diagonal elements of the precision matrix of the second condition
# rho1_init: an initial value for the partial correlation matrix of the first condition
# rho2_init: an initial value for the partial correlation matrix of the second condition
# niter: a total number of iterations in the block-wise coordinate descent
# tol: a tolerance for the convergence
# Values
# rho1: an estimated partial correlatioin matrix of the first condition
# rho2: an estimated partial correlatioin matrix of the first condition
dpcid_l1_core <- function(A,B,lambda1,lambda2,wd1,wd2,rho1_init,rho2_init,niter=1000,tol=1e-6)
{
   n1 = nrow(A)
   n2 = nrow(B)
   p = ncol(A)

   out = .C('dpcid_l1', as.integer(n1), as.integer(n2), as.integer(p),
   as.double(A), as.double(B), as.double(lambda1), as.double(lambda2),
   as.double(wd1),as.double(wd2),
   as.double(rho1_init), as.double(rho2_init), as.integer(niter), as.double(tol),
   rho1 = as.double(matrix(0,p,p)), rho2 = as.double(matrix(0,p,p)),
   resid1 = as.double(rep(0,n1*p)),resid2=as.double(rep(0,n2*p)))

   out1 = matrix(out$rho1,p,p)
   out2 = matrix(out$rho2,p,p)
   resid1 = matrix(out$resid1,n1,p)
   resid2 = matrix(out$resid2,n2,p)

   return(list(rho1 = out1,rho2=out2,resid1=resid1,resid2=resid2))
}


# Differential Partial Correlation Estimation for two conditions
# Arguments
# A: observed dataset from the first condition
# B: observed dataset from the second condition
# lambda1: a tuning parameter for the ridge penalty
# lambda2: a tuning parameter for the fusion penalty between two precision matrices
# niter: a total number of iterations in the block-wise coordinate descent
# tol: a tolerance for the convergence
# scaling: a logical flag for scaling variable to have unit variance.
#          Default is scaling=FALSE.
# Values
# rho1: an estimated partial correlatioin matrix of the first condition
# rho2: an estimated partial correlatioin matrix of the second condition
# wd1: A vector of estimated diagonal elements of the first precision matrices
# wd2: A vector of estimated diagonal elements of the second precision matrices
# diff_edge: an index matrix of different edges between two conditions
# n_diff: the number of different edges between two conditions
dpcid <- function(A,B,lambda1,lambda2,niter=1000,tol=1e-6,scaling=FALSE)
{
     n1 = nrow(A)
     n2 = nrow(B)
     p = ncol(A)

     if(scaling==T)
     {
        A = scale(A,center=T,scale=T)
        B = scale(B,center=T,scale=T)
     } else
     {
        A = scale(A,center=T,scale=F)
        B = scale(B,center=T,scale=F)
     }

     
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

     rho1 = rho1_init
     rho2 = rho1_init
    
     dpc = dpcid_core(A,B,lambda1,lambda2,wd1,wd2,rho1,rho2,niter,tol);

     rho1 = dpc$rho1;
     rho2 = dpc$rho2;
     est_mat = (abs(rho1-rho2)>tol)
     diag(est_mat) = 0

     diff_edge = which(est_mat==1,T)
     diff_edge = matrix(diff_edge[diff_edge[,1]<diff_edge[,2],],ncol=2)
     n_diff = nrow(diff_edge)
     wd = list()
     wd[[1]] = diag(PM1)
     wd[[2]] = diag(PM2)
     return(list(rho1=rho1,rho2=rho2,wd1=wd1,wd2=wd2,diff_edge=diff_edge,n_diff=n_diff))
}


# K-fold Crossvalidation for the lambda1
# Arguments
# A: observed dataset from the first condition
# B: observed dataset from the second condition
# nfold: the number of folds in the crossvalidation(i.e., K in K-fold cross validation)
# seq_lambda1: a sequence of tuning parameters for the ridge penalty
# niter: a total number of iterations in the block-wise coordinate descent
# tol: a tolerance for the convergence
# scaling: a logical flag for scaling variable to have unit variance.
#          Default is scaling=FALSE.
# Values
# cv_err: a vector of crossvalidated errors corresponding to a given sequence of tuning paramters
cv.lambda1 <- function(A,B, nfold,seq_lambda1,niter=1000,tol=1e-6,scaling=FALSE)
{
    n1 = nrow(A)
    n2 = nrow(B)
    p = ncol(A)


     if(scaling==T)
     {
        A = scale(A,center=T,scale=T)
        B = scale(B,center=T,scale=T)
     } else
     {
        A = scale(A,center=T,scale=F)
        B = scale(B,center=T,scale=F)
     }

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


     m = length(seq_lambda1)
     cv_err = rep(0,m)
     names(cv_err) = seq_lambda1

     len1 = floor(n1/nfold);
     indx1 = list()
     tag1 = sample(n1)
     
     for(i in 1:(nfold-1))
     {
        fr = len1*(i-1)+1
        ba = len1*i
        indx1[[i]] = tag1[fr:ba]
     }
     fr = len1*(nfold-1)
     indx1[[nfold]] = tag1[fr:n1]
     
     len2 = floor(n2/nfold);
     indx2 = list()
     tag2 = sample(n2)
     
     for(i in 1:(nfold-1))
     {
        fr = len2*(i-1)+1
        ba = len2*i
        indx2[[i]] = tag2[fr:ba]
     }
     fr = len2*(nfold-1)
     indx2[[nfold]] = tag2[fr:n1]
     
     
     for(i in 1:nfold)
     {
        tr_A = A[unlist(indx1[-i]),]
        test_A = A[indx1[[i]],]
        tr_B = B[unlist(indx2[-i]),]
        test_B = B[indx2[[i]],]
         
        rho1 = rho1_init
        rho2 = rho2_init
        for(j in 1:m)
        {
            lambda1 = seq_lambda1[j]   
            dpc= dpcid_core(tr_A,tr_B,lambda1,0,wd1,wd2,rho1,rho2,niter=1000,tol=1e-6)
            rho1 = dpc$rho1
            rho2 = dpc$rho2
            
            for(k in 1:p)
            {
                beta1 = rho1[-k,k]*sqrt(wd1[-k]/wd1[k])
                cv_err[j] = cv_err[j] + sum((test_A[,k]-test_A[,-k]%*%beta1)^2)
                beta2 = rho2[-k,k]*sqrt(wd2[-k]/wd2[k])
                cv_err[j] = cv_err[j] + sum((test_B[,k]-test_B[,-k]%*%beta2)^2)
            }
        }
     }
    return(list(cv=cv_err,pm1=PM1,pm2=PM2))
}




# AIC and BIC for the DPCID
# Arguments
# A: observed dataset from the first condition
# B: observed dataset from the second condition
# lambda1: the chosen tuning parameter lambda1 by cv.lambda1
# seq_lambda2: a sequence of tuning parameters for the fusion penalty
# wd1: the estimate of diagonal elements of the precision matrix of the first condition
# wd2: the estimate of diagonal elements of the precision matrix of the second condition
# rho1_init: an initial value for the partial correlation matrix of the first condition
# rho2_init: an initial value for the partial correlation matrix of the second condition
# niter: a total number of iterations in the block-wise coordinate descent
# tol: a tolerance for the convergence
# scaling: a logical flag for scaling variable to have unit variance.
#          Default is scaling=FALSE.
# Values
# aic: a vector of aic values corresponding to a given sequence of tuning paramters
# bic: a vector of bic values corresponding to a given sequence of tuning paramters
crit.dpcid <- function(A,B,l1,seq_l2,wd1,wd2,rho1_init,rho2_init,niter=1000,tol=1e-6,scaling=FALSE)
{
    n1 = nrow(A)
    n2 = nrow(B)
    p = ncol(A)

     if(scaling==T)
     {
        A = scale(A,center=T,scale=T)
        B = scale(B,center=T,scale=T)
     } else
     {
        A = scale(A,center=T,scale=F)
        B = scale(B,center=T,scale=F)
     }

     
     m = length(seq_l2)
     aic_value = rep(0,m)
     names(aic_value) = seq_l2

     bic_value = rep(0,m)
     names(bic_value) = seq_l2
     rho1 = rho1_init
     rho2 = rho2_init
     for(j in 1:m)
     {
         lambda2 = seq_l2[j]
         
         dpc = dpcid_core(A,B,l1,lambda2,wd1,wd2,rho1,rho2,niter,tol);

         rho1 = dpc$rho1;
         rho2 = dpc$rho2;
         est_mat = (abs(rho1-rho2)>tol)
         diag(est_mat) = 0

         diff_edge = which(est_mat==1,T)
         #diff_edge = matrix( diff_edge[ diff_edge[,1]< diff_edge[,2],],ncol=2)
         n_diff = nrow(diff_edge)

         rss1 = apply(dpc$resid1,2,function(x) return(sum(x^2)))
         rss2 = apply(dpc$resid2,2,function(x) return(sum(x^2)))
         
         temp = sum(rss1*wd1)+sum(rss2*wd2) 
         aic_value[j] = temp +2*n_diff
         bic_value[j] = temp +log(n1+n2)*n_diff
    }
   return(list(aic=aic_value,bic=bic_value))
}

