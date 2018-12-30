
#include<R.h>
#include<Rmath.h>
#include<R_ext/BLAS.h>
#include<memory.h>

#define max(A,B) (((A)>(B))?(A):(B))

// sign_ftn
double sign_ftn(double x)
{
	double res = ((double)(x>0) - (double)(x<0));
    return res;
}

// Soft-thresholding 
double soft(double x, double lam1)
{
  double s = (sign_ftn(x)*max(fabs(x)-lam1,0.0));
  return s;
}

// Arguments : X1, X2, lambda1, lambda2, CM1, CM2, niter, tol
void dpcid_l1(int *N1, int *N2, int *Dim,
double *X1, double *X2, double *Lam1, double *Lam2, 
double *wd1, double *wd2, double *rho1_init, double *rho2_init, int *Niter, double *Tol,
double *rho1, double *rho2, double *resid1, double *resid2)
{
	double *rho1_old, *rho2_old;
	double lambda1=*Lam1, lambda2=*Lam2,tol=*Tol;
	double x1e_pij,x2e_pij,*x1_pii,*x2_pii, x1e_pji, x2e_pji;
	//double *resid1, *resid2;
	double *temp_coef, temp=0,temp11,temp12,temp21,temp22;
	double XX1,XX2,Xe1,Xe2;
	double max_diff = 0;
	int n1=*N1,n2=*N2,p=*Dim;
	
	int niter=*Niter;
	int i,j,k;
	
	char NoTrans = 'N';
    double One = 1.0, MinusOne = -1.0;
    int IntOne = 1;

	rho1_old = (double *) malloc(p*p*sizeof(double));
	rho2_old = (double *) malloc(p*p*sizeof(double));
	//resid1 = (double *) malloc(n1*p*sizeof(double));
	//resid2 = (double *) malloc(n2*p*sizeof(double));

	// Initializatioin of rho
	memcpy(rho1_old,rho1_init,p*p*sizeof(double));
	memcpy(rho2_old,rho2_init, p*p*sizeof(double));
	memcpy(rho1, rho1_init, p*p*sizeof(double));
	memcpy(rho2, rho2_init, p*p*sizeof(double));

	memset(resid1,0,n1*p*sizeof(double));
	memset(resid2,0,n2*p*sizeof(double));

	temp_coef = (double *) malloc(p*sizeof(double));

	memset(temp_coef,0,p*sizeof(double));

	x1_pii = (double*)malloc(p*sizeof(double));
	x2_pii = (double*)malloc(p*sizeof(double));
	
	
	for(i=0;i<p;i++)
	{               
		x1_pii[i] = F77_NAME(ddot)(&n1,&X1[n1*i],&IntOne, &X1[n1*i],&IntOne);
    	x2_pii[i] = F77_NAME(ddot)(&n2,&X2[n2*i],&IntOne, &X2[n2*i],&IntOne);
	}

	
	// Update residuals

    memcpy(resid1,X1,n1*p*sizeof(double));
	memcpy(resid2,X2,n2*p*sizeof(double));

	for(i=0;i<p;i++)
	{
		//memcpy(temp_coef,rho_old[i*p],p*sizeof(double));
		for(j=0;j<p;j++)
			temp_coef[j] = sqrt(wd1[j]/wd1[i])*rho1_old[i+p*j];
		temp_coef[i] = 0;
        F77_NAME(dgemv)(&NoTrans, &n1, &p, &MinusOne, X1, &n1, temp_coef, &IntOne, &One, &resid1[n1*i], &IntOne);
		
		for(j=0;j<p;j++)
			temp_coef[j] = sqrt(wd2[j]/wd2[i])*rho2_old[i+p*j];
		temp_coef[i] = 0;
		F77_NAME(dgemv)(&NoTrans, &n2, &p, &MinusOne, X2, &n2, temp_coef, &IntOne, &One, &resid2[n2*i], &IntOne);
	}

    
	// Block Coordinate descent
	for(k=0;k<niter;k++)
	{

		memcpy(rho1_old,rho1,p*p*sizeof(double));
		memcpy(rho2_old,rho2,p*p*sizeof(double));

		max_diff = 0;
		// update rho
		for(i=0;i<(p-1);i++)
		{
			for(j=i+1;j<p;j++)
			{
					
			// Y_i^T e_j
			x1e_pij = sqrt(wd1[i]/wd1[j])*F77_NAME(ddot)(&n1,&X1[n1*i],&IntOne,&resid1[n1*j],&IntOne);
			x2e_pij = sqrt(wd2[i]/wd2[j])*F77_NAME(ddot)(&n2, &X2[n2*i], &IntOne, &resid2[n2*j], &IntOne);

			// Y_j^T e_i
			x1e_pji = sqrt(wd1[j]/wd1[i])*F77_NAME(ddot)(&n1, &X1[n1*j], &IntOne, &resid1[n1*i], &IntOne);
			x2e_pji = sqrt(wd1[j]/wd1[i])*F77_NAME(ddot)(&n2, &X2[n2*j], &IntOne, &resid2[n2*i], &IntOne);

  			Xe1 = (x1e_pji + x1e_pij);
			XX1 = x1_pii[j]*(wd1[j]/wd1[i]) + x1_pii[i]*(wd1[i]/wd1[j]);

			Xe2 = (x2e_pji  + x2e_pij);
			XX2 =  x2_pii[j]*(wd2[j]/wd2[i]) + x2_pii[i]*(wd2[i]/wd2[j]);

			temp11 = soft(Xe1+rho1_old[i+p*j]*XX1 + lambda2*n1,n1*lambda1)/XX1;
			temp12 = soft(Xe1+rho1_old[i+p*j]*XX1 - lambda2*n1,n1*lambda1)/XX1;
		
			temp21 = soft(Xe2+rho2_old[i+p*j]*XX2 - lambda2*n2,n2*lambda1)/XX2;
			temp22 = soft(Xe2+rho2_old[i+p*j]*XX2 + lambda2*n2,n2*lambda1)/XX2;
			
		
            // rho1
			if(temp11<temp21)
			{ 
				rho1[i+p*j] = temp11;
				rho2[i+p*j] = temp21;
				//printf("\n T11 = %f \t T12 = %f\n",temp11,temp12);
			} else if(temp12>temp22)
			{
				rho1[i+p*j] = temp12;
				rho2[i+p*j] = temp22;
				//printf(" T12 = %f \t T22 = %f\n",temp12,temp22);
			} else
			{
				rho1[i+p*j] = rho2[i+p*j] = soft(Xe1+rho1_old[i+p*j]*XX1+Xe2+rho2_old[i+p*j]*XX2,2*(n1+n2)*lambda1)/(XX1+XX2);
    			//printf(" T31 = %f \t T32 = %f \n",rho1[i+p*j],rho2[i+p*j]);
			}

			if(rho1[i+p*j]>1) rho1[i+p*j]=1;
			if(rho1[i+p*j]<-1) rho1[i+p*j]=-1;
		    if(rho2[i+p*j]>1) rho2[i+p*j]=1;
			if(rho2[i+p*j]<-1) rho2[i+p*j]=-1;
		
			rho1[j+p*i] = rho1[i+p*j];
			rho2[j+p*i] = rho2[i+p*j];


			temp = rho1_old[i+p*j]-rho1[i+p*j];
			if(max_diff<fabs(temp))
			     max_diff = fabs(temp);

			temp = rho2_old[i+p*j]-rho2[i+p*j];
			if(max_diff<fabs(temp))
			     max_diff = fabs(temp);

			//Update residual

			// X_i <- X_j
			temp = sqrt(wd1[j]/wd1[i])*(rho1_old[i+p*j]-rho1[i+p*j]);
			F77_NAME(daxpy)(&n1, &temp, &X1[n1*j],&IntOne, &resid1[n1*i],&IntOne);
				   
			temp = sqrt(wd2[j]/wd1[i])*(rho2_old[i+p*j]-rho2[i+p*j]);
			F77_NAME(daxpy)(&n2, &temp, &X2[n2*j], &IntOne,&resid2[n2*i],&IntOne);

			// X_j <- X_i
			temp = sqrt(wd1[i]/wd1[j])*(rho1_old[i+p*j]-rho1[i+p*j]);
			F77_NAME(daxpy)(&n1, &temp, &X1[n1*i], &IntOne,&resid1[n1*j],&IntOne);
				   
			temp = sqrt(wd2[i]/wd2[j])*(rho2_old[i+p*j]-rho2[i+p*j]);
			F77_NAME(daxpy)(&n2, &temp, &X2[n2*i], &IntOne,&resid2[n2*j],&IntOne);
	   
			}
		}

		if(max_diff<=tol)
				break;
	}
	
	for(i=0;i<p;i++)
	{
		rho1[i+p*i] = 1.0;
		rho2[i+p*i] = 1.0;
	}


	//Rprintf("\n All set - Max difference = %.4f \n",max_diff);

	free(rho1_old);
	free(rho2_old);
	//free(resid1);
	//free(resid2);
	free(temp_coef);
	free(x1_pii);
	free(x2_pii);
}



// Arguments : X1, X2, lambda1, lambda2, CM1, CM2, niter, tol
void dpcid_c(int *N1, int *N2, int *Dim,
	double *X1, double *X2, double *Lam1, double *Lam2,
	double *wd1, double *wd2, double *rho1_init, double *rho2_init, int *Niter, double *Tol,
	double *rho1, double *rho2, double *resid1, double *resid2)
{
	double *rho1_old, *rho2_old;
	double lambda1=*Lam1, lambda2=*Lam2, tol=*Tol;
	double x1e_pij, x2e_pij, *x1_pii, *x2_pii, x1e_pji, x2e_pji;
	//double *resid1, *resid2;
	double *temp_coef, temp=0, temp11, temp12, temp21, temp22;
	double XX1, XX2, Xe1, Xe2;
	double max_diff = 0;
	int n1=*N1, n2=*N2, p=*Dim;

	int niter=*Niter;
	int i, j, k;

	char NoTrans = 'N';
	double One = 1.0, MinusOne = -1.0;
	int IntOne = 1;

	rho1_old = (double *)malloc(p*p*sizeof(double));
	rho2_old = (double *)malloc(p*p*sizeof(double));
	//resid1 = (double *) malloc(n1*p*sizeof(double));
	//resid2 = (double *) malloc(n2*p*sizeof(double));

	// Initializatioin of rho
	memcpy(rho1_old, rho1_init, p*p*sizeof(double));
	memcpy(rho2_old, rho2_init, p*p*sizeof(double));
	memcpy(rho1, rho1_init, p*p*sizeof(double));
	memcpy(rho2, rho2_init, p*p*sizeof(double));

	memset(resid1, 0, n1*p*sizeof(double));
	memset(resid2, 0, n2*p*sizeof(double));

	temp_coef = (double *)malloc(p*sizeof(double));

	memset(temp_coef, 0, p*sizeof(double));

	x1_pii = (double*)malloc(p*sizeof(double));
	x2_pii = (double*)malloc(p*sizeof(double));


	for(i=0;i<p;i++)
	{
		x1_pii[i] = F77_NAME(ddot)(&n1, &X1[n1*i], &IntOne, &X1[n1*i], &IntOne);
		x2_pii[i] = F77_NAME(ddot)(&n2, &X2[n2*i], &IntOne, &X2[n2*i], &IntOne);
	}


	// Update residuals

	memcpy(resid1, X1, n1*p*sizeof(double));
	memcpy(resid2, X2, n2*p*sizeof(double));

	for(i=0;i<p;i++)
	{
		//memcpy(temp_coef,rho_old[i*p],p*sizeof(double));
		for(j=0;j<p;j++)
			temp_coef[j] = sqrt(wd1[j]/wd1[i])*rho1_old[i+p*j];
		temp_coef[i] = 0;
		F77_NAME(dgemv)(&NoTrans, &n1, &p, &MinusOne, X1, &n1, temp_coef, &IntOne, &One, &resid1[n1*i], &IntOne);

		for(j=0;j<p;j++)
			temp_coef[j] = sqrt(wd2[j]/wd2[i])*rho2_old[i+p*j];
		temp_coef[i] = 0;
		F77_NAME(dgemv)(&NoTrans, &n2, &p, &MinusOne, X2, &n2, temp_coef, &IntOne, &One, &resid2[n2*i], &IntOne);
	}


	// Block Coordinate descent
	for(k=0;k<niter;k++)
	{

		memcpy(rho1_old, rho1, p*p*sizeof(double));
		memcpy(rho2_old, rho2, p*p*sizeof(double));

		max_diff = 0;
		// update rho
		for(i=0;i<(p-1);i++)
		{
			for(j=i+1;j<p;j++)
			{

				// Y_i^T e_j
				x1e_pij = sqrt(wd1[i]/wd1[j])*F77_NAME(ddot)(&n1, &X1[n1*i], &IntOne, &resid1[n1*j], &IntOne);
				x2e_pij = sqrt(wd2[i]/wd2[j])*F77_NAME(ddot)(&n2, &X2[n2*i], &IntOne, &resid2[n2*j], &IntOne);

				// Y_j^T e_i
				x1e_pji = sqrt(wd1[j]/wd1[i])*F77_NAME(ddot)(&n1, &X1[n1*j], &IntOne, &resid1[n1*i], &IntOne);
				x2e_pji = sqrt(wd1[j]/wd1[i])*F77_NAME(ddot)(&n2, &X2[n2*j], &IntOne, &resid2[n2*i], &IntOne);

				Xe1 = (x1e_pji + x1e_pij);
				XX1 = x1_pii[j]*(wd1[j]/wd1[i]) + x1_pii[i]*(wd1[i]/wd1[j]);

				Xe2 = (x2e_pji  + x2e_pij);
				XX2 =  x2_pii[j]*(wd2[j]/wd2[i]) + x2_pii[i]*(wd2[i]/wd2[j]);

				temp11 = (Xe1+rho1_old[i+p*j]*XX1 + lambda2*n1)/(XX1+2*lambda1*n1);
				temp12 = (Xe1+rho1_old[i+p*j]*XX1 - lambda2*n1)/(XX1+2*lambda1*n1);

				temp21 = (Xe2+rho2_old[i+p*j]*XX2 - lambda2*n2)/(XX2+2*lambda1*n2);
				temp22 = (Xe2+rho2_old[i+p*j]*XX2 + lambda2*n2)/(XX2+2*lambda1*n2);


				// rho1
				if(temp11<temp21)
				{
					rho1[i+p*j] = temp11;
					rho2[i+p*j] = temp21;
					//printf("\n T11 = %f \t T12 = %f\n",temp11,temp12);
				}
				else if(temp12>temp22)
				{
					rho1[i+p*j] = temp12;
					rho2[i+p*j] = temp22;
					//printf(" T12 = %f \t T22 = %f\n",temp12,temp22);
				}
				else
				{
					rho1[i+p*j] = rho2[i+p*j] = (Xe1+rho1_old[i+p*j]*XX1+Xe2+rho2_old[i+p*j]*XX2)/(XX1+XX2+2*lambda1*n1+2*lambda1*n2);
					//printf(" T31 = %f \t T32 = %f \n",rho1[i+p*j],rho2[i+p*j]);
				}

				if(rho1[i+p*j]>1) rho1[i+p*j]=1;
				if(rho1[i+p*j]<-1) rho1[i+p*j]=-1;
				if(rho2[i+p*j]>1) rho2[i+p*j]=1;
				if(rho2[i+p*j]<-1) rho2[i+p*j]=-1;

				rho1[j+p*i] = rho1[i+p*j];
				rho2[j+p*i] = rho2[i+p*j];


				temp = rho1_old[i+p*j]-rho1[i+p*j];
				if(max_diff<fabs(temp))
					max_diff = fabs(temp);

				temp = rho2_old[i+p*j]-rho2[i+p*j];
				if(max_diff<fabs(temp))
					max_diff = fabs(temp);

				//Update residual

				// X_i <- X_j
				temp = sqrt(wd1[j]/wd1[i])*(rho1_old[i+p*j]-rho1[i+p*j]);
				F77_NAME(daxpy)(&n1, &temp, &X1[n1*j], &IntOne, &resid1[n1*i], &IntOne);

				temp = sqrt(wd2[j]/wd1[i])*(rho2_old[i+p*j]-rho2[i+p*j]);
				F77_NAME(daxpy)(&n2, &temp, &X2[n2*j], &IntOne, &resid2[n2*i], &IntOne);

				// X_j <- X_i
				temp = sqrt(wd1[i]/wd1[j])*(rho1_old[i+p*j]-rho1[i+p*j]);
				F77_NAME(daxpy)(&n1, &temp, &X1[n1*i], &IntOne, &resid1[n1*j], &IntOne);

				temp = sqrt(wd2[i]/wd2[j])*(rho2_old[i+p*j]-rho2[i+p*j]);
				F77_NAME(daxpy)(&n2, &temp, &X2[n2*i], &IntOne, &resid2[n2*j], &IntOne);

			}
		}

		if(max_diff<=tol)
			break;
	}

	for(i=0;i<p;i++)
	{
		rho1[i+p*i] = 1.0;
		rho2[i+p*i] = 1.0;
	}


	//Rprintf("\n All set - Max difference = %.4f \n",max_diff);

	free(rho1_old);
	free(rho2_old);
	//free(resid1);
	//free(resid2);
	free(temp_coef);
	free(x1_pii);
	free(x2_pii);
}
