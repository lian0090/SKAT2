#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"


SEXP C_getDL(SEXP R_var_e, SEXP R_taud, SEXP R_d1, SEXP R_n, SEXP R_tU1y, SEXP R_tU1X, SEXP R_tXX, SEXP R_tXy, SEXP R_tyy, SEXP R_tauw, SEXP R_kw, SEXP R_nw,  SEXP R_tU1W, SEXP R_tXW, SEXP R_tWW, SEXP R_tWy, SEXP R_tZtZt, SEXP R_tU1Zt, SEXP R_tXZt, SEXP R_tyZt, SEXP R_tWZt, SEXP R_getQ, SEXP R_getS, SEXP R_getNeg2Log, SEXP R_REML, SEXP R_getAlphaHat)
{
	int nPROTECT=0;
	double var_e=NUMERIC_VALUE(R_var_e);
	double taud=NUMERIC_VALUE(R_taud);
	PROTECT(R_d1=AS_NUMERIC(R_d1));nPROTECT+=1;
	double *d1=NUMERIC_POINTER(R_d1);
	int n=INTEGER_VALUE(R_n);
	PROTECT(R_tU1y=AS_NUMERIC(R_tU1y));nPROTECT+=1;
	double *tU1y=NUMERIC_POINTER(R_tU1y);
	PROTECT(R_tU1X=AS_NUMERIC(R_tU1X)); nPROTECT+=1;
	double *tU1X=NUMERIC_POINTER(R_tU1X);
	PROTECT(R_tXX=AS_NUMERIC(R_tXX));nPROTECT+=1;
	double *tXX=NUMERIC_POINTER(R_tXX);
	PROTECT(R_tXy=AS_NUMERIC(R_tXy));nPROTECT+=1;
	double *tXy=NUMERIC_POINTER(R_tXy);
	double tyy=NUMERIC_VALUE(R_tyy);
	PROTECT(R_tauw=AS_NUMERIC(R_tauw));nPROTECT+=1;
	double *tauw=NUMERIC_POINTER(R_tauw);
	PROTECT(R_kw=AS_NUMERIC(R_kw));nPROTECT+=1;
	double *kw=NUMERIC_POINTER(R_kw);
	PROTECT(R_tU1W=AS_NUMERIC(R_tU1W));nPROTECT+=1;
	double *tU1W=NUMERIC_POINTER(R_tU1W);
	PROTECT(R_tXW=AS_NUMERIC(R_tXW));nPROTECT+=1;
	double *tXW=NUMERIC_POINTER(R_tXW);
	PROTECT(R_tWW=AS_NUMERIC(R_tWW));nPROTECT+=1;
	double *tWW=NUMERIC_POINTER(R_tWW);
	PROTECT(R_tWy=AS_NUMERIC(R_tWy));nPROTECT+=1;
	double *tWy=NUMERIC_POINTER(R_tWy);
	PROTECT(R_tZtZt=AS_NUMERIC(R_tZtZt));nPROTECT+=1;
	double *tZtZt=NUMERIC_POINTER(R_tZtZt);
	PROTECT(R_tU1Zt=AS_NUMERIC(R_tU1Zt));nPROTECT+=1;
	double *tU1Zt=NUMERIC_POINTER(R_tU1Zt);
	PROTECT(R_tXZt=AS_NUMERIC(R_tXZt));nPROTECT+=1;
	double *tXZt=NUMERIC_POINTER(R_tXZt);
	PROTECT(R_tyZt=AS_NUMERIC(R_tyZt));nPROTECT+=1;
	double *tyZt=NUMERIC_POINTER(R_tyZt);
	PROTECT(R_tWZt=AS_NUMERIC(R_tWZt));nPROTECT+=1;
	double *tWZt=NUMERIC_POINTER(R_tWZt);
	int getQ=INTEGER_VALUE(R_getQ);
	int getS=INTEGER_VALUE(R_getS);
	int REML=INTEGER_VALUE(R_REML);
	int getNeg2Log=INTEGER_VALUE(R_getNeg2Log);
    int getAlphaHat=INTEGER_VALUE(R_getAlphaHat);
    int nw=INTEGER_VALUE(R_nw);

	
	//output
	int nlistRout=7;
	SEXP list,listnames; 
	PROTECT(list = allocVector(VECSXP, nlistRout));nPROTECT+=1;
  	PROTECT(listnames = allocVector(VECSXP, nlistRout));nPROTECT+=1;
  	SET_VECTOR_ELT(listnames, 0, PROTECT(mkString("neg2logLik")));nPROTECT+=1;
  	SET_VECTOR_ELT(listnames, 1, PROTECT(mkString("Q")));nPROTECT+=1;
  	SET_VECTOR_ELT(listnames, 2,PROTECT(mkString("lambda")));nPROTECT+=1;
  	SET_VECTOR_ELT(listnames, 3, PROTECT(mkString("S")));nPROTECT+=1;
	SET_VECTOR_ELT(listnames, 4, PROTECT(mkString("sdS")));nPROTECT+=1;
    SET_VECTOR_ELT(listnames, 5, PROTECT(mkString("hat_alpha")));nPROTECT+=1;
    SET_VECTOR_ELT(listnames, 6, PROTECT(mkString("invtXVinvX")));nPROTECT+=1;

	//dimensions
	int kd=length(R_d1);
    int kx= INTEGER(getAttrib(R_tU1X,R_DimSymbol))[1];
		
	 
	int i,j,k;
	double *d_sharp=(double *)R_alloc(kd,sizeof(double));
    double *d_tau=(double *)R_alloc(kd,sizeof(double));
	double *tXVinvX=(double *)R_alloc(kx*kx,sizeof(double));
	double *tXVinvy=(double *)R_alloc(kx,sizeof(double));
	double *tXU1d = (double *)R_alloc(kx*kd,sizeof(double));//tXU1dsharp or tXU1dtau
	double *tWU1d, *tXVdW,*tWVdy,*tWVdW,*Gamma,*Vgamma,*Cgamma,*tehatVdW,*tWVdZt ;//declare these upfront,otherwise, there might be error
   
	for(i=0;i<kd;i++){
		d_sharp[i]=1/(d1[i]*taud+var_e);
	}
	
	if(kd<n){
		for(i=0;i<kd;i++){
			d_tau[i]=d_sharp[i]-1/var_e;
		}
		//tXU1d
		sweep_prod(tU1X,d_tau,tXU1d,kd,kx,'T',2);
		//tXVdX
		memcpy(tXVinvX,tXX,kx*kx*sizeof(double));
		matprod(tXU1d,kx,kd,tU1X,kd,kx,tXVinvX,'N','N',1,1/var_e);
		//tXVdy;
		memcpy(tXVinvy,tXy,kx*sizeof(double));
		matprod(tXU1d,kx,kd,tU1y,kd,1,tXVinvy,'N','N',1,1/var_e);

	}else{
		//tXU1d
		sweep_prod(tU1X,d_sharp,tXU1d,kd,kx,'T',2);
		//tXVdX
		matprod(tXU1d,kx,kd,tU1X,kd,kx,tXVinvX,'N','N',1,0);
		//tXVdy
		matprod(tXU1d,kx,kd,tU1y,kd,1,tXVinvy,'N','N',1,0);
	}
    
    int kwT;//total number of columns for W.
	if(nw>0){
        kwT=0;
        for(i=0;i<nw;i++){
            kwT+=kw[i];
        }
		tWU1d = (double *)R_alloc(kwT*kd,sizeof(double));//tWU1dsharp or tWU1dtau
    	tXVdW = (double *)R_alloc(kx*kwT,sizeof(double));
        tWVdy = (double *)R_alloc(kwT,sizeof(double));
        tWVdW = (double *)R_alloc(kwT*kwT,sizeof(double));

		if(kd<n){
			//tWU1d
			sweep_prod(tU1W,d_tau,tWU1d,kd,kwT,'T',2);
			//tXVdW
			memcpy(tXVdW,tXW,kx*kwT*sizeof(double));
			matprod(tXU1d,kx,kd,tU1W,kd,kwT,tXVdW,'N','N',1,1/var_e);
			//tWVdW
			memcpy(tWVdW,tWW,kwT*kwT*sizeof(double));
			matprod(tWU1d,kwT,kd,tU1W,kd,kwT,tWVdW,'N','N',1,1/var_e);
			//tWVdy
			memcpy(tWVdy,tWy,kwT*sizeof(double));
			matprod(tWU1d,kwT,kd,tU1y,kd,1,tWVdy,'N','N',1,1/var_e);
		}else{
			//tWU1
			sweep_prod(tU1W,d_sharp,tWU1d,kd,kwT,'T',2);
			//tXVdW
			matprod(tXU1d,kx,kd,tU1W,kd,kwT,tXVdW,'N','N',1,0);
			//tWVdW
			matprod(tWU1d,kwT,kd,tU1W,kd,kwT,tWVdW,'N','N',1,0);
			//tWVdy
			matprod(tWU1d,kwT,kd,tU1y,kd,1,tWVdy,'N','N',1,0);
		}

		Gamma =(double *) R_alloc(kwT,sizeof(double));
		int countGamma;
		countGamma=0;
        for(i=0;i<nw;i++){
    		for(j=0;j<kw[i];j++){
    			Gamma[countGamma]=tauw[i];
    			countGamma+=1;
    		}
    	}

    	Vgamma=(double *) R_alloc(kwT*kwT,sizeof(double));
    	sweep_prod(tWVdW,Gamma,Vgamma,kwT,kwT,'N',1);
        for(i=0;i<kwT;i++){
    	Vgamma[i+i*kwT]+=1;
    	}
    	Cgamma = (double *) R_alloc(kwT * kwT, sizeof(double));
    	memcpy(Cgamma, Vgamma, sizeof(double)*kwT * kwT);
        matinv(Cgamma,kwT);
    	sweep_prod(Cgamma,Gamma,Cgamma,kwT,kwT,'N',2);
    	
    	double *tXVdW_Cgamma= (double *) calloc(kx*kwT,sizeof(double));
    	if(tXVdW_Cgamma==NULL)error("not enough memory for tXVdW_Cgamma");
    	double *tXVdW_Cgamma_tWVdX= (double *) calloc(kx*kx,sizeof(double));
    	if(tXVdW_Cgamma_tWVdX==NULL)error("not enough memory for tXVdW_Cgamma_tWVdX");
    	double *tXVdW_Cgamma_tWVdy= (double *) calloc(kx,sizeof(double));
    	if(tXVdW_Cgamma_tWVdy==NULL)error("not enough memory for tXVdW_Cgamma_tWVdy");
    	matprod(tXVdW,kx,kwT,Cgamma,kwT,kwT,tXVdW_Cgamma,'N','N',1,0);
    	matprod(tXVdW_Cgamma,kx,kwT,tXVdW,kx,kwT,tXVdW_Cgamma_tWVdX,'N','T',1,0);
    	matprod(tXVdW_Cgamma,kx,kwT,tWVdy,kwT,1,tXVdW_Cgamma_tWVdy,'N','N',1,0);
    	free(tXVdW_Cgamma);
		for(j=0;j<kx*kx;j++) tXVinvX[j]-=tXVdW_Cgamma_tWVdX[j] ;
    	for(j=0;j<kx;j++)tXVinvy[j]-=tXVdW_Cgamma_tWVdy[j];
    	free(tXVdW_Cgamma_tWVdX);
    	free(tXVdW_Cgamma_tWVdy);    	    	    	    	
	}
    SEXP R_invtXVinvX;
    PROTECT(R_invtXVinvX=allocMatrix(REALSXP,kx,kx)); nPROTECT+=1;
	double *invtXVinvX;
	invtXVinvX=NUMERIC_POINTER(R_invtXVinvX);
	memcpy(invtXVinvX, tXVinvX, sizeof(double)*kx * kx);
    matinv(invtXVinvX,kx);
    //hat_alpha
    SEXP R_hat_alpha;
    PROTECT(R_hat_alpha=allocVector(REALSXP,kx));nPROTECT+=1;
   	double *hat_alpha;
   	hat_alpha=NUMERIC_POINTER(R_hat_alpha);
    matprod(invtXVinvX,kx,kx,tXVinvy,kx,1,hat_alpha,'N','N',1,0);
    if(getAlphaHat){
    SET_VECTOR_ELT(list,5,R_hat_alpha);
    SET_VECTOR_ELT(list,6,R_invtXVinvX);
    }
    
    if(nw>0){
    tehatVdW =(double *)R_alloc(kwT,sizeof(double));
	memcpy(tehatVdW,tWVdy,sizeof(double)*kwT);
	matprod(hat_alpha,kx,1,tXVdW,kx,kwT,tehatVdW,'T','N',-1,1);	
    }
     
    if(getNeg2Log==1){
    	//tU1_ehat
    	double *tU1_ehat;
    	tU1_ehat = (double *)R_alloc(kd,sizeof(double));
    	memcpy(tU1_ehat,tU1y,sizeof(double)*kd);
    	matprod(tU1X,kd,kx,hat_alpha,kx,1,tU1_ehat,'N','N',-1,1);

    	double logDetVd=0;
   		double tehat_Vd_ehat=0;
    	if(kd<n){
   			//logDetVd
   			for(i=0;i<kd;i++){
    			logDetVd+=log(d1[i]*taud+var_e);   
    		}
    		logDetVd+=(n-kd)*log(var_e);
    		
    		//tehat_Vd_ehat
    		for(i=0;i<kx;i++){
    			for(j=0;j<kx;j++){
    				tehat_Vd_ehat+=hat_alpha[i]*tXX[i+j*kx]*hat_alpha[j];
    			}
    			tehat_Vd_ehat-=2*hat_alpha[i]*tXy[i];
    		}
    		tehat_Vd_ehat+=tyy;
    		tehat_Vd_ehat=tehat_Vd_ehat/var_e;
    		for(i=0;i<kd;i++){
    			tehat_Vd_ehat+=pow(tU1_ehat[i],2)*d_tau[i];
    		}
    	}else{
    		//kd=n
    		//logDetVd
   			for(i=0;i<kd;i++){
    			logDetVd+=log(d1[i]*taud+var_e);   
    		}
   			for(i=0;i<kd;i++)tehat_Vd_ehat+=(pow(tU1_ehat[i],2))*d_sharp[i];
    	}
        		

    	double neg2logLik1;
    	double neg2logLik2;
    	double neg2logLik3;
   		if(nw==0){
    		neg2logLik1=logDetVd;
    		neg2logLik3=tehat_Vd_ehat;
    	}else{
    		neg2logLik1=logDetVd+determinant(Vgamma,kwT,1);
    		double tehat_Vd_W_Cgamma_tW_Vd_ehat=0;
    		for(i=0;i<kwT;i++){
    			for(j=0;j<kwT;j++){
    				tehat_Vd_W_Cgamma_tW_Vd_ehat+=tehatVdW[i]*Cgamma[i+j*kwT]*tehatVdW[j];
    			}
    		}
    		neg2logLik3=tehat_Vd_ehat-tehat_Vd_W_Cgamma_tW_Vd_ehat;
		}
		
		neg2logLik2=determinant(tXVinvX,kx,1);
		double neg2logLik;
		if(REML==1){
			neg2logLik=neg2logLik1+neg2logLik2+neg2logLik3;
			//compatible with emma
   			//out<- sum(neg2logLik1,neg2logLik2,neg2logLik3)+(n-kx)*log(2*pi)-log(det(tXX))
		}else{
			neg2logLik=neg2logLik1+neg2logLik3+n*log(2*M_PI);
		}
		
		SET_VECTOR_ELT(list,0,PROTECT(ScalarReal(neg2logLik)));nPROTECT+=1;
		
		}
    
    
	if(getQ==1|getS==1){
        
		if(ISNA(tU1Zt[0]))error("tU1Zt is NULL, cannot get Q and S");
        int kt= INTEGER(getAttrib(R_tU1Zt,R_DimSymbol))[1];

        double *tXVinvZt,*tyVdZt,*tZtVinvZt;
        tXVinvZt=(double *)R_alloc(kx*kt,sizeof(double));
		tyVdZt =(double *)R_alloc(kt,sizeof(double));
		tZtVinvZt =(double *)R_alloc(kt*kt,sizeof(double));
		
		if(kd<n){
			//tXVinvZt
           	 memcpy(tXVinvZt,tXZt,sizeof(double)*kx*kt);
			matprod(tXU1d,kx,kd,tU1Zt,kd,kt,tXVinvZt,'N','N',1,1/var_e);
			//tyVdZt
			for(j=0;j<kt;j++){
				tyVdZt[j]=0;
				for(k=0;k<kd;k++){
					tyVdZt[j]+=tU1y[k]*d_tau[k]*tU1Zt[k+kd*j];
				}
				tyVdZt[j]+=tyZt[j]/var_e;
			}
			//tZtVinvZt
			for(i=0;i<kt;i++){
				for(j=0;j<kt;j++){
					tZtVinvZt[i+kt*j]=0;
					for(k=0;k<kd;k++){
						tZtVinvZt[i+kt*j]+= tU1Zt[k+i*kd]*d_tau[k]*tU1Zt[k+j*kd];
					}
					tZtVinvZt[i+kt*j]+=tZtZt[i+kt*j]/var_e;
				}
			}
					

			//tWVdZt
			if(nw>0){
				tWVdZt =(double *)R_alloc(kwT*kt,sizeof(double));
				memcpy(tWVdZt,tWZt,sizeof(double)*kwT*kt);
				matprod(tWU1d,kwT,kd,tU1Zt,kd,kt,tWVdZt,'N','N',1,1/var_e);	
			}					
		}else{

			//tXVinvZt
			matprod(tXU1d,kx,kd,tU1Zt,kd,kt,tXVinvZt,'N','N',1,0);
			//tyVdZt		
			for(j=0;j<kt;j++){
				tyVdZt[j]=0;
				for(k=0;k<kd;k++){
					tyVdZt[j]+=tU1y[k]*d_sharp[k]*tU1Zt[k+kd*j];
				}
			}
			//tZtVinvZt
			for(i=0;i<kt;i++){
				for(j=0;j<kt;j++){
					tZtVinvZt[i+kt*j]=0;
					for(k=0;k<kd;k++){
						tZtVinvZt[i+kt*j]+= tU1Zt[k+i*kd]*d_sharp[k]*tU1Zt[k+j*kd];
					}
				}
			}
			//tWVdZt
			if(nw>0){
				tWVdZt =(double *)R_alloc(kwT*kt,sizeof(double));
				matprod(tWU1d,kwT,kd,tU1Zt,kd,kt,tWVdZt,'N','N',1,0);	
			}			
		}
       
		
	    double *LQ =(double *)R_alloc(kt,sizeof(double));
	    memcpy(LQ,tyVdZt,sizeof(double)*kt);
		matprod(hat_alpha,kx,1,tXVinvZt,kx,kt,LQ,'T','N',-1,1);
		if(nw>0){
			double *Cgamma_tWVdZt=(double *)calloc(kwT*kt,sizeof(double));
            if(Cgamma_tWVdZt==NULL)error("NULL pointer, no memory available for Cgamma_tWVdZt\n");
			matprod(Cgamma,kwT,kwT,tWVdZt,kwT,kt,Cgamma_tWVdZt,'N','N',1,0);
			matprod(tehatVdW,1,kwT,Cgamma_tWVdZt,kwT,kt,LQ,'N','N',-1,1);
			matprod(tWVdZt,kwT,kt,Cgamma_tWVdZt,kwT,kt,tZtVinvZt,'T','N',-1,1);
			matprod(tXVdW,kx,kwT,Cgamma_tWVdZt,kwT,kt,tXVinvZt,'N','N',-1,1);		
			free(Cgamma_tWVdZt);
		}
       

		//tZtPZt
		double *invtXVinvX_tXVinvZt;
		invtXVinvX_tXVinvZt=(double*) calloc(kx*kt,sizeof(double));
		matprod(invtXVinvX,kx,kx,tXVinvZt,kx,kt,invtXVinvX_tXVinvZt,'N','N',1,0);
		
		double *tZtPZt;
		tZtPZt =(double *) calloc(kt*kt,sizeof(double));
		if(tZtPZt==NULL)error("NULL pointer, no memory available for tZtPZt");
		//tZtPZt =(double *) R_alloc(kt*kt,sizeof(double));

		memcpy(tZtPZt,tZtVinvZt,sizeof(double)*kt*kt);
		
		matprod(tXVinvZt,kx,kt,invtXVinvX_tXVinvZt,kx,kt,tZtPZt,'T','N',-1,1);
		free(invtXVinvX_tXVinvZt);
				

		double Q=0;
		for(i=0;i<kt;i++){Q+=pow(LQ[i],2);}
		Q=Q/2;
		if(getQ){
			SEXP R_lambda;
			
    		PROTECT(R_lambda=allocVector(REALSXP,kt));nPROTECT+=1;
    		
			double *lambda=REAL(R_lambda);
			sEigenValue(tZtPZt,kt,lambda);
			for(i=0;i<kt;i++)lambda[i]=lambda[i]/2;
			//extremly small values of lambda may not be the same with R. for example,  1.192657e-13 versus 1.484921e-14, actually, both of them should be 0.
			SET_VECTOR_ELT(list,1,PROTECT(ScalarReal(Q)));nPROTECT+=1;
			SET_VECTOR_ELT(list,2,R_lambda);
			
		}
		if(getS){
			double S=0;
			for(i=0;i<kt;i++){
				S+=tZtPZt[i+kt*i];
			}
			S=Q-S/2;
			double sdS=0;
			for(i=0;i<kt*kt;i++){
				sdS+=pow(tZtPZt[i],2);
			}
			sdS=sqrtf(sdS/2);
			
			nlistRout+=2;
			SET_VECTOR_ELT(list,3,PROTECT(ScalarReal(S)));nPROTECT+=1;
			SET_VECTOR_ELT(list,4,PROTECT(ScalarReal(sdS)));nPROTECT+=1;
			
		}	
		
	}
		
    setAttrib(list, R_NamesSymbol,listnames);
	UNPROTECT(nPROTECT);
	return(list);
}