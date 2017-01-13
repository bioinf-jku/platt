/* Copyright (C) 2014 Andreas Mayr

/*
target=c(rep(FALSE, 10), rep(TRUE, 20))
out=target+rnorm(length(target), 0, 0.3)
dyn.load("plattScaling.so")
.Call("plattScaling", out, target, as.numeric(sum(!target)), as.numeric(sum(target)))
dyn.unload("plattScaling.so")
 */

#include <R.h>
#include <Rinternals.h>
#include <limits>

extern "C" SEXP plattScaling(SEXP outS, SEXP targetS, SEXP prior0S, SEXP prior1S) {
	int len=length(outS);
	double *deci=REAL(outS);
	int *label=INTEGER(targetS);
	double prior0=REAL(prior0S)[0];
	double prior1=REAL(prior1S)[0];

	int maxiter=100; //Maximum number of iterations
	double minstep=1e-10; //Minimum step taken in line search
	double sigma=1e-12; //Set to any value > 0
	//Construct initial values: target support in array t,
	// initial function value in fval
	double hiTarget=(prior1+1.0)/(prior1+2.0);
	double loTarget=1/(prior0+2.0);
	len=prior1+prior0; // Total number of data
	double* t=(double*)malloc(sizeof(double)*len);
	for(int i=0; i<len; i++) {
		if (label[i] > 0)
			t[i]=hiTarget;
		else
			t[i]=loTarget;
	}
	double A=0.0;
	double B=log((prior0+1.0)/(prior1+1.0));
	double fval=0.0;
	for(int i=0; i<len; i++) {
		double fApB=deci[i]*A+B;
		if (fApB >= 0)
			fval += t[i]*fApB+log(1+exp(-fApB));
		else
			fval += (t[i]-1)*fApB+log(1+exp(fApB));
	}
	for(int it=1; it<maxiter; it++) {
		//Update Gradient and Hessian (use H’ = H + sigma I)
		double g1=0.0;
		double g2=0.0;
		double h11=sigma;
		double h22=sigma;
		double h21=0.0;
		for(int i=0; i<len; i++) {
			double fApB=deci[i]*A+B;
			double p;
			double q;
			if (fApB >= 0) {
				p=exp(-fApB)/(1.0+exp(-fApB));
				q=1.0/(1.0+exp(-fApB));
			}
			else {
				p=1.0/(1.0+exp(fApB));
				q=exp(fApB)/(1.0+exp(fApB));
			}
			double d2=p*q;
			h11 += deci[i]*deci[i]*d2;
			h22 += d2;
			h21 += deci[i]*d2;
			double d1=t[i]-p;
			g1 += deci[i]*d1;
			g2 += d1;
		}
		if (abs(g1)<1e-5 && abs(g2)<1e-5) //Stopping criteria
				break;
		//Compute modified Newton directions
		double det=h11*h22-h21*h21;
		double dA=-(h22*g1-h21*g2)/det;
		double dB=-(-h21*g1+h11*g2)/det;
		double gd=g1*dA+g2*dB;
		double stepsize=1;
		while (stepsize >= minstep){ //Line search
			double newA=A+stepsize*dA;
			double newB=B+stepsize*dB;
			double newf=0.0;
			for(int i=0; i<len; i++) {
				double fApB=deci[i]*newA+newB;
				if (fApB >= 0)
					newf += t[i]*fApB+log(1+exp(-fApB));
				else
					newf += (t[i]-1)*fApB+log(1+exp(fApB));
			}
			if (newf<fval+0.0001*stepsize*gd){
				A=newA;
				B=newB;
				fval=newf;
				break;
			}
			else
				stepsize /= 2.0;
		}
		if (stepsize < minstep){
			break;
		}

		if (A > 0) {
			A=0.0;
		}
	}

	SEXP ARET;
	PROTECT(ARET = allocVector(REALSXP, 1));
	REAL(ARET)[0]=A;
	SEXP BRET;
	PROTECT(BRET = allocVector(REALSXP, 1));
	REAL(BRET)[0]=B;

	SEXP namesRET;
	PROTECT(namesRET = allocVector(STRSXP, 2));
	SET_STRING_ELT(namesRET, 0, mkChar("A"));
	SET_STRING_ELT(namesRET, 1, mkChar("B"));

	SEXP RET;
	PROTECT(RET = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(RET, 0, ARET);
	SET_VECTOR_ELT(RET, 1, BRET);
	setAttrib(RET, R_NamesSymbol, namesRET);

	free(t);

	UNPROTECT(4);
	return(RET);
}
