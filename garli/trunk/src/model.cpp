// GARLI version 0.93 source code
// Copyright  2005 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	Integrative Biology, UT
//	1 University Station, C0930
//	Austin, TX  78712
//  email: zwickl@mail.utexas.edu
//
//	Note: In 2006  moving to NESCENT (The National
//	Evolutionary Synthesis Center) for a postdoc

#include <iostream>

using namespace std;

#include "memchk.h"
#include "utility.h"
#include "parameters.h"
#include "linalg.h"
#include "model.h"
#include "rng.h"

extern rng rnd;
double Model::mutationShape;
double Model::maxPropInvar;

double PointNormal (double prob);
double IncompleteGamma (double x, double alpha, double LnGamma_alpha);
double PointChi2 (double prob, double v);

Model::~Model(){
	if(rates!=NULL) delete []rates;
	if(nst==6){
		delete []eigvals;
		delete []eigvalsimag;
		delete []iwork;
		delete []work;
		delete []col;
		delete []indx;
		delete []c_ijk;
		delete []EigValexp;
		delete []*eigvecs;
		delete []eigvecs;
		delete []*teigvecs;
		delete []teigvecs;
		delete []*inveigvecs;
		delete []inveigvecs;
//		delete []pmat;
		Delete3DArray(pmat);
//		delete []p;
		delete []*qmat;
		delete []qmat;
		delete []*tempqmat;
		delete []tempqmat;
		Delete3DArray(deriv1);
		Delete3DArray(deriv2);
		}
	}

void Model::AllocateEigenVariables(){
	//a bunch of allocation here for all of the qmatrix->eigenvector->pmatrix related variables
	eigvals=new double[nstates];//eigenvalues
	eigvalsimag=new double[nstates];
	iwork=new int[nstates];
	work=new double[nstates];
	col=new double[nstates];
	indx=new int[nstates];
	c_ijk=new double[nstates*nstates*nstates];	
	EigValexp=new double[nstates*nRateCats];	

	//create the matrix for the eigenvectors
	temp=new double[nstates*nstates];
	eigvecs=new double*[nstates];
	for(int x=0;x<nstates;x++){
		eigvecs[x]=&temp[nstates*x];
		}

	//create a temporary matrix to hold the eigenvectors that will be destroyed during the invertization
	temp=new double[nstates*nstates];
	teigvecs=new double*[nstates];
	for(int x=0;x<nstates;x++){
		teigvecs[x]=&temp[nstates*x];
		}

	//create the matrix for the inverse eigenvectors
	temp=new double[nstates*nstates];
	inveigvecs=new double*[nstates];
	for(int x=0;x<nstates;x++){
		inveigvecs[x]=&temp[nstates*x];
		}
	
	//allocate the pmat and p
	pmat=New3DArray<double>(nRateCats, nstates, nstates);
/*	pmat=new double[nstates*nstates*nRateCats];
	p=new double*[nstates];
	for(int i=0;i<nstates;i++){
		p[i]=&pmat[i*nstates];
		}
*/	
	//allocate qmat and tempqmat
	temp=new double[nstates*nstates];
	qmat=new double*[nstates];
	for(int x=0;x<nstates;x++){
		qmat[x]=&temp[nstates*x];
		}
		
	temp=new double[nstates*nstates];
	tempqmat=new double*[nstates];
	for(int x=0;x<nstates;x++){
		tempqmat[x]=&temp[nstates*x];
		}	

	deriv1=New3DArray<double>(nRateCats, nstates, nstates);
	deriv2=New3DArray<double>(nRateCats, nstates, nstates);
	}

void Model::UpdateQMat(){
	//recalculate the qmat from the basefreqs and rates
	
	if(nstates==4){
		qmat[0][1]=rates[0]*pi[1];  //a * piC
		qmat[0][2]=rates[1]*pi[2];  //b * piG
		qmat[0][3]=rates[2]*pi[3];  //c * piT
		qmat[1][2]=rates[3]*pi[2];  //d * piG
		qmat[1][3]=rates[4]*pi[3];  //e * piT
		qmat[2][3]=pi[3];  			//f(=1) * piT 
		qmat[1][0]=rates[0]*pi[0];  //a * piA
		qmat[2][0]=rates[1]*pi[0];  //b * piA
		qmat[2][1]=rates[3]*pi[1];  //d * piC
		qmat[3][0]=rates[2]*pi[0];  //c * piA
		qmat[3][1]=rates[4]*pi[1];  //e * piC
		qmat[3][2]=pi[2]; 			//f(=1) * piG
		}
	else {
		//general nstate x nstate method	
		int rnum=0;
		for(int i=0;i<nstates;i++){
			for(int j=i+1;j<nstates;j++){
				qmat[i][j]=rates[rnum] * pi[j];
				qmat[j][i]=rates[rnum] * pi[i];
				rnum++;
				}
			}
		}
		
	//set diags to sum rows to 0
	double sum;
	for(int x=0;x<nstates;x++){
		sum=0.0;
		for(int y=0;y<nstates;y++){
			if(x!=y) sum+=qmat[x][y];
			}
		qmat[x][x]=-sum;
		}
	}

void Model::CalcEigenStuff(){
				//if rate params or basefreqs have been altered, requiring the recalculation of the eigenvectors and c_ijk
				UpdateQMat();

				memcpy(*tempqmat, *qmat, nstates*nstates*sizeof(double));
				EigenRealGeneral(nstates, tempqmat, eigvals, eigvalsimag, eigvecs, iwork, work);

				memcpy(*teigvecs, *eigvecs, nstates*nstates*sizeof(double));
				InvertMatrix(teigvecs, nstates, col, indx, inveigvecs);
				CalcCijk(c_ijk, nstates, (const double**) eigvecs, (const double**) inveigvecs);
				blen_multiplier=(.5/((qmat[0][1]*pi[0])+(qmat[0][2]*pi[0])+(qmat[0][3]*pi[0])+(qmat[1][2]*pi[1])+(qmat[1][3]*pi[1])+(qmat[2][3]*pi[2])));
				dirty=false;
				}

void Model::CalcPmat(double blen, double *metaPmat, bool flip /*=false*/){
	//this will be a wacky pmat calculation that combines the pmats for all of the rates
	//assert(blen>0.0 && blen<10);
		if(nst==6){
			if(dirty==true){
				//if rate params or basefreqs have been altered, requiring the recalculation of the eigenvectors and c_ijk
				UpdateQMat();

				memcpy(*tempqmat, *qmat, nstates*nstates*sizeof(double));
				EigenRealGeneral(nstates, tempqmat, eigvals, eigvalsimag, eigvecs, iwork, work);

				memcpy(*teigvecs, *eigvecs, nstates*nstates*sizeof(double));
				InvertMatrix(teigvecs, nstates, col, indx, inveigvecs);
				CalcCijk(c_ijk, nstates, (const double**) eigvecs, (const double**) inveigvecs);
				blen_multiplier=(.5/((qmat[0][1]*pi[0])+(qmat[0][2]*pi[0])+(qmat[0][3]*pi[0])+(qmat[1][2]*pi[1])+(qmat[1][3]*pi[1])+(qmat[2][3]*pi[2])));
				dirty=false;
				}
			double tempblen=(blen * blen_multiplier) / (1.0-propInvar);
			CalcPij(c_ijk, nstates, eigvals, 1, tempblen, pmat[0], EigValexp);
			}
			
		else if(nst==2){
/*			//remember that rates[0] is kappa for nst=2 models
			double PI, A;
			if(dirty){//the blen multiplier needs to be recalced if the qmat changes
				R=pi[0]+pi[2];
				Y=1.0 - R;
				blen_multiplier=(-.5/((R*Y)+rates[0]*(*(pi+0)*(*(pi+2))+*(pi+1)*(*(pi+3)))));
				}
			blen*=blen_multiplier;	
			double expblen=exp(blen);
			dirty=false;	
			for(register int f=0;f<4;f++){
				for(register int t=0;t<4;t++){	
					if(f==t){
						if(t==0||t==2) PI = R;
						else PI = Y;
						A=1.0 + PI * (rates[0] - 1.0);
						pmat[f*4+t]=*(pi+t)+*(pi+t)*((1.0/PI)-1.0)*expblen+((PI-*(pi+t))/PI)*exp(A*blen);
						}
					else if((f+t)%2){
						pmat[f*4+t]=(*(pi+t))*(1.0-expblen);//tranversion
						}
					else{
						if(t==0||t==2) PI=R;
						else PI = Y;
						A=1.0 + PI * (rates[0] - 1.0);
						pmat[f*4+t]=*(pi+t)+*(pi+t)*((1.0/PI)-1.0)*expblen-(*(pi+t)/PI)*exp(A*blen);//transition
						}
						
					}
				}
*/			}
		
		if(flip==true){
			assert(0);
			//copy and flip the calculated pmat into the metaPmat
/*			for(int i=0;i<4;i++)
				for(int j=0;j<4;j++){
					metaPmat[i*16+j+r*4]=pmat[i+j*4];
					if(metaPmat[i*16+j+r*4]<0.0){
						ofstream p("pmat.log", ios::app);
						p << "flip" << "\t" << blen << "\t" << metaPmat[i*16+j+r*4] << endl;
						metaPmat[i*16+j+r*4]=0.0;
						}
					//assert(metaPmat[i*16+j+r*4]>=0.0);
					}
*/			}
		else{
			//Copy the pmats into the metaPmat in order
			for(int i=0;i<4;i++)
				for(int j=0;j<4;j++){
					metaPmat[i*4+j]=pmat[0][0][i*4+j];
					if(metaPmat[i*4+j]<0.0){
						ofstream p("pmat.log", ios::app);
						p << "noflip" << "\t" << blen << "\t" << metaPmat[i*4+j] << endl;
						metaPmat[i*4+j]=0.0;
						}
					
					//assert(metaPmat[i*4+j+r*16]>=0.0);
					}
			}
	}	

void Model::CalcPmatRateHet(double blen, double *metaPmat, bool flip /*=false*/){
	//this will be a wacky pmat calculation that combines the pmats for all of the rates
	//assert(blen>0.0 && blen<10);
	for(int r=0;r<nRateCats;r++){
		if(nst==6){
			if(dirty==true)
				CalcEigenStuff();
			double tempblen=(blen * blen_multiplier * gammaRates[r]) / (1.0-propInvar);
			CalcPij(c_ijk, nstates, eigvals, 1, tempblen, pmat[0], EigValexp);
			}
			
/*		else if(nst==2){
			//remember that rates[0] is kappa for nst=2 models
			double PI, A;
			if(dirty){//the blen multiplier needs to be recalced if the qmat changes
				R=pi[0]+pi[2];
				Y=1.0 - R;
				blen_multiplier=(-.5/((R*Y)+rates[0]*(*(pi+0)*(*(pi+2))+*(pi+1)*(*(pi+3)))));
				}
			blen*=blen_multiplier;	
			double expblen=exp(blen);
			dirty=false;	
			for(register int f=0;f<4;f++){
				for(register int t=0;t<4;t++){	
					if(f==t){
						if(t==0||t==2) PI = R;
						else PI = Y;
						A=1.0 + PI * (rates[0] - 1.0);
						pmat[f*4+t]=*(pi+t)+*(pi+t)*((1.0/PI)-1.0)*expblen+((PI-*(pi+t))/PI)*exp(A*blen);
						}
					else if((f+t)%2){
						pmat[f*4+t]=(*(pi+t))*(1.0-expblen);//tranversion
						}
					else{
						if(t==0||t==2) PI=R;
						else PI = Y;
						A=1.0 + PI * (rates[0] - 1.0);
						pmat[f*4+t]=*(pi+t)+*(pi+t)*((1.0/PI)-1.0)*expblen-(*(pi+t)/PI)*exp(A*blen);//transition
						}
						
					}
				}
			}
*/		
		if(flip==true){
			//copy and flip the calculated pmat into the metaPmat
			for(int i=0;i<4;i++)
				for(int j=0;j<4;j++){
					metaPmat[i*16+j+r*4]=pmat[0][0][i+j*4];
					if(metaPmat[i*16+j+r*4]<0.0){
						ofstream p("pmat.log", ios::app);
						p << "flip\t" << blen << "\t" << metaPmat[i*16+j+r*4] << endl;
						metaPmat[i*16+j+r*4]=0.0;
						}
					//assert(metaPmat[i*16+j+r*4]>=0.0);
					}
			}
		else{
			//Copy the pmats into the metaPmat in order
			for(int i=0;i<4;i++)
				for(int j=0;j<4;j++){
					metaPmat[i*4+j+r*16]=pmat[0][0][i*4+j];
					if(metaPmat[i*4+j+r*16]<0.0){
						ofstream p("pmat.log", ios::app);
						p << "noflip\t" << blen << "\t" << metaPmat[i*4+j+r*16] << endl;
						metaPmat[i*4+j+r*16]=0.0;
						}
					
					//assert(metaPmat[i*4+j+r*16]>=0.0);
					}
			}
	}	
	}	

/*
	CalcCijk(c_ijk, nstates, (const double**) eigvecs, (const double**) inveigvecs);
void CalcCijk (double *c_ijk, int n, const double **u, const double **v)

{
		double *pc = c_ijk;
		for (int i=0; i<n; i++)
			for (int j=0; j<n; j++)
				for (int k=0; k<n; k++)
				 	*pc++ = u[i][k] * v[k][j];	// (note: pc = &c[i][j][k]) 
}
*/

void Model::CalcDerivatives(double dlen, double ***&pr, double ***&one, double ***&two){

	if(dirty==true)
		CalcEigenStuff();

	double *EigValderiv=new double[nstates*nRateCats];
	double *EigValderiv2=new double[nstates*nRateCats];

	for(int rate=0;rate<nRateCats;rate++){
		const unsigned rateOffset = nstates*rate; 
		for(int k=0; k<nstates; k++){
			const double scaledEigVal = eigvals[k]*gammaRates[rate]*blen_multiplier/(1.0-propInvar);
			EigValexp[k+rateOffset] = exp(scaledEigVal * dlen);
			EigValderiv[k+rateOffset] = scaledEigVal*EigValexp[k+rateOffset];
			EigValderiv2[k+rateOffset] = scaledEigVal*EigValderiv[k+rateOffset];
			}
		}

	for(int rate=0;rate<nRateCats;rate++)
		{
		const unsigned rateOffset = nstates*rate;
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < 4; j++){
				double sum_p=0.0;
				double sum_d1p=0.0;
				double sum_d2p = 0.0;
				for (int k = 0; k < 4; k++){ 
					//const double x = eigvecs[i][k]*inveigvecs[j][k];
					const double x = eigvecs[i][k]*inveigvecs[k][j];
					sum_p   += x*EigValexp[k+rateOffset];
					sum_d1p += x*EigValderiv[k+rateOffset];
					sum_d2p += x*EigValderiv2[k+rateOffset];
					}

				pmat[rate][i][j] = (sum_p > 0.0 ? sum_p : 0.0);
				deriv1[rate][i][j] = sum_d1p;
				deriv2[rate][i][j] = sum_d2p;
				}
			}
		}
	//try root transformation on fully calced matrices
/*	double tD1=1.0/(4.0*dlen*dlen*dlen);
	double tD2=1.0/(12.0*dlen*dlen);
	for(int i=0;i<64;i++){
		deriv1[0][0][i] /= tD1;
		deriv2[0][0][i] /= tD2;		
		}		
*/	one=deriv1;
	two=deriv2;
	pr=pmat;
	
	delete []EigValderiv;
	delete []EigValderiv2;
	}

/*	

			CalcPij(c_ijk, nstates, eigvals, 1, tempblen, p, EigValexp);
	{

	register int		nsq = n * n;
	double				sum;
	const double *ptr;
	double *pMat = p[0];
	double vr = v * r;
	double *g = EigValexp;
	for (int k=0; k<n; k++)
		*g++ = exp(*eigenValues++ * vr);

	ptr = c_ijk;
#if 1
	for(int i=0; i<nsq; i++){
		g = EigValexp;
		sum = 0.0;
		for(int k=0; k<n; k++)
			sum += (*ptr++) * (*g++);
		*pMat++ = (sum < 0.0) ? 0.0 : sum;
		}
#else
	for(i=0; i<n; i++)
		{
		for(j=0; j<n; j++)
			{
			g = EigValexp;
			sum = 0.0;
			for(k=0; k<n; k++)
				sum += (*ptr++) * (*g++);
			//p[i][j] = (sum < 0.0) ? 0.0 : sum;
			*pMat++ = (sum < 0.0) ? 0.0 : sum;
					}
				}
#endif

			}
		
	
						}
					
*/

void Model::MutateRates(){
	int rateToChange=int(rnd.uniform()*(nst));
	
	if(rateToChange<nst-1){
		rates[rateToChange] *= rnd.gamma( Model::mutationShape );
//		rates[rateToChange] *= exp(MODEL_CHANGE_SCALER * (params->rnd.uniform()-.5));
		if(rates[rateToChange]>99.9) rates[rateToChange]=99.9;
		}

	else{//if we alter the reference rate GT (fixed to 1.0)
		//scale all of the other rates
		//double scaler=exp(MODEL_CHANGE_SCALER * (params->rnd.uniform()-.5));
		double scaler= rnd.gamma( Model::mutationShape );
		for(int i=0;i<nst-1;i++){
			rates[i] /= scaler;
			}
		}

	// don't let rates[0] become greater than 99.0
	// if this upper limit is changed, be sure to check consequences
	// in AllocMigrantStrings function in gamlmain.C (i.e., currently
	// only allow 3 characters plus number needed for precision)
	dirty=true;
	}

void Model::MutateAlpha(){
	
//	alpha *= exp(MODEL_CHANGE_SCALER * (params->rnd.uniform()-.5));
	alpha *=rnd.gamma( Model::mutationShape );
	dirty=true;
	DiscreteGamma();
	}

void Model::MutatePis(){
	//dirichlet type proposal to all the pi's
/*	double dirParams[4];
	double multiplier=1000;
		
	for(int i=0;i<4;i++)
		dirParams[i]=pi[i]*multiplier;
		
	params->rnd.DirichletRandomVariable(dirParams, pi, 4);
*/	
	//alternative:change one pi with a multiplier and rescale the rest
	int piToChange=int(rnd.uniform()*4.0);
	
//	double newPi=pi[piToChange] * exp(MODEL_CHANGE_SCALER * (params->rnd.uniform()-.5));
	double newPi=pi[piToChange] * rnd.gamma( Model::mutationShape );
	for(int b=0;b<4;b++)
		if(b!=piToChange) pi[b] *= (1.0-newPi)/(1.0-pi[piToChange]);
	pi[piToChange]=newPi;
	dirty=true;
	}
	
void Model::MutatePropInvar(){
//	propInvar *= exp(MODEL_CHANGE_SCALER * (params->rnd.uniform()-.5));
	double mult=rnd.gamma( Model::mutationShape );
	if(propInvar == maxPropInvar && (mult > 1.0)) mult=1.0/mult;
	propInvar *= mult;
	propInvar = (propInvar > maxPropInvar ? maxPropInvar : propInvar);
	dirty=true;
	}
	
void Model::CopyModel(const Model *from){
	for(int i=0;i<nst-1;i++)
		rates[i]=from->rates[i];
		
	memcpy(pi, from->pi, sizeof(double)*4);
	alpha=from->alpha;
	memcpy(gammaRates, from->gammaRates, sizeof(double)*4);
	propInvar=from->propInvar;
	dirty=true;
	}	

void Model::SetModel(double *model_string){
	int slot=0;
	for(int i=0;i<nst-1;i++)
		rates[i]=model_string[slot++];
	for(int j=0;j<4;j++)
		pi[j]=model_string[slot++];
	if(nRateCats>1) alpha=model_string[slot++];
	DiscreteGamma();
	//using whether or not this individual had a PI of >0 in the first
	//place to decide whether we should expect one in the string.
	//Seems safe.
	if(propInvar!=0.0) propInvar=model_string[slot++];
	dirty=true;
	}

double Model::TRatio() const{
	double numerator = rates[0] * ( pi[0]*pi[2] + pi[1]*pi[3] );
	double denominator = ( pi[0] + pi[2] ) * ( pi[1] + pi[3] );
	return ( numerator / denominator );
	}

bool Model::IsModelEqual(const Model *other) const {
	//this will need to be generalized if other models are introduced
	for(int i=0;i<nst-1;i++)
		if(rates[i]!=other->rates[i]) return false;
	if(pi[0]!=other->pi[0]) return false;
	if(pi[1]!=other->pi[1]) return false;
	if(pi[2]!=other->pi[2]) return false;
	if(alpha!=other->alpha) return false;
	if(propInvar!=other->propInvar) return false;
	return true;
	}

//a bunch of the gamma rate het machinery
//from MrBayes

/*-------------------------------------------------------------------------------
|                                                                               |
|  Discretization of gamma distribution with equal proportions in each          |
|  category.                                                                    |
|                                                                               |
-------------------------------------------------------------------------------*/  
#define POINTGAMMA(prob,alpha,beta) 		PointChi2(prob,2.0*(alpha))/(2.0*(beta))

/* ------------------------------------------------------------------------------
|                                                                               |
|  Returns z so That Prob{x<z} = prob where x ~ N(0,1) and                      |
|  (1e-12) < prob < 1-(1e-12).  Returns (-9999) if in error.                    |
|                                                                               |
|  Odeh, R. E. and J. O. Evans.  1974.  The percentage points of the normal     |
|     distribution.  Applied Statistics, 22:96-97 (AS70)                        |
|                                                                               |
|  Newer methods:                                                               |
|                                                                               |
|  Wichura, M. J.  1988.  Algorithm AS 241: The percentage points of the        |
|     normal distribution.  37:477-484.                                         |
|  Beasley, JD & S. G. Springer.  1977.  Algorithm AS 111: The percentage       |
|     points of the normal distribution.  26:118-121.                           |
|                                                                               |
-------------------------------------------------------------------------------*/   
double PointNormal (double prob){

	double 		a0 = -0.322232431088, a1 = -1.0, a2 = -0.342242088547, a3 = -0.0204231210245,
 					a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495,
 					b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634,
 					y, z = 0, p = prob, p1;

	p1 = (p<0.5 ? p : 1-p);
	if (p1<1e-20) 
	   return (-9999);
	y = sqrt (log(1/(p1*p1)));   
	z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	return (p<0.5 ? -z : z);

}

/*-------------------------------------------------------------------------------
|                                                                               |
|  Returns the incomplete gamma ratio I(x,alpha) where x is the upper           |
|  limit of the integration and alpha is the shape parameter.  Returns (-1)     |
|  if in error.                                                                 |
|  LnGamma_alpha = ln(Gamma(alpha)), is almost redundant.                      |
|  (1) series expansion     if (alpha>x || x<=1)                                |
|  (2) continued fraction   otherwise                                           |
|                                                                               |
|  RATNEST FORTRAN by                                                           |
|  Bhattacharjee, G. P.  1970.  The incomplete gamma integral.  Applied         |
|     Statistics, 19:285-287 (AS32)                                             |
|                                                                               |
-------------------------------------------------------------------------------*/   
double IncompleteGamma (double x, double alpha, double LnGamma_alpha){
	int 			i;
	double 		p = alpha, g = LnGamma_alpha,
					accurate = 1e-8, overflow = 1e30,
					factor, gin = 0.0, rn = 0.0, a = 0.0, b = 0.0, an = 0.0, 
					dif = 0.0, term = 0.0, pn[6];

	if (x == 0.0) 
		return (0.0);
	if (x < 0 || p <= 0) 
		return (-1.0);

	factor = exp(p*log(x)-x-g);   
	if (x>1 && x>=p) 
		goto l30;
	gin = 1.0;  
	term = 1.0;  
	rn = p;
	l20:
		rn++;
		term *= x/rn;   
		gin += term;
		if (term > accurate) 
			goto l20;
		gin *= factor/p;
		goto l50;
	l30:
		a = 1.0-p;   
		b = a+x+1.0;  
		term = 0.0;
		pn[0] = 1.0;  
		pn[1] = x;  
		pn[2] = x+1;  
		pn[3] = x*b;
		gin = pn[2]/pn[3];
	l32:
		a++;  
		b += 2.0;  
		term++;   
		an = a*term;
		for (i=0; i<2; i++) 
			pn[i+4] = b*pn[i+2]-an*pn[i];
		if (pn[5] == 0) 
			goto l35;
		rn = pn[4]/pn[5];   
		dif = fabs(gin-rn);
		if (dif>accurate) 
			goto l34;
		if (dif<=accurate*rn) 
			goto l42;
	l34:
		gin = rn;
	l35:
		for (i=0; i<4; i++) 
			pn[i] = pn[i+2];
		if (fabs(pn[4]) < overflow) 
			goto l32;
		for (i=0; i<4; i++) 
			pn[i] /= overflow;
		goto l32;
	l42:
		gin = 1.0-factor*gin;
	l50:
		return (gin);

}

inline double LnGamma (double alp){
/*	double cof[6];
	cof[0]=76.18009172947146;
    cof[1]=-86.50532032941677;
    cof[2]=24.01409824083091;
    cof[3]=-1.231739572450155;
    cof[4]=0.1208650973866179e-2;
    cof[5]=-0.5395239384953e-5;	
	double xx=alp;
	double yy=alp;
	double tmp=xx + 5.5 - (xx + 0.5) * log(xx + 5.5);
	double ser = 1.000000000190015;
	for(int j=0;j<5;j++){
		ser += (cof[j] / ++yy);
		}
	return log(2.5066282746310005*ser/xx)-tmp;
	}
*/
	double x = alp, f=0.0, z;
	
	if (x < 7) 
		{
		f = 1.0;  
		z = x-1.0;
		while (++z < 7.0)  
			f *= z;
		x = z;   
		f = -log(f);
		}
	z = 1.0/(x*x);
	return  (f + (x-0.5)*log(x) - x + 0.918938533204673 + 
			(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
			0.083333333333333)/x); 
	}

double PointChi2 (double prob, double v){
	double 		e = 0.5e-6, aa = 0.6931471805, p = prob, g,
					xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0, 
					x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6;

	if (p < 0.000002 || p > 0.999998 || v <= 0.0) 
		return (-1.0);
	g = LnGamma (v/2.0);
	xx = v/2.0;   
	c = xx - 1.0;
	if (v >= -1.24*log(p)) 
		goto l1;
	ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
	if (ch-e<0) 
		return (ch);
	goto l4;
	l1:
		if (v > 0.32) 
			goto l3;
		ch = 0.4;   
		a = log(1.0-p);
	l2:
		q = ch;  
		p1 = 1.0+ch*(4.67+ch);  
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-1.0)-0.01 <= 0.0) 
			goto l4;
		else                       
			goto l2;
	l3: 
		x = PointNormal (p);
		p1 = 0.222222/v;   
		ch = v*pow((x*sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*v+6.0)  
			ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
	l4:
		q = ch;   
		p1 = 0.5*ch;
		if ((t = IncompleteGamma (p1, xx, g)) < 0.0) 
			{
			printf ("\nerr IncompleteGamma");
			return (-1.0);
			}
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));   
		b = t/ch;  
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a))/2520.0;
		s6 = (120.0+c*(346.0+127.0*c))/5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (fabs(q/ch-1.0) > e) 
			goto l4;
		return (ch);

}

//function taken from MB and hard wired for use here	
void Model::DiscreteGamma(){
	bool median=false;
	int 	i;
	double 	gap05 = 1.0/(2.0*nRateCats), t, factor = alpha/alpha*nRateCats, lnga1;

	if (median){
		for (i=0; i<nRateCats; i++) 
			gammaRates[i] = POINTGAMMA((i*2.0+1)*gap05, alpha, alpha);
		for (i=0,t=0; i<nRateCats; i++) 
			t += gammaRates[i];
		for (i=0; i<nRateCats; i++)     
			gammaRates[i] *= factor/t;
		}
	else {
		lnga1 = LnGamma(alpha+1);
		
		//DZ HACK
		//I don't think that these lines are needed, since the frequencies are fixed at .25 anyway.
		for (i=0; i<nRateCats-1; i++) 
			gammaProps[i] = POINTGAMMA((i+1.0)/nRateCats, alpha, alpha);
		for (i=0; i<nRateCats-1; i++) 
			gammaProps[i] = IncompleteGamma(gammaProps[i]*alpha, alpha+1, lnga1);
		
		
		gammaRates[0] = gammaProps[0]*factor;
		gammaRates[nRateCats-1] = (1-gammaProps[nRateCats-2])*factor;
		for (i=1; i<nRateCats-1; i++)  
			gammaRates[i] = (gammaProps[i]-gammaProps[i-1])*factor;
		}
	for (i=0; i<nRateCats; i++) 
		gammaProps[i]=1.0/nRateCats;
	}	
	
	
void Model::OutputPaupBlockForModel(ofstream &outf, const char *treefname){
	outf << "begin paup;\nclear;\ngett file=" << treefname << " storebr;\nlset userbr ";
	if(Nst() == 2) outf << "nst=2 k= " << Rates(0);
	else outf << "nst=6 rmat=(" << Rates(0) << " " << Rates(1) << " " << Rates(2) << " " << Rates(3) << " " << Rates(4);
	outf << ") base=(" << Pi(0) << " " << Pi(1) << " " << Pi(2);
	if(NRateCats()>1) outf << ") rates=gamma shape= " << Alpha();
	if(ProportionInvariant()!=0.0) outf << " pinv= " << ProportionInvariant();
	outf << ";\nend;\n";
	}

void Model::OutputGamlFormattedModel(ostream &outf){
	if(Nst() == 2) outf << "k " << Rates(0);
	else outf << " r " << Rates(0) << " " << Rates(1) << " " << Rates(2) << " " << Rates(3) << " " << Rates(4);
	outf << " b " << Pi(0) << " " << Pi(1) << " " << Pi(2) << " " << Pi(3);
	if(NRateCats()>1) outf << " a " << Alpha();
	if(ProportionInvariant()!=0.0) outf << " p " << ProportionInvariant();
	outf << " ";
	}
/*	
void Model::ReadModelFromFile(NexusToken &token){
	token.GetNextToken();

	do{
		if(token.Equals("r")){//rate parameters
			token.GetNextToken();
			double r[5];
			for(int i=0;i<5;i++){
				r[i]=atof(token.GetToken().c_str());
				token.GetNextToken();
				}
			SetRmat(r);
			if(token.IsNumericalToken()) token.GetNextToken();//this is necessary incase GT is included
			}
		else if(token.Equals("b")){
			token.GetNextToken();
			double b[3];
			for(int i=0;i<3;i++){
				b[i]=atof(token.GetToken().c_str());
				token.GetNextToken();
				}
			SetPis(b);
			if(token.IsNumericalToken()) token.GetNextToken();//this is necessary incase T is included						
			}
		else if(token.Equals("a")){
			token.GetNextToken();
			SetAlpha(atof(token.GetToken().c_str()));
			token.GetNextToken();
			}				
		else if(token.Equals("p")){
			token.GetNextToken();
			SetPinv(atof(token.GetToken().c_str()));
			token.GetNextToken();
			}
		else if(token.Begins("(") == false){
			token.GetNextToken();
			}
		}while(token.Begins("(") == false);
	UpdateQMat();
}
*/	
	
