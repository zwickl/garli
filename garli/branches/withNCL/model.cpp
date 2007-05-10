// GARLI version 0.951 source code
// Copyright  2005-2006 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	National Evolutionary Synthesis Center
//	2024 W. Main Street, Suite A200
//	Durham, NC 27705
//  email: zwickl@nescent.org
//


#include <iostream>
#include <string>

using namespace std;

#include "defs.h"
#include "memchk.h"
#include "utility.h"
#include "linalg.h"
#include "model.h"
#include "individual.h"
#include "mlhky.h"
#include "rng.h"

#undef ALIGN_MODEL

extern rng rnd;
double Model::mutationShape;

double Model::maxPropInvar;
//bool Model::useFlexRates;
//bool Model::noPinvInModel;
//int Model::nRateCats;

double PointNormal (double prob);
double IncompleteGamma (double x, double alpha, double LnGamma_alpha);
double PointChi2 (double prob, double v);

Model::~Model(){
	if(stateFreqs.empty() == false){
		for(int i=0;i<(int)stateFreqs.size();i++)
			delete stateFreqs[i];
		}

	if(relRates.empty() == false){
		if(nst==6){
			for(int i=0;i<(int)relRates.size();i++)
				delete relRates[i];
			}
		else if(nst==2){
			delete relRates[0];
			delete relRates[1];
			}
		else if(nst==1) delete relRates[0];
		}

	if(propInvar != NULL) delete propInvar;

	if(alpha != NULL) delete alpha;

	for(vector<BaseParameter*>::iterator delit=paramsToMutate.begin();delit!=paramsToMutate.end();delit++)
		delete *(delit);

	if(nst==6){
		delete []eigvals;
		delete []eigvalsimag;
		delete []iwork;
		delete []work;
		delete []col;
		delete []indx;
		delete []c_ijk;
		delete []EigValexp;

#ifndef ALIGN_MODEL
		Delete2DArray(eigvecs);
		Delete2DArray(teigvecs);
		Delete2DArray(inveigvecs);
		Delete3DArray(pmat);
		Delete2DArray(qmat);
		Delete2DArray(tempqmat);
		Delete3DArray(deriv1);
		Delete3DArray(deriv2);
#else
		Delete2DAlignedArray(eigvecs);
		Delete2DAlignedArray(teigvecs);
		Delete2DAlignedArray(inveigvecs);
		Delete3DAlignedArray(pmat);
		Delete2DAlignedArray(qmat);
		Delete2DAlignedArray(tempqmat);
		Delete3DAlignedArray(deriv1);
		Delete3DAlignedArray(deriv2);
#endif
		}
	}

void Model::AllocateEigenVariables(){
#ifndef ALIGN_MODEL
	//a bunch of allocation here for all of the qmatrix->eigenvector->pmatrix related variables
	eigvals=new double[nstates];//eigenvalues
	eigvalsimag=new double[nstates];
	iwork=new int[nstates];
	work=new double[nstates];
	col=new double[nstates];
	indx=new int[nstates];
	c_ijk=new double[nstates*nstates*nstates];	
	EigValexp=new double[nstates*NRateCats()];	

	//create the matrix for the eigenvectors
	eigvecs=New2DArray<double>(nstates,nstates);

	//create a temporary matrix to hold the eigenvectors that will be destroyed during the invertization
	teigvecs=New2DArray<double>(nstates,nstates);

	//create the matrix for the inverse eigenvectors
	inveigvecs=New2DArray<double>(nstates,nstates);	

	//allocate the pmat
	pmat=New3DArray<double>(NRateCats(), nstates, nstates);

	//allocate qmat and tempqmat
	qmat=New2DArray<double>(nstates,nstates);
	tempqmat=New2DArray<double>(nstates,nstates);

	deriv1=New3DArray<double>(NRateCats(), nstates, nstates);
	deriv2=New3DArray<double>(NRateCats(), nstates, nstates);
#else

	//a bunch of allocation here for all of the qmatrix->eigenvector->pmatrix related variables
	eigvals=new double[nstates];//eigenvalues
	eigvalsimag=new double[nstates];
	iwork=new int[nstates];
	work=new double[nstates];
	col=new double[nstates];
	indx=new int[nstates];
	c_ijk=new double[nstates*nstates*nstates];	
	EigValexp=new double[nstates*NRateCats()];	

	//create the matrix for the eigenvectors
	eigvecs=New2DAlignedArray<double>(nstates,nstates);

	//create a temporary matrix to hold the eigenvectors that will be destroyed during the invertization
	teigvecs=New2DAlignedArray<double>(nstates,nstates);

	//create the matrix for the inverse eigenvectors
	inveigvecs=New2DAlignedArray<double>(nstates,nstates);	

	//allocate the pmat
	pmat=New3DAlignedArray<double>(NRateCats(), nstates, nstates);

	//allocate qmat and tempqmat
	qmat=New2DAlignedArray<double>(nstates,nstates);
	tempqmat=New2DAlignedArray<double>(nstates,nstates);

	deriv1=New3DAlignedArray<double>(NRateCats(), nstates, nstates);
	deriv2=New3DAlignedArray<double>(NRateCats(), nstates, nstates);
#endif
	}

void Model::UpdateQMat(){
	//recalculate the qmat from the basefreqs and rates
	
	if(nstates==4){
		qmat[0][1]=*relRates[0] * *stateFreqs[1];  //a * piC
		qmat[0][2]=*relRates[1] * *stateFreqs[2];  //b * piG
		qmat[0][3]=*relRates[2] * *stateFreqs[3];  //c * piT
		qmat[1][2]=*relRates[3] * *stateFreqs[2];  //d * piG
		qmat[1][3]=*relRates[4] * *stateFreqs[3];  //e * piT
		qmat[2][3]=*stateFreqs[3];  			//f(=1) * piT 
		qmat[1][0]=*relRates[0] * *stateFreqs[0];  //a * piA
		qmat[2][0]=*relRates[1] * *stateFreqs[0];  //b * piA
		qmat[2][1]=*relRates[3] * *stateFreqs[1];  //d * piC
		qmat[3][0]=*relRates[2] * *stateFreqs[0];  //c * piA
		qmat[3][1]=*relRates[4] * *stateFreqs[1];  //e * piC
		qmat[3][2]=*stateFreqs[2]; 			//f(=1) * piG
		}
	else {
		//general nstate x nstate method	
		int rnum=0;
		for(int i=0;i<nstates;i++){
			for(int j=i+1;j<nstates;j++){
				qmat[i][j]=*relRates[rnum] * *stateFreqs[j];
				qmat[j][i]=*relRates[rnum] * *stateFreqs[i];
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
	blen_multiplier=(.5/((qmat[0][1]**stateFreqs[0])+(qmat[0][2]**stateFreqs[0])+(qmat[0][3]**stateFreqs[0])+(qmat[1][2]**stateFreqs[1])+(qmat[1][3]**stateFreqs[1])+(qmat[2][3]**stateFreqs[2])));
	eigenDirty=false;
	}

void Model::CalcPmat(double blen, double *metaPmat, bool flip /*=false*/){
	//this will be a wacky pmat calculation that combines the pmats for all of the rates
	//assert(blen>0.0 && blen<10);

	//assuming 4 state for now
	double tmpFreqs[4];
	for(int i=0;i<nstates;i++) tmpFreqs[i] = *stateFreqs[i];

	for(int r=0;r<NRateCats();r++){
		if(nst==6){
			if(eigenDirty==true)
				CalcEigenStuff();

			double tempblen;
			if(NoPinvInModel()==true || modSpec.flexRates==true)//if we're using flex rates, pinv should already be included
				//in the rate normalization, and doesn't need to be figured in here
				tempblen=(blen * blen_multiplier * rateMults[r]);
			else
				tempblen=(blen * blen_multiplier * rateMults[r]) / (1.0-*propInvar);
			
			CalcPij(c_ijk, nstates, eigvals, 1, tempblen, pmat[0], EigValexp);
			}
		else if(nst==2 || modSpec.equalStateFreqs == false){
			//remember that relRates[1] is kappa for nst=2 models
			double PI, A, K=*relRates[1];
			double R=tmpFreqs[0]+tmpFreqs[2];
			double Y=1.0 - R;
			blen_multiplier=(.5/((R*Y)+K*((tmpFreqs[0])*((tmpFreqs[2]))+(tmpFreqs[1])*((tmpFreqs[3])))));
			double tempblen ;
			if(NoPinvInModel()==true || modSpec.flexRates==true)//if we're using flex rates, pinv should already be included
				//in the rate normalization, and doesn't need to be figured in here
				tempblen=(blen * blen_multiplier * rateMults[r]);
			else
				tempblen=(blen * blen_multiplier * rateMults[r]) / (1.0-*propInvar);
			double expblen=exp(-tempblen);

			for(register int f=0;f<4;f++){
				for(register int t=0;t<4;t++){	
					if(f==t){
						if(t==0||t==2) PI = R;
						else PI = Y;
						A=1.0 + PI * (K - 1.0);
						(**pmat)[f*4+t]=(tmpFreqs[t])+(tmpFreqs[t])*((1.0/PI)-1.0)*expblen+((PI-(tmpFreqs[t]))/PI)*exp(-A*tempblen);
						assert((**pmat)[f*4+t] > 0.0);
						assert((**pmat)[f*4+t] < 1.0);
						}
					else if((f+t)%2){
						(**pmat)[f*4+t]=((tmpFreqs[t]))*(1.0-expblen);//tranversion
						assert((**pmat)[f*4+t] > 0.0);
						assert((**pmat)[f*4+t] < 1.0);
						}
					else{
						if(t==0||t==2) PI=R;
						else PI = Y;
						A=1.0 + PI * (K - 1.0);
						(**pmat)[f*4+t]=(tmpFreqs[t])+(tmpFreqs[t])*((1.0/PI)-1.0)*expblen-((tmpFreqs[t])/PI)*exp(-A*tempblen);//transition
						assert((**pmat)[f*4+t] > 0.0);
						assert((**pmat)[f*4+t] < 1.0);
						}
					}
				}
			}
		else if(nst==1){
			blen_multiplier=4.0/3.0;
			//	}
			double tempblen ;
			if(NoPinvInModel()==true || modSpec.flexRates==true)//if we're using flex rates, pinv should already be included
				//in the rate normalization, and doesn't need to be figured in here
				tempblen=(blen * blen_multiplier * rateMults[r]);
			else
				tempblen=(blen * blen_multiplier * rateMults[r]) / (1.0-*propInvar);
			double expblen=exp(-tempblen);			
			for(register int f=0;f<4;f++){
				for(register int t=0;t<4;t++){
					if(f==t)
						(**pmat)[f*4+t]=expblen+0.25*(1.0-expblen);
					else 
						(**pmat)[f*4+t]=0.25*(1.0-expblen);
					}
				}
			}
	
		if(flip==true){
			//copy and flip the calculated pmat into the metaPmat
			for(int i=0;i<4;i++)
				for(int j=0;j<4;j++){
					metaPmat[i*16+j+r*4]=pmat[0][0][i+j*4];
					}
			}
		else{
			//Copy the pmats into the metaPmat in order
			for(int i=0;i<4;i++)
				for(int j=0;j<4;j++){
					metaPmat[i*4+j+r*16]=pmat[0][0][i*4+j];
					}
			}
	}	
}	

void Model::CalcDerivatives(double dlen, double ***&pr, double ***&one, double ***&two){

	if(eigenDirty==true)
		CalcEigenStuff();

	double * const EigValderiv=new double[nstates*NRateCats()];
	double * const EigValderiv2=new double[nstates*NRateCats()];

	for(int rate=0;rate<NRateCats();rate++){
		const unsigned rateOffset = nstates*rate; 
		for(int k=0; k<nstates; k++){
			double scaledEigVal;
			if(NoPinvInModel()==true || modSpec.flexRates==true)//if we're using flex rates, pinv should already be included
				//in the rate normalization, and doesn't need to be figured in here
				scaledEigVal = eigvals[k]*rateMults[rate]*blen_multiplier;	
			else
				scaledEigVal = eigvals[k]*rateMults[rate]*blen_multiplier/(1.0-*propInvar);

			EigValexp[k+rateOffset] = exp(scaledEigVal * dlen);
			EigValderiv[k+rateOffset] = scaledEigVal*EigValexp[k+rateOffset];
			EigValderiv2[k+rateOffset] = scaledEigVal*EigValderiv[k+rateOffset];
			}
		}

	for(int rate=0;rate<NRateCats();rate++)
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
	one=deriv1;
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

void Model::SetDefaultModelParameters(const HKYData *data){
	//some of these depend on having read the data already
	//also note that this resets the values in the case of 
	//bootstrapping.  Any of this could be overridden by
	//values specified in a start file

	for(vector<BaseParameter*>::iterator pit=paramsToMutate.begin();pit != paramsToMutate.end();pit++){
		(*pit)->SetToDefaultValues();
		}
	if(modSpec.numRateCats > 1) DiscreteGamma(rateMults, rateProbs, *alpha);

	if(modSpec.equalStateFreqs == false){
		double f[4];
		data->GetEmpiricalFreqs(f);
		SetPis(f, false);
		}

	if(modSpec.includeInvariantSites==false){
		SetPinv(0.0, false);
		SetMaxPinv(0.0);
		}
	else{
		SetPinv(0.25 * ((double)data->NConstant()/(data->NConstant()+data->NInformative()+data->NAutapomorphic())), false);
		SetMaxPinv((double)data->NConstant()/(data->NConstant()+data->NInformative()+data->NAutapomorphic()));
		if(modSpec.flexRates == true) NormalizeRates();
		else AdjustRateProportions();
		}
	}

void Model::MutateRates(){
	//paramsToMutate[1]->Mutator(Model::mutationShape);
	//assert(Rates(0) == Rates(2));

/*	int rateToChange=int(rnd.uniform()*(nst));
	
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
*/	eigenDirty=true;
	}

void Model::MutatePis(){
	assert(0);
	//basetest->Mutator(Model::mutationShape);
//	paramsToMutate[0]->Mutator(Model::mutationShape);
//	dirty=true;

	//alternative:change one pi with a multiplier and rescale the rest
/*	int piToChange=int(rnd.uniform()*4.0);
	
	double newPi=pi[piToChange] * rnd.gamma( Model::mutationShape );
	for(int b=0;b<4;b++)
		if(b!=piToChange) pi[b] *= (1.0-newPi)/(1.0-pi[piToChange]);
	pi[piToChange]=newPi;
	dirty=true;
*/	}
/*
void Model::MutateRateProbs(){
	int ProbToChange=int(rnd.uniform()*(double) NRateCats());
	
	double newProb=rateProbs[ProbToChange] * rnd.gamma( Model::mutationShape / 10.0 );
	for(int b=0;b<NRateCats();b++)
		if(b!=ProbToChange) rateProbs[b] *= (1.0-newProb)/(1.0-rateProbs[ProbToChange]);
	rateProbs[ProbToChange]=newProb;
	NormalizeRates();
	}

void Model::MutateRateMults(){
	int rateToChange=int(rnd.uniform()*NRateCats());
	rateMults[rateToChange] *= rnd.gamma( Model::mutationShape / 10.0);
	NormalizeRates();
	}
	
void Model::MutateAlpha(){
//	alpha *= exp(MODEL_CHANGE_SCALER * (params->rnd.uniform()-.5));
	*alpha *=rnd.gamma( Model::mutationShape );
	DiscreteGamma(rateMults, rateProbs, *alpha);
	//change the proportion of rates in each gamma cat
	}
	
void Model::MutatePropInvar(){
//	propInvar *= exp(MODEL_CHANGE_SCALER * (params->rnd.uniform()-.5));
	double mult=rnd.gamma( Model::mutationShape );
	if(*propInvar == maxPropInvar && (mult > 1.0)) mult=1.0/mult;
	*propInvar *= mult;
	*propInvar = (*propInvar > maxPropInvar ? maxPropInvar : *propInvar);
	//change the proportion of rates in each gamma cat
	for(int i=0;i<NRateCats();i++){
		rateProbs[i]=(1.0-*propInvar)/NRateCats();
		}
	}
*/
void Model::CopyModel(const Model *from){
	for(int i=0;i<6;i++)
		*relRates[i]=*(from->relRates[i]);
	
	for(int i=0;i<nstates;i++)
		*stateFreqs[i]=*(from->stateFreqs[i]);

	//memcpy(pi, from->pi, sizeof(double)*4);

	memcpy(rateMults, from->rateMults, sizeof(double)*NRateCats());
	memcpy(rateProbs, from->rateProbs, sizeof(double)*NRateCats());

	if(modSpec.flexRates == false && modSpec.numRateCats > 1)
		*alpha=*(from->alpha);
	*propInvar=*(from->propInvar);

	eigenDirty=true;
	}	

void Model::SetModel(double *model_string){
	int slot=0;
	for(int i=0;i<nst-1;i++)
		*relRates[i]=model_string[slot++];
	for(int j=0;j<4;j++)
		*stateFreqs[j]=model_string[slot++];
		
	if(NRateCats()>1) *alpha=model_string[slot++];
	DiscreteGamma(rateMults, rateProbs, *alpha);
	//using whether or not this individual had a PI of >0 in the first
	//place to decide whether we should expect one in the string.
	//Seems safe.
	if(*propInvar!=0.0) *propInvar=model_string[slot++];
	eigenDirty=true;
	}

double Model::TRatio() const{
	double numerator = *relRates[1] * ( *stateFreqs[0]**stateFreqs[2] + *stateFreqs[1]**stateFreqs[3] );
	double denominator = ( *stateFreqs[0] + *stateFreqs[2] ) * ( *stateFreqs[1] + *stateFreqs[3] );
	return ( numerator / denominator );
	}

bool Model::IsModelEqual(const Model *other) const {
	//this will need to be generalized if other models are introduced
	for(int i=0;i<6;i++)
		if(*relRates[i]!=*(other->relRates[i])) return false;
	
	for(int i=0;i<nstates;i++)
		if(*stateFreqs[i]!=*(other->stateFreqs[i])) return false;

	if(rateMults[0] != other->rateMults[0]) return false;
	if(rateMults[1] != other->rateMults[1]) return false;
	if(rateMults[2] != other->rateMults[2]) return false;
	if(rateMults[3] != other->rateMults[3]) return false;
	
	if(rateProbs[0] != other->rateProbs[0]) return false;
	if(rateProbs[1] != other->rateProbs[1]) return false;
	if(rateProbs[2] != other->rateProbs[2]) return false;
	if(rateProbs[3] != other->rateProbs[3]) return false;

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
void Model::DiscreteGamma(double *rates, double *props, double shape){
	bool median=false;
	int 	i;
	double 	gap05 = 1.0/(2.0*NRateCats()), t, factor = shape/shape*NRateCats(), lnga1;

	if (median){
		for (i=0; i<NRateCats(); i++) 
			rates[i] = POINTGAMMA((i*2.0+1)*gap05, shape, shape);
		for (i=0,t=0; i<NRateCats(); i++) 
			t += rates[i];
		for (i=0; i<NRateCats(); i++)     
			rates[i] *= factor/t;
		}
	else {
		lnga1 = LnGamma(shape+1);
		
		//DZ HACK
		//I don't think that these lines are needed, since the frequencies are fixed at .25 anyway.
		for (i=0; i<NRateCats()-1; i++) 
			props[i] = POINTGAMMA((i+1.0)/NRateCats(), shape, shape);
		for (i=0; i<NRateCats()-1; i++) 
			props[i] = IncompleteGamma(props[i]*shape, shape+1, lnga1);
		
		
		rates[0] = props[0]*factor;
		rates[NRateCats()-1] = (1-props[NRateCats()-2])*factor;
		for (i=1; i<NRateCats()-1; i++)  
			rates[i] = (props[i]-props[i-1])*factor;
		}
	for (i=0; i<NRateCats(); i++) 
		props[i]=(1.0-*propInvar)/NRateCats();
	}	
	
void Model::OutputPaupBlockForModel(ofstream &outf, const char *treefname) const{
	outf << "begin paup;\nclear;\ngett file=" << treefname << " storebr;\nlset userbr ";
	if(Nst() == 2) outf << "nst=2 trat= " << TRatio();
	else if(Nst() == 1) outf << "nst=1 ";
	else outf << "nst=6 rmat=(" << Rates(0) << " " << Rates(1) << " " << Rates(2) << " " << Rates(3) << " " << Rates(4) << ")";
	
	if(modSpec.equalStateFreqs == true) outf << " base=eq ";
	else if(modSpec.empiricalStateFreqs == true) outf << " base=emp ";
	else outf << " base=(" << StateFreq(0) << " " << StateFreq(1) << " " << StateFreq(2) << ")";
	
	if(modSpec.flexRates==false){
		if(NRateCats()>1) outf << " rates=gamma shape= " << Alpha() << " ncat=" << NRateCats();
		else outf << " rates=equal";
		outf << " pinv= " << PropInvar();
		outf << ";\nend;\n";
		}
	else{
		outf << " pinv= " << PropInvar();
		outf << " [FLEX RATES:\t";
		for(int i=0;i<NRateCats();i++){
			outf << rateMults[i] << "\t";
			outf << rateProbs[i] << "\t";
			}
		outf << "];\nend;\n";
		outf << "[!THIS TREE INFERRED UNDER FLEX RATE MODEL WITH GARLI.\nNO COMPARABLE MODEL IS AVAILABLE IN PAUP!]" << endl;
		}
	}

//this should ONLY be used for screen output, since the precision is pretty low
void Model::CreateGarliFormattedModel(string &s) const{
	char buf[200];
	sprintf(buf, " r %.2f %.2f %.2f %.2f %.2f b %.2f %.2f %.2f %.2f", Rates(0), Rates(1), Rates(2), Rates(3), Rates(4), StateFreq(0), StateFreq(1), StateFreq(2), StateFreq(3));
	s += buf;
	
	if(modSpec.flexRates==true){
		s += " f ";
		for(int i=0;i<NRateCats();i++){
			sprintf(buf, " %.2f %.2f", rateMults[i], rateProbs[i]);
			s += buf;
			}
		}
	else if(NRateCats()>1){
		sprintf(buf, " a %.2f", Alpha());
		s += buf;
		}
	if(PropInvar()!=0.0){
		sprintf(buf, " p %.2f", PropInvar());
		s += buf;		
		}
	s += " ";
	}

void Model::OutputGarliFormattedModel(ostream &outf) const{
	outf << " r " << Rates(0) << " " << Rates(1) << " " << Rates(2) << " " << Rates(3) << " " << Rates(4);
	outf << " b " << StateFreq(0) << " " << StateFreq(1) << " " << StateFreq(2) << " " << StateFreq(3);

	
	if(modSpec.flexRates==true){
		outf << " f ";
		for(int i=0;i<NRateCats();i++){
			outf << " " << rateMults[i] << "\t";
			outf << rateProbs[i] << "\t";
			}
		}
	else{
		if(NRateCats()>1) outf << " a " << Alpha();
		}
	if(PropInvar()!=0.0) outf << " p " << PropInvar();
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
	

void Model::CreateModelFromSpecification(){
	nstates = modSpec.nstates;
	nst = modSpec.nst;

	//deal with rate het models
	propInvar = new double;
	if(modSpec.includeInvariantSites){
		*propInvar=0.2;
		if(modSpec.fixInvariantSites == false){
			ProportionInvariant *pi = new ProportionInvariant("proportion invariant", (double **) &propInvar);
			pi->SetWeight(1);
			paramsToMutate.push_back(pi);
			}			
		}
	else *propInvar=0.0;

	if(NRateCats() > 1){
		alpha = new double;
		*alpha = 0.5;
		
		if(modSpec.flexRates == false){
			DiscreteGamma(rateMults, rateProbs, *alpha);
			if(modSpec.fixAlpha == false){
				AlphaShape *a= new AlphaShape("alpha", &alpha);
				a->SetWeight(1);
				paramsToMutate.push_back(a);
				}
			}
		else{
			//start the flex rates out being equivalent to
			//a gamma with alpha=.5
			DiscreteGamma(rateMults, rateProbs, 0.5);
			if(modSpec.includeInvariantSites == true) NormalizeRates();

			vector<double*> dummy;
			dummy.reserve(NRateCats());
			
			for(int i=0;i<NRateCats();i++)
				dummy.push_back(&rateProbs[i]);
			RateProportions *rateP=new RateProportions(&dummy[0], NRateCats());
			rateP->SetWeight(NRateCats());
			paramsToMutate.push_back(rateP);

			dummy.clear();
			for(int i=0;i<NRateCats();i++)
				dummy.push_back(&rateMults[i]);
			RateMultipliers *rateM=new RateMultipliers(&dummy[0], NRateCats());
			rateM->SetWeight(NRateCats());
			paramsToMutate.push_back(rateM);
			}
		}
	else{
		rateMults[0]=1.0;
		rateProbs[0]=1.0;
		alpha=NULL;
		}

	//deal with the state frequencies
	for(int i=0;i<nstates;i++){
		double *f=new double;
		*f=(1.0/(double) nstates);
		stateFreqs.push_back(f);
		}
	if(modSpec.fixStateFreqs == false){
		StateFrequencies *s=new StateFrequencies(&stateFreqs[0], nstates);
		s->SetWeight(4);
		paramsToMutate.push_back(s);
		}

	//deal with the relative rates
	if(nstates!=4) throw(ErrorException("Sorry, GARLI currently only deals with 4 state models"));
	else{
		if(modSpec.nst==6){
			//make the transitions higher to begin with
			for(int i=0;i<6;i++){
				double *d=new double;
				relRates.push_back(d);
				}
			*relRates[0]=*relRates[2]=*relRates[3]=*relRates[5] = 1.0;
			*relRates[1]=*relRates[4] = 4.0;
			if(modSpec.fixRelativeRates == false){
				RelativeRates *r=new RelativeRates("Rate matrix", &relRates[0], 6);
				r->SetWeight(6);
				paramsToMutate.push_back(r);
				}
			}
		else if(modSpec.nst==2){
			double *a=new double;
			double *b=new double;
			*a=1.0;
			*b=4.0;
			relRates.push_back(a);
			relRates.push_back(b);
			relRates.push_back(a);
			relRates.push_back(a);
			relRates.push_back(b);
			relRates.push_back(a);
			if(modSpec.fixRelativeRates == false){
				RelativeRates *r=new RelativeRates("Rate matrix", &b, 1);
				r->SetWeight(2);
				paramsToMutate.push_back(r);
				}
			}
		else if(modSpec.nst==1){
			double *a=new double;
			*a=1.0;
			for(int i=0;i<6;i++)
				relRates.push_back(a);
			}
		
		AllocateEigenVariables();//these need to be allocated regardless of
		//nst because I don't feel like simplifying the deriv calcs for simpler
		//models.  Pmat calcs for simpler models are simplified, and don't
		//require the Eigen stuff
		UpdateQMat();
		}
	eigenDirty=true;
	}

int Model::PerformModelMutation(){
	if(paramsToMutate.empty()) return 0;
	BaseParameter *mut = SelectModelMutation();
	assert(mut != NULL);
	mut->Mutator(mutationShape);
	int retType;

	if(mut->Type() == RELATIVERATES){
		UpdateQMat();
		retType=Individual::rates;
		eigenDirty=true;
		}
	else if(mut->Type() == STATEFREQS){
		UpdateQMat();
		retType=Individual::pi;
		eigenDirty=true;
		}
	
	else if(mut->Type() == PROPORTIONINVARIANT){
		//this max checking should really be rolled into the parameter class
		*propInvar = (*propInvar > maxPropInvar ? maxPropInvar : *propInvar);
		//the non invariant rates need to be rescaled even if there is only 1
		if(modSpec.flexRates == false) AdjustRateProportions();
		else NormalizeRates();
		retType=Individual::pinv;
		}
	else if(mut->Type() == ALPHASHAPE){
		DiscreteGamma(rateMults, rateProbs, *alpha);
		retType=Individual::alpha;
		}
	else if(mut->Type() == RATEPROPS || mut->Type() == RATEMULTS){
		assert(modSpec.flexRates == true);
		NormalizeRates();
		retType=Individual::alpha;
		}
	
	return retType;
	}


BaseParameter *Model::SelectModelMutation(){
	CalcMutationProbsFromWeights();
	if(paramsToMutate.empty() == true) return NULL;
	double r=rnd.uniform();
	vector<BaseParameter*>::iterator it;
	for(it=paramsToMutate.begin();it!=paramsToMutate.end();it++){
		if((*it)->GetProb() > r) return *it;
		}
	it--;
	return *it;
	}

void Model::CalcMutationProbsFromWeights(){
	double tot=0.0, running=0.0;
	for(vector<BaseParameter*>::iterator it=paramsToMutate.begin();it!=paramsToMutate.end();it++){
		tot += (*it)->GetWeight();
		}
	for(vector<BaseParameter*>::iterator it=paramsToMutate.begin();it!=paramsToMutate.end();it++){
		running += (*it)->GetWeight() / tot;
		(*it)->SetProb(running);
		}
	}
