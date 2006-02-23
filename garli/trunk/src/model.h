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

#ifndef _MODEL_
#define _MODEL_

#include "memchk.h"
#include <iostream>

class Parameters;

class Model{
	int nst;
	int nstates;
	double pi[4];
	double R, Y;  
	double *rates; //this will be kappa in the case of nst=2 or 5 rel rates in the case of nst=6
	double dirty;
	double blen_multiplier;
	
	//rate het stuff
	const int nRateCats;
	double alpha;
	double gammaRates[4];
	double gammaProps[4];
	double propInvar;
	static double maxPropInvar;
	
	//variables used for the eigen process if nst=6
	int *iwork, *indx;
	double *eigvals, *eigvalsimag, **eigvecs, **inveigvecs, **teigvecs, *work, *temp, *col, *c_ijk, *EigValexp;//, **p;
	double **qmat, ***pmat;
	double **tempqmat;
	
	//Newton Raphson crap
	double ***deriv1, ***deriv2;
	
	Parameters *params;
	
	public:
	static double mutationShape;
	Model(int _nst)
		:nRateCats(4)
		{
		nstates=4;
		nst=_nst;
		if(nst==2){
			rates=new double;	
			*rates=2;
			qmat=NULL;
		//	pmat=new double[16];
			}
		else if(nst==6){
			rates=new double[5];
			//make the transitions higher to begin with
			rates[0]=rates[2]=rates[3]=1.0;
			rates[1]=rates[4]=4.0;

			AllocateEigenVariables();
			UpdateQMat();
			}
		else cout << _nst << ": unknown value for nst!" << endl;
		
		//the default base freq starting point is equal, but this
		//will generally be overridden by the empirical freqs in Randomize()
		pi[0]=pi[1]=pi[2]=pi[3]=0.25;
		dirty=true;
		if(nRateCats>1){
			alpha=.5;
			DiscreteGamma();
			}
		propInvar=0.2;
		}
	void SetRmat(double *r){
		for(int i=0;i<5;i++) rates[i]=r[i];
		}
	void SetPis(double *b){
		for(int i=0;i<3;i++) pi[i]=b[i];
		pi[3]=1.0 - pi[0] - pi[1] - pi[2];
		}
	void SetAlpha(double a){
		alpha=a;
		DiscreteGamma();
		}
	void SetPinv(double p){
		propInvar=p;
		}
	void SetMaxPinv(double p){
		Model::maxPropInvar=p;
		}
	double MaxPinv(){
		return maxPropInvar;
		}
	private:
	void AllocateEigenVariables();
	void CalcEigenStuff();
	
	public:
	~Model();
		
	void CalcPmat(double t, double *metaPmat, bool flip=false);
	void CalcPmatRateHet(double blen, double *metaPmat, bool flip =false);
	void CalcDerivatives(double, double ***&, double ***&, double ***&);
	void UpdateQMat();
	void DiscreteGamma();
	
	void CopyModel(const Model *from);
	void SetModel(double *model_string);
	bool IsModelEqual(const Model *other) const ;
	void MutateRates();
	void MutatePis();
	void MutateAlpha();
	void MutatePropInvar();
	double TRatio() const;
	inline double Pi(int p) const{ return pi[p];}
	inline double Alpha() const {return alpha;}
	inline int Nst() const {return nst;}
	inline double Rates(int r) const { return rates[r];}
	inline double ProportionInvariant() const { return propInvar;}
	void SetParams(Parameters *p) {params=p;}
	int NRateCats() const {return nRateCats;}
#ifdef GANESH
    int NStates() const {return nstates;}
#endif

	void SetDirty(bool tf){
		if(tf) dirty=true;
		else dirty=false;
		}
	void OutputPaupBlockForModel(ofstream &outf, const char *treefname);
	void OutputGamlFormattedModel(ostream &outf);
//	void ReadModelFromFile(NexusToken &);
	};
	
#endif
