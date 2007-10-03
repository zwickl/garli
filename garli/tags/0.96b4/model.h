// GARLI version 0.96b4 source code
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

#ifndef _MODEL_
#define _MODEL_

#if !defined(_MSC_VER)
#define _stricmp strcasecmp
#endif

#include <iostream>
#include <cassert>
#include <math.h>
#include <vector>

using namespace std;

#include "rng.h"
#include "mlhky.h"
#include "configoptions.h"
#include "errorexception.h"

class ModelSpecification;
class MFILE;

extern rng rnd;
extern ModelSpecification modSpec;

	enum{//the types
		STATEFREQS = 1,
		RELATIVERATES = 2,
		ALPHASHAPE = 3,
		RATEMULTS = 4,
		RATEPROPS = 5,
		PROPORTIONINVARIANT = 6
		};

class BaseParameter {
protected:
	NxsString name;
	int type;
	int numElements;
	bool fixed;
	FLOAT_TYPE maxv,minv;
	FLOAT_TYPE mutationWeight;
	FLOAT_TYPE mutationProb;
	vector<FLOAT_TYPE*> vals; //this will be aliased to the actual 
						//parameter value within the model class (sneaky!)
	vector<FLOAT_TYPE> default_vals;

public:
	BaseParameter()	{
		numElements=1;
		maxv=minv=0.0;
		}

	BaseParameter(const char *n, FLOAT_TYPE **dv, int t, int numE, FLOAT_TYPE mn, FLOAT_TYPE mx) {
		vals.reserve(6);
		default_vals.reserve(6);
		name=n;
		type=t;
		numElements=numE;
		for(int i=0;i<numElements;i++){
			default_vals.push_back(*(dv[i]));
			vals.push_back(dv[i]);
			}
		minv=mn;
		maxv=mx;
		fixed=false;
		}
	~BaseParameter(){};
/*	void Report(ostream &out){
		out << "Type:\t" << name << "\n";
		if(numElements > 1)
			out << "Current value:\t";
		else 
			out << "Current values:\t";
		for(int i=0;i<numElements;i++)
			out << vals[

		}
*/
	void SetFixed(bool tf) {fixed=tf;}
	bool IsFixed() const {return fixed;}
	int Type() const {return type;}
	void SetWeight(FLOAT_TYPE w){mutationWeight=w;}
	FLOAT_TYPE GetWeight(){return mutationWeight;}
	void SetProb(FLOAT_TYPE p){mutationProb=p;}
	FLOAT_TYPE GetProb(){return mutationProb;}
	virtual void Mutator(FLOAT_TYPE) = 0;
	void SetToDefaultValues(){
		for(int e=0;e<numElements;e++) *vals[e] = default_vals[e];
		}
	};

class StateFrequencies:public BaseParameter{

public:
	StateFrequencies(FLOAT_TYPE **dv, int numE):BaseParameter("Base frequencies", dv, STATEFREQS, numE, 0.0, 1.0){};

	void Mutator(FLOAT_TYPE mutationShape){
		int freqToChange=int(rnd.uniform()*numElements);
		FLOAT_TYPE newFreq=*vals[freqToChange] * rnd.gamma( mutationShape );
		for(int b=0;b<numElements;b++)
			if(b!=freqToChange) *vals[b] *= (FLOAT_TYPE)((1.0-newFreq)/(1.0-*vals[freqToChange]));
		*vals[freqToChange]=newFreq;
		}
	};

class RelativeRates:public BaseParameter{
public:
	// 5/9/06 now enforcing non-zero minimum relative rate to avoid problems in the linear algebra functions
	RelativeRates(const char *c, FLOAT_TYPE **dv, int numE):BaseParameter(c, dv, RELATIVERATES, numE, (FLOAT_TYPE)1.0e-6, (FLOAT_TYPE)999.9){};

	void Mutator(FLOAT_TYPE mutationShape){
		if(numElements > 1){
			int rateToChange=int(rnd.uniform()*(numElements));
			
			if(rateToChange<numElements-1){
				*vals[rateToChange] *= rnd.gamma( mutationShape );
				if(*vals[rateToChange]>maxv) *vals[rateToChange]=maxv;
				if(*vals[rateToChange]<minv) *vals[rateToChange]=minv;
				}

			else{//if we alter the reference rate, which we are assuming
				//is the last one (GT for DNA models, fixed to 1.0)
				//scale all of the other rates
				FLOAT_TYPE scaler= rnd.gamma( mutationShape );
				for(int i=0;i<numElements-1;i++){
					*vals[i] /= scaler;
					if(*vals[i]>maxv) *vals[i]=maxv;
					if(*vals[i]<minv) *vals[i]=minv;
					}
				}
			}
		else {
			*vals[0] *= rnd.gamma( mutationShape );
			if(*vals[0]>maxv) *vals[0]=maxv;
			if(*vals[0]<minv) *vals[0]=minv;
			}
		}
	};

class RateProportions:public BaseParameter{
public:
	RateProportions(FLOAT_TYPE **dv, int numE):BaseParameter("Rate props", dv, RATEPROPS, numE, 0.0, 1.0){};
	void Mutator(FLOAT_TYPE mutationShape){
		int rateToChange=int(rnd.uniform()*(numElements));
		*vals[rateToChange] *= rnd.gamma( mutationShape );
		if(*vals[rateToChange]>maxv) *vals[rateToChange]=maxv;		
		}
	};

class RateMultipliers:public BaseParameter{
public:
	RateMultipliers(FLOAT_TYPE **dv, int numE):BaseParameter("Rate mults", dv, RATEMULTS, numE, (FLOAT_TYPE)0.0, (FLOAT_TYPE)999.9){};
	void Mutator(FLOAT_TYPE mutationShape){
		int rateToChange=int(rnd.uniform()*(numElements));
		*vals[rateToChange] *= rnd.gamma( mutationShape );
		if(*vals[rateToChange]>maxv) *vals[rateToChange]=maxv;
		}
	};

class AlphaShape:public BaseParameter{
public:
	AlphaShape(const char *c, FLOAT_TYPE **dv):BaseParameter(c, dv, ALPHASHAPE, 1, (FLOAT_TYPE)0.0, (FLOAT_TYPE)999.9){};
	void Mutator(FLOAT_TYPE mutationShape){
		*vals[0] *=rnd.gamma( mutationShape );
		}
	};

class ProportionInvariant:public BaseParameter{
public:
	ProportionInvariant(const char *c, FLOAT_TYPE **dv):BaseParameter(c, dv, PROPORTIONINVARIANT, 1, (FLOAT_TYPE)0.0, (FLOAT_TYPE)1.0){};
	void Mutator(FLOAT_TYPE mutationShape){
		*vals[0] *=rnd.gamma( mutationShape );
		}
	};

class ModelSpecification{
	//this will hold the model specification as a global variable
	//so that any models allocated will immediately know what they are
public:
	int nstates;
	int nst;
	int numRateCats;
	int numOmegaCats;

	bool equalStateFreqs;
	bool empiricalStateFreqs;
	bool fixStateFreqs;
	bool fixRelativeRates;

//	bool fixSubstitutionRates;
	bool flexRates;

	bool fixInvariantSites;
	bool fixAlpha;
	bool includeInvariantSites;
	
	bool gotRmatFromFile;
	bool gotStateFreqsFromFile;
	bool gotAlphaFromFile;
	bool gotFlexFromFile;
	bool gotPinvFromFile;

	enum{
		DNA = 0,
		RNA = 1,
		CODON = 2,
		AMINOACID = 3,
		CODONAMINOACID = 4
		}datatype;
	
	enum{
		NONE = 0,
		JONES = 1,
		DAYHOFF = 2,
		POISSON = 3,
		WAG = 4
		}AAmatrix, AAfreqs;
		
	enum{
		STANDARD = 0,
		VERTMITO = 1
		}geneticCode;	

	ModelSpecification(){
		nstates=4;
		//this is the default model
		SetGTR();
		SetGammaRates();
		SetNumRateCats(4, false);
		SetInvariantSites();
		datatype=DNA;
		gotRmatFromFile = gotStateFreqsFromFile = gotAlphaFromFile = gotFlexFromFile = gotPinvFromFile = false;
		geneticCode=STANDARD;
		AAmatrix = NONE;
		AAfreqs = NONE;
		}

	bool IsCodon() {return datatype == CODON;}
	bool IsNucleotide() {return (datatype == DNA || datatype == RNA);}
	bool IsAminoAcid() {return datatype == AMINOACID;}
	bool IsCodonAminoAcid() {return datatype == CODONAMINOACID;}
	//A number of canned model setups
	

	void SetJC(){
		nstates=4;
		SetNst(1);
		SetEqualStateFreqs();
		fixRelativeRates=true;
		}

	void K2P(){
		nstates=4;
		SetNst(2);
		SetEqualStateFreqs();
		fixRelativeRates=false;
		}

	void SetF81(){
		nstates=4;
		SetNst(1);
		SetEstimateStateFreqs();
		fixRelativeRates=true;	
		}

	void SetHKY(){
		nstates=4;
		SetNst(2);
		SetEstimateStateFreqs();
		fixRelativeRates=false;
		}

	void SetGTR(){
		nstates=4;
		SetNst(6);
		SetEstimateStateFreqs();
		fixRelativeRates=false;
		}

	void SetNst(int n){
		nst=n;
		}

	void SetCodon(){
		datatype = CODON;
		nstates = 61;
		numOmegaCats = 1;
		}

	void SetAminoAcid(){
		datatype = AMINOACID;
		nstates = 20;
		}

	void SetCodonAminoAcid(){
		datatype = CODONAMINOACID;
		nstates = 20;
		}

	void SetGammaRates(){
		flexRates=false;
		fixAlpha=false;
		}


	void SetFlexRates(){
		if(includeInvariantSites==true) throw(ErrorException("Sorry, invariant sites models cannot be used with the \"flex\" model of rate heterogeneity"));
		flexRates=true;		
		}

	void SetNumRateCats(int nrates, bool test){//correct behavior here depends on the fact that the default 
		//model includes gamma with 4 rate cats
		if(test ==true){
			if(numRateCats == 1 && nrates > 1)
				throw(ErrorException("ratehetmodel set to \"none\", but numratecats is equal to %d!", nrates));
			if(numRateCats > 1 && nrates == 1){
				if(flexRates == false && fixAlpha == false)
					throw(ErrorException("ratehetmodel set to \"gamma\", but numratecats is equal to 1!"));
				else if(flexRates == false && fixAlpha == true)
					throw(ErrorException("ratehetmodel set to \"gammafixed\", but numratecats is equal to 1!"));
				else if(flexRates == true)
					throw(ErrorException("ratehetmodel set to \"flex\", but numratecats is equal to 1!"));
				}
			}
		
		if(nrates < 1) throw(ErrorException("1 is the minimum value for numratecats."));
		numRateCats=nrates;
		}

	void SetInvariantSites(){
		includeInvariantSites=true;
		fixInvariantSites=false;
		}

	void RemoveInvariantSites(){
		includeInvariantSites=false;
		fixInvariantSites=false;
		}

	void SetEmpiricalStateFreqs(){
		empiricalStateFreqs=fixStateFreqs=true;
		equalStateFreqs=false;
		}

	void SetFixedAlpha(){
		fixAlpha=true;
		}
	void SetFixedInvariantSites(){
		fixInvariantSites=true;
		includeInvariantSites=true;
		}
	void SetEqualStateFreqs(){
		equalStateFreqs=fixStateFreqs=true;
		empiricalStateFreqs=false;
		}
	void SetFixedStateFreqs(){
		fixStateFreqs=true;
		equalStateFreqs=empiricalStateFreqs=false;
		}
	void SetEstimateStateFreqs(){
		equalStateFreqs=fixStateFreqs=empiricalStateFreqs=false;
		}
	void SetFixedRateMatrix(){
		fixRelativeRates=true;
		}
	void SetJonesAAMatrix(){
		AAmatrix = JONES;
		}
	void SetPoissonAAMatrix(){
		AAmatrix = POISSON;
		}	
	void SetDayhoffAAMatrix(){
		AAmatrix = DAYHOFF;
		}
	void SetWAGAAMatrix(){
		AAmatrix = WAG;
		}	
	
	void SetJonesAAFreqs(){
		fixStateFreqs=true;
		equalStateFreqs=empiricalStateFreqs=false;
		AAfreqs = JONES;
		}
	void SetDayhoffAAFreqs(){
		fixStateFreqs=true;
		equalStateFreqs=empiricalStateFreqs=false;
		AAfreqs = DAYHOFF;	
		}
	void SetWAGAAFreqs(){
		fixStateFreqs=true;
		equalStateFreqs=empiricalStateFreqs=false;
		AAfreqs = WAG;	
		}

	bool IsJonesAAFreqs() {return (AAfreqs == JONES);}
	bool IsJonesAAMatrix() {return (AAmatrix == JONES);}
	bool IsDayhoffAAFreqs() {return (AAfreqs == DAYHOFF);}
	bool IsDayhoffAAMatrix() {return (AAmatrix == DAYHOFF);}
	bool IsWAGAAFreqs() {return (AAfreqs == WAG);}
	bool IsWAGAAMatrix() {return (AAmatrix == WAG);}
	bool IsVertMitoCode() {return (geneticCode == VERTMITO);}
	bool IsPoissonAAMatrix() {return (AAmatrix == POISSON);}

	void SetStateFrequencies(const char *str){
		if(_stricmp(str, "equal") == 0) SetEqualStateFreqs();
		else if(_stricmp(str, "estimate") == 0){
			//if((datatype == CODON || datatype == AMINOACID)) throw(ErrorException("Sorry, estimation of equilibrium state frequencies not yet supported with Codon data"));
			SetEstimateStateFreqs();
			}
		else if(_stricmp(str, "empirical") == 0) SetEmpiricalStateFreqs();
		else if(_stricmp(str, "fixed") == 0) SetFixedStateFreqs();
		else if(_stricmp(str, "jones") == 0) SetJonesAAFreqs();
		else if(_stricmp(str, "dayhoff") == 0) SetDayhoffAAFreqs();
		else if(_stricmp(str, "wag") == 0) SetWAGAAFreqs();
		else throw(ErrorException("Unknown setting for statefrequencies: %s\n\t(options are equal, estimate, empirical, fixed, dayhoff, jones, wag)", str));
		}
	void SetRateMatrix(const char *str){
		if(datatype == AMINOACID || datatype == CODONAMINOACID){
			if(_stricmp(str, "jones") == 0) SetJonesAAMatrix();
			else if(_stricmp(str, "dayhoff") == 0) SetDayhoffAAMatrix();
			else if(_stricmp(str, "poisson") == 0) SetPoissonAAMatrix();
			else if(_stricmp(str, "wag") == 0) SetWAGAAMatrix();
			else throw(ErrorException("Sorry, %s is not a valid AminoAcid rate matrix. \n\t(Options are: dayhoff, jones, poisson, wag)", str));
			}
		else{
			if(_stricmp(str, "6rate") == 0) SetNst(6);
			else if(_stricmp(str, "2rate") == 0) SetNst(2);
			else if(_stricmp(str, "1rate") == 0) SetNst(1);
			else if(_stricmp(str, "fixed") == 0) SetFixedRateMatrix();
			else throw(ErrorException("Unknown setting for ratematrix: %s\n\t(options are: 6rate, 2rate, 1rate, fixed)", str));
			}
		}
	void SetProportionInvariant(const char *str){
		if(_stricmp(str, "none") == 0) RemoveInvariantSites();
		else if(datatype == CODON || datatype == AMINOACID) throw(ErrorException("Sorry, invariant sites not yet supported with Codon/Aminoacid data"));
		else if(_stricmp(str, "fixed") == 0) SetFixedInvariantSites();
		else if(_stricmp(str, "estimate") == 0) SetInvariantSites();
		else throw(ErrorException("Unknown setting for proportioninvariant: %s\n\t(options are: none, fixed, estimate)", str));
		}
	void SetRateHetModel(const char *str){
	//	if((datatype != DNA) && (datatype != AMINOACID) && _stricmp(str, "none")) throw(ErrorException("Sorry, rate heterogeneity not yet supported with Codon/Aminoacid data"));
		if(_stricmp(str, "gamma") == 0) SetGammaRates();
		else if(_stricmp(str, "gammafixed") == 0){
			SetGammaRates();
			SetFixedAlpha();		
			}
		else if(_stricmp(str, "flex") == 0) SetFlexRates();
		else if(_stricmp(str, "none") == 0) SetNumRateCats(1, false);
		else throw(ErrorException("Unknown setting for ratehetmodel: %s\n\t(options are: gamma, gammafixed, flex, none)", str));
		}
	void SetDataType(const char *str){
		if(_stricmp(str, "codon") == 0) SetCodon();
		else if(_stricmp(str, "codon-aminoacid") == 0) SetCodonAminoAcid();
		else if(_stricmp(str, "aminoacid") == 0) SetAminoAcid();
		else if(_stricmp(str, "protein") == 0) SetAminoAcid();
		else if(_stricmp(str, "dna") == 0) str;
		else if(_stricmp(str, "rna") == 0) str;
		else if(_stricmp(str, "nucleotide") == 0) str;
		else throw(ErrorException("Unknown setting for datatype: %s\n\t(options are: codon, codon-aminoacid, aminoacid, dna, rna)", str));
		}
	void SetGeneticCode(const char *str){
		if(_stricmp(str, "standard") == 0) geneticCode = STANDARD;
		else if(_stricmp(str, "vertmito") == 0) geneticCode = VERTMITO;
		else throw(ErrorException("Unknown genetic code: %s\n\t(options are: standard, vertmito)", str));
		}

	void SetupModSpec(const GeneralGamlConfig &conf){
		SetDataType(conf.datatype.c_str());
		SetGeneticCode(conf.geneticCode.c_str());
		SetStateFrequencies(conf.stateFrequencies.c_str());
		SetRateMatrix(conf.rateMatrix.c_str());
		SetProportionInvariant(conf.proportionInvariant.c_str());
		SetRateHetModel(conf.rateHetModel.c_str());
		SetNumRateCats(conf.numRateCats, true);
		}
	};

class Model{
	int nst;
	int nstates;

	vector<FLOAT_TYPE*> stateFreqs;
	vector<FLOAT_TYPE*> relNucRates;
	vector<FLOAT_TYPE*> omegas;
	vector<FLOAT_TYPE*> omegaProbs;

	bool eigenDirty;
	FLOAT_TYPE blen_multiplier;
	
	FLOAT_TYPE rateMults[20];
	FLOAT_TYPE rateProbs[20];
	
	FLOAT_TYPE *alpha;
	FLOAT_TYPE *propInvar;

	//variables used for the eigen process if nst=6
	int *iwork, *indx;
	FLOAT_TYPE *eigvals, *eigvalsimag, **eigvecs, **inveigvecs, **teigvecs, *work, *temp, *col, *c_ijk, *EigValexp, *EigValderiv, *EigValderiv2;
	FLOAT_TYPE **qmat, ***pmat;
	FLOAT_TYPE **tempqmat;
	
	//Newton Raphson crap
	FLOAT_TYPE ***deriv1, ***deriv2;

	int *qmatLookup;

	public:
//	static bool noPinvInModel;
//	static bool useFlexRates;
//	static int nRateCats;
	static FLOAT_TYPE mutationShape;
	static FLOAT_TYPE maxPropInvar;

	vector<BaseParameter*> paramsToMutate;

	~Model();

	Model(){
		stateFreqs.reserve(4);
		relNucRates.reserve(6);
		paramsToMutate.reserve(5);
		CreateModelFromSpecification();
		}

	void CalcMutationProbsFromWeights();
	BaseParameter *SelectModelMutation();
	int PerformModelMutation();
	void CreateModelFromSpecification();

	private:
	void AllocateEigenVariables();
	void CalcEigenStuff();

	public:
	void CalcPmat(FLOAT_TYPE blen, FLOAT_TYPE *metaPmat, bool flip =false);
	void CalcPmatNState(FLOAT_TYPE blen, FLOAT_TYPE *metaPmat);
	void CalcDerivatives(FLOAT_TYPE, FLOAT_TYPE ***&, FLOAT_TYPE ***&, FLOAT_TYPE ***&);
	void UpdateQMat();
	void UpdateQMatCodon();
	void UpdateQMatAminoAcid();
	void DiscreteGamma(FLOAT_TYPE *, FLOAT_TYPE *, FLOAT_TYPE);
	bool IsModelEqual(const Model *other) const ;	
	void CopyModel(const Model *from);
	void SetModel(FLOAT_TYPE *model_string);
	void OutputPaupBlockForModel(ofstream &, const char *) const;
	void FillPaupBlockStringForModel(string &str, const char *treefname) const;
	void OutputGarliFormattedModel(ostream &) const;
	void FillGarliFormattedModelString(string &s) const;
	void OutputBinaryFormattedModel(OUTPUT_CLASS &) const;

	void ReadBinaryFormattedModel(FILE *);
	void FillQMatLookup();
	void SetJonesAAFreqs();
	void SetDayhoffAAFreqs();
	void SetWAGAAFreqs();
	void MultiplyByJonesAAMatrix();
	void MultiplyByDayhoffAAMatrix();
	void MultiplyByWAGAAMatrix();

	//model mutations
	void MutateRates();
	void MutatePis();
	void MutateAlpha();
	void MutatePropInvar();
	void MutateRateProbs();
	void MutateRateMults();
	
	//Accessor functions
	FLOAT_TYPE StateFreq(int p) const{ return *stateFreqs[p];}
	FLOAT_TYPE TRatio() const;
	int Nst() const {return nst;}
	FLOAT_TYPE Rates(int r) const { return *relNucRates[r];}
	int NRateCats() const {return modSpec.numRateCats;}
	int NOmegaCats() const {return modSpec.numOmegaCats;}
	FLOAT_TYPE *GetRateMults() {return rateMults;}
	FLOAT_TYPE Alpha() const {return *alpha;}
	FLOAT_TYPE PropInvar() const { return *propInvar;}
	bool NoPinvInModel() const { return ! (modSpec.includeInvariantSites);}
	FLOAT_TYPE MaxPinv() const{return maxPropInvar;}
	int NStates() const {return nstates;}
	int NumMutatableParams() const {return (int) paramsToMutate.size();}

	//Setting things
	void SetDefaultModelParameters(const HKYData *data);
	void SetRmat(FLOAT_TYPE *r, bool checkValidity){
		assert(modSpec.IsAminoAcid() == false);
		if(checkValidity == true){
			if(nst==1){
				if((r[0]==r[1] && r[1]==r[2] &&
					r[2]==r[3] && r[3]==r[4])==false)
					throw(ErrorException("Config file specifies ratematrix = 1rate, but starting model has nonequal rates!\n"));
				}
			if(nst==2){
				if(((r[0]==r[2]  && r[2]==r[3] && r[1]==r[4]))==false)
					throw(ErrorException("Config file specifies ratematrix = 2rate, but starting model does not match!\n"));
				}
			}
		for(int i=0;i<5;i++) *relNucRates[i]=r[i];
		*relNucRates[5]=1.0;
		eigenDirty=true;
		}
	void SetPis(FLOAT_TYPE *b, bool checkValidity){
		//7/12/07 we'll now assume that all freqs have been passed in, rather than calcing the last
		//from the others
		if(checkValidity == true){
			if(modSpec.IsNucleotide()){
				if(modSpec.equalStateFreqs==true && (b[0]==b[1] && b[1]==b[2]) == false) 
					throw(ErrorException("Config file specifies equal statefrequencies,\nbut starting model has nonequal frequencies!\n"));
				if(modSpec.empiricalStateFreqs==true) 
					throw(ErrorException("Config file specifies empirical statefrequencies,\nbut starting model specifies frequencies!\n"));
				}
			}
		FLOAT_TYPE freqTot = 0.0;
		for(int i=0;i<nstates;i++){
			*stateFreqs[i]=b[i];
			freqTot += *stateFreqs[i];
			}
		if(fabs(ONE_POINT_ZERO - freqTot) > (FLOAT_TYPE) 1.0e-5)
			throw(ErrorException("State frequencies do not appear to add up to 1.0!\n"));
		eigenDirty=true;
		}

	void SetFlexRates(FLOAT_TYPE *rates, FLOAT_TYPE *probs){
		if(modSpec.flexRates == false) throw ErrorException("Error: Flex rate values specified in start file,\nbut ratehetmodel is = flex in conf file.");
		for(int r=0;r<NRateCats();r++){
			rateMults[r]=rates[r];
			rateProbs[r]=probs[r];
			}		
		}
	void SetAlpha(FLOAT_TYPE a, bool checkValidity){
		if(checkValidity == true)
			if(modSpec.numRateCats==1) throw(ErrorException("Config file specifies ratehetmodel = none, but starting model contains alpha!\n"));
		*alpha=a;
		DiscreteGamma(rateMults, rateProbs, *alpha);
		//This is odd, but we need to call normalize rates here if we are just using a gamma distrib to get starting rates for 
		//flex.  Flex expects that the rates will be normalized including pinv elsewhere
		if(modSpec.flexRates == true) NormalizeRates();
		}
	void SetPinv(FLOAT_TYPE p, bool checkValidity){
		if(checkValidity == true)
			if(modSpec.includeInvariantSites==false && p!=0.0) throw(ErrorException("Config file specifies invariantsites = none, but starting model contains it!\n"));
		*propInvar=p;
		//change the proportion of rates in each gamma cat
		for(int i=0;i<NRateCats();i++){
			rateProbs[i]=(FLOAT_TYPE)(1.0-*propInvar)/NRateCats();
			}
		}
	void SetMaxPinv(FLOAT_TYPE p){
		Model::maxPropInvar=p;
		}
	void SetDirty(bool tf){
		if(tf) eigenDirty=true;
		else eigenDirty=false;
		}

	void AdjustRateProportions(){
		//this will change the gamma class probs when pinv changes
		for(int i=0;i<NRateCats();i++) rateProbs[i]=(FLOAT_TYPE)(1.0-*propInvar)/NRateCats();
#ifndef NDEBUG
		FLOAT_TYPE sum=0.0;
		for(int i=0;i<NRateCats();i++){
			sum += rateProbs[i];
			}		
		sum += *propInvar;
		assert(fabs(sum - 1.0) < 0.0001);
#endif
		}

	void NormalizeRates(){
		assert(modSpec.flexRates == true);

		FLOAT_TYPE sum=0.0;

		for(int i=0;i<NRateCats();i++){
			sum += rateProbs[i];
			}

		//pinv is figured into the normalization here, but it won't be changed itself.
		if(NoPinvInModel()==false){
			sum = sum / (FLOAT_TYPE)(1.0-*propInvar);
			}	


		for(int i=0;i<NRateCats();i++)	{
			rateProbs[i] /= sum;
			}

		sum=0.0;
		for(int i=0;i<NRateCats();i++){
			sum += rateMults[i]*rateProbs[i];
			}
		for(int i=0;i<NRateCats();i++)	{
			rateMults[i] /= sum;
			}

#ifndef NDEBUG
		sum=0.0;
		for(int i=0;i<NRateCats();i++){
			sum += rateProbs[i];
			}		
		sum += *propInvar;
		assert(fabs(sum - 1.0) < 0.0001);
		sum=0.0;
		for(int i=0;i<NRateCats();i++){
			sum += rateMults[i]*rateProbs[i];
			}
		assert(fabs(sum - 1.0) < 0.0001);

#endif
		}

	const FLOAT_TYPE *GetRateProbs() {
		
/*		FLOAT_TYPE sum=0.0;
		for(int i=0;i<NRateCats();i++){
			sum += rateProbs[i];
			}
		sum+=*propInvar;
		assert(fabs(1.0-sum) < .001);
*/		
		return rateProbs;
		}
	};
	
#endif
