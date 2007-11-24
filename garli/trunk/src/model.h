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
#include "sequencedata.h"
#include "configoptions.h"
#include "errorexception.h"

class ModelSpecification;
class MFILE;

extern rng rnd;
extern ModelSpecification modSpec;
extern bool FloatingPointEquals(const FLOAT_TYPE first, const FLOAT_TYPE sec, const FLOAT_TYPE epsilon);

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

		FLOAT_TYPE newTot = 1.0 - *vals[rateToChange];
		FLOAT_TYPE oldTot = 0.0;
		for(int i=0;i<numElements;i++)
			if(i != rateToChange) oldTot += *vals[i];
		for(int i=0;i<numElements;i++)
			if(i != rateToChange) *vals[i] *= newTot / oldTot;
		newTot = 0.0;
		for(int i=0;i<numElements;i++) newTot += *vals[i];
		assert(fabs(newTot - 1.0) < 0.0001);
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
	int numRateCats;

	bool fixStateFreqs;
	bool fixRelativeRates;

//	bool fixSubstitutionRates;
	//bool flexRates;

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
		EQUAL = 0,
		EMPIRICAL = 1,
		ESTIMATE = 2,
		F1X4 = 3,
		F3X4 = 4,
		JONES = 5,
		DAYHOFF = 6,
		WAG = 8,
		USERSPECIFIED = 20
		}stateFrequencies;

	enum{
		NST1 = 0,
		NST2 = 1,
		NST6 = 2,
		JONESMAT = 5,
		DAYHOFFMAT = 6,
		POISSON = 7,
		WAGMAT = 8,
		USERSPECIFIEDMAT = 20
		}rateMatrix;
	
	enum{
		NONE = 0,
		GAMMA = 1,
		FLEX  = 2,
		NONSYN = 3
		}rateHetType;

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
		}

	bool IsCodon() {return datatype == CODON;}
	bool IsNucleotide() {return (datatype == DNA || datatype == RNA);}
	bool IsAminoAcid() {return (datatype == AMINOACID || datatype == CODONAMINOACID);}//for most purposes codon-aminoacid should be considered AA
	bool IsCodonAminoAcid() {return datatype == CODONAMINOACID;}
	//A number of canned model setups
	
	void SetJC(){
		nstates=4;
		rateMatrix = NST1;
		SetEqualStateFreqs();
		fixRelativeRates=true;
		}

	void K2P(){
		nstates=4;
		rateMatrix = NST2;
		SetEqualStateFreqs();
		fixRelativeRates=false;
		}

	void SetF81(){
		nstates=4;
		rateMatrix = NST1;
		SetEstimateStateFreqs();
		fixRelativeRates=true;	
		}

	void SetHKY(){
		nstates=4;
		rateMatrix = NST2;
		SetEstimateStateFreqs();
		fixRelativeRates=false;
		}

	void SetGTR(){
		nstates=4;
		rateMatrix = NST6;
		SetEstimateStateFreqs();
		fixRelativeRates=false;
		}

	void SetCodon(){
		datatype = CODON;
		rateMatrix = NST2;
		stateFrequencies = EQUAL;
		nstates = 61;
		numRateCats = 1;
		fixRelativeRates=false;
		}

	void SetAminoAcid(){
		datatype = AMINOACID;
		rateMatrix = WAGMAT;
		stateFrequencies = WAG;
		nstates = 20;
		fixRelativeRates=true;
		}

	void SetCodonAminoAcid(){
		datatype = CODONAMINOACID;
		rateMatrix = WAGMAT;
		stateFrequencies = WAG;
		nstates = 20;
		fixRelativeRates=true;
		}

	void SetGammaRates(){
		rateHetType = GAMMA;
		fixAlpha=false;
		}

	void SetFlexRates(){
		if(includeInvariantSites==true) throw(ErrorException("Invariant sites models should not be (and don't need to be) used\n     with the \"flex\" model of rate heterogeneity, since flex is able to\n     incorporate a class of sites with a rate of effectively zero"));
		rateHetType = FLEX;
		}

	void SetNumRateCats(int nrates, bool test){//correct behavior here depends on the fact that the default 
		//model includes gamma with 4 rate cats
		if(test ==true){
			if(rateHetType == NONE && nrates > 1)
				throw(ErrorException("ratehetmodel set to \"none\", but numratecats is equal to %d!", nrates));
			if(rateHetType != NONE && nrates == 1){
				if(rateHetType == GAMMA && fixAlpha == false)
					throw(ErrorException("ratehetmodel set to \"gamma\", but numratecats is equal to 1!"));
				else if(rateHetType == GAMMA && fixAlpha == true)
					throw(ErrorException("ratehetmodel set to \"gammafixed\", but numratecats is equal to 1!"));
				else if(rateHetType == FLEX)
					throw(ErrorException("ratehetmodel set to \"flex\", but numratecats is equal to 1!"));
				else if(rateHetType == NONSYN)
					throw(ErrorException("ratehetmodel set to \"nonsynonymous\", but numratecats is equal to 1!"));				
				}
			}
		
		if(nrates < 1) throw(ErrorException("1 is the minimum value for numratecats."));
		if(nrates > 20) throw(ErrorException("20 is the maximum value for numratecats."));
		numRateCats=nrates;
		}

	void SetInvariantSites(){
		if(rateHetType == NONSYN)
			throw(ErrorException("invariant sites cannot be used with nonsynonymous rate variation"));			
		includeInvariantSites=true;
		fixInvariantSites=false;
		}

	void RemoveInvariantSites(){
		includeInvariantSites=false;
		fixInvariantSites=false;
		}

	void SetEmpiricalStateFreqs(){
		stateFrequencies = EMPIRICAL;
		fixStateFreqs=true;
		}

	void SetEqualStateFreqs(){
		stateFrequencies = EQUAL;
		fixStateFreqs=true;
		}
	void SetUserSpecifiedStateFreqs(){
		stateFrequencies = USERSPECIFIED;
		fixStateFreqs=true;
		}
	void SetEstimateStateFreqs(){
		stateFrequencies = ESTIMATE;
		fixStateFreqs=false;
		}
	void SetF1X4StateFreqs(){
		stateFrequencies = F1X4;
		fixStateFreqs = true;
		}
	void SetF3X4StateFreqs(){
		stateFrequencies = F3X4;
		fixStateFreqs = true;
		}
	void SetFixedAlpha(){
		fixAlpha=true;
		}
	void SetFixedInvariantSites(){
		fixInvariantSites=true;
		includeInvariantSites=true;
		}
	void SetUserSpecifiedRateMatrix(){
		rateMatrix = USERSPECIFIEDMAT;
		fixRelativeRates=true;
		}
	void SetJonesAAMatrix(){
		rateMatrix = JONESMAT;
		fixRelativeRates=true;
		}
	void SetPoissonAAMatrix(){
		rateMatrix = POISSON;
		fixRelativeRates=true;
		}	
	void SetDayhoffAAMatrix(){
		rateMatrix = DAYHOFFMAT;
		fixRelativeRates=true;
		}
	void SetWAGAAMatrix(){
		rateMatrix = WAGMAT;
		fixRelativeRates=true;
		}	
	
	void SetJonesAAFreqs(){
		stateFrequencies = JONES;
		fixStateFreqs=true;
		}
	void SetDayhoffAAFreqs(){
		stateFrequencies = DAYHOFF;
		fixStateFreqs=true;	
		}
	void SetWAGAAFreqs(){
		stateFrequencies = WAG;
		fixStateFreqs=true;
		}

	void SetM3Model(){
		rateHetType = NONSYN;
		numRateCats = 3;
		}

	int Nst() const {
		if(rateMatrix == NST1) return 1;
		else if(rateMatrix == NST2) return 2;
		else if(rateMatrix == NST6) return 6;
		else if(rateMatrix == USERSPECIFIEDMAT && datatype != AMINOACID && datatype != CODONAMINOACID) return 6;
		else assert(0);
		return -1;
		}

	bool IsJonesAAFreqs() {return (stateFrequencies == JONES);}
	bool IsJonesAAMatrix() {return (rateMatrix == JONES);}
	bool IsDayhoffAAFreqs() {return (stateFrequencies == DAYHOFF);}
	bool IsDayhoffAAMatrix() {return (rateMatrix == DAYHOFFMAT);}
	bool IsWAGAAFreqs() {return (stateFrequencies == WAG);}
	bool IsWAGAAMatrix() {return (rateMatrix == WAGMAT);}
	bool IsVertMitoCode() {return (geneticCode == VERTMITO);}
	bool IsPoissonAAMatrix() {return (rateMatrix == POISSON);}
	bool IsUserSpecifiedRateMatrix(){return rateMatrix == USERSPECIFIEDMAT;}

	bool IsEqualStateFrequencies(){return stateFrequencies == EQUAL;}
	bool IsEmpiricalStateFrequencies(){return stateFrequencies == EMPIRICAL;}
	bool IsUserSpecifiedStateFrequencies(){return stateFrequencies == USERSPECIFIED;}
	bool IsF3x4StateFrequencies(){return stateFrequencies == F3X4;}
	bool IsF1x4StateFrequencies(){return stateFrequencies == F1X4;}

	bool IsFlexRateHet(){return rateHetType == FLEX;}
	bool IsGammaRateHet(){return rateHetType == GAMMA;}
	bool IsNonsynonymousRateHet(){return rateHetType == NONSYN;}

	void SetStateFrequencies(const char *str){
		if(_stricmp(str, "equal") == 0) SetEqualStateFreqs();
		else if(_stricmp(str, "estimate") == 0) SetEstimateStateFreqs();
		else if(_stricmp(str, "empirical") == 0) SetEmpiricalStateFreqs();
		else if(_stricmp(str, "fixed") == 0) SetUserSpecifiedStateFreqs();
		else if(datatype == CODON && _stricmp(str, "f1x4") == 0) SetF1X4StateFreqs();
		else if(datatype == CODON && _stricmp(str, "f3x4") == 0) SetF3X4StateFreqs();
		else if((datatype == AMINOACID || datatype == CODONAMINOACID) && _stricmp(str, "jones") == 0) SetJonesAAFreqs();
		else if((datatype == AMINOACID || datatype == CODONAMINOACID) && _stricmp(str, "dayhoff") == 0) SetDayhoffAAFreqs();
		else if((datatype == AMINOACID || datatype == CODONAMINOACID) && _stricmp(str, "wag") == 0) SetWAGAAFreqs();
		else throw(ErrorException("Invalid setting for statefrequencies: %s\n\tOptions for all datatypes: equal, estimate, empirical, fixed\n\tFor aminoacid datatype only: dayhoff, jones, wag\n\tFor codon datatype only: F1X4, F3X4", str));
		}
	void SetRateMatrix(const char *str){
		if(datatype == AMINOACID || datatype == CODONAMINOACID){
			if(_stricmp(str, "jones") == 0) SetJonesAAMatrix();
			else if(_stricmp(str, "dayhoff") == 0) SetDayhoffAAMatrix();
			else if(_stricmp(str, "poisson") == 0) SetPoissonAAMatrix();
			else if(_stricmp(str, "wag") == 0) SetWAGAAMatrix();
			else throw(ErrorException("Sorry, %s is not a valid aminoacid rate matrix. \n\t(Options are: dayhoff, jones, poisson, wag)", str));
			}
		else{
			if(_stricmp(str, "6rate") == 0) rateMatrix = NST6;
			else if(_stricmp(str, "2rate") == 0) rateMatrix = NST2;
			else if(_stricmp(str, "1rate") == 0) rateMatrix = NST1;
			else if(_stricmp(str, "fixed") == 0) SetUserSpecifiedRateMatrix();
			else throw(ErrorException("Unknown setting for ratematrix: %s\n\t(options are: 6rate, 2rate, 1rate, fixed)", str));
			}
		}
	void SetProportionInvariant(const char *str){
		if(_stricmp(str, "none") == 0) RemoveInvariantSites();
		//else if(datatype == CODON || datatype == AMINOACID) throw(ErrorException("Sorry, invariant sites not yet supported with Codon/Aminoacid data"));
		else if(datatype == CODON){
			if(_stricmp(str, "fixed") == 0 || _stricmp(str, "estimate") == 0) 
				throw ErrorException("Invariant sites cannot be used with codon models.\n     Try ratehetmodel = nonsynonymous to allow dN/dS variation across sites");
			else throw(ErrorException("Unknown setting for proportioninvariant: %s\n\t(only valid option for codon models is none)", str));
			}
		else if(_stricmp(str, "fixed") == 0) SetFixedInvariantSites();
		else if(_stricmp(str, "estimate") == 0) SetInvariantSites();
		else throw(ErrorException("Unknown setting for proportioninvariant: %s\n\t(options are: none, fixed, estimate)", str));
		}
	void SetRateHetModel(const char *str){
	//	if((datatype != DNA) && (datatype != AMINOACID) && _stricmp(str, "none")) throw(ErrorException("Sorry, rate heterogeneity not yet supported with Codon/Aminoacid data"));
		if(datatype == CODON){
			if(_stricmp(str, "nonsynonymous") == 0) SetM3Model();
			else if(_stricmp(str, "none") == 0){
				SetNumRateCats(1, false);
				rateHetType = NONE;
				}
			else if(_stricmp(str, "gamma") == 0) throw ErrorException("Gamma rate heterogeneity cannot be used with codon models.\n     Try ratehetmodel = nonsynonymous to allow dN/dS variation across sites");
			else throw(ErrorException("Unknown setting for ratehetmodel: %s\n\t(options for codon datatype are: nonsynonymous, none)", str));
			}
		else{			
			if(_stricmp(str, "gamma") == 0) SetGammaRates();
			else if(_stricmp(str, "gammafixed") == 0){
				SetGammaRates();
				SetFixedAlpha();
				}
			else if(_stricmp(str, "flex") == 0) SetFlexRates();
			else if(_stricmp(str, "none") == 0){
				SetNumRateCats(1, false);
				rateHetType = NONE;
				}
			else throw(ErrorException("Unknown setting for ratehetmodel: %s\n\t(options are for nucleotide or aminoacid data are: gamma, gammafixed, flex, none)", str));
			}
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
	int effectiveModels;//this is the number of models with different Q matrices
						//it does not include things like gamma or flex rates,
						//in which it is only the overall rate that varies

	vector<FLOAT_TYPE*> stateFreqs;
	vector<FLOAT_TYPE*> relNucRates;
	vector<FLOAT_TYPE*> omegas;
	vector<FLOAT_TYPE*> omegaProbs;

	bool eigenDirty;
	FLOAT_TYPE *blen_multiplier;
	
	FLOAT_TYPE rateMults[20];
	FLOAT_TYPE rateProbs[20];
	
	FLOAT_TYPE *alpha;
	FLOAT_TYPE *propInvar;

	//variables used for the eigen process if nst=6
	int *iwork, *indx;
	FLOAT_TYPE **eigvals, *eigvalsimag, ***eigvecs, ***inveigvecs, **teigvecs, *work, *temp, *col, **c_ijk, *EigValexp, *EigValderiv, *EigValderiv2;
	FLOAT_TYPE ***qmat, ***pmat;
	FLOAT_TYPE ***tempqmat;
	
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
		//DEBUG - we should probably move this out of here
		CreateModelFromSpecification(0);
		}

	void CalcMutationProbsFromWeights();
	BaseParameter *SelectModelMutation();
	int PerformModelMutation();
	void CreateModelFromSpecification(int);

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
	void CopyEigenVariables(const Model *from);
	void SetModel(FLOAT_TYPE *model_string);
	void OutputPaupBlockForModel(ofstream &, const char *) const;
	void FillPaupBlockStringForModel(string &str, const char *treefname) const;
	void OutputGarliFormattedModel(ostream &) const;
	void FillGarliFormattedModelString(string &s) const;
	void OutputBinaryFormattedModel(OUTPUT_CLASS &) const;
	void ReadGarliFormattedModelString(string &);

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
	FLOAT_TYPE Rates(int r) const { return *relNucRates[r];}
	int NRateCats() const {return modSpec.numRateCats;}
	FLOAT_TYPE *GetRateMults() {return rateMults;}
	FLOAT_TYPE Alpha() const {return *alpha;}
	FLOAT_TYPE PropInvar() const { return *propInvar;}
	bool NoPinvInModel() const { return ! (modSpec.includeInvariantSites);}
	FLOAT_TYPE MaxPinv() const{return maxPropInvar;}
	int NStates() const {return nstates;}
	int NumMutatableParams() const {return (int) paramsToMutate.size();}

	//Setting things
	void SetDefaultModelParameters(const SequenceData *data);
	void SetRmat(FLOAT_TYPE *r, bool checkValidity){
		assert(modSpec.IsAminoAcid() == false);
		if(checkValidity == true){
			if(nst==1){
				if((FloatingPointEquals(r[0], r[1], 1.0e-5) &&
					FloatingPointEquals(r[1], r[2], 1.0e-5) &&
					FloatingPointEquals(r[2], r[3], 1.0e-5) &&
					FloatingPointEquals(r[3], r[4], 1.0e-5) &&
					FloatingPointEquals(r[4], r[5], 1.0e-5)) == false)
					throw(ErrorException("Config file specifies ratematrix = 1rate, but starting model has nonequal rates!\n"));
				}
			if(nst==2){
				if((FloatingPointEquals(r[0], r[2], 1.0e-5) &&
					FloatingPointEquals(r[2], r[3], 1.0e-5) &&
					FloatingPointEquals(r[1], r[4], 1.0e-5) &&
					FloatingPointEquals(r[3], r[5], 1.0e-5)) == false)
					throw(ErrorException("Config file specifies ratematrix = 2rate, but starting model parameters do not match!\n"));
				}
			}
		if(FloatingPointEquals(r[5], ONE_POINT_ZERO, 1.0e-5) == false){
			//if an alternate GTR paramterization is used in which GT != 1.0, rescale the rates
			for(int i=0;i<5;i++)
				r[i] /= r[5];
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
				if(modSpec.IsEqualStateFrequencies() && (FloatingPointEquals(b[0], b[1], 1.0e-5) && FloatingPointEquals(b[1], b[2], 1.0e-5)) == false) 
					throw(ErrorException("Config file specifies equal statefrequencies,\nbut starting model has nonequal frequencies!\n"));
				if(modSpec.IsEmpiricalStateFrequencies()) 
					throw(ErrorException("Config file specifies empirical statefrequencies,\nbut starting model specifies frequencies!\n"));
				}
			}
		FLOAT_TYPE freqTot = 0.0;
		for(int i=0;i<nstates;i++){
			*stateFreqs[i]=b[i];
			freqTot += *stateFreqs[i];
			}
		if(FloatingPointEquals(freqTot, ONE_POINT_ZERO, 1.0e-5) == false)
			throw(ErrorException("State frequencies do not appear to add up to 1.0!\n"));
		eigenDirty=true;
		}

	void SetFlexRates(FLOAT_TYPE *rates, FLOAT_TYPE *probs){
		if(modSpec.IsFlexRateHet() == false) throw ErrorException("Error: Flex rate values specified in start file,\nbut ratehetmodel is not flex in conf file.");
		for(int r=0;r<NRateCats();r++){
			rateMults[r]=rates[r];
			rateProbs[r]=probs[r];
			}		
		}

	FLOAT_TYPE FlexRate(int which){
		assert(which < NRateCats());
		return rateMults[which];
		}

	FLOAT_TYPE FlexProb(int which){
		assert(which < NRateCats());
		return rateProbs[which];
		}

	void SetFlexRate(int which, FLOAT_TYPE val){
		assert(which < NRateCats());
		rateMults[which] = val;
		NormalizeRates(which);
		eigenDirty = true;
		}

	void SetFlexProb(int which, FLOAT_TYPE val){
		assert(which < NRateCats());
		rateProbs[which] = val;
		NormalizeRates(which);
		eigenDirty = true;
		}

	void SetOmegas(const FLOAT_TYPE *rates, const FLOAT_TYPE *probs){
		FLOAT_TYPE tot=0.0;
		for(int r=0;r<NRateCats();r++){
			*omegas[r]=rates[r];
			*omegaProbs[r]=probs[r];
			tot += *omegaProbs[r];
			}
		if(FloatingPointEquals(tot, ONE_POINT_ZERO, 1.0e-5) == false) throw ErrorException("omega category proportions add up to %f, not 1.0.", tot);
		eigenDirty = true;
		}

	void SetOmega(int which, FLOAT_TYPE val){
		assert(which < NRateCats());
		*omegas[which] = val;
		eigenDirty = true;
		}

	FLOAT_TYPE Omega(int which){
		assert(which < NRateCats());
		return *omegas[which];
		}

	void SetOmegaProb(int which, FLOAT_TYPE val){
		assert(which < NRateCats());
		*omegaProbs[which] = val;

		FLOAT_TYPE newTot = 1.0 - *omegaProbs[which];
		FLOAT_TYPE oldTot = 0.0;
		for(int i=0;i<NRateCats();i++)
			if(i != which) oldTot += *omegaProbs[i];
		for(int i=0;i<NRateCats();i++)
			if(i != which) *omegaProbs[i] *= newTot / oldTot;
		newTot = 0.0;
		for(int i=0;i<NRateCats();i++) newTot += *omegaProbs[i];
		assert(FloatingPointEquals(newTot, ONE_POINT_ZERO, 1.0e-5));
		eigenDirty = true;
		}

	FLOAT_TYPE OmegaProb(int which){
		assert(which < NRateCats());
		return *omegaProbs[which];
		}

	void SetAlpha(FLOAT_TYPE a, bool checkValidity){
		if(checkValidity == true)
			if(modSpec.numRateCats==1) throw(ErrorException("Config file specifies ratehetmodel = none, but starting model contains alpha!\n"));
		*alpha=a;
		DiscreteGamma(rateMults, rateProbs, *alpha);
		//This is odd, but we need to call normalize rates here if we are just using a gamma distrib to get starting rates for 
		//flex.  Flex expects that the rates will be normalized including pinv elsewhere
		if(modSpec.IsFlexRateHet()) NormalizeRates();
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

	void CheckAndCorrectRateOrdering(){
		assert(NRateCats() > 1);
		if(modSpec.IsNonsynonymousRateHet()){
			for(int f=0;f<NRateCats()-1;f++){
				if(*omegas[f] > *omegas[f+1]){
					//outman.UserMessage("prevented: %f %f", *omegas[f], *omegas[f+1]); 
					FLOAT_TYPE dum = *omegas[f+1];
					*omegas[f+1] = *omegas[f];
					*omegas[f] = dum;
					dum = *omegaProbs[f+1];
					*omegaProbs[f+1] = *omegaProbs[f];
					*omegaProbs[f] = dum;
					}
				}
			}
		else if(modSpec.IsFlexRateHet()){
			for(int f=0;f<NRateCats()-1;f++){
				if(rateMults[f] > rateMults[f+1]){
					FLOAT_TYPE dum = rateMults[f+1];
					rateMults[f+1] = rateMults[f];
					rateMults[f] = dum;
					dum = rateProbs[f+1];
					rateProbs[f+1] = rateProbs[f];
					rateProbs[f] = dum;
					}
				}
			}
		else assert(0);
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
		assert(FloatingPointEquals(sum, ONE_POINT_ZERO, 1.0e-5));
#endif
		}

	void NormalizeRates(int toRemainConstant = -1){
		//optionally, pass the number of one of the rate/prob pairs to hold constant

		FLOAT_TYPE sum=0.0;

		for(int i=0;i<NRateCats();i++){
			if(i != toRemainConstant) sum += rateProbs[i];
			}

		//pinv is figured into the normalization here, but it won't be changed itself.
		if(NoPinvInModel()==false){
			sum = sum / (FLOAT_TYPE)(1.0-*propInvar);
			}

		if(toRemainConstant > -1) sum /= (ONE_POINT_ZERO - rateProbs[toRemainConstant]);
		for(int i=0;i<NRateCats();i++)	{
			if(i != toRemainConstant) rateProbs[i] /= sum;
			}

		sum=0.0;
		for(int i=0;i<NRateCats();i++){
			if(i != toRemainConstant) sum += rateMults[i]*rateProbs[i];
			}
		if(toRemainConstant > -1) sum /= (ONE_POINT_ZERO - (rateMults[toRemainConstant] * rateProbs[toRemainConstant]));
		for(int i=0;i<NRateCats();i++){
			if(i != toRemainConstant) rateMults[i] /= sum;
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
		assert(FloatingPointEquals(sum, 1.0, 1.0e-5));

#endif
		}

	const FLOAT_TYPE *GetRateProbs() {
		//this is silly, but use the rateProbs as a holder to return the omegaProbs, which are in a vector of double pointers
		if(modSpec.IsNonsynonymousRateHet())
			for(int i=0;i<NRateCats();i++)
				rateProbs[i] = *omegaProbs[i];

#ifndef NDEBUG
		FLOAT_TYPE sum=0.0;
		for(int i=0;i<NRateCats();i++){
			sum += rateProbs[i];
			}
		sum+=*propInvar;
		assert(fabs(1.0-sum) < .001);
#endif		
		return rateProbs;
		}
	};
	
	class ModelSet{//this is a set of models that are applied to a _single_ set of sites
					//i.e., model mixtures, although only ones in which the models have separate
					//Q matrices, eigen variables, etc (not just rates, so gamma and flex rates don't count)
		int numModels;
		vector<Model*> mods;
		vector<FLOAT_TYPE> modelProbs;
		ModelSet(){
			numModels = 1;
			for(int i=0;i<numModels;i++){
				Model *mod = new Model();
				mod->CreateModelFromSpecification(i);
				}
			}
		ModelSet(const ModelSet &m){
			CopyModelSet(m);
			}
		void CopyModelSet(const ModelSet &m){
			int num = 0;
			for(vector<Model*>::const_iterator modit = m.mods.begin();modit != m.mods.end();modit++){
				Model *mod;
				if(mods.empty()){
					mod = new Model();
					mods.push_back(mod);
					}
				else mod = mods[num];
				mod->CopyModel(*modit);
				num++;
				}
			}
		};

typedef void (Model::*SetParamFunc) (int, FLOAT_TYPE);
#define CALL_SET_PARAM_FUNCTION(object, ptrToMember) ((object).*(ptrToMember))

#endif
