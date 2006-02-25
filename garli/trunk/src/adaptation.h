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

#ifndef ADAPTION_H
#define ADAPTION_H

#include <fstream>

using namespace std;

#include "configoptions.h"
#include "hashdefines.h"

class Adaptation{
	public:
	int intervalLength;
	int intervalsToStore;
	double lastgenscore;
	double laststepscore;

	double improveOverStoredIntervals;
	double *improvetotal;
	
	bool reset;

	double *randNNI;
	int *randNNInum;
	
	double *exNNI;
	int *exNNInum;

	double *randSPR;
	int *randSPRnum;
	
	double *limSPR;
	int *limSPRnum;

	double *exlimSPR;
	int *exlimSPRnum;

#ifdef GANESH
	double *randPECR;
	int *randPECRnum;
#endif 


	double *taxonSwap;
	int *taxonSwapnum;
	
	double *randRecom;
	int *randRecomnum;
	
	double *bipartRecom;
	int *bipartRecomnum;

	double *onlyBrlen;
	int *onlyBrlennum;

	double *anyModel;
	int *anyModelnum;
	
	double *fromRemoteSubtree;
	double *fromRemoteNonSubtree;
	
	double *slopes;
	
	double startOptPrecision;
	double branchOptPrecision;
	double minOptPrecision;
	double precReductionFactor;
	int numPrecReductions;
	
	double topoWeight;
	double modWeight;
	double brlenWeight;
	
	double randNNIweight;
	double origRandNNIweight;
	double randSPRweight;
	double limSPRweight;
#ifdef GANESH
	double randPECRweight; 
#endif
	
	double recTopImproveSize;
	double exNNIprob;
	double exlimSPRprob;

	double topoMutateProb;
	double randNNIprob;
	double randSPRprob;
#ifdef GANESH
	double randPECRprob;
#endif
	double limSPRprob;
	double taxonSwapprob;
	
	double modelMutateProb;

	int limSPRrange;

#ifdef MPI_VERSION
	int *bestFromRemoteNum;
	double *bestFromRemote;
#endif

	Adaptation(const GeneralGamlConfig *gc);

	public:
	void PrepareForNextInterval();
	void UpdateProbs();
	void OutputProbs(ofstream &plog, int gen);
	void BeginProbLog(ofstream &plot);
	bool ReducePrecision(){
		if(branchOptPrecision==minOptPrecision) return false;
		if(topoMutateProb > .1 || topoWeight==0.0){
			branchOptPrecision*=precReductionFactor;
			if(branchOptPrecision < minOptPrecision) branchOptPrecision=minOptPrecision;
			return true;
			}
		else return false;
		}

};

#endif
