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


#ifndef ADAPTION_H
#define ADAPTION_H

#include <fstream>

using namespace std;

#include "configoptions.h"
#include "hashdefines.h"

class Adaptation{
	public:
	//here are all of the scalars:
	//4 ints, 1 bool, 23 doubles
	int intervalLength; 
	int intervalsToStore;
	double lastgenscore;
	double laststepscore;

	double improveOverStoredIntervals;
	
	bool reset;

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
	
	double recTopImproveSize;
	double exNNIprob;
	double exlimSPRprob;

	double topoMutateProb;
	double randNNIprob;
	double randSPRprob;
	double limSPRprob;
	
	double modelMutateProb;

	int limSPRrange;

	//the arrays. All will be of length intervalsToStore
	double *improvetotal;

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

	double *randRecom;
	int *randRecomnum;
	
	double *bipartRecom;
	int *bipartRecomnum;

	double *onlyBrlen;
	int *onlyBrlennum;

	double *anyModel;
	int *anyModelnum;
	
#ifdef MPI_VERSION
	double *fromRemoteSubtree;
	double *fromRemoteNonSubtree;
	int *bestFromRemoteNum;
	double *bestFromRemote;
#endif

	Adaptation(const GeneralGamlConfig *gc);
	~Adaptation();
	public:
	void SetChangeableVariablesFromConfAfterReadingCheckpoint(const GeneralGamlConfig *gc);
	void PrepareForNextInterval();
	void UpdateProbs();
	void OutputProbs(ofstream &plog, int gen);
	void BeginProbLog(ofstream &plot);
	bool ReducePrecision(){
		if(branchOptPrecision==minOptPrecision || numPrecReductions == 0) return false;
		if(topoMutateProb > .01 || topoWeight==0.0){
			//changing this to a linear reduction in prec.  Geometric was too fast
			//branchOptPrecision*=precReductionFactor;
			branchOptPrecision -= precReductionFactor;
			if((branchOptPrecision - 0.001) < minOptPrecision) branchOptPrecision=minOptPrecision;
			return true;
			}
		else return false;
		}
	void WriteToCheckpoint(ofstream &out);
	void ReadFromCheckpoint(ifstream &in);
	
};

#endif
