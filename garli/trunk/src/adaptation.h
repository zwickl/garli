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


#ifndef ADAPTION_H
#define ADAPTION_H

#include <fstream>

using namespace std;

#include "configoptions.h"
#include "hashdefines.h"

class MFILE;

class Adaptation{
	public:
	//here are all of the scalars:
	//4 ints, 1 bool, 23 FLOAT_TYPEs
	unsigned intervalLength; 
	unsigned intervalsToStore;
	FLOAT_TYPE lastgenscore;
	FLOAT_TYPE laststepscore;

	FLOAT_TYPE improveOverStoredIntervals;
	
	bool reset;

	FLOAT_TYPE startOptPrecision;
	FLOAT_TYPE branchOptPrecision;
	FLOAT_TYPE minOptPrecision;
	FLOAT_TYPE precReductionFactor;
	int numPrecReductions;
	
	FLOAT_TYPE topoWeight;
	FLOAT_TYPE modWeight;
	FLOAT_TYPE brlenWeight;
	
	FLOAT_TYPE randNNIweight;
	FLOAT_TYPE origRandNNIweight;
	FLOAT_TYPE randSPRweight;
	FLOAT_TYPE limSPRweight;
	
	FLOAT_TYPE recTopImproveSize;
	FLOAT_TYPE exNNIprob;
	FLOAT_TYPE exlimSPRprob;

	FLOAT_TYPE topoMutateProb;
	FLOAT_TYPE randNNIprob;
	FLOAT_TYPE randSPRprob;
	FLOAT_TYPE limSPRprob;
	
	FLOAT_TYPE modelMutateProb;

	unsigned limSPRrange;

	//the arrays. All will be of length intervalsToStore
	FLOAT_TYPE *improvetotal;

	FLOAT_TYPE *randNNI;
	int *randNNInum;
	
	FLOAT_TYPE *exNNI;
	int *exNNInum;

	FLOAT_TYPE *randSPR;
	int *randSPRnum;
	
	FLOAT_TYPE *limSPR;
	int *limSPRnum;

	FLOAT_TYPE *exlimSPR;
	int *exlimSPRnum;

	FLOAT_TYPE *randRecom;
	int *randRecomnum;
	
	FLOAT_TYPE *bipartRecom;
	int *bipartRecomnum;

	FLOAT_TYPE *onlyBrlen;
	int *onlyBrlennum;

	FLOAT_TYPE *anyModel;
	int *anyModelnum;
	
#ifdef MPI_VERSION
	FLOAT_TYPE *fromRemoteSubtree;
	FLOAT_TYPE *fromRemoteNonSubtree;
	int *bestFromRemoteNum;
	FLOAT_TYPE *bestFromRemote;
#endif

	Adaptation(const GeneralGamlConfig *gc);
	~Adaptation();
	public:
	void SetChangeableVariablesFromConfAfterReadingCheckpoint(const GeneralGamlConfig *gc);
	void PrepareForNextInterval();
	void UpdateProbs();
	void OutputProbs(ofstream &plog, int gen);
	void BeginProbLog(ofstream &plot, int gen);
	bool ReducePrecision(){
		if(branchOptPrecision==minOptPrecision || numPrecReductions == 0) return false;
		if(topoMutateProb > .01 || topoWeight==0.0){
			//changing this to a linear reduction in prec.  Geometric was too fast
			//branchOptPrecision*=precReductionFactor;
			branchOptPrecision -= precReductionFactor;
			if((branchOptPrecision < minOptPrecision) || (branchOptPrecision - minOptPrecision) < (precReductionFactor / 2.0)) branchOptPrecision=minOptPrecision;
			return true;
			}
		else return false;
		}

	void WriteToCheckpoint(OUTPUT_CLASS &) const;
	void ReadFromCheckpoint(FILE *);

};

#endif
