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



#ifndef CONFIGOPTIONS_H
#define CONFIGOPTIONS_H

#include <string>

using std::string;

#include "hashdefines.h"

class GeneralGamlConfig{
	public:
	//these options will be the same regardless of whether a population is master or remote
	unsigned logevery;
	unsigned saveevery;
	double megsClaMemory;
	double availableMemory;
	int randseed;

	string datafname;
	string method;
	string ofprefix;
	string streefname;
	string constraintfile;
	bool refineStart;
	bool outputTreelog;
	bool outputMostlyUselessFiles;
	bool outputPhylipTree;
	bool restart;
	bool checkpoint;

	bool enforceTermConditions;
	unsigned lastTopoImproveThresh;
	double improveOverStoredIntervalsThresh;
	double significantTopoChange;

	//model settings
	string stateFrequencies; //equal, estimate, emprical, fixed
	string rateMatrix;		 //6rate, 2rate, 1rate, fixed, custom(
	string proportionInvariant; //none, fixed, estimate
	string rateHetModel;			//gamma, gammafixed, flex, none

//	bool dontInferProportionInvariant;
//	bool useflexrates;
	unsigned numRateCats;	

	//all of the following options can vary between master and remote
	//general population stuff
	unsigned nindivs;
	unsigned holdover;
	double selectionIntensity;
	double holdoverPenalty;
	unsigned stopgen;
	//unsigned stopgen;
	unsigned stoptime;

	double startOptPrec;
	double minOptPrec;
	double precReductionFactor;
	int numPrecReductions;

	double treeRejectionThreshold;

	//parameters affecting proportion of mutations
	double topoWeight;
		double randNNIweight;
		double randSPRweight;
		double limSPRweight;
        double randPECRweight;
	double modWeight;
	double brlenWeight;

	unsigned intervalLength;
	unsigned intervalsToStore;

	unsigned searchReps;
	unsigned bootstrapReps;
	bool inferInternalStateProbs;

	//parameters affecting other details of mutations				
	double meanBrlenMuts;
	unsigned gammaShapeBrlen;
	unsigned gammaShapeModel;
	unsigned limSPRrange;		
	double uniqueSwapBias;
	double distanceSwapBias;

	//perturbation parameters
	int pertType;			
	double pertThresh;
	int minPertInterval;
	int maxPertsNoImprove;
	bool restartAfterAbandon;
	int gensBeforeRestart;
	
	double ratchetProportion;
	double ratchetOffThresh;
	int ratchetMaxGen;
	
	int nniTargetAccepts;
	int nniMaxAttempts;
	
	int numSprCycles;
	int sprPertRange;

	//the number of seconds between remote tree sends
	double sendInterval;
	
 	// methods
	GeneralGamlConfig();
	int Read(const char*, bool isMaster=false);
	int Serialize(char**, int*) const;
	int Deserialize(char*, int);
	bool operator==(const GeneralGamlConfig&) const;
	};
		
class MasterGamlConfig: public GeneralGamlConfig{
	public:
	//parallel behavior parameters-stored in pop->paraMan on master only
	double startUpdateThresh;
	double minUpdateThresh;
	double updateReductionFactor;
		
	int subtreeInterval;
	double subtreeStartThresh;
	int minSubtreeSize;
	int targetSubtreeSize;
	double orphanFactor;
	
	int maxRecomIndivs;
/*	
	int pertType;			
	double pertThresh;
	double pertAmount;
	int maxPertsNoImprove;
	
	double ratchetProportion;
	double ratchetOffThresh;
	int ratchetMaxGen;
	
	double nniAcceptThresh;
	int numSprCycles;
	int sprPertRange;
*/
	int bootstrapReps;
	double bootTermThresh;

 	// methods
	MasterGamlConfig();
	int Read(const char*, bool isMaster=false);
	};

#endif

